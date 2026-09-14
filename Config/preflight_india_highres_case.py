#!/usr/bin/env python3
"""Fail-fast checks for a staged India high-resolution HydroPol2D case."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform_bounds


STATIC = (
    "DEM.tif", "LULC.tif", "SOIL.tif", "DTB.tif", "GW_table.tif", "LAI.tif", "Albedo.tif",
    "Initial_Soil_Moisture.tif", "RiverWidths.tif", "RiverDepths.tif", "D4_flow_direction.tif", "D4_idx_facc.tif",
)


def active_mask(dataset: rasterio.DatasetReader, values: np.ndarray) -> np.ndarray:
    if dataset.nodata is None or (isinstance(dataset.nodata, float) and np.isnan(dataset.nodata)):
        return np.isfinite(values)
    return np.isfinite(values) & (values != dataset.nodata)


def float32_raster_tolerance(reference: np.ndarray) -> np.ndarray:
    """Return one float32 unit in the last place, with a 0.1-mm floor."""
    return np.maximum(1e-4, np.abs(np.spacing(reference.astype("float32"))).astype(float))


def source_tree_sha256(model_root: Path) -> str:
    digest = hashlib.sha256()
    roots = [model_root / "HydroPol2D_V115.m", model_root / "HydroPol2D_Functions", model_root / "Config"]
    files = []
    for root in roots:
        if root.is_file():
            files.append(root)
        elif root.is_dir():
            files.extend(
                path for path in root.rglob("*")
                if path.is_file() and path.suffix.lower() in {".m", ".py", ".json", ".sbatch"}
            )
    for path in sorted(files):
        digest.update(str(path.relative_to(model_root)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_state(model_root: Path) -> tuple[str | None, bool | None]:
    if not (model_root / ".git").exists():
        return None, None
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=model_root, text=True).strip()
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=model_root, text=True)
    return commit, not bool(status.strip())


def check_rainfall(case: Path, start: datetime, end: datetime, dem_crs, dem_bounds) -> dict:
    rainfall = case / "Forcing" / "Rainfall"
    time, paths = start, []
    while time <= end:
        path = rainfall / f"IMERG_30min_mmhr_IndiaRegion_{time:%Y_%m_%d_%H_%M}.tif"
        if not path.is_file():
            raise FileNotFoundError(path)
        paths.append(path)
        time += timedelta(minutes=30)
    for path in (paths[0], paths[-1]):
        with rasterio.open(path) as source:
            if source.count != 1:
                raise RuntimeError(f"Rainfall raster must have one band: {path}")
            rain_bounds = transform_bounds(source.crs, dem_crs, *source.bounds, densify_pts=21)
            if not (rain_bounds[0] <= dem_bounds.left and rain_bounds[1] <= dem_bounds.bottom
                    and rain_bounds[2] >= dem_bounds.right and rain_bounds[3] >= dem_bounds.top):
                raise RuntimeError(f"Rainfall raster does not cover the full model domain: {path}")
    return {"maps": len(paths), "first": str(paths[0]), "last": str(paths[-1])}


def check_etp(case: Path, start: datetime, end: datetime) -> dict:
    path = case / "Forcing" / "Evapotranspiration" / "ETP_input_data.xlsx"
    if not path.is_file():
        raise FileNotFoundError(path)
    table = pd.read_excel(path, header=None)
    observed = pd.DatetimeIndex(
        value for value in table.iloc[:, 1] if isinstance(value, (datetime, pd.Timestamp))
    ).normalize()
    expected = pd.date_range(start.date(), end.date(), freq="D")
    missing = expected.difference(observed.unique())
    if len(missing):
        raise RuntimeError(f"ETP forcing is missing {len(missing)} day(s), beginning with {missing[0].date()}")
    return {"days": len(expected), "first": str(expected[0].date()), "last": str(expected[-1].date())}


def check_d4(direction: np.ndarray, active: np.ndarray) -> dict:
    valid_codes = {0, 1, 2, 3, 4}
    codes = set(np.unique(direction[active & np.isfinite(direction)]).astype(int).tolist())
    if not codes.issubset(valid_codes):
        raise RuntimeError(f"Unexpected D4 direction codes: {sorted(codes - valid_codes)}")
    invalid_receivers = 0
    for code, (dr, dc) in {1: (-1, 0), 2: (0, 1), 3: (1, 0), 4: (0, -1)}.items():
        rows, cols = np.where(active & (direction == code))
        rr, cc = rows + dr, cols + dc
        inside = (rr >= 0) & (rr < active.shape[0]) & (cc >= 0) & (cc < active.shape[1])
        valid = np.zeros_like(inside)
        valid[inside] = active[rr[inside], cc[inside]]
        invalid_receivers += int((~valid).sum())
    if invalid_receivers:
        raise RuntimeError(f"D4 has {invalid_receivers:,} links leaving active terrain without an outlet code")
    return {"codes": sorted(codes), "outlet_cells": int(np.count_nonzero(active & (direction == 0)))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--model-root", required=True, type=Path)
    parser.add_argument("--allow-dirty-model", action="store_true")
    args = parser.parse_args()
    manifest_path = args.case / "case_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "staged":
        raise RuntimeError(f"Case status is {manifest.get('status')!r}; run the stage step first")
    if manifest["sources"]["water_table"]["product"].lower().find("fan") < 0:
        raise RuntimeError("The staged case does not document Fan water-table depth")
    source_manifest = args.case / manifest["source_download_manifest"]
    if file_sha256(source_manifest) != manifest["source_download_manifest_sha256"]:
        raise RuntimeError("Archived source-download manifest checksum does not match the case manifest")
    for relative, expected_hash in manifest.get("staged_static_sha256", {}).items():
        if file_sha256(args.case / relative) != expected_hash:
            raise RuntimeError(f"Staged static checksum changed: {relative}")
    commit, clean = git_state(args.model_root)
    locked = manifest["model_source"]
    if commit is not None and commit != locked["commit"]:
        raise RuntimeError(f"Model commit changed: case={locked['commit']}, current={commit}")
    current_tree_hash = source_tree_sha256(args.model_root)
    if current_tree_hash != locked["source_tree_sha256"]:
        raise RuntimeError("Executable model source differs from the source-tree hash stored in the case manifest")
    if clean is False and not args.allow_dirty_model:
        raise RuntimeError("Model working tree is dirty; commit the mass-balance correction or pass --allow-dirty-model")

    static = args.case / "Static"
    for name in STATIC:
        if not (static / name).is_file():
            raise FileNotFoundError(static / name)
    with rasterio.open(static / "DEM.tif") as ref:
        dem = ref.read(1).astype(float)
        active = active_mask(ref, dem)
        resolution = abs(ref.transform.a)
        if not np.isclose(resolution, manifest["resolution_m"], atol=1e-5):
            raise RuntimeError(f"DEM resolution {resolution} differs from manifest {manifest['resolution_m']}")
        arrays = {}
        for name in STATIC[1:]:
            with rasterio.open(static / name) as source:
                if source.shape != ref.shape or source.crs != ref.crs or not source.transform.almost_equals(ref.transform):
                    raise RuntimeError(f"Static grid mismatch: {name}")
                arrays[name] = source.read(1).astype(float)
        rain_summary = check_rainfall(
            args.case, datetime.fromisoformat(manifest["simulation_start"]),
            datetime.fromisoformat(manifest["simulation_end"]), ref.crs, ref.bounds,
        )
        etp_summary = check_etp(
            args.case, datetime.fromisoformat(manifest["simulation_start"]),
            datetime.fromisoformat(manifest["simulation_end"]),
        )
        transform = ref.transform

    minimum_dtb = float(manifest["dtb_conditioning"]["minimum_effective_depth_m"])
    if not np.isfinite(arrays["DTB.tif"][active]).all() or (arrays["DTB.tif"][active] < minimum_dtb - 1e-6).any():
        raise RuntimeError(f"Depth to bedrock must be finite and at least {minimum_dtb:g} m in every active cell")
    for name in ("GW_table.tif", "Initial_Soil_Moisture.tif"):
        if not np.isfinite(arrays[name][active]).all():
            raise RuntimeError(f"{name} contains missing values in active cells")
    if (arrays["Initial_Soil_Moisture.tif"][active] < 0).any():
        raise RuntimeError("Initial soil storage contains negative values")
    floor = dem - arrays["DTB.tif"]
    rounding_tolerance = float32_raster_tolerance(dem)
    if (
        (arrays["GW_table.tif"][active] < floor[active] - rounding_tolerance[active]).any()
        or (arrays["GW_table.tif"][active] > dem[active] + rounding_tolerance[active]).any()
    ):
        raise RuntimeError("Initial groundwater head falls outside the soil/bedrock column")
    river = np.isfinite(arrays["RiverWidths.tif"]) & (arrays["RiverWidths.tif"] > 0)
    if not river.any() or (arrays["RiverDepths.tif"][river] <= 0).any():
        raise RuntimeError("River width/depth rasters do not define a valid network")
    if (arrays["RiverWidths.tif"][river] > resolution + 1e-6).any():
        raise RuntimeError("Neal channel width exceeds the model-cell resolution")
    d4_summary = check_d4(arrays["D4_flow_direction.tif"], active)

    gauge_path = args.case / "Forcing" / "Observed_Gauges" / "camels_india_gauges_2km.csv"
    gauges = pd.read_csv(gauge_path, dtype={"Gauge": str})
    rows, cols = rasterio.transform.rowcol(transform, gauges["Easting_m"], gauges["Northing_m"])
    rows, cols = np.asarray(rows), np.asarray(cols)
    if not ((rows >= 0) & (rows < active.shape[0]) & (cols >= 0) & (cols < active.shape[1])).all():
        raise RuntimeError("A CAMELS gauge falls outside the model grid")
    if not river[rows, cols].all():
        raise RuntimeError("A CAMELS gauge was not snapped to an active river cell")

    observed = pd.read_csv(args.case / "Validation" / "CAMELS_IND_observed_streamflow.csv")
    coverage = {gauge: float(observed[gauge].notna().mean()) for gauge in gauges["Gauge"]}
    summary = {
        "case": str(args.case),
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "model_commit": commit,
        "model_working_tree_clean": clean,
        "model_source_tree_sha256": current_tree_hash,
        "resolution_m": resolution,
        "shape": list(active.shape),
        "active_cells": int(active.sum()),
        "river_cells": int(river.sum()),
        "d4": d4_summary,
        "hydrobathydem_warnings": [
            item for item in manifest.get("hydrobathydem_qa", []) if item.get("status") != "ok"
        ],
        "rainfall": rain_summary,
        "meteorology": etp_summary,
        "gauges": {"count": len(gauges), "all_on_river": True, "observed_daily_coverage_fraction": coverage},
        "validation_note": "CAMELS streamflow is validation-only and is not used as forcing or calibration.",
        "status": "passed",
    }
    output = args.case / "preflight_india_highres.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
