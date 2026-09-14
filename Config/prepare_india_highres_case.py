#!/usr/bin/env python3
"""Prepare reproducible 2 km, 90 m, and nested 30 m India HydroPol2D cases.

The script has two intentionally separate steps:

``domain`` selects CAMELS-IND catchments and writes HydroBathyDEM configs.
``stage`` assembles a runnable HydroPol2D case only after source data exist.

No national 2-km raster is accepted as a high-resolution static input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timedelta, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.warp import reproject


HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "india_highres_cases.json"
SOIL_THETA_R = np.array(
    [0.068, 0.089, 0.075, 0.095, 0.089, 0.065, 0.067, 0.078, 0.100, 0.034, 0.049, 0.045],
    dtype="float32",
)
SOIL_LAYER_M = np.array([0.07, 0.21, 0.72, 1.89], dtype="float32")


def load_catalog(path: Path) -> dict:
    catalog = json.loads(path.read_text())
    if catalog.get("schema_version") != 1:
        raise ValueError(f"Unsupported catalog schema: {catalog.get('schema_version')}")
    return catalog


def parse_time(value: str) -> datetime:
    return datetime.fromisoformat(value)


def simulation_start(case: dict, defaults: dict) -> datetime:
    return parse_time(case["event_start"]) - timedelta(days=int(defaults["warmup_days"]))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def git_lock(model_root: Path) -> dict:
    lock = {"source_tree_sha256": source_tree_sha256(model_root)}
    if not (model_root / ".git").exists():
        return {
            **lock,
            "commit": None,
            "working_tree_clean": None,
            "note": "Deployment copy has no Git metadata; source_tree_sha256 fixes the executable source exactly.",
        }
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=model_root, text=True).strip()
    diff = subprocess.check_output(["git", "diff", "--binary", "HEAD"], cwd=model_root)
    untracked = subprocess.check_output(
        ["git", "ls-files", "--others", "--exclude-standard"], cwd=model_root, text=True
    ).splitlines()
    return {
        **lock,
        "commit": commit,
        "tracked_patch_sha256": hashlib.sha256(diff).hexdigest(),
        "working_tree_clean": not diff and not untracked,
        "untracked_files": sorted(untracked),
        "note": "Commit, tracked patch hash, and source-tree hash identify this candidate until the mass-balance fix is committed.",
    }


def request_memory_gib(active_cells: int, defaults: dict) -> int:
    memory = defaults["memory"]
    estimate = float(memory["base_gib"]) + active_cells * float(memory["bytes_per_active_cell"]) / 1024**3
    for bucket in memory["request_buckets_gib"]:
        if estimate <= bucket:
            return int(bucket)
    raise RuntimeError(f"Estimated memory {estimate:.1f} GiB exceeds the configured 128 GiB ceiling")


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def select_camels(case: dict, camels_root: Path) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, Path, Path]:
    root = camels_root / "shapefiles_catchment" / case["camels_subdir"]
    catchment_path = root / case["catchment_file"]
    station_path = root / case["station_file"]
    if not catchment_path.is_file():
        catchment_path = camels_root / "shapefiles_catchment" / "catchments.shp"
    if not station_path.is_file():
        station_path = camels_root / "shapefiles_catchment" / "gauge_stations.shp"
    catchments = gpd.read_file(catchment_path)
    stations = gpd.read_file(station_path)
    for frame in (catchments, stations):
        if "gauge_id" not in frame:
            raise ValueError("CAMELS-IND vector is missing gauge_id")
        frame["gauge_id"] = frame["gauge_id"].astype(str)
    gauge_ids = set(case["gauge_ids"])
    catchments = catchments[catchments["gauge_id"].isin(gauge_ids)].copy()
    stations = stations[stations["gauge_id"].isin(gauge_ids)].copy()
    if set(catchments["gauge_id"]) != gauge_ids or set(stations["gauge_id"]) != gauge_ids:
        raise RuntimeError(f"Not all requested CAMELS gauges were found: {sorted(gauge_ids)}")
    return catchments, stations, catchment_path, station_path


def source_file_records(path: Path) -> list[dict]:
    files = sorted(path.parent.glob(f"{path.stem}.*")) if path.suffix.lower() == ".shp" else [path]
    return [{"path": str(file.resolve()), "bytes": file.stat().st_size, "sha256": sha256(file)} for file in files]


def hbd_config(dem: Path, out_dir: Path, resolution: int, defaults: dict, spatial: bool) -> dict:
    hbd = defaults["hydrobathydem"]
    config = {
        "dem": str(dem),
        "out-dir": str(out_dir),
        "resample-dem": False,
        "auto-rivers-d4": True,
        "min-area": defaults["d4_min_river_area_km2"],
        "river-geometry-source": "spatial_coefficients_or_power_law" if spatial else "power_law",
        "beta-1": hbd["beta_1"],
        "beta-2": hbd["beta_2"],
        "alfa-1": hbd["alfa_1"],
        "alfa-2": hbd["alfa_2"],
        "carve-mode": "wide",
        "channel-cell-width-m": float(resolution),
        "river-width-cap-m": hbd["river_width_cap_m"],
        "river-depth-cap-m": hbd["river_depth_cap_m"],
        "max-H-abg-m": hbd["max_h_abg_m"],
        "max-nodata-fill-pixels": hbd["max_nodata_fill_pixels"],
        "slope-percentile": hbd["slope_percentile"],
        "smooth-filter-cells": hbd["smooth_filter_cells"],
        "protect-stream-buffer-m": float(3 * resolution),
        "breach-dist-cells": max(5, int(math.ceil(9000 / resolution))),
        "breach-flat-increment": hbd["breach_flat_increment"],
        "fill-max-depth-m": hbd["fill_max_depth_m"],
    }
    if spatial:
        threshold = int(hbd["lin_calibration_threshold_km2"])
        calibration = out_dir.parent / "Lin2020" / "calibration"
        config.update(
            {
                "spatial-beta-1-raster": str(calibration / f"D4_beta_1_width_{threshold}km2.tif"),
                "spatial-beta-2-raster": str(calibration / f"D4_beta_2_width_{threshold}km2.tif"),
                "spatial-alfa-1-raster": str(calibration / f"D4_alfa_1_depth_{threshold}km2.tif"),
                "spatial-alfa-2-raster": str(calibration / f"D4_alfa_2_depth_{threshold}km2.tif"),
            }
        )
    return config


def resumable_command(command: str, required: list[Path]) -> str:
    checks = " && ".join(f"test -s {path}" for path in required)
    return f"if ! ( {checks} ); then {command}; fi; {checks}"


def domain_command(args: argparse.Namespace) -> None:
    catalog = load_catalog(args.config)
    case = catalog["cases"][args.case]
    defaults = catalog["defaults"]
    if args.resolution == defaults["nested_resolution_m"] and args.aoi is None:
        raise ValueError("A 30 m run requires --aoi with an explicit urban/nested boundary; no city boundary is invented")

    catchments, stations, catchment_path, station_path = select_camels(case, args.camels_root)
    selected = catchments if args.aoi is None else gpd.read_file(args.aoi)
    if selected.crs is None:
        raise ValueError("AOI has no CRS")
    selected = gpd.GeoDataFrame(
        {"case_id": [args.case], "resolution_m": [args.resolution]},
        geometry=[selected.to_crs(defaults["target_crs"]).geometry.union_all()],
        crs=defaults["target_crs"],
    )
    if not selected.geometry.iloc[0].is_valid:
        selected.geometry = selected.geometry.make_valid()

    case_dir = args.workspace / args.case / f"{args.resolution}m"
    domain_dir = case_dir / "Domain"
    build_dir = case_dir / "Build" / "HydroBathyDEM"
    domain_dir.mkdir(parents=True, exist_ok=True)
    aoi_path = domain_dir / "model_aoi.gpkg"
    selected.to_file(aoi_path, layer="model_aoi", driver="GPKG")
    catchments.to_crs(defaults["target_crs"]).to_file(
        domain_dir / "camels_reference_catchments.gpkg", layer="catchments", driver="GPKG"
    )
    station_out = stations.to_crs(defaults["target_crs"])
    station_out.to_file(domain_dir / "camels_reference_stations.gpkg", layer="stations", driver="GPKG")
    station_table = pd.DataFrame(
        {
            "Gauge": station_out["gauge_id"].astype(str),
            "Label_Name": [next(item["name"] for item in case["gauges"] if item["id"] == gid) for gid in station_out["gauge_id"]],
            "Easting_m_original": station_out.geometry.x,
            "Northing_m_original": station_out.geometry.y,
            "Longitude": stations.to_crs(4326).geometry.x,
            "Latitude": stations.to_crs(4326).geometry.y,
        }
    ).sort_values("Gauge")
    station_table.to_csv(domain_dir / "camels_gauges_original.csv", index=False)

    area_km2 = float(selected.geometry.area.sum() / 1e6)
    active_cells = int(math.ceil(area_km2 * 1e6 / args.resolution**2))
    if args.resolution == 30 and active_cells > int(defaults["full_basin_30m_cell_limit"]):
        raise RuntimeError(
            f"The supplied 30 m AOI has about {active_cells:,} active cells, above the "
            f"{defaults['full_basin_30m_cell_limit']:,}-cell safety limit. Use a smaller nested urban AOI."
        )

    fabdem_dir = build_dir / "FABDEM"
    prefix = f"DEM_fabdem_{args.case}"
    dem = fabdem_dir / f"{prefix}_{args.resolution}m_filled.tif"
    first_dir, final_dir = build_dir / "FirstPass", build_dir / "Final"
    first_config = hbd_config(dem, first_dir, args.resolution, defaults, spatial=False)
    final_config = hbd_config(dem, final_dir, args.resolution, defaults, spatial=True)
    write_json(build_dir / "hydrobathydem_first_pass.json", first_config)
    write_json(build_dir / "hydrobathydem_final.json", final_config)

    lin_dir = build_dir / "Lin2020"
    first_dem = first_dir / "dem/DEM_hydrologically_conditioned_pre_bathymetry.tif"
    first_mask = first_dir / "d4/D4_idx_facc.tif"
    first_direction = first_dir / "d4/D4_flow_direction.tif"
    first_area = first_dir / "d4/D4_Wshed_Properties_fac_area_km2.tif"
    lin_gpkg = lin_dir / "processed/lin2020_dem_domain_width_depth.gpkg"
    calibration = lin_dir / "calibration"
    hydraulic_decision = calibration / "hydraulic_geometry_decision.json"
    final_conditioned_dem = final_dir / "dem/DEM_hydraulic_conditioned.tif"
    final_exported_dem = final_dir / "dem/DEM_conditioned_no_bathymetry.tif"
    final_width = final_dir / "d4/D4_River_Width_river_cells_m.tif"
    final_depth = final_dir / "d4/D4_River_Depth_river_cells_m.tif"
    commands = [
        resumable_command(
            f"hydrobathydem-build-dem --download --aoi {aoi_path} --output-dir {fabdem_dir} --output-prefix {prefix} --target-crs {defaults['target_crs']} --target-resolution {args.resolution} --download-cache {args.fabdem_cache} --fill-gaps",
            [dem],
        ),
        resumable_command(
            f"hydrobathydem-condition --config {build_dir / 'hydrobathydem_first_pass.json'}",
            [first_dem, first_mask, first_direction, first_area],
        ),
        resumable_command(
            f"hydrobathydem-prepare-lin2020 --download --raw-dir {args.lin_cache} --processed-dir {lin_dir / 'processed'} --dem-grid {first_dem} --d4-mask {first_mask} --manning-n {defaults['manning_n']} --max-H-abg-m {defaults['hydrobathydem']['max_h_abg_m']} --carve-mode wide",
            [lin_gpkg],
        ),
        f"python {HERE / 'calibrate_india_highres_hydraulics.py'} --case-dir {case_dir} && test -s {hydraulic_decision}",
        f"if python -c \"import json,sys; sys.exit(0 if json.load(open('{hydraulic_decision}'))['mode'] == 'spatial_coefficients_or_power_law' else 1)\"; then hydrobathydem-preflight --config {build_dir / 'hydrobathydem_final.json'} --require-lin --require-spatial-coefficients; else hydrobathydem-preflight --config {build_dir / 'hydrobathydem_final.json'}; fi",
        resumable_command(
            f"hydrobathydem-condition --config {build_dir / 'hydrobathydem_final.json'}",
            [final_conditioned_dem],
        ),
        resumable_command(
            f"hydrobathydem-export-geometry --out-dir {final_dir} --overwrite",
            [final_exported_dem, final_width, final_depth],
        ),
    ]
    (build_dir / "workflow_commands.txt").write_text("\n".join(commands) + "\n")

    resource = {
        "case": args.case,
        "resolution_m": args.resolution,
        "aoi_source": "CAMELS-IND catchment union" if args.aoi is None else str(args.aoi),
        "area_km2": area_km2,
        "estimated_active_cells": active_cells,
        "requested_memory_gib": request_memory_gib(active_cells, defaults),
        "simulation_start": simulation_start(case, defaults).isoformat(),
        "event_start": case["event_start"],
        "simulation_end": case["simulation_end"],
        "coastal_free_outflow": case["coastal_free_outflow"],
        "model_source": git_lock(args.model_root),
    }
    write_json(case_dir / "resource_plan.json", resource)
    manifest = {
        **resource,
        "status": "domain_planned",
        "sources": catalog["sources"],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "original_input_files": {
            "camels_catchments": source_file_records(catchment_path),
            "camels_gauge_stations": source_file_records(station_path),
            "camels_observed_streamflow": source_file_records(
                args.camels_root / "streamflow_timeseries" / "streamflow_observed.csv"
            ),
        },
        "files": {
            "model_aoi": {"path": str(aoi_path), "sha256": sha256(aoi_path)},
            "gauges": {"path": str(domain_dir / "camels_gauges_original.csv"), "sha256": sha256(domain_dir / "camels_gauges_original.csv")},
            "hydrobathydem_first_pass": str(build_dir / "hydrobathydem_first_pass.json"),
            "hydrobathydem_final": str(build_dir / "hydrobathydem_final.json"),
        },
    }
    write_json(case_dir / "case_manifest.json", manifest)
    print(json.dumps(resource, indent=2))


def active_mask(dataset: rasterio.DatasetReader, values: np.ndarray) -> np.ndarray:
    if dataset.nodata is None or (isinstance(dataset.nodata, float) and np.isnan(dataset.nodata)):
        return np.isfinite(values)
    return np.isfinite(values) & (values != dataset.nodata)


def align_raster(src: Path, dst: Path, template: Path, method: Resampling, dtype: str | None = None) -> None:
    with rasterio.open(template) as ref, rasterio.open(src) as source:
        out_dtype = dtype or source.dtypes[0]
        nodata = np.nan if np.issubdtype(np.dtype(out_dtype), np.floating) else 0
        profile = ref.profile.copy()
        profile.update(dtype=out_dtype, count=1, nodata=nodata, compress="deflate", BIGTIFF="IF_SAFER")
        data = np.full(ref.shape, nodata, dtype=out_dtype)
        reproject(
            rasterio.band(source, 1), data,
            src_transform=source.transform, src_crs=source.crs, src_nodata=source.nodata,
            dst_transform=ref.transform, dst_crs=ref.crs, dst_nodata=nodata, resampling=method,
        )
        ref_values = ref.read(1)
        data[~active_mask(ref, ref_values)] = nodata
        dst.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(dst, "w", **profile) as output:
            output.write(data, 1)


def write_like(path: Path, values: np.ndarray, template: Path, dtype: str = "float32") -> None:
    with rasterio.open(template) as ref:
        profile = ref.profile.copy()
        profile.update(dtype=dtype, count=1, nodata=np.nan, compress="deflate", BIGTIFF="IF_SAFER")
        with rasterio.open(path, "w", **profile) as output:
            output.write(values.astype(dtype), 1)


def condition_neal_channel_geometry(
    static: Path, hydrobathydem: Path, sources: Path, resolution_m: float
) -> dict:
    """Represent rivers wider than one cell with equivalent full-cell geometry."""
    width_path = static / "RiverWidths.tif"
    depth_path = static / "RiverDepths.tif"
    equivalent_depth_path = hydrobathydem / "d4" / "D4_H_abg_m.tif"
    with rasterio.open(static / "DEM.tif") as ref:
        dem = ref.read(1)
        active = active_mask(ref, dem)
    with rasterio.open(width_path) as source:
        width = source.read(1).astype("float64")
    with rasterio.open(depth_path) as source:
        depth = source.read(1).astype("float64")
    with rasterio.open(equivalent_depth_path) as source:
        if source.shape != width.shape:
            raise RuntimeError("HydroBathyDEM equivalent channel depth does not match the staged grid")
        equivalent_depth = source.read(1).astype("float64")

    channel = active & np.isfinite(width) & np.isfinite(depth) & (width > 0) & (depth > 0)
    wide = channel & (width > resolution_m)
    if np.any(wide):
        if not np.all(np.isfinite(equivalent_depth[wide]) & (equivalent_depth[wide] > 0)):
            raise RuntimeError("Wide channels lack valid HydroBathyDEM equivalent depths")
        sources.mkdir(parents=True, exist_ok=True)
        shutil.copy2(width_path, sources / "RiverWidths_source_hydrobathydem.tif")
        shutil.copy2(depth_path, sources / "RiverDepths_source_hydrobathydem.tif")
        shutil.copy2(equivalent_depth_path, sources / "RiverEquivalentDepth_source_hydrobathydem.tif")

        maximum_source_width_m = float(np.max(width[wide]))
        source_conveyance = width[wide] * np.power(depth[wide], 5.0 / 3.0)
        equivalent_conveyance = resolution_m * np.power(equivalent_depth[wide], 5.0 / 3.0)
        relative_error = np.abs(equivalent_conveyance - source_conveyance) / source_conveyance
        width[wide] = resolution_m
        depth[wide] = equivalent_depth[wide]
        width[~active] = np.nan
        depth[~active] = np.nan
        write_like(width_path, width, static / "DEM.tif")
        write_like(depth_path, depth, static / "DEM.tif")
    else:
        maximum_source_width_m = None
        relative_error = np.array([0.0])

    return {
        "method": "HydroBathyDEM equivalent full-cell channel geometry",
        "reason": (
            "The Neal subgrid represents one channel per centerline cell. Channels wider than the "
            "model cell are set to the full cell width and assigned HydroBathyDEM's equivalent depth, "
            "which preserves wide-channel Manning conveyance."
        ),
        "resolution_m": resolution_m,
        "river_cells": int(channel.sum()),
        "converted_wide_channel_cells": int(wide.sum()),
        "converted_river_fraction": float(wide.sum() / channel.sum()) if np.any(channel) else 0.0,
        "maximum_source_width_m": maximum_source_width_m,
        "maximum_equivalent_depth_m": float(np.max(equivalent_depth[wide])) if np.any(wide) else None,
        "maximum_relative_wide_manning_conveyance_error": float(np.max(relative_error)),
        "unaltered_width_raster": "Sources/RiverWidths_source_hydrobathydem.tif" if np.any(wide) else None,
        "unaltered_depth_raster": "Sources/RiverDepths_source_hydrobathydem.tif" if np.any(wide) else None,
        "equivalent_depth_raster": "Sources/RiverEquivalentDepth_source_hydrobathydem.tif" if np.any(wide) else None,
    }


def condition_dtb(static: Path, sources: Path, minimum_depth_m: float) -> dict:
    path = static / "DTB.tif"
    with rasterio.open(static / "DEM.tif") as ref:
        dem = ref.read(1)
        active = active_mask(ref, dem)
    with rasterio.open(path) as source:
        dtb = source.read(1).astype("float32")
    if not np.isfinite(dtb[active]).all():
        raise RuntimeError("Depth-to-bedrock source has missing values within the active domain")
    adjusted = active & (dtb < minimum_depth_m)
    sources.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, sources / "DTB_source_aligned.tif")
    dtb[adjusted] = minimum_depth_m
    dtb[~active] = np.nan
    write_like(path, dtb, static / "DEM.tif")
    return {
        "source": "projects/ee-marcusep2025/assets/Depth_to_bedrock",
        "minimum_effective_depth_m": minimum_depth_m,
        "adjusted_cells": int(adjusted.sum()),
        "adjusted_active_fraction": float(adjusted.sum() / active.sum()),
        "reason": "Nonpositive source values were set to the numerical soil-column floor already used by the India workflow.",
        "unaltered_aligned_raster": "Sources/DTB_source_aligned.tif",
    }


def build_initial_states(static: Path, source_static: Path) -> dict:
    align_raster(source_static / "Fan_WTD.tif", static / "Fan_WTD_aligned.tif", static / "DEM.tif", Resampling.bilinear, "float32")
    with rasterio.open(static / "DEM.tif") as ref:
        dem = ref.read(1).astype("float32")
        active = active_mask(ref, dem)
    with rasterio.open(static / "DTB.tif") as src:
        dtb = src.read(1).astype("float32")
    with rasterio.open(static / "Fan_WTD_aligned.tif") as src:
        wtd = src.read(1).astype("float32")
    wtd = np.where(np.isfinite(wtd), np.maximum(wtd, 0), dtb)
    lower = dem - dtb
    upper = dem
    head = np.minimum(np.maximum(dem - wtd, lower), upper)
    head[~active] = np.nan
    write_like(static / "GW_table.tif", head, static / "DEM.tif")

    with rasterio.open(source_static / "ERA5Land_initial_soil_water.tif") as src:
        if src.count != 4:
            raise RuntimeError("ERA5Land_initial_soil_water.tif must contain four soil-water bands")
        layers = []
        with rasterio.open(static / "DEM.tif") as ref:
            for band in range(1, 5):
                layer = np.full(ref.shape, np.nan, dtype="float32")
                reproject(
                    rasterio.band(src, band), layer,
                    src_transform=src.transform, src_crs=src.crs, src_nodata=src.nodata,
                    dst_transform=ref.transform, dst_crs=ref.crs, dst_nodata=np.nan,
                    resampling=Resampling.bilinear,
                )
                layers.append(layer)
    with rasterio.open(static / "SOIL.tif") as src:
        soil = src.read(1).astype("int16")
    zwt = np.clip(dem - head, 0, np.maximum(dtb, 0))
    remaining, water = zwt.copy(), np.zeros_like(zwt)
    for values, thickness in zip(layers, SOIL_LAYER_M):
        used = np.minimum(remaining, thickness)
        water += np.where(np.isfinite(values), values, 0) * used
        remaining -= used
    water += np.where(np.isfinite(layers[-1]), layers[-1], 0) * np.maximum(remaining, 0)
    theta_r = np.full(soil.shape, SOIL_THETA_R[0], dtype="float32")
    valid_soil = (soil >= 1) & (soil <= len(SOIL_THETA_R))
    theta_r[valid_soil] = SOIL_THETA_R[soil[valid_soil] - 1]
    initial_mm = np.maximum(water / np.maximum(zwt, 1e-6) - theta_r, 0) * zwt * 1000
    initial_mm[~active] = np.nan
    write_like(static / "Initial_Soil_Moisture.tif", initial_mm, static / "DEM.tif")
    return {
        "water_table_source": "Fan et al. annual-mean depth",
        "groundwater_head_clipped_cells": int(np.count_nonzero(active & ~np.isclose(head, dem - wtd, equal_nan=True))),
        "median_initial_soil_storage_mm": float(np.nanmedian(initial_mm)),
    }


def snap_gauges(case_dir: Path, static: Path, max_distance_m: float = 5000) -> dict:
    gauges = pd.read_csv(case_dir / "Domain" / "camels_gauges_original.csv", dtype={"Gauge": str})
    with rasterio.open(static / "RiverWidths.tif") as ref:
        width = ref.read(1)
        river_rows, river_cols = np.where(np.isfinite(width) & (width > 0))
        if river_rows.size == 0:
            raise RuntimeError("No active river cells are available for gauge snapping")
        river_x, river_y = rasterio.transform.xy(ref.transform, river_rows, river_cols, offset="center")
        transformer = Transformer.from_crs(4326, ref.crs, always_xy=True)
    river_x, river_y = np.asarray(river_x), np.asarray(river_y)
    snapped = []
    for record in gauges.itertuples(index=False):
        x0, y0 = transformer.transform(record.Longitude, record.Latitude)
        distance2 = (river_x - x0) ** 2 + (river_y - y0) ** 2
        index = int(np.argmin(distance2))
        distance = float(np.sqrt(distance2[index]))
        if distance > max_distance_m:
            raise RuntimeError(f"Gauge {record.Gauge} is {distance:.0f} m from the nearest active river cell")
        snapped.append(
            {
                "Gauge": record.Gauge,
                "Easting_m": float(river_x[index]),
                "Northing_m": float(river_y[index]),
                "Label_Name": record.Label_Name,
                "Original_Easting_m": float(record.Easting_m_original),
                "Original_Northing_m": float(record.Northing_m_original),
                "Snap_Distance_m": distance,
            }
        )
    output = case_dir / "Forcing" / "Observed_Gauges" / "camels_india_gauges_2km.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(snapped).to_csv(output, index=False)
    return {"count": len(snapped), "maximum_snap_distance_m": max(item["Snap_Distance_m"] for item in snapped)}


def stage_rainfall(case_dir: Path, rainfall_source: Path, start: datetime, end: datetime) -> int:
    output = case_dir / "Forcing" / "Rainfall"
    output.mkdir(parents=True, exist_ok=True)
    time, count = start, 0
    while time <= end:
        name = f"IMERG_30min_mmhr_IndiaRegion_{time:%Y_%m_%d_%H_%M}.tif"
        source = rainfall_source / name
        if not source.is_file():
            raise FileNotFoundError(f"Missing IMERG interval: {source}")
        destination = output / name
        if not destination.exists():
            destination.symlink_to(source.resolve())
        count += 1
        time += timedelta(minutes=30)
    return count


def subset_streamflow(source: Path, destination: Path, gauges: list[str], start: datetime, end: datetime) -> dict:
    observed = pd.read_csv(source, usecols=["year", "month", "day", *gauges])
    dates = pd.to_datetime(observed[["year", "month", "day"]])
    selected = observed[(dates >= start) & (dates < end)].copy()
    destination.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(destination, index=False)
    coverage = {gauge: float(selected[gauge].notna().mean()) for gauge in gauges}
    return {"rows": len(selected), "coverage_fraction": coverage, "used_for_forcing": False}


def stage_command(args: argparse.Namespace) -> None:
    catalog = load_catalog(args.config)
    case = catalog["cases"][args.case]
    defaults = catalog["defaults"]
    case_dir = args.workspace / args.case / f"{args.resolution}m"
    static = case_dir / "Static"
    static.mkdir(parents=True, exist_ok=True)
    hbd = case_dir / "Build" / "HydroBathyDEM" / "Final"
    source_static = args.source_data / args.case / f"Static_{args.resolution}m"
    source_forcing = args.source_data / args.case / "Forcing"
    source_manifest = args.source_data / args.case / "source_download_manifest.json"
    if not source_manifest.is_file():
        raise FileNotFoundError(source_manifest)
    source_record = json.loads(source_manifest.read_text())
    required_components = {f"static_{args.resolution}m", "forcing"}
    missing_components = required_components - set(source_record.get("components", {}))
    if missing_components:
        raise RuntimeError(f"Source manifest is missing components: {sorted(missing_components)}")
    archived_manifest = case_dir / "Sources" / "source_download_manifest.json"
    archived_manifest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_manifest, archived_manifest)
    hbd_files = {
        "DEM.tif": hbd / "dem" / "DEM_conditioned_no_bathymetry.tif",
        "RiverWidths.tif": hbd / "d4" / "D4_River_Width_river_cells_m.tif",
        "RiverDepths.tif": hbd / "d4" / "D4_River_Depth_river_cells_m.tif",
        "D4_flow_direction.tif": hbd / "d4" / "D4_flow_direction.tif",
        "D4_idx_facc.tif": hbd / "d4" / "D4_idx_facc.tif",
    }
    for name, source in hbd_files.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, static / name)
    channel_geometry = condition_neal_channel_geometry(
        static, hbd, case_dir / "Sources", float(args.resolution)
    )
    with rasterio.open(static / "DEM.tif") as dem:
        resolution = abs(dem.transform.a)
        if not np.isclose(resolution, args.resolution, atol=1e-5):
            raise RuntimeError(f"HydroBathyDEM resolution is {resolution}, expected {args.resolution} m")

    static_sources = {
        "LULC.tif": ("WorldCover_LULC.tif", Resampling.nearest, None),
        "SOIL.tif": ("OpenLandMap_SOIL.tif", Resampling.nearest, None),
        "DTB.tif": ("Depth_to_bedrock.tif", Resampling.bilinear, "float32"),
        "LAI.tif": ("MODIS_LAI_median.tif", Resampling.bilinear, "float32"),
    }
    for destination, (source_name, method, dtype) in static_sources.items():
        source = source_static / source_name
        if not source.is_file():
            raise FileNotFoundError(source)
        align_raster(source, static / destination, static / "DEM.tif", method, dtype)
    dtb_summary = condition_dtb(
        static, case_dir / "Sources", float(defaults["minimum_effective_soil_depth_m"])
    )
    with rasterio.open(static / "DEM.tif") as ref:
        dem = ref.read(1)
        albedo = np.full(ref.shape, float(defaults["albedo"]), dtype="float32")
        albedo[~active_mask(ref, dem)] = np.nan
    write_like(static / "Albedo.tif", albedo, static / "DEM.tif")
    if (source_static / "MODIS_Albedo_median.tif").is_file():
        archive = case_dir / "Sources" / "MODIS_Albedo_median.tif"
        archive.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_static / "MODIS_Albedo_median.tif", archive)
    state_summary = build_initial_states(static, source_static)
    gauge_summary = snap_gauges(case_dir, static)
    start = simulation_start(case, defaults)
    end = parse_time(case["simulation_end"])
    rain_count = stage_rainfall(case_dir, source_forcing / "Rainfall", start, end)
    etp_source = source_forcing / "Evapotranspiration" / "ETP_input_data.xlsx"
    if not etp_source.is_file():
        raise FileNotFoundError(etp_source)
    etp_destination = case_dir / "Forcing" / "Evapotranspiration" / "ETP_input_data.xlsx"
    etp_destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(etp_source, etp_destination)
    streamflow = subset_streamflow(
        args.camels_root / "streamflow_timeseries" / "streamflow_observed.csv",
        case_dir / "Validation" / "CAMELS_IND_observed_streamflow.csv",
        case["gauge_ids"], start, end,
    )
    manifest_path = case_dir / "case_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    hydraulic_decision_path = (
        case_dir / "Build" / "HydroBathyDEM" / "Lin2020" / "calibration" / "hydraulic_geometry_decision.json"
    )
    if not hydraulic_decision_path.is_file():
        raise FileNotFoundError(hydraulic_decision_path)
    hbd_qa_path = case_dir / "Build" / "HydroBathyDEM" / "Final" / "reports" / "qa_scorecard.json"
    if not hbd_qa_path.is_file():
        raise FileNotFoundError(hbd_qa_path)
    manifest.update(
        {
            "status": "staged",
            "staged_utc": datetime.now(timezone.utc).isoformat(),
            "model_source": git_lock(args.model_root),
            "rainfall_maps": rain_count,
            "initial_states": state_summary,
            "dtb_conditioning": dtb_summary,
            "neal_channel_geometry": channel_geometry,
            "observed_gauges": gauge_summary,
            "observed_streamflow": streamflow,
            "hydraulic_geometry": json.loads(hydraulic_decision_path.read_text()),
            "hydrobathydem_qa": json.loads(hbd_qa_path.read_text()),
            "validation_is_not_forcing": True,
            "source_download_manifest": str(archived_manifest.relative_to(case_dir)),
            "source_download_manifest_sha256": sha256(archived_manifest),
        }
    )
    checksums = {}
    for path in sorted(static.glob("*.tif")):
        checksums[str(path.relative_to(case_dir))] = sha256(path)
    manifest["staged_static_sha256"] = checksums
    write_json(manifest_path, manifest)
    print(json.dumps({"case": args.case, "case_dir": str(case_dir), "rainfall_maps": rain_count, **gauge_summary}, indent=2))


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    subparsers = command.add_subparsers(dest="command", required=True)
    domain = subparsers.add_parser("domain", help="Select CAMELS catchments and write HydroBathyDEM configs")
    domain.add_argument("--case", required=True)
    domain.add_argument("--resolution", required=True, type=int, choices=[30, 90, 2000])
    domain.add_argument("--workspace", required=True, type=Path)
    domain.add_argument("--camels-root", required=True, type=Path)
    domain.add_argument("--model-root", type=Path, default=HERE.parent)
    domain.add_argument("--aoi", type=Path, help="Required explicit nested/urban AOI for 30 m")
    domain.add_argument("--fabdem-cache", type=Path, required=True)
    domain.add_argument("--lin-cache", type=Path, required=True)
    domain.set_defaults(func=domain_command)
    stage = subparsers.add_parser("stage", help="Assemble a HydroPol2D case from original-source derivatives")
    stage.add_argument("--case", required=True)
    stage.add_argument("--resolution", required=True, type=int, choices=[30, 90, 2000])
    stage.add_argument("--workspace", required=True, type=Path)
    stage.add_argument("--source-data", required=True, type=Path)
    stage.add_argument("--camels-root", required=True, type=Path)
    stage.add_argument("--model-root", type=Path, default=HERE.parent)
    stage.set_defaults(func=stage_command)
    return command


def main() -> None:
    args = parser().parse_args()
    catalog = load_catalog(args.config)
    if args.case not in catalog["cases"]:
        raise ValueError(f"Unknown case {args.case}; choose from {', '.join(catalog['cases'])}")
    args.func(args)


if __name__ == "__main__":
    main()
