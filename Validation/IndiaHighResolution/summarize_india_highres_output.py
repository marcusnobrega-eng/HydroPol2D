#!/usr/bin/env python3
"""Summarize one completed India high-resolution HydroPol2D simulation."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timedelta
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask


DEPTH_PATTERN = re.compile(
    r"Flood_Depths_(\d{4})_(\d{2})_(\d{2})_(\d{2})_(\d{2})_(\d{2})\.tif$"
)


def model_dates(values: pd.Series, simulation_start: pd.Timestamp) -> pd.DatetimeIndex:
    """Interpret the two HydroPol2D time formats written by post_processing.m."""
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.notna().mean() < 0.9:
        return pd.DatetimeIndex(pd.to_datetime(values, errors="coerce"))
    if numeric.median() > 100_000:
        # MATLAB serial day: 719529 is 1970-01-01.
        return pd.DatetimeIndex(pd.to_datetime(numeric - 719_529, unit="D", origin="unix"))
    return pd.DatetimeIndex(simulation_start + pd.to_timedelta(numeric, unit="min"))


def flow_columns(table: pd.DataFrame) -> list[str]:
    columns = [
        name for name in table.columns
        if "flow" in name.lower() and "discharge" in name.lower()
    ]
    if not columns:
        raise RuntimeError("Rating_Curve_Gauges.csv has no modeled-discharge columns")
    return columns


def hydrograph_metrics(case: Path, manifest: dict, minimum_coverage: float) -> list[dict]:
    model_path = case / "Outputs" / "Modeling_Results" / "Tables_CSV" / "Rating_Curve_Gauges.csv"
    if not model_path.is_file():
        raise FileNotFoundError(
            f"Missing {model_path}; run the complete simulation with post-processing enabled"
        )
    modeled = pd.read_csv(model_path)
    gauges = pd.read_csv(
        case / "Forcing" / "Observed_Gauges" / "camels_india_gauges_2km.csv",
        dtype={"Gauge": str},
    )
    flows = flow_columns(modeled)
    if len(flows) != len(gauges):
        raise RuntimeError(f"Found {len(flows)} modeled series for {len(gauges)} gauges")

    start = pd.Timestamp(manifest["simulation_start"])
    event_start = pd.Timestamp(manifest["event_start"])
    event_end = pd.Timestamp(manifest["simulation_end"])
    modeled_dates = model_dates(modeled.iloc[:, 0], start).normalize()
    observed = pd.read_csv(case / "Validation" / "CAMELS_IND_observed_streamflow.csv")
    observed_dates = pd.to_datetime(observed[["year", "month", "day"]]).dt.normalize()
    rows: list[dict] = []

    for index, gauge in gauges.reset_index(drop=True).iterrows():
        model_daily = pd.Series(
            pd.to_numeric(modeled[flows[index]], errors="coerce").to_numpy(),
            index=modeled_dates,
        ).groupby(level=0).mean()
        observed_daily = pd.Series(
            pd.to_numeric(observed[gauge["Gauge"]], errors="coerce").to_numpy(),
            index=observed_dates,
        )
        event_days = pd.date_range(event_start.normalize(), event_end.normalize() - pd.Timedelta(days=1), freq="D")
        observed_event = observed_daily.reindex(event_days)
        coverage = float(observed_event.notna().mean())
        record = {
            "gauge_id": gauge["Gauge"],
            "gauge_name": gauge["Label_Name"],
            "event_observation_coverage": coverage,
            "quantitative_comparison": coverage >= minimum_coverage,
        }
        if coverage >= minimum_coverage:
            common = event_days[observed_event.notna() & model_daily.reindex(event_days).notna()]
            if len(common) == 0:
                raise RuntimeError(f"Gauge {gauge['Gauge']} has no common modeled and observed event days")
            sim = model_daily.reindex(common)
            obs = observed_daily.reindex(common)
            sim_peak_day, obs_peak_day = sim.idxmax(), obs.idxmax()
            sim_volume, obs_volume = float(sim.sum() * 86_400), float(obs.sum() * 86_400)
            record.update(
                {
                    "comparison_days": len(common),
                    "simulated_peak_m3s": float(sim.max()),
                    "observed_peak_m3s": float(obs.max()),
                    "peak_timing_error_days": int((sim_peak_day - obs_peak_day).days),
                    "simulated_volume_m3_on_common_days": sim_volume,
                    "observed_volume_m3_on_common_days": obs_volume,
                    "volume_bias_percent": 100 * (sim_volume - obs_volume) / obs_volume if obs_volume else None,
                }
            )
        rows.append(record)
    return rows


def depth_timestamp(path: Path) -> datetime:
    match = DEPTH_PATTERN.search(path.name)
    if not match:
        raise ValueError(path.name)
    return datetime(*map(int, match.groups()))


def inundation_metrics(
    case: Path, manifest: dict, thresholds: list[float], comparison_aoi: Path | None
) -> list[dict]:
    raster_dir = case / "Outputs" / "Modeling_Results" / "Rasters_Water_Depths"
    paths = sorted(raster_dir.glob("Flood_Depths_*.tif"), key=depth_timestamp)
    if not paths:
        raise FileNotFoundError(
            f"No daily depth rasters in {raster_dir}; run the complete simulation with post-processing enabled"
        )
    event_start = datetime.fromisoformat(manifest["event_start"])
    event_end = datetime.fromisoformat(manifest["simulation_end"])
    paths = [path for path in paths if event_start <= depth_timestamp(path) < event_end]
    if not paths:
        raise RuntimeError("No exported depth raster falls inside the event period")

    rows: list[dict] = []
    aoi = gpd.read_file(comparison_aoi) if comparison_aoi else None
    for path in paths:
        with rasterio.open(path) as source:
            depth = source.read(1, masked=True)
            valid = ~np.ma.getmaskarray(depth)
            if aoi is not None:
                projected = aoi.to_crs(source.crs)
                inside = geometry_mask(
                    projected.geometry, source.shape, source.transform, invert=True, all_touched=False
                )
                valid &= inside
            cell_area_km2 = abs(source.transform.a * source.transform.e) / 1e6
            valid_cells = int(valid.sum())
            values = np.asarray(depth.filled(np.nan), dtype=float)
            for threshold in thresholds:
                flooded = valid & np.isfinite(values) & (values >= threshold)
                rows.append(
                    {
                        "date": depth_timestamp(path).isoformat(),
                        "threshold_m": threshold,
                        "flooded_area_km2": float(flooded.sum() * cell_area_km2),
                        "fraction_of_comparison_area": float(flooded.sum() / valid_cells) if valid_cells else None,
                        "comparison_area_km2": float(valid_cells * cell_area_km2),
                    }
                )
    return rows


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    parser.add_argument("--comparison-aoi", type=Path)
    parser.add_argument("--minimum-observation-coverage", type=float, default=0.80)
    parser.add_argument("--thresholds", type=float, nargs="+", default=[0.01, 0.05, 0.15, 0.30])
    args = parser.parse_args()

    manifest = json.loads((args.case / "case_manifest.json").read_text())
    output = args.case / "Validation" / "HighResolutionSummary"
    output.mkdir(parents=True, exist_ok=True)
    hydro = hydrograph_metrics(args.case, manifest, args.minimum_observation_coverage)
    inundation = inundation_metrics(args.case, manifest, args.thresholds, args.comparison_aoi)
    pd.DataFrame(hydro).to_csv(output / "hydrograph_metrics.csv", index=False)
    pd.DataFrame(inundation).to_csv(output / "daily_inundated_area.csv", index=False)
    summary = {
        "case": manifest["case"],
        "resolution_m": manifest["resolution_m"],
        "event_start": manifest["event_start"],
        "simulation_end": manifest["simulation_end"],
        "comparison_aoi": str(args.comparison_aoi.resolve()) if args.comparison_aoi else "full model domain",
        "minimum_observation_coverage": args.minimum_observation_coverage,
        "flood_depth_thresholds_m": args.thresholds,
        "hydrographs": hydro,
        "maximum_daily_flooded_area_km2": {
            str(threshold): max(
                item["flooded_area_km2"] for item in inundation if item["threshold_m"] == threshold
            )
            for threshold in args.thresholds
        },
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
