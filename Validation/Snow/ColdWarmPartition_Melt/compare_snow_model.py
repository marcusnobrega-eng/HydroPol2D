#!/usr/bin/env python3
"""Compare P1-SNOW-001 HydroPol2D snow output against the reference solution."""

from __future__ import annotations

import argparse
import csv
from math import isfinite, sqrt
from pathlib import Path
import sys


CASE_ID = "P1-SNOW-001"
TOLERANCE_M3 = 1.0e-6
SERIES_TOLERANCE = 1.0e-8


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"no rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def keyed(rows: list[dict[str, str]]) -> dict[tuple[str, float], dict[str, str]]:
    return {(row["scenario"], float(row["time_day"])): row for row in rows}


def rmse(values: list[float]) -> float:
    return sqrt(sum(value * value for value in values) / len(values)) if values else 0.0


def compare(reference_rows: list[dict[str, str]], model_rows: list[dict[str, str]], output_dir: Path) -> int:
    ref_by_key = keyed(reference_rows)
    model_by_key = keyed(model_rows)

    missing = sorted(set(ref_by_key) - set(model_by_key))
    extra = sorted(set(model_by_key) - set(ref_by_key))
    if missing or extra:
        raise ValueError(f"row mismatch: missing={missing}, extra={extra}")

    series = {
        "swe_error_mm": ("swe_mm", "mm"),
        "snow_depth_error_mm": ("snow_depth_mm", "mm"),
        "melt_error_mm": ("melt_mm", "mm"),
        "snowfall_error_mm": ("snowfall_mm", "mm"),
        "rain_error_mm": ("rain_mm", "mm"),
        "sublimation_error_mm": ("sublimation_mm", "mm"),
        "rho_snow_error_kg_m3": ("rho_snow_kg_m3", "kg/m3"),
    }
    diffs_by_metric: dict[str, list[float]] = {metric: [] for metric in series}
    max_abs_by_metric = {metric: 0.0 for metric in series}
    max_abs_residual_m3 = 0.0
    detail_rows: list[dict[str, object]] = []
    finite_values = True

    for key in sorted(ref_by_key):
        ref = ref_by_key[key]
        model = model_by_key[key]
        detail: dict[str, object] = {
            "case_id": CASE_ID,
            "scenario": key[0],
            "time_day": key[1],
        }
        for metric, (column, _) in series.items():
            error = f(model, column) - f(ref, column)
            diffs_by_metric[metric].append(error)
            max_abs_by_metric[metric] = max(max_abs_by_metric[metric], abs(error))
            detail[metric] = error
            finite_values = finite_values and isfinite(error)

        residual_m3 = f(model, "mass_residual_mm") / 1000.0 * f(ref, "cell_area_m2")
        detail["mass_residual_m3"] = residual_m3
        finite_values = finite_values and isfinite(residual_m3)
        max_abs_residual_m3 = max(max_abs_residual_m3, abs(residual_m3))
        detail_rows.append(detail)

    metric_rows: list[dict[str, object]] = []
    for metric, (_, units) in series.items():
        rmse_metric = metric.replace("_error_", "_rmse_")
        metric_rows.append(
            {
                "case_id": CASE_ID,
                "module": "Snow module",
                "benchmark_type": "independent_reference",
                "metric": rmse_metric,
                "value": rmse(diffs_by_metric[metric]),
                "units": units,
                "threshold": SERIES_TOLERANCE,
                "pass": rmse(diffs_by_metric[metric]) <= SERIES_TOLERANCE,
                "notes": "Compared against independent Python mirror of Snow_Model_Function formulas.",
            }
        )
        metric_rows.append(
            {
                "case_id": CASE_ID,
                "module": "Snow module",
                "benchmark_type": "independent_reference",
                "metric": f"max_abs_{metric}",
                "value": max_abs_by_metric[metric],
                "units": units,
                "threshold": SERIES_TOLERANCE,
                "pass": max_abs_by_metric[metric] <= SERIES_TOLERANCE,
                "notes": "Compared against independent Python mirror of Snow_Model_Function formulas.",
            }
        )

    metric_rows.append(
        {
            "case_id": CASE_ID,
            "module": "Snow module",
            "benchmark_type": "mass_balance",
            "metric": "max_abs_mass_residual_m3",
            "value": max_abs_residual_m3,
            "units": "m3",
            "threshold": TOLERANCE_M3,
            "pass": max_abs_residual_m3 <= TOLERANCE_M3,
            "notes": "Maximum absolute per-step snowpack residual converted using representative 20 m cell area.",
        }
    )

    benchmark_pass = all(row["pass"] for row in metric_rows if row["benchmark_type"] == "independent_reference")
    mass_balance_pass = max_abs_residual_m3 <= TOLERANCE_M3
    diagnostics_pass = bool(benchmark_pass and mass_balance_pass and finite_values)

    mass_rows = [
        {
            "case_id": CASE_ID,
            "control_volume": "snowpack",
            "time_start_day": 0,
            "time_end_day": max(f(row, "time_day") for row in reference_rows),
            "initial_storage": 0,
            "inflow": "see row-level detail",
            "outflow": "see row-level detail",
            "source": "snowfall",
            "sink": "melt_and_sublimation",
            "final_storage": "see row-level detail",
            "residual": max_abs_residual_m3,
            "relative_error_pct": "not_normalized",
            "units": "m3",
            "notes": "Maximum absolute per-step residual across v-tilted snow scenarios.",
        }
    ]
    pass_rows = [
        {
            "case_id": CASE_ID,
            "module": "Snow module",
            "status": "pass" if diagnostics_pass else "fail",
            "primary_reason": (
                "HydroPol2D MATLAB snow output passes formula and mass-balance diagnostics."
                if diagnostics_pass
                else "One or more snow diagnostics exceeded tolerance."
            ),
            "mass_balance_pass": str(mass_balance_pass).lower(),
            "benchmark_pass": str(benchmark_pass).lower(),
            "nan_check_pass": str(finite_values).lower(),
            "report_ready": str(diagnostics_pass).lower(),
            "reviewer": "compare_snow_model.py",
            "date_utc": "2026-06-15",
        }
    ]

    write_rows(output_dir / "Tables" / "P1-SNOW-001_detail_errors.csv", detail_rows)
    write_rows(output_dir / "Metric_Summary.csv", metric_rows)
    write_rows(output_dir / "Mass_Balance.csv", mass_rows)
    write_rows(output_dir / "Pass_Fail.csv", pass_rows)
    return 0 if diagnostics_pass else 1


def parse_args() -> argparse.Namespace:
    case_dir = Path(__file__).resolve().parent
    model_root = case_dir.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        type=Path,
        default=model_root / "Validation" / "Reference_Outputs" / "Phase1" / "snow_degree_day" / "P1-SNOW-001_reference.csv",
    )
    parser.add_argument(
        "--model-output",
        type=Path,
        default=case_dir / "Outputs" / "Validation" / "HydroPol2D_Snow_Model.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=case_dir / "Outputs" / "Validation")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.model_output.exists():
        print(
            f"ERROR: model output not found: {args.model_output}\n"
            "Run run_snow_model.m in MATLAB first.",
            file=sys.stderr,
        )
        return 2
    return compare(read_rows(args.reference), read_rows(args.model_output), args.output_dir)


if __name__ == "__main__":
    raise SystemExit(main())
