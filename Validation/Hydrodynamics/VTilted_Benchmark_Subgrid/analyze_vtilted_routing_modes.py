#!/usr/bin/env python3
"""Analyze exact-config V-tilted routing-mode HydroPol2D runs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CASE_DIR = Path(__file__).resolve().parent
ROUTING_DIR = CASE_DIR / "Outputs" / "FullConfigRouting"
REF_PATH = CASE_DIR / "Reference" / "VTilted_Benchmark_Hydrograph.csv"
FIG_DIR = ROUTING_DIR / "Figures"

MODES = [
    ("full_momentum", "Full momentum"),
    ("local_inertial", "Local inertial"),
    ("cellular_automata", "Cellular automata"),
    ("kinematic", "Kinematic"),
    ("diffusive", "Diffusive"),
]


def interp_at(t_source: np.ndarray, q_source: np.ndarray, t_eval: np.ndarray) -> np.ndarray:
    return np.interp(t_eval, t_source, q_source, left=np.nan, right=np.nan)


def cumulative_trapezoid_m3(t_min: np.ndarray, q_m3s: np.ndarray) -> np.ndarray:
    out = np.zeros_like(t_min, dtype=float)
    if len(t_min) < 2:
        return out
    dt_s = np.diff(t_min) * 60.0
    area = 0.5 * (q_m3s[:-1] + q_m3s[1:]) * dt_s
    out[1:] = np.cumsum(area)
    return out


def score_mode(mode: str, label: str, ref: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    mode_dir = ROUTING_DIR / mode
    hydro = pd.read_csv(mode_dir / "FullConfig_Hydrograph.csv")
    mass = pd.read_csv(mode_dir / "FullConfig_Mass_Summary.csv").iloc[0]

    t_model = hydro["Time_min"].to_numpy(float)
    q_model = hydro["Q_m3s"].to_numpy(float)
    t_end = float(np.nanmax(t_model))

    ref_window = ref.loc[ref["Time_min"] <= t_end + 1e-9].copy()
    t_ref = ref_window["Time_min"].to_numpy(float)
    q_ref = ref_window["Benchmark_m3s"].to_numpy(float)
    q_model_i = interp_at(t_model, q_model, t_ref)
    valid = np.isfinite(q_model_i) & np.isfinite(q_ref)

    err = q_model_i[valid] - q_ref[valid]
    rmse = float(np.sqrt(np.mean(err**2)))
    mae = float(np.mean(np.abs(err)))
    max_error = float(np.max(np.abs(err)))
    rel_l2 = float(np.sqrt(np.sum(err**2)) / max(np.sqrt(np.sum(q_ref[valid] ** 2)), np.finfo(float).eps))
    nse = float(1.0 - np.sum(err**2) / max(np.sum((q_ref[valid] - np.mean(q_ref[valid])) ** 2), np.finfo(float).eps))

    i_peak_model = int(np.nanargmax(q_model))
    i_peak_ref = int(np.nanargmax(q_ref))
    peak_model = float(q_model[i_peak_model])
    peak_ref = float(q_ref[i_peak_ref])
    peak_time_model = float(t_model[i_peak_model])
    peak_time_ref = float(t_ref[i_peak_ref])
    peak_time_error = abs(peak_time_model - peak_time_ref)
    peak_mag_error_pct = 100.0 * abs(peak_model - peak_ref) / max(abs(peak_ref), np.finfo(float).eps)

    t_dense = t_model[(t_model >= 0) & (t_model <= t_end)]
    q_ref_dense = interp_at(ref["Time_min"].to_numpy(float), ref["Benchmark_m3s"].to_numpy(float), t_dense)
    q_model_dense = interp_at(t_model, q_model, t_dense)
    volume_ref_m3 = float(np.trapezoid(q_ref_dense, t_dense * 60.0))
    volume_model_m3 = float(np.trapezoid(q_model_dense, t_dense * 60.0))
    outlet_volume_error_pct = 100.0 * abs(volume_model_m3 - volume_ref_m3) / max(abs(volume_ref_m3), np.finfo(float).eps)

    passed_hydrograph_shape = nse > 0.95 and peak_time_error <= 20.0
    passed_mass = abs(float(mass["mass_error_pct"])) < 0.1
    passed_volume = outlet_volume_error_pct < 5.0

    metrics = {
        "case_id": f"P1-HYDRO-VTILT-{mode.upper().replace('_', '-')}-001",
        "routing_mode": mode,
        "routing_label": label,
        "rmse_m3s": rmse,
        "mae_m3s": mae,
        "max_error_m3s": max_error,
        "relative_l2": rel_l2,
        "nse": nse,
        "peak_time_model_min": peak_time_model,
        "peak_time_benchmark_min": peak_time_ref,
        "peak_time_error_min": peak_time_error,
        "peak_model_m3s": peak_model,
        "peak_benchmark_m3s": peak_ref,
        "peak_magnitude_error_pct": peak_mag_error_pct,
        "benchmark_volume_m3_0_to_model_end": volume_ref_m3,
        "model_volume_m3_0_to_model_end": volume_model_m3,
        "outlet_volume_error_pct": outlet_volume_error_pct,
        "rain_volume_m3": float(mass["rain_volume_m3"]),
        "outlet_volume_m3": float(mass["outlet_volume_m3"]),
        "final_storage_m3": float(mass["final_storage_m3"]),
        "mass_residual_m3": float(mass["mass_residual_m3"]),
        "mass_error_pct": float(mass["mass_error_pct"]),
        "passed_hydrograph_shape": passed_hydrograph_shape,
        "passed_mass_balance": passed_mass,
        "passed_volume_error_5pct": passed_volume,
        "passed_overall_exact_config": passed_hydrograph_shape and passed_mass,
    }

    comp = pd.DataFrame(
        {
            "Time_min": t_ref,
            "Benchmark_m3s": q_ref,
            f"{mode}_m3s": q_model_i,
            f"{mode}_error_m3s": q_model_i - q_ref,
        }
    )
    return metrics, comp


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ref = pd.read_csv(REF_PATH)

    metrics_rows = []
    comparison: pd.DataFrame | None = None
    dense_series = []

    for mode, label in MODES:
        metrics, comp = score_mode(mode, label, ref)
        metrics_rows.append(metrics)
        if comparison is None:
            comparison = comp[["Time_min", "Benchmark_m3s"]].copy()
        comparison[f"{mode}_m3s"] = comp[f"{mode}_m3s"]
        comparison[f"{mode}_error_m3s"] = comp[f"{mode}_error_m3s"]

        hydro = pd.read_csv(ROUTING_DIR / mode / "FullConfig_Hydrograph.csv")
        hydro["routing_mode"] = mode
        hydro["routing_label"] = label
        hydro["cumulative_outlet_volume_m3"] = cumulative_trapezoid_m3(
            hydro["Time_min"].to_numpy(float), hydro["Q_m3s"].to_numpy(float)
        )
        dense_series.append(hydro)

    metrics_df = pd.DataFrame(metrics_rows)
    mass_df = metrics_df[
        [
            "routing_mode",
            "rain_volume_m3",
            "outlet_volume_m3",
            "final_storage_m3",
            "mass_residual_m3",
            "mass_error_pct",
        ]
    ].copy()
    pass_fail = metrics_df[
        [
            "case_id",
            "routing_mode",
            "passed_hydrograph_shape",
            "passed_mass_balance",
            "passed_volume_error_5pct",
            "passed_overall_exact_config",
        ]
    ].copy()

    metrics_df.to_csv(ROUTING_DIR / "VTilted_Routing_Benchmark_Metrics.csv", index=False)
    mass_df.to_csv(ROUTING_DIR / "VTilted_Routing_Mass_Summary.csv", index=False)
    pass_fail.to_csv(ROUTING_DIR / "VTilted_Routing_Pass_Fail.csv", index=False)
    assert comparison is not None
    comparison.to_csv(ROUTING_DIR / "VTilted_Routing_Hydrograph_Comparison.csv", index=False)

    series_df = pd.concat(dense_series, ignore_index=True)
    series_df.to_csv(ROUTING_DIR / "VTilted_Routing_Model_TimeSeries.csv", index=False)

    make_hydrograph_figure(ref, dense_series)
    make_volume_figure(ref, dense_series)
    print(metrics_df.to_string(index=False))


def make_hydrograph_figure(ref: pd.DataFrame, dense_series: list[pd.DataFrame]) -> None:
    plt.figure(figsize=(8.2, 4.8))
    ref_180 = ref.loc[ref["Time_min"] <= 180]
    plt.plot(ref_180["Time_min"], ref_180["Benchmark_m3s"], "ko-", lw=2.2, ms=4, label="Benchmark")
    for hydro in dense_series:
        label = str(hydro["routing_label"].iloc[0])
        plt.plot(hydro["Time_min"], hydro["Q_m3s"], lw=1.7, label=label)
    plt.xlabel("Time (min)")
    plt.ylabel("Outlet discharge (m$^3$ s$^{-1}$)")
    plt.xlim(0, 180)
    plt.grid(True, alpha=0.25)
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "P1_HYDRO_VTILT_ROUTING_HYDROGRAPH.png", dpi=220)
    plt.close()


def make_volume_figure(ref: pd.DataFrame, dense_series: list[pd.DataFrame]) -> None:
    plt.figure(figsize=(8.2, 4.8))
    ref_180 = ref.loc[ref["Time_min"] <= 180].copy()
    ref_180["cumulative_volume_m3"] = cumulative_trapezoid_m3(
        ref_180["Time_min"].to_numpy(float), ref_180["Benchmark_m3s"].to_numpy(float)
    )
    plt.plot(ref_180["Time_min"], ref_180["cumulative_volume_m3"], "ko-", lw=2.2, ms=4, label="Benchmark")
    for hydro in dense_series:
        label = str(hydro["routing_label"].iloc[0])
        plt.plot(hydro["Time_min"], hydro["cumulative_outlet_volume_m3"], lw=1.7, label=label)
    plt.xlabel("Time (min)")
    plt.ylabel("Cumulative outlet volume (m$^3$)")
    plt.xlim(0, 180)
    plt.grid(True, alpha=0.25)
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "P1_HYDRO_VTILT_ROUTING_VOLUME.png", dpi=220)
    plt.close()


if __name__ == "__main__":
    main()
