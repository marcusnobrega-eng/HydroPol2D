#!/usr/bin/env python3
"""Analyze V-tilted 20 m vs 40 m ordinary vs 40 m subgrid runs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CASE_DIR = Path(__file__).resolve().parent
OUT = CASE_DIR / "Outputs" / "VTilted_40m_Comparison"
MAPS = OUT / "Maps"
FIGS = OUT / "Figures"
FIGS.mkdir(parents=True, exist_ok=True)


def read_hydro(scenario: str) -> pd.DataFrame:
    return pd.read_csv(OUT / scenario / "Hydrograph.csv")


def cumtrapz(t_min: np.ndarray, q: np.ndarray) -> float:
    if len(t_min) < 2:
        return 0.0
    return float(np.trapezoid(q, t_min * 60.0))


def score_hydro(test: pd.DataFrame, ref: pd.DataFrame, name: str) -> dict:
    t = ref["Time_min"].to_numpy(float)
    q_ref = ref["Q_m3s"].to_numpy(float)
    q_test = np.interp(t, test["Time_min"].to_numpy(float), test["Q_m3s"].to_numpy(float))
    err = q_test - q_ref
    rmse = float(np.sqrt(np.mean(err**2)))
    mae = float(np.mean(np.abs(err)))
    nse = float(1.0 - np.sum(err**2) / max(np.sum((q_ref - q_ref.mean()) ** 2), np.finfo(float).eps))
    vol_ref = cumtrapz(t, q_ref)
    vol_test = cumtrapz(t, q_test)
    peak_ref = float(np.max(q_ref))
    peak_test = float(np.max(q_test))
    t_peak_ref = float(t[np.argmax(q_ref)])
    t_peak_test = float(t[np.argmax(q_test)])
    return {
        "case": name,
        "hydrograph_rmse_m3s": rmse,
        "hydrograph_mae_m3s": mae,
        "hydrograph_nse": nse,
        "reference_outlet_volume_m3": vol_ref,
        "test_outlet_volume_m3": vol_test,
        "outlet_volume_error_pct": 100.0 * (vol_test - vol_ref) / max(vol_ref, np.finfo(float).eps),
        "reference_peak_m3s": peak_ref,
        "test_peak_m3s": peak_test,
        "peak_error_pct": 100.0 * (peak_test - peak_ref) / max(peak_ref, np.finfo(float).eps),
        "peak_time_error_min": abs(t_peak_test - t_peak_ref),
    }


def wet_metrics(ref: np.ndarray, pred: np.ndarray, threshold: float) -> tuple[float, float]:
    r = ref > threshold
    p = pred > threshold
    wet_area_error = 100.0 * (p.sum() - r.sum()) / max(r.sum(), 1)
    csi = (r & p).sum() / max((r | p).sum(), 1)
    return float(wet_area_error), float(csi)


def score_map(pred: np.ndarray, ref: np.ndarray, name: str) -> dict:
    err = pred - ref
    wet001, csi001 = wet_metrics(ref, pred, 0.01)
    wet010, csi010 = wet_metrics(ref, pred, 0.10)
    return {
        "case": name,
        "max_depth_rmse_m": float(np.sqrt(np.nanmean(err**2))),
        "max_depth_mae_m": float(np.nanmean(np.abs(err))),
        "max_depth_bias_m": float(np.nanmean(err)),
        "max_depth_max_error_m": float(np.nanmax(np.abs(err))),
        "reference_depth_sum_m3_like": float(np.nansum(ref) * 20.0 * 20.0),
        "test_depth_sum_m3_like": float(np.nansum(pred) * 20.0 * 20.0),
        "depth_sum_error_pct": float(100.0 * (np.nansum(pred) - np.nansum(ref)) / max(np.nansum(ref), np.finfo(float).eps)),
        "wet_area_error_pct_001m": wet001,
        "csi_001m": csi001,
        "wet_area_error_pct_010m": wet010,
        "csi_010m": csi010,
    }


def main() -> None:
    h_ref = read_hydro("reference20m")
    h_base = read_hydro("baseline40m")
    h_sub = read_hydro("subgrid40m")

    hydro_rows = [
        score_hydro(h_base, h_ref, "baseline40m"),
        score_hydro(h_sub, h_ref, "subgrid40m"),
    ]
    pd.DataFrame(hydro_rows).to_csv(OUT / "Hydrograph_Metrics.csv", index=False)

    comp = h_ref.rename(columns={"Q_m3s": "reference20m_m3s"}).copy()
    comp["baseline40m_m3s"] = np.interp(comp["Time_min"], h_base["Time_min"], h_base["Q_m3s"])
    comp["subgrid40m_m3s"] = np.interp(comp["Time_min"], h_sub["Time_min"], h_sub["Q_m3s"])
    comp.to_csv(OUT / "Hydrograph_Comparison.csv", index=False)

    ref = np.loadtxt(MAPS / "reference20m_max_depth_20m_m.csv", delimiter=",")
    base = np.loadtxt(MAPS / "baseline40m_projected_max_depth_20m_m.csv", delimiter=",")
    sub = np.loadtxt(MAPS / "subgrid40m_projected_max_depth_20m_m.csv", delimiter=",")
    map_rows = [
        score_map(base, ref, "baseline40m"),
        score_map(sub, ref, "subgrid40m"),
    ]
    pd.DataFrame(map_rows).to_csv(OUT / "MaxDepth_MapMetrics.csv", index=False)

    mass_rows = []
    for scenario in ["reference20m", "baseline40m", "subgrid40m"]:
        row = pd.read_csv(OUT / scenario / "Mass_Summary.csv").iloc[0].to_dict()
        row["case"] = scenario
        if scenario == "subgrid40m" and abs(row["rain_volume_m3"] - 25920.0) > 1e-6:
            row["rain_volume_m3_note"] = "reported by runner using subgrid wet area; corrected rainfall volume is 25920"
            row["rain_volume_m3_corrected"] = 25920.0
            row["mass_residual_m3_corrected"] = row["outlet_volume_m3"] + row["final_storage_m3"] - 25920.0
            row["mass_error_pct_corrected"] = 100.0 * row["mass_residual_m3_corrected"] / 25920.0
        mass_rows.append(row)
    pd.DataFrame(mass_rows).to_csv(OUT / "Mass_Summary_Comparison.csv", index=False)

    make_hydrograph_figure(comp)
    make_map_figure(ref, base, sub)
    make_scatter_figure(ref, base, sub)

    print(pd.DataFrame(hydro_rows).to_string(index=False))
    print(pd.DataFrame(map_rows).to_string(index=False))


def make_hydrograph_figure(comp: pd.DataFrame) -> None:
    plt.figure(figsize=(8.0, 4.8))
    plt.plot(comp["Time_min"], comp["reference20m_m3s"], "k-", lw=2.4, label="20 m reference")
    plt.plot(comp["Time_min"], comp["baseline40m_m3s"], lw=2.0, label="40 m ordinary")
    plt.plot(comp["Time_min"], comp["subgrid40m_m3s"], lw=2.0, label="40 m lookup subgrid")
    plt.xlabel("Time (min)")
    plt.ylabel("Outlet discharge (m$^3$ s$^{-1}$)")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGS / "P1_SUBGRID_VTILT_40M_HYDROGRAPH.png", dpi=220)
    plt.close()


def make_map_figure(ref: np.ndarray, base: np.ndarray, sub: np.ndarray) -> None:
    vmax = np.nanpercentile(ref, 99.5)
    vmax = max(vmax, 0.2)
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.2), constrained_layout=True)
    imgs = [
        (ref, "20 m reference", 0, vmax),
        (base, "40 m ordinary projected", 0, vmax),
        (sub, "40 m subgrid projected", 0, vmax),
        (base - ref, "ordinary - reference", -vmax, vmax),
        (sub - ref, "subgrid - reference", -vmax, vmax),
        (sub - base, "subgrid - ordinary", -vmax, vmax),
    ]
    for ax, (arr, title, vmin, vmax_i) in zip(axes.flat, imgs):
        cmap = "viridis" if vmin == 0 else "coolwarm"
        im = ax.imshow(arr, vmin=vmin, vmax=vmax_i, cmap=cmap)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, shrink=0.75)
    plt.savefig(FIGS / "P1_SUBGRID_VTILT_40M_MAX_DEPTH_MAPS.png", dpi=220)
    plt.close()


def make_scatter_figure(ref: np.ndarray, base: np.ndarray, sub: np.ndarray) -> None:
    mask = np.isfinite(ref) & np.isfinite(base) & np.isfinite(sub)
    lim = max(float(np.nanpercentile(ref[mask], 99.9)), float(np.nanpercentile(base[mask], 99.9)), 0.2)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), constrained_layout=True)
    for ax, pred, title in [(axes[0], base, "40 m ordinary"), (axes[1], sub, "40 m subgrid")]:
        ax.scatter(ref[mask], pred[mask], s=6, alpha=0.35)
        ax.plot([0, lim], [0, lim], "k--", lw=1)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel("20 m reference max depth (m)")
        ax.set_ylabel("Projected max depth (m)")
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    plt.savefig(FIGS / "P1_SUBGRID_VTILT_40M_DEPTH_SCATTER.png", dpi=220)
    plt.close()


if __name__ == "__main__":
    main()
