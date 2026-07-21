#!/usr/bin/env python3
"""Generate simple the analytical/reference solutions for HydroPol2D.

These outputs are not HydroPol2D model results. They are reference truths used
to compare formula implementations, storage ledgers, and simple dynamics.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import math


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"no rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def canopy_bucket(outdir: Path) -> None:
    cell_area_m2 = 400.0
    coefficient_mm_per_lai = 0.2
    rainfall_mm = [0.0, 0.08, 0.50, 0.0, 0.0, 0.0]
    potential_evap_mm = [0.0, 0.0, 0.0, 0.10, 0.20, 0.20]
    zones = [
        ("vtilted_left_hillslope", 1, 0.5),
        ("vtilted_channel_strip", 2, 0.0),
        ("vtilted_right_hillslope", 3, 2.0),
    ]
    rows = []
    for scenario, zone_id, lai in zones:
        storage_mm = 0.0
        smax_mm = coefficient_mm_per_lai * lai
        for i, (rain, ep) in enumerate(zip(rainfall_mm, potential_evap_mm)):
            beta = storage_mm / smax_mm if smax_mm > 0 else 0.0
            evaporation_mm = min(beta * ep, rain + storage_mm)
            if smax_mm <= 0.0:
                evaporation_mm = 0.0
            provisional_storage_mm = storage_mm + rain - evaporation_mm
            throughfall_mm = max(provisional_storage_mm - smax_mm, 0.0)
            storage_next_mm = min(provisional_storage_mm, smax_mm)
            residual_mm = (
                storage_next_mm
                - storage_mm
                - (rain - evaporation_mm - throughfall_mm)
            )
            rows.append(
                {
                    "scenario": scenario,
                    "zone_id": zone_id,
                    "time_s": i * 3600.0,
                    "cell_area_m2": cell_area_m2,
                    "lai": lai,
                    "coefficient_mm_per_lai": coefficient_mm_per_lai,
                    "smax_mm": smax_mm,
                    "gross_rainfall_mm": rain,
                    "potential_evaporation_mm": ep,
                    "evaporation_mm": evaporation_mm,
                    "canopy_storage_mm": storage_next_mm,
                    "throughfall_mm": throughfall_mm,
                    "mass_residual_mm": residual_mm,
                    "gross_rainfall_m3": rain / 1000.0 * cell_area_m2,
                    "evaporation_m3": evaporation_mm / 1000.0 * cell_area_m2,
                    "canopy_storage_m3": storage_next_mm / 1000.0 * cell_area_m2,
                    "throughfall_m3": throughfall_mm / 1000.0 * cell_area_m2,
                    "mass_residual_m3": residual_mm / 1000.0 * cell_area_m2,
                }
            )
            storage_mm = storage_next_mm
    write_csv(outdir / "VAL-CANOPY-001_reference.csv", rows)


def green_ampt_cumulative(t_s: float, ks_m_s: float, psi_dtheta_m: float) -> float:
    if t_s <= 0:
        return 0.0
    # Solve F - psi*dtheta*ln(1 + F/(psi*dtheta)) = Ks*t by bisection.
    target = ks_m_s * t_s
    lo = 0.0
    hi = max(ks_m_s * t_s + psi_dtheta_m, 1e-9)
    def f(value: float) -> float:
        return value - psi_dtheta_m * math.log1p(value / psi_dtheta_m) - target
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def green_ampt(outdir: Path) -> None:
    ks = 1.0e-6
    psi_dtheta = 0.055
    rows = []
    previous_f = 0.0
    previous_t = 0.0
    for step in range(0, 49):
        t = step * 1800.0
        cumulative = green_ampt_cumulative(t, ks, psi_dtheta)
        rate = 0.0 if step == 0 else (cumulative - previous_f) / (t - previous_t)
        rows.append(
            {
                "time_s": t,
                "cumulative_infiltration_m": cumulative,
                "infiltration_rate_m_s": rate,
                "Ks_m_s": ks,
                "psi_delta_theta_m": psi_dtheta,
            }
        )
        previous_f = cumulative
        previous_t = t
    write_csv(outdir / "VAL-INFIL-GA-001_reference.csv", rows)


def philip_infiltration(outdir: Path) -> None:
    sorptivity = 2.5e-4
    k_term = 5.0e-7
    rows = []
    for step in range(0, 49):
        t = step * 1800.0
        cumulative = sorptivity * math.sqrt(t) + k_term * t if t > 0 else 0.0
        rate = 0.5 * sorptivity / math.sqrt(t) + k_term if t > 0 else 0.0
        rows.append(
            {
                "time_s": t,
                "cumulative_infiltration_m": cumulative,
                "infiltration_rate_m_s": rate,
                "sorptivity_m_sqrt_s": sorptivity,
                "K_term_m_s": k_term,
            }
        )
    write_csv(outdir / "VAL-INFIL-PHILIP-001_reference.csv", rows)


def et_availability(outdir: Path) -> None:
    storage = 0.006
    wilting_storage = 0.001
    demands = [0.0010, 0.0025, 0.0025, 0.0010, 0.0010]
    rows = []
    for day, demand in enumerate(demands):
        available = max(0.0, storage - wilting_storage)
        actual = min(demand, available)
        next_storage = storage - actual
        rows.append(
            {
                "time_day": day,
                "et_demand_m": demand,
                "actual_et_m": actual,
                "soil_storage_m": next_storage,
                "wilting_storage_m": wilting_storage,
                "storage_residual_m": storage - actual - next_storage,
            }
        )
        storage = next_storage
    write_csv(outdir / "VAL-ET-001_reference.csv", rows)


def linear_reservoir(outdir: Path) -> None:
    storage0 = 1000.0
    k = 2.5e-5
    rows = []
    for step in range(0, 49):
        t = step * 3600.0
        storage = storage0 * math.exp(-k * t)
        q = k * storage
        rows.append(
            {
                "time_s": t,
                "storage_m3": storage,
                "discharge_m3_s": q,
                "recession_constant_s_1": k,
            }
        )
    write_csv(outdir / "VAL-GW-LINRES-001_reference.csv", rows)


def manning_rectangular(outdir: Path) -> None:
    width = 5.0
    slope = 0.001
    n = 0.035
    rows = []
    for i in range(1, 21):
        depth = 0.05 * i
        area = width * depth
        perimeter = width + 2.0 * depth
        radius = area / perimeter
        discharge = (1.0 / n) * area * (radius ** (2.0 / 3.0)) * math.sqrt(slope)
        rows.append(
            {
                "depth_m": depth,
                "width_m": width,
                "slope_m_m": slope,
                "manning_n": n,
                "area_m2": area,
                "hydraulic_radius_m": radius,
                "discharge_m3_s": discharge,
            }
        )
    write_csv(outdir / "VAL-MANNING-RECTANGULAR_reference.csv", rows)
    write_csv(outdir / "VAL-HYDRO-STEADY-001_reference.csv", rows)
    write_csv(outdir / "VAL-SUBGRID-001_reference.csv", rows)


def reservoir_rating_curve(outdir: Path) -> None:
    area = 1000.0
    h0 = 0.0
    h_initial = 2.0
    a = 0.4
    exponent = 1.0
    dt_s = 300.0
    h = h_initial
    rows = []
    for step in range(0, 49):
        q = a * max(0.0, h - h0) ** exponent
        storage = area * h
        rows.append(
            {
                "time_s": step * dt_s,
                "stage_m": h,
                "storage_m3": storage,
                "outflow_m3_s": q,
                "rating_a": a,
                "rating_exponent": exponent,
            }
        )
        h = h * math.exp(-(a / area) * dt_s)
    write_csv(outdir / "VAL-RES-001_reference.csv", rows)


def hydrograph_volume(outdir: Path) -> None:
    dt_s = 600.0
    q_values = [0.0, 1.0, 3.0, 5.0, 3.0, 1.0, 0.0]
    cumulative = 0.0
    rows = []
    for i, q in enumerate(q_values):
        volume = q * dt_s
        cumulative += volume
        rows.append(
            {
                "time_s": i * dt_s,
                "inflow_m3_s": q,
                "interval_volume_m3": volume,
                "cumulative_volume_m3": cumulative,
            }
        )
    write_csv(outdir / "VAL-BC-INFLOW-001_reference.csv", rows)


def stage_volume(outdir: Path) -> None:
    length = 1000.0
    width = 10.0
    bed = 0.0
    rows = []
    for i in range(0, 13):
        stage = 0.05 * i
        depth = max(0.0, stage - bed)
        volume = length * width * depth
        rows.append(
            {
                "stage_m": stage,
                "depth_m": depth,
                "length_m": length,
                "width_m": width,
                "volume_m3": volume,
            }
        )
    write_csv(outdir / "VAL-BC-STAGE-001_reference.csv", rows)


def raster_rainfall_totals(outdir: Path) -> None:
    cell_area = 100.0
    rasters_m = [
        [[0.001, 0.002], [0.003, 0.004]],
        [[0.000, 0.001], [0.002, 0.003]],
    ]
    rows = []
    total_volume = 0.0
    for t, raster in enumerate(rasters_m):
        for row_i, row in enumerate(raster):
            for col_i, depth in enumerate(row):
                volume = depth * cell_area
                total_volume += volume
                rows.append(
                    {
                        "time_index": t,
                        "row": row_i,
                        "col": col_i,
                        "rainfall_depth_m": depth,
                        "cell_area_m2": cell_area,
                        "rainfall_volume_m3": volume,
                        "cumulative_domain_volume_m3": total_volume,
                    }
                )
    write_csv(outdir / "VAL-RAIN-MAP-001_reference.csv", rows)


def washoff(outdir: Path) -> None:
    mass0 = 10.0
    build_rate = 0.02
    runoff = 1.5
    k = 0.08
    decay = k * runoff
    dt_s = 600.0
    mass = mass0
    rows = []
    for step in range(0, 49):
        t = step * dt_s
        washoff_rate = decay * mass
        rows.append(
            {
                "time_s": t,
                "surface_mass_kg": mass,
                "build_rate_kg_s": build_rate,
                "washoff_rate_kg_s": washoff_rate,
                "runoff_reference": runoff,
                "washoff_coefficient": k,
            }
        )
        steady = build_rate / decay
        mass = steady + (mass - steady) * math.exp(-decay * dt_s)
    write_csv(outdir / "VAL-WQ-001_reference.csv", rows)


def human_risk_thresholds(outdir: Path) -> None:
    rows = []
    depths = [0.0, 0.1, 0.3, 0.6, 1.0]
    velocities = [0.0, 0.5, 1.0, 2.0]
    for depth in depths:
        for velocity in velocities:
            dv = depth * velocity
            if depth < 0.1 and dv < 0.05:
                risk = "low"
            elif depth < 0.5 and dv < 0.4:
                risk = "moderate"
            elif depth < 1.0 and dv < 1.0:
                risk = "high"
            else:
                risk = "extreme"
            rows.append(
                {
                    "depth_m": depth,
                    "velocity_m_s": velocity,
                    "depth_velocity_m2_s": dv,
                    "expected_class": risk,
                }
            )
    write_csv(outdir / "VAL-HR-001_reference.csv", rows)


CASES = {
    "canopy_bucket": canopy_bucket,
    "green_ampt": green_ampt,
    "philip_infiltration": philip_infiltration,
    "et_availability": et_availability,
    "linear_reservoir": linear_reservoir,
    "manning_rectangular": manning_rectangular,
    "reservoir_rating_curve": reservoir_rating_curve,
    "hydrograph_volume": hydrograph_volume,
    "stage_volume": stage_volume,
    "raster_rainfall_totals": raster_rainfall_totals,
    "washoff": washoff,
    "human_risk_thresholds": human_risk_thresholds,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=sorted(CASES) + ["all"],
        default="all",
        help="Reference case to generate.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("HydroPol2D_Model/Validation/Reference_Outputs/Validation"),
        help="Output directory for reference CSVs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    selected = CASES if args.case == "all" else {args.case: CASES[args.case]}
    for name, generator in selected.items():
        case_outdir = args.output / name
        generator(case_outdir)
        print(f"generated {name}: {case_outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
