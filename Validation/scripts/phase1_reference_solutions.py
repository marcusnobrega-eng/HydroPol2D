#!/usr/bin/env python3
"""Generate simple Phase 1 analytical/reference solutions for HydroPol2D.

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
    write_csv(outdir / "P1-CANOPY-001_reference.csv", rows)


def snow_degree_day(outdir: Path) -> None:
    cell_area_m2 = 400.0
    alpha = 0.8
    epsilon = 0.98
    c_e = 0.001
    ddf_mm_c_day = 2.0
    t_thresh_c = 0.0
    rho_snow_init = 100.0
    rho_max = 400.0
    k_t = 0.1
    k_swe = 0.001
    k_d = 0.02
    doy = 30.0
    lat_deg = 39.0
    wind_m_s = 2.0
    precip_mm = [10.0, 15.0, 0.0, 20.0, 0.0, 0.0]
    t_min_offset_c = 4.0
    zones = [
        ("vtilted_left_hillslope", 1, 2000, 125.0, [-5.0, 0.0, 2.0, 5.0, 8.0, 10.0]),
        ("vtilted_channel_strip", 2, 100, 100.0, [1.0, 5.0, 8.0, 2.0, 6.5, 9.0]),
        ("vtilted_right_hillslope", 3, 1950, 150.0, [-3.0, 4.5, 7.0, 1.0, 5.5, 8.0]),
    ]

    def radiation(t_air_c: float) -> tuple[float, float, float]:
        sigma = 5.67e-8
        t_kelvin = t_air_c + 273.15
        decl = 23.45 * math.sin(math.radians(360.0 * (doy - 81.0) / 365.0))
        theta_z = math.acos(
            math.sin(math.radians(lat_deg)) * math.sin(math.radians(decl))
            + math.cos(math.radians(lat_deg)) * math.cos(math.radians(decl))
        )
        s_down = 0.75 * 1367.0 * math.cos(theta_z)
        s_down = max(s_down, 0.0)
        epsilon_a = 0.642 + 0.035 * math.sqrt(max(t_air_c, 0.0))
        l_down = epsilon_a * sigma * t_kelvin**4
        t_snow_kelvin = 273.15
        l_up = epsilon * sigma * t_snow_kelvin**4
        q_net = (
            s_down * (1.0 - alpha)
            - l_up
            + l_down
            - epsilon * sigma * (t_snow_kelvin**4 - t_kelvin**4)
        )
        return s_down, l_down, q_net

    def snow_partition(precip: float, t_air_c: float) -> float:
        # Match Snow_Model_Function.m: hard-coded linear transition from 4 to 7 C.
        fraction = (7.0 - t_air_c) / (7.0 - 4.0)
        fraction = min(max(fraction, 0.0), 1.0)
        return precip * fraction

    def sublimation(wind: float, t_air_c: float, t_min_c: float, swe_prev: float, swe_after_melt: float) -> float:
        e_s = 6.112 * math.exp((17.67 * t_air_c) / (t_air_c + 243.5))
        e_a = 10.0 * (0.61 * math.exp((17.27 * t_min_c) / (t_min_c + 237.3)))
        q_s = 0.622 * e_s / (1013.25 - e_s)
        q_a = 0.622 * e_a / (1013.25 - e_a)
        e_sub = c_e * wind * (q_s - q_a) * swe_prev
        e_sub = max(e_sub, 0.0)
        return min(e_sub, swe_after_melt)

    rows = []
    for scenario, zone_id, zone_cells, h_snow_t0_mm, temp_c in zones:
        swe_prev = 0.0
        h_snow_prev = 0.0
        for i, (p_mm, t_air_c) in enumerate(zip(precip_mm, temp_c)):
            t_min_c = t_air_c - t_min_offset_c
            p_snow_mm = snow_partition(p_mm, t_air_c)
            p_rain_mm = p_mm - p_snow_mm
            swe_with_snow_mm = swe_prev + p_snow_mm
            s_down_w_m2, l_down_w_m2, q_net_w_m2 = radiation(t_air_c)
            rho_snow = rho_snow_init + k_t * t_air_c + k_swe * swe_prev + k_d * h_snow_t0_mm
            rho_snow = min(rho_snow, rho_max)
            h_snow_mm = swe_prev / rho_snow if swe_prev != 0.0 else 0.0
            m_temp_mm = max(ddf_mm_c_day * t_air_c, 0.0)
            m_rad_mm = (1.0 - alpha) * q_net_w_m2 / 334.0
            melt_mm = max(m_temp_mm + m_rad_mm, 0.0)
            melt_mm = min(melt_mm, swe_with_snow_mm)
            swe_after_melt_mm = swe_with_snow_mm - melt_mm
            e_sub_mm = sublimation(wind_m_s, t_air_c, t_min_c, swe_prev, swe_after_melt_mm)
            swe_next_mm = swe_after_melt_mm - e_sub_mm
            residual_mm = swe_prev + p_snow_mm - melt_mm - e_sub_mm - swe_next_mm
            rows.append(
                {
                    "scenario": scenario,
                    "zone_id": zone_id,
                    "zone_cells": zone_cells,
                    "time_day": i,
                    "cell_area_m2": cell_area_m2,
                    "temperature_c": t_air_c,
                    "t_min_c": t_min_c,
                    "precipitation_mm": p_mm,
                    "wind_m_s": wind_m_s,
                    "lat_deg": lat_deg,
                    "doy": doy,
                    "h_snow_t0_mm": h_snow_t0_mm,
                    "alpha": alpha,
                    "epsilon": epsilon,
                    "c_e": c_e,
                    "ddf_mm_c_day": ddf_mm_c_day,
                    "t_thresh_c": t_thresh_c,
                    "rho_snow_init_kg_m3": rho_snow_init,
                    "rho_max_kg_m3": rho_max,
                    "k_t": k_t,
                    "k_swe": k_swe,
                    "k_d": k_d,
                    "s_down_w_m2": s_down_w_m2,
                    "l_down_w_m2": l_down_w_m2,
                    "q_net_w_m2": q_net_w_m2,
                    "rain_mm": p_rain_mm,
                    "snowfall_mm": p_snow_mm,
                    "melt_mm": melt_mm,
                    "sublimation_mm": e_sub_mm,
                    "swe_mm": swe_next_mm,
                    "snow_depth_mm": h_snow_mm,
                    "rho_snow_kg_m3": rho_snow,
                    "mass_residual_mm": residual_mm,
                    "rain_m3": p_rain_mm / 1000.0 * cell_area_m2,
                    "snowfall_m3": p_snow_mm / 1000.0 * cell_area_m2,
                    "melt_m3": melt_mm / 1000.0 * cell_area_m2,
                    "sublimation_m3": e_sub_mm / 1000.0 * cell_area_m2,
                    "swe_m3": swe_next_mm / 1000.0 * cell_area_m2,
                    "mass_residual_m3": residual_mm / 1000.0 * cell_area_m2,
                }
            )
            swe_prev = swe_next_mm
            h_snow_prev = h_snow_mm
    write_csv(outdir / "P1-SNOW-001_reference.csv", rows)


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
    write_csv(outdir / "P1-INFIL-GA-001_reference.csv", rows)


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
    write_csv(outdir / "P1-INFIL-PHILIP-001_reference.csv", rows)


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
    write_csv(outdir / "P1-ET-001_reference.csv", rows)


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
    write_csv(outdir / "P1-GW-LINRES-001_reference.csv", rows)


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
    write_csv(outdir / "P1-MANNING-RECTANGULAR_reference.csv", rows)
    write_csv(outdir / "P1-HYDRO-STEADY-001_reference.csv", rows)
    write_csv(outdir / "P1-SUBGRID-001_reference.csv", rows)


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
    write_csv(outdir / "P1-RES-001_reference.csv", rows)


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
    write_csv(outdir / "P1-BC-INFLOW-001_reference.csv", rows)


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
    write_csv(outdir / "P1-BC-STAGE-001_reference.csv", rows)


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
    write_csv(outdir / "P1-RAIN-MAP-001_reference.csv", rows)


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
    write_csv(outdir / "P1-WQ-001_reference.csv", rows)


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
    write_csv(outdir / "P1-HR-001_reference.csv", rows)


CASES = {
    "canopy_bucket": canopy_bucket,
    "snow_degree_day": snow_degree_day,
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
        default=Path("HydroPol2D_Model/Validation/Reference_Outputs/Phase1"),
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
