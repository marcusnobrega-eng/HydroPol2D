"""Python reference for VAL-POND-001: the HydroPol2D-Python detention-pond kernel.

Python repository: HydroPol2D-Python, branch release/studio-0.9.0rc2, commit
3cbd8ddb6c3767416279c32cb8b5d3da196413fe (src/hydropol2d/routing/detention.py).

Writes, next to this script:
  kernel_case.json              the inputs (grid, DEM, ponds, devices, forcing series),
                                which run_detention_ponds_check.m reads, so both codes
                                see the same case
  python_reference_geometry.csv footprints, faces and initial storage of each pond
  python_reference_steps.csv    one row per step: pond storage, stage, rain,
                                evaporation, seepage, recharge, critical inflow, device
                                discharge and pump switch, and the surface volume

The loop calls the public kernel functions in the engine's order: pond_forcing_step,
add_groundwater_exfiltration_to_ponds, then pond_exchange_step. Between them the
surface gets the rain of the ordinary cells, the exfiltration left outside the ponds
and, at one step, a volume poured on a destination cell (to reverse an orifice).
After the exchange every ordinary cell loses DRAIN_FRACTION of its volume, a stand-in
for routing: without it the destination cells fill until they are level with the
pond, and there the explicit orifice law (Q ~ sqrt(dH)) overshoots every step and
amplifies round-off by ~1.4 per step in both codes.

Run (from any folder outside the Python repository):
  PYTHONPATH=<HydroPol2D-Python>/src <HydroPol2D-Python>/.venv/bin/python make_python_reference.py
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402

import hydropol2d  # noqa: E402
from hydropol2d.config.models import (  # noqa: E402
    CaseConfig,
    DetentionPond,
    GridConfig,
    PondOrifice,
    PondPump,
    PondSpillway,
    StaticFields,
    TimeConfig,
)
from hydropol2d.routing.detention import (  # noqa: E402
    add_groundwater_exfiltration_to_ponds,
    initial_pond_storage_m3,
    pond_exchange_step,
    pond_forcing_step,
    pond_stages_m,
    prepare_detention_ponds,
)

HERE = Path(__file__).resolve().parent
NROWS, NCOLS, DX = 8, 10, 10.0
DT_S = 10.0
N_STEPS = 600
CELL_AREA = DX * DX
DRAIN_FRACTION = DT_S / 600.0


def python_commit() -> str:
    repo = Path(hydropol2d.__file__).resolve().parents[2]
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def elevation() -> np.ndarray:
    rows, cols = np.mgrid[0:NROWS, 0:NCOLS]
    z = 100.5 - 0.04 * cols - 0.02 * rows
    # Pond A: a 3 x 3 dug basin (rows 2-4, columns 1-3, 0-based).
    z[2:5, 1:4] = np.array(
        [[100.20, 100.10, 100.30], [100.00, 99.90, 100.15], [100.25, 100.05, 100.35]]
    )
    return z


def cell_square(row0: int, row1: int, col0: int, col1: int) -> list[list[float]]:
    """Polygon around rows row0..row1 and columns col0..col1 (0-based, inclusive)."""
    west, east = col0 * DX, (col1 + 1) * DX
    north, south = (NROWS - row0) * DX, (NROWS - row1 - 1) * DX
    # Inset by a quarter cell so no cell centre lies on an edge.
    return [
        [west + 2.5, south + 2.5],
        [east - 2.5, south + 2.5],
        [east - 2.5, north - 2.5],
        [west + 2.5, north - 2.5],
    ]


def pond_definitions() -> list[dict]:
    return [
        {
            "identifier": "pond_a",
            "vertices_m": cell_square(2, 4, 1, 3),
            "storage_method": "dtm",
            "initial_stage_m": 100.12,
            "lined": False,
            "bottom_ksat_mm_h": 20.0,
            "devices": [
                {"kind": "orifice", "identifier": "low_orifice", "destination_row": 3,
                 "destination_column": 4, "invert_elevation_m": 100.0,
                 "opening_area_m2": 0.05, "discharge_coefficient": 0.62},
                {"kind": "spillway", "identifier": "overflow", "destination_row": 5,
                 "destination_column": 2, "crest_elevation_m": 100.32,
                 "crest_length_m": 2.0, "discharge_coefficient": 1.7},
                {"kind": "pump", "identifier": "lift_pump", "destination_row": 0,
                 "destination_column": 2, "inlet_row": 3, "inlet_column": 2,
                 "maximum_flow_m3_s": 0.08, "shutoff_head_m": 3.0,
                 "additional_head_loss_m": 0.2, "start_stage_m": 100.33,
                 "stop_stage_m": 100.18},
            ],
        },
        {
            "identifier": "pond_b",
            "vertices_m": cell_square(1, 2, 6, 8),
            "storage_method": "power_law",
            "initial_stage_m": 100.60,
            "invert_elevation_m": 100.30,
            "storage_coefficient": 150.0,
            "storage_exponent": 1.5,
            "lined": True,
            "bottom_ksat_mm_h": 0.0,
            "devices": [
                {"kind": "orifice", "identifier": "b_orifice", "destination_row": 3,
                 "destination_column": 7, "invert_elevation_m": 100.40,
                 "opening_area_m2": 0.04, "discharge_coefficient": 0.62},
                {"kind": "spillway", "identifier": "b_spill", "destination_row": 1,
                 "destination_column": 9, "crest_elevation_m": 100.66,
                 "crest_length_m": 1.5, "discharge_coefficient": 1.7},
                {"kind": "orifice", "identifier": "b_closed", "destination_row": 0,
                 "destination_column": 5, "invert_elevation_m": 100.30,
                 "opening_area_m2": 0.10, "enabled": False},
            ],
        },
    ]


def build_device(definition: dict):
    kind = definition["kind"]
    fields = {key: value for key, value in definition.items() if key != "kind"}
    return {"orifice": PondOrifice, "spillway": PondSpillway, "pump": PondPump}[kind](**fields)


def forcing_series() -> dict:
    step = np.arange(N_STEPS)
    minutes = step * DT_S / 60.0
    rain = np.where(minutes < 5.0, 0.0, np.where(minutes < 25.0, 150.0, 0.0))
    evaporation = np.full(N_STEPS, 60.0)
    exfiltration_mm = np.where((minutes >= 35.0) & (minutes < 40.0), 0.2, 0.0)
    pour_m3 = np.where(step == int(47.5 * 60.0 / DT_S), 60.0, 0.0)
    return {
        "rain_mm_h": rain.tolist(),
        "evaporation_mm_day": evaporation.tolist(),
        "exfiltration_mm": exfiltration_mm.tolist(),
        "pour_m3": pour_m3.tolist(),
        "pour_row": 3,
        "pour_column": 4,
    }


def spatial_patterns() -> tuple[np.ndarray, np.ndarray]:
    rows, cols = np.mgrid[0:NROWS, 0:NCOLS]
    return 1.0 + 0.1 * cols, 1.0 + 0.05 * rows


def main() -> None:
    z = elevation()
    ponds = pond_definitions()
    forcing = forcing_series()
    rain_pattern, evaporation_pattern = spatial_patterns()
    initial_depth_m = 0.02
    case_inputs = {
        "python_commit": python_commit(),
        "nrows": NROWS,
        "ncols": NCOLS,
        "dx_m": DX,
        "dt_s": DT_S,
        "n_steps": N_STEPS,
        "x_origin_m": 0.0,
        "y_top_m": NROWS * DX,
        "elevation_m": z.tolist(),
        "rain_pattern": rain_pattern.tolist(),
        "evaporation_pattern": evaporation_pattern.tolist(),
        "initial_surface_depth_m": initial_depth_m,
        "drain_fraction": DRAIN_FRACTION,
        "ponds": ponds,
        "forcing": forcing,
        "note": "Rows and columns are 0-based here (Python); MATLAB adds 1.",
    }
    (HERE / "kernel_case.json").write_text(json.dumps(case_inputs, indent=1))

    config_ponds = tuple(
        DetentionPond(
            **{key: value for key, value in pond.items() if key != "devices"},
            devices=tuple(build_device(device) for device in pond["devices"]),
        )
        for pond in ponds
    )
    grid = GridConfig(nrows=NROWS, ncols=NCOLS, dx_m=DX)
    case = CaseConfig(
        name="val-pond-001",
        grid=grid,
        time=TimeConfig(end_s=N_STEPS * DT_S, dt_s=DT_S),
        detention_ponds=config_ponds,
    )
    static = StaticFields(
        grid=grid,
        transform=from_origin(0.0, NROWS * DX, DX, DX),
        valid_mask=np.ones((NROWS, NCOLS), dtype=bool),
        elevation_m=z,
    )
    prepared = prepare_detention_ponds(case, static)
    assert prepared is not None

    with (HERE / "python_reference_geometry.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["pond", "cells", "faces", "face_width_sum_m", "face_sill_sum_m",
             "initial_storage_m3", "initial_stage_m"]
        )
        storage0 = np.asarray(initial_pond_storage_m3(prepared))
        stage0 = np.asarray(pond_stages_m(storage0, prepared))
        for index, pond in enumerate(prepared.ponds):
            writer.writerow([
                pond.config.identifier,
                int(pond.rows.size),
                len(pond.critical_faces),
                repr(sum(face.width_m for face in pond.critical_faces)),
                repr(sum(face.sill_elevation_m for face in pond.critical_faces)),
                repr(float(storage0[index])),
                repr(float(stage0[index])),
            ])

    ordinary = ~prepared.pond_mask
    surface = np.where(ordinary, initial_depth_m * CELL_AREA, 0.0)
    storage = initial_pond_storage_m3(prepared)
    device_count = prepared.device_count
    pump_active = jnp.zeros((device_count,), dtype=bool)
    cumulative = jnp.zeros((device_count,), dtype=jnp.float64)
    terrain = jnp.asarray(z)
    names = [name.replace("/", "__") for name in prepared.device_identifiers]
    pond_names = [pond.config.identifier for pond in prepared.ponds]

    header = ["step", "time_s"]
    for name in pond_names:
        header += [f"{name}_{field}" for field in (
            "storage_m3", "stage_m", "rainfall_m3", "evaporation_m3", "seepage_m3",
            "recharge_m3", "exfiltration_m3", "lateral_m3")]
    for name in names:
        header += [f"{name}_q_m3_s", f"{name}_pump_on", f"{name}_cumulative_m3"]
    header += ["surface_m3", "ordinary_rain_m3", "pour_m3", "drained_m3"]

    rows = []
    for step in range(N_STEPS):
        rain = forcing["rain_mm_h"][step] * rain_pattern
        evaporation = forcing["evaporation_mm_day"][step] * evaporation_pattern
        result = pond_forcing_step(
            storage_m3=storage,
            rainfall_rate_mm_h=jnp.asarray(rain),
            open_water_evaporation_mm_day=jnp.asarray(evaporation),
            prepared=prepared,
            cell_area_m2=CELL_AREA,
            dt_s=DT_S,
        )
        storage = result.storage_m3
        recharge_m3 = np.zeros(len(pond_names))
        recharge_map = np.asarray(result.groundwater_recharge_rate_m_s)
        for index in range(len(pond_names)):
            recharge_m3[index] = float(
                np.sum(recharge_map[prepared.pond_index == index]) * CELL_AREA * DT_S
            )
        ordinary_rain = np.where(ordinary, rain / 1000.0 / 3600.0 * CELL_AREA * DT_S, 0.0)
        surface = surface + ordinary_rain

        exfiltration = np.full((NROWS, NCOLS), forcing["exfiltration_mm"][step])
        storage_before = np.asarray(storage)
        storage, remaining = add_groundwater_exfiltration_to_ponds(
            storage_m3=storage,
            exfiltration_depth_mm=jnp.asarray(exfiltration),
            prepared=prepared,
            cell_area_m2=CELL_AREA,
        )
        exfiltration_m3 = np.asarray(storage) - storage_before
        surface = surface + np.where(ordinary, np.asarray(remaining) / 1000.0 * CELL_AREA, 0.0)
        pour = forcing["pour_m3"][step]
        surface[forcing["pour_row"], forcing["pour_column"]] += pour

        exchange = pond_exchange_step(
            surface_volume_m3=jnp.asarray(surface),
            storage_m3=storage,
            pump_active=pump_active,
            cumulative_device_volume_m3=cumulative,
            terrain_elevation_m=terrain,
            prepared=prepared,
            cell_area_m2=CELL_AREA,
            dt_s=DT_S,
        )
        surface = np.asarray(exchange.surface_volume_m3)
        drained = np.where(ordinary, surface * DRAIN_FRACTION, 0.0)
        surface = surface - drained
        storage = exchange.storage_m3
        pump_active = exchange.pump_active
        cumulative = exchange.cumulative_device_volume_m3
        stage = np.asarray(pond_stages_m(storage, prepared))

        row = [step + 1, (step + 1) * DT_S]
        for index in range(len(pond_names)):
            row += [
                float(storage[index]), float(stage[index]),
                float(result.rainfall_m3[index]), float(result.evaporation_m3[index]),
                float(result.seepage_m3[index]), float(recharge_m3[index]),
                float(exfiltration_m3[index]), float(exchange.lateral_inflow_m3[index]),
            ]
        for index in range(device_count):
            row += [
                float(exchange.device_discharge_m3_s[index]),
                int(bool(pump_active[index])),
                float(cumulative[index]),
            ]
        row += [float(np.sum(surface)), float(np.sum(ordinary_rain)), pour, float(np.sum(drained))]
        rows.append(row)

    with (HERE / "python_reference_steps.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow([repr(value) if isinstance(value, float) else value for value in row])


if __name__ == "__main__":
    main()
