"""Python reference for VAL-GI-001 (MATLAB green infrastructure vs HydroPol2D-Python).

Made with HydroPol2D-Python, branch release/studio-0.9.0rc2; the commit used is
written to python_reference_info.csv (3cbd8dd when this reference was made).

Inputs (written by run_green_infrastructure_check.m with HYDROPOL_GI_WRITE_REPLAY=1):
  gi_case.json                 grid, terrain, designs, polygons, underdrain targets
  replay_inputs_<variant>.csv  per accepted MATLAB step: dt, host-cell volumes before
                               the GI step, gross rain and ET rates per compartment,
                               and the fraction of the native-exfiltration request
                               the MATLAB soil accepted
Outputs:
  python_geometry.csv                         compartments and connections of
                                              prepare_green_infrastructure
  python_reference_<variant>_compartments.csv states and cumulative fluxes every
                                              10th step and the last
  python_reference_<variant>_hosts.csv        host-cell volumes after the GI step
  python_limiter_as_is.csv                    what the unpatched kernel does
  limiter_case.json, python_limiter_case.csv  a 4 x 5 case (default designs, fixed
                                              8 m3 host volumes, 20 steps of 30 s)
                                              where the unpatched limiter overdraws
  python_reference_info.csv                   commit and settings

The replay is green_infrastructure_step plus the engine's coupling
(_advance_green_infrastructure): the soil takes the recorded fraction of the native
request, the rest returns to the gravel. Python's _limited_edge_transfer caps each
edge by its donor's storage on its own, so a donor with several lateral edges can go
negative and a gravel store can be filled above capacity, after which the clip at
the next step creates or deletes water (code review C14). MATLAB caps all outflows
of a donor together and gravel inflow by the free capacity. For the reference that
function, and only that function, is replaced by the MATLAB rule; the unpatched
kernel is run as well and its own mass error is reported.

Run from anywhere:
  /path/to/HydroPol2D-Python-gui-polish/.venv/bin/python make_python_reference.py \
      [--python-repo /path/to/HydroPol2D-Python-gui-polish]
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_REPO = Path("/Users/mngomes/Documents/GitHub/HydroPol2D-Python-gui-polish")
VARIANTS = ("legacy", "layered")
NAMES = (
    "surface_m3", "media_m3", "gravel_m3", "cum_rain_m3", "cum_capture_m3",
    "cum_infiltration_m3", "cum_overflow_m3", "cum_et_m3", "cum_native_request_m3",
    "cum_native_m3", "cum_underdrain_m3",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-repo", type=Path, default=DEFAULT_REPO)
    args = parser.parse_args()
    sys.path.insert(0, str(args.python_repo / "src"))
    import jax

    jax.config.update("jax_enable_x64", True)
    import hydropol2d.modules.green_infrastructure as gim

    commit = subprocess.run(
        ["git", "-C", str(args.python_repo), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    case = json.loads((HERE / "gi_case.json").read_text())
    prepared = build_prepared(case)
    write_geometry(prepared)

    host_rows = np.atleast_1d(np.asarray(case["host_rows"], dtype=int))
    host_cols = np.atleast_1d(np.asarray(case["host_columns"], dtype=int))
    original = gim._limited_edge_transfer
    as_is_rows = []
    for variant in VARIANTS:
        inputs = np.genfromtxt(HERE / f"replay_inputs_{variant}.csv", delimiter=",", names=True)
        gim._limited_edge_transfer = matlab_limiter(prepared)
        fixed = replay(gim, prepared, inputs, host_rows, host_cols)
        gim._limited_edge_transfer = original
        as_is = replay(gim, prepared, inputs, host_rows, host_cols)
        write_reference(variant, fixed)
        as_is_rows.append({
            "variant": variant,
            "steps": len(inputs),
            "matlab_rule_gi_mass_residual_m3": fixed["residual"],
            "python_as_is_gi_mass_residual_m3": as_is["residual"],
            "python_as_is_negative_store_steps": as_is["negative_steps"],
            "python_as_is_max_state_difference_m3": float(np.max(np.abs(
                as_is["table"][:, 2:5] - fixed["table"][:, 2:5]))),
        })
        print(variant, as_is_rows[-1])
    as_is_rows.append(limiter_case(gim, original))
    print("limiter_case", as_is_rows[-1])
    write_rows(HERE / "python_limiter_as_is.csv", as_is_rows)
    write_rows(HERE / "python_reference_info.csv", [{
        "python_repo": str(args.python_repo), "python_commit": commit,
        "kernel": "hydropol2d.modules.green_infrastructure.green_infrastructure_step",
        "limiter": "donor-total and gravel-capacity cap (MATLAB rule) in place of _limited_edge_transfer",
        "rows": "every 10th step and the last",
    }])
    return 0


def limiter_case(gim, original) -> dict:
    """Small case where Python's per-edge cap lets a compartment go negative."""
    from hydropol2d.modules.green_infrastructure import default_green_infrastructure_designs

    x = np.arange(5.0, 50.0, 10.0)
    y = np.arange(35.0, 0.0, -10.0)
    elevation = 10.0 - 0.01 * np.tile(x, (4, 1)) - 0.02 * np.tile((40.0 - y)[:, None], (1, 5))
    case = {
        "note": "Written by make_python_reference.py. Rows and columns are 0-based.",
        "nrows": 4, "ncols": 5, "dx_m": 10.0, "x_left_m": 0.0, "y_top_m": 40.0,
        "elevation_m": elevation.tolist(), "valid": np.ones((4, 5), bool).tolist(),
        "native_ksat_mm_h": np.full((4, 5), 5.0).tolist(),
        "designs": [d.model_dump() for d in default_green_infrastructure_designs()],
        "placements": [
            {"identifier": "bio", "design_identifier": "bioretention_default",
             "vertices_m": [[12, 33], [38, 31], [37, 24], [13, 26]],
             "underdrain_destination": {"kind": "surface_cell", "row": 3, "column": 4}},
            {"identifier": "pp", "design_identifier": "permeable_pavement_default",
             "vertices_m": [[3, 3], [18, 3], [18, 12], [3, 12]],
             "underdrain_destination": {"kind": "outside_domain", "row": -1, "column": -1}},
        ],
        "courant": 0.5, "minimum_substep_s": 0.01, "maximum_substeps": 10000,
        "steps": 20, "dt_s": 30.0, "host_m3": 8.0, "rain_mm_h": 60.0, "et_mm_day": 4.0,
    }
    (HERE / "limiter_case.json").write_text(json.dumps(case, indent=1))
    prepared = build_prepared(case)
    n = prepared.compartment_count
    # Every cell is a host column of the replay table, so the host volume resets to
    # 8 m3 in each cell at each step.
    host_rows = np.arange(4).repeat(5)
    host_cols = np.tile(np.arange(5), 4)
    dtype = [("t_min", float), ("dt_s", float)] + [(f"host_{i}", float) for i in range(20)] + \
        [(f"rain_{i}", float) for i in range(n)] + [(f"et_{i}", float) for i in range(n)] + \
        [(f"f_{i}", float) for i in range(n)]
    inputs = np.zeros(case["steps"], dtype=dtype)
    for name in inputs.dtype.names:
        inputs[name] = 1.0 if name.startswith("f_") else 0.0
    inputs["dt_s"] = case["dt_s"]
    for i in range(20):
        inputs[f"host_{i}"] = case["host_m3"]
    for i in range(n):
        inputs[f"rain_{i}"] = case["rain_mm_h"]
        inputs[f"et_{i}"] = case["et_mm_day"]
    gim._limited_edge_transfer = matlab_limiter(prepared)
    fixed = replay(gim, prepared, inputs, host_rows, host_cols, every_step=True)
    gim._limited_edge_transfer = original
    as_is = replay(gim, prepared, inputs, host_rows, host_cols, every_step=True)
    with (HERE / "python_limiter_case.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "compartment", "surface_m3", "media_m3", "gravel_m3",
                         "as_is_surface_m3", "as_is_media_m3", "as_is_gravel_m3"])
        for a, b in zip(fixed["table"], as_is["table"]):
            writer.writerow([int(a[0]), int(a[1]), *(f"{v:.17g}" for v in a[2:5]),
                             *(f"{v:.17g}" for v in b[2:5])])
    return {
        "variant": "limiter_case", "steps": case["steps"],
        "matlab_rule_gi_mass_residual_m3": fixed["residual"],
        "python_as_is_gi_mass_residual_m3": as_is["residual"],
        "python_as_is_negative_store_steps": as_is["negative_steps"],
        "python_as_is_max_state_difference_m3": float(np.max(np.abs(
            as_is["table"][:, 2:5] - fixed["table"][:, 2:5]))),
    }


def build_prepared(case: dict):
    from rasterio.transform import from_origin

    from hydropol2d.config.models import (
        CaseConfig, GreenInfrastructureConfig, GreenInfrastructureDesign,
        GreenInfrastructurePlacement, GreenInfrastructureUnderdrainDestination,
        GridConfig, ModelFlags, StaticFields, TimeConfig,
    )
    from hydropol2d.modules.green_infrastructure import prepare_green_infrastructure

    designs = case["designs"] if isinstance(case["designs"], list) else [case["designs"]]
    design_objects = tuple(GreenInfrastructureDesign(**d) for d in designs)
    placements = []
    for item in case["placements"]:
        dest = item["underdrain_destination"]
        if dest["kind"] == "none":
            destination = None
        elif dest["kind"] == "surface_cell":
            destination = GreenInfrastructureUnderdrainDestination(
                kind="surface_cell", row=int(dest["row"]), column=int(dest["column"]))
        else:
            destination = GreenInfrastructureUnderdrainDestination(kind=dest["kind"])
        placements.append(GreenInfrastructurePlacement(
            identifier=item["identifier"], design_identifier=item["design_identifier"],
            vertices_m=tuple(tuple(float(v) for v in row) for row in item["vertices_m"]),
            underdrain_destination=destination))
    grid = GridConfig(nrows=case["nrows"], ncols=case["ncols"], dx_m=case["dx_m"])
    config = CaseConfig(
        name="val-gi-001", grid=grid, time=TimeConfig(end_s=10.0, dt_s=1.0),
        flags=ModelFlags(green_infrastructure=True),
        green_infrastructure=GreenInfrastructureConfig(
            designs=design_objects, placements=tuple(placements), courant=case["courant"],
            minimum_substep_s=case["minimum_substep_s"],
            maximum_substeps=int(case["maximum_substeps"])),
    )
    static = StaticFields(
        grid=grid,
        transform=from_origin(case["x_left_m"], case["y_top_m"], case["dx_m"], case["dx_m"]),
        valid_mask=np.asarray(case["valid"], dtype=bool),
        elevation_m=np.asarray(case["elevation_m"], dtype=np.float64),
    )
    return prepare_green_infrastructure(
        config, static, native_ksat_mm_h=np.asarray(case["native_ksat_mm_h"], dtype=np.float64))


def matlab_limiter(prepared):
    """_limited_edge_transfer with MATLAB's caps; calls alternate surface, gravel."""
    import jax.numpy as jnp

    gravel_capacity = jnp.asarray(
        prepared.area_m2 * prepared.gravel_depth_m * prepared.gravel_void_fraction)
    calls = {"n": 0}

    def limited(volume, start, end, proposed_q_m3_s, dt_s):
        capacity = gravel_capacity if calls["n"] % 2 == 1 else None
        calls["n"] += 1
        donor = jnp.where(proposed_q_m3_s >= 0.0, start, end)
        receiver = jnp.where(proposed_q_m3_s >= 0.0, end, start)
        amount = jnp.abs(proposed_q_m3_s) * dt_s
        outflow = jnp.zeros_like(volume).at[donor].add(amount)
        scale = jnp.where(outflow > 0, jnp.minimum(
            1.0, jnp.maximum(volume, 0.0) / jnp.where(outflow > 0, outflow, 1.0)), 1.0)
        amount = amount * scale[donor]
        if capacity is not None:
            inflow = jnp.zeros_like(volume).at[receiver].add(amount)
            room = jnp.maximum(capacity - volume, 0.0)
            scale = jnp.where(inflow > 0, jnp.minimum(
                1.0, room / jnp.where(inflow > 0, inflow, 1.0)), 1.0)
            amount = amount * scale[receiver]
        transfer = jnp.where(proposed_q_m3_s < 0.0, -amount, amount)
        return volume.at[start].add(-transfer).at[end].add(transfer), transfer / dt_s

    return limited


def replay(gim, prepared, inputs, host_rows, host_cols, every_step: bool = False) -> dict:
    import jax.numpy as jnp

    n = prepared.compartment_count
    nh = host_rows.size
    rows, cols = prepared.rows, prepared.columns
    surface = np.zeros(n)
    media = np.asarray(prepared.initial_media_water_volume_m3, dtype=np.float64)
    gravel = np.zeros(n)
    edge_q = np.zeros(prepared.connection_count)
    cumulative = np.zeros((n, 8))
    external = 0.0
    residual = 0.0
    negative_steps = 0
    keep = set(range(10, len(inputs) + 1, 10)) | {len(inputs)}
    if every_step:
        keep = set(range(1, len(inputs) + 1))
    table, hosts = [], []
    fields = inputs.dtype.names
    for step, record in enumerate(inputs, start=1):
        values = np.asarray([record[name] for name in fields], dtype=np.float64)
        dt = values[1]
        host_in = values[2:2 + nh]
        rain = values[2 + nh:2 + nh + n]
        et = values[2 + nh + n:2 + nh + 2 * n]
        fraction = values[2 + nh + 2 * n:2 + nh + 3 * n]
        host = np.zeros(prepared.shape)
        host[host_rows, host_cols] = host_in
        rain_grid = np.zeros(prepared.shape)
        rain_grid[rows, cols] = rain
        et_grid = np.zeros(prepared.shape)
        et_grid[rows, cols] = et
        before = surface.sum() + media.sum() + gravel.sum()
        result = gim.green_infrastructure_step(
            prepared=prepared, host_surface_volume_m3=jnp.asarray(host),
            surface_volume_m3=jnp.asarray(surface), media_water_volume_m3=jnp.asarray(media),
            gravel_volume_m3=jnp.asarray(gravel),
            surface_connection_discharge_m3_s=jnp.asarray(edge_q),
            rainfall_rate_mm_h=jnp.asarray(rain_grid), et_rate_mm_day=jnp.asarray(et_grid),
            apply_et_crop_coefficient=True, drainage_node_count=0, dt_s=float(dt))
        request = np.asarray(result.native_exfiltration_m3)
        native = request * fraction
        surface = np.asarray(result.surface_volume_m3)
        media = np.asarray(result.media_water_volume_m3)
        gravel = np.asarray(result.gravel_volume_m3) + np.maximum(request - native, 0.0)
        edge_q = np.asarray(result.surface_connection_discharge_m3_s)
        fluxes = np.column_stack([
            result.rainfall_m3, result.capture_m3, result.media_infiltration_m3,
            result.overflow_m3, result.evapotranspiration_m3, request, native,
            result.underdrain_m3])
        cumulative += fluxes
        external += float(result.external_underdrain_m3)
        after = surface.sum() + media.sum() + gravel.sum()
        gain = fluxes[:, 0].sum() + fluxes[:, 1].sum() - fluxes[:, 3].sum() - \
            fluxes[:, 4].sum() - fluxes[:, 6].sum() - fluxes[:, 7].sum()
        residual += after - before - gain
        if min(surface.min(), media.min(), gravel.min()) < -1e-12:
            negative_steps += 1
        if step in keep:
            host_after = np.asarray(result.host_surface_volume_m3)[host_rows, host_cols]
            for c in range(n):
                table.append([step, c + 1, surface[c], media[c], gravel[c], *cumulative[c]])
            for h in range(nh):
                hosts.append([step, h + 1, host_after[h], external])
    return {"table": np.asarray(table), "hosts": np.asarray(hosts), "residual": residual,
            "negative_steps": negative_steps}


def write_reference(variant: str, result: dict) -> None:
    with (HERE / f"python_reference_{variant}_compartments.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "compartment", *NAMES])
        for row in result["table"]:
            writer.writerow([int(row[0]), int(row[1]), *(f"{v:.17g}" for v in row[2:])])
    with (HERE / f"python_reference_{variant}_hosts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "host", "host_after_m3", "cumulative_external_m3"])
        for row in result["hosts"]:
            writer.writerow([int(row[0]), int(row[1]), f"{row[2]:.17g}", f"{row[3]:.17g}"])


def write_geometry(prepared) -> None:
    rows = []
    for c in range(prepared.compartment_count):
        rows.append({"kind": "compartment", "row": int(prepared.rows[c]) + 1,
                     "column": int(prepared.columns[c]) + 1,
                     "facility": int(prepared.facility_index[c]) + 1,
                     "value": f"{prepared.area_m2[c]:.17g}"})
    for e in range(prepared.connection_count):
        rows.append({"kind": "edge", "row": int(prepared.edge_start[e]) + 1,
                     "column": int(prepared.edge_end[e]) + 1, "facility": "",
                     "value": f"{prepared.edge_width_m[e]:.17g}"})
    write_rows(HERE / "python_geometry.csv", rows)


def write_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
