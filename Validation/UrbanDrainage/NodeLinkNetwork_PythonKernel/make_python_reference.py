"""Write the drainage network and the HydroPol2D-Python reference for VAL-UD-001.

Run with the HydroPol2D-Python virtual environment, for example:

    PY=/Users/mngomes/Documents/GitHub/HydroPol2D-Python-gui-polish
    $PY/.venv/bin/python make_python_reference.py --python-repo $PY

The script only imports HydroPol2D-Python; it writes nothing there.

Outputs (in this folder):
    network_kernel.gpkg        the network, written by write_drainage_geopackage
    python_reference_kernel.csv one row per 5 s surface step: the sources that
                               were added, then node heads and volumes, pipe and
                               inlet discharges, inlet-cell volumes and outfall
                               volumes after the drainage advance
    python_reference_inputs.csv the 4 x 5 terrain grid
    network_grid.gpkg, Grid/DEM.tif, Grid/LULC.tif, Grid/SOIL.tif
                               the MATLAB grid run's inputs (16 x 24 cells of 5 m);
                               the network passes Python's audit and
                               prepare_urban_drainage on this DEM

The advance is the engine's own sub-cycle (_compiled_urban_drainage_advance),
with prepare_urban_drainage(maximum_velocity_m_s=10) and the default
hydraulics (courant 0.5, min_dt 0.01 s, max_dt 5 s, 10 000 substeps).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
DT_S = 5.0
STEPS = 400
DX_M = 10.0
TERRAIN_M = np.array(
    [
        [10.0, 9.9, 9.8, 9.7, 9.6],
        [9.8, 9.7, 9.5, 9.4, 9.3],
        [9.4, 9.3, 9.2, 9.0, 8.9],
        [7.6, 7.5, 7.7, 7.75, 7.8],
    ]
)
SOURCE_CELLS = ((0, 0), (1, 2))
INITIAL_POND = ((3, 4), 30.0)  # cell, m3


def source_rates_m3_s(t_mid_s: float) -> tuple[float, float]:
    """Triangular 0.25 m3/s peak at 300 s on cell (0, 0); 0.08 m3/s to 600 s on (1, 2)."""
    if t_mid_s <= 0.0 or t_mid_s >= 900.0:
        first = 0.0
    elif t_mid_s <= 300.0:
        first = 0.25 * t_mid_s / 300.0
    else:
        first = 0.25 * (900.0 - t_mid_s) / 600.0
    second = 0.08 if t_mid_s < 600.0 else 0.0
    return first, second


def cell_point(row: int, column: int):
    from shapely.geometry import Point

    return Point(DX_M * (column + 0.5), DX_M * (TERRAIN_M.shape[0] - row - 0.5))


def build_network(repo: Path):
    import geopandas as gpd
    from shapely.geometry import LineString

    from hydropol2d.urban_drainage import UrbanDrainageNetwork

    crs = "EPSG:32633"
    node_rows = [
        # id, type, rim, invert, shaft, initial head, cell
        ("J1", "junction", 10.0, 8.0, 1.2, np.nan, (0, 0)),
        ("J2", "junction", 9.5, 7.6, 1.0, np.nan, (1, 2)),
        ("J3", "junction", 9.0, 7.2, 1.0, np.nan, (2, 3)),
        ("OF1", "outfall", 8.8, 7.0, 1.0, np.nan, (2, 4)),
        ("OF2", "outfall", 9.3, 7.3, 1.0, np.nan, (1, 4)),
        ("J4", "junction", 7.5, 6.0, 1.0, 6.5, (3, 1)),
        ("J5", "junction", 7.8, 6.3, 0.8, np.nan, (3, 4)),
        ("OF3", "outfall", 7.6, 5.8, 1.0, np.nan, (3, 0)),
    ]
    nodes = gpd.GeoDataFrame(
        [
            {
                "node_id": node_id,
                "node_type": kind,
                "rim_elevation_m": rim,
                "invert_elevation_m": invert,
                "shaft_area_m2": shaft,
                "initial_head_m": head,
                "design_status": "valid",
                "geometry": cell_point(*cell),
            }
            for node_id, kind, rim, invert, shaft, head, cell in node_rows
        ],
        geometry="geometry",
        crs=crs,
    )
    where = {row[0]: cell_point(*row[6]) for row in node_rows}
    pipe_rows = [
        # id, upstream, downstream, length, diameter, upstream invert, downstream invert
        ("P1", "J1", "J2", 40.0, 0.45, 8.0, 7.6),
        ("P2", "J2", "J3", 40.0, 0.45, 7.6, 7.2),
        ("P3", "J3", "OF1", 30.0, 0.15, 7.2, 7.05),
        ("P4", "J2", "OF2", 35.0, 0.20, 7.6, 7.4),
        ("P5", "J5", "J4", 25.0, 0.30, 6.3, 6.05),
        ("P6", "J4", "OF3", 20.0, 0.30, 6.0, 5.85),
    ]
    pipes = gpd.GeoDataFrame(
        [
            {
                "pipe_id": pipe_id,
                "upstream_node": start,
                "downstream_node": end,
                "length_m": length,
                "diameter_m": diameter,
                "slope": (up - down) / length,
                "upstream_invert_m": up,
                "downstream_invert_m": down,
                "manning_n": 0.013,
                "design_flow_m3_s": 0.0,
                "design_status": "valid",
                "geometry": LineString([where[start], where[end]]),
            }
            for pipe_id, start, end, length, diameter, up, down in pipe_rows
        ],
        geometry="geometry",
        crs=crs,
    )
    inlet_rows = [
        # id, node, cell, opening area, discharge coefficient, reverse flow allowed
        ("I1", "J1", (0, 0), 0.10, 0.67, True),
        ("I6", "J2", (0, 0), 0.05, 0.67, True),
        ("I2", "J2", (1, 2), 0.10, 0.67, True),
        ("I3", "J3", (2, 3), 0.10, 0.60, True),
        ("I4", "J4", (3, 1), 0.08, 0.67, True),
        ("I5", "J5", (3, 4), 0.10, 0.67, False),
    ]
    inlets = gpd.GeoDataFrame(
        [
            {
                "inlet_id": number,
                "identifier": identifier,
                "node_id": node_id,
                "row": cell[0],
                "column": cell[1],
                "placement_reason": "spacing",
                "surface_elevation_m": float(TERRAIN_M[cell]),
                "crest_elevation_m": float(TERRAIN_M[cell]),
                "opening_area_m2": area,
                "discharge_coefficient": coefficient,
                "reverse_flow_allowed": reverse,
                "geometry": cell_point(*cell),
            }
            for number, (identifier, node_id, cell, area, coefficient, reverse) in enumerate(
                inlet_rows, 1
            )
        ],
        geometry="geometry",
        crs=crs,
    )
    outfalls = gpd.GeoDataFrame(
        [
            {"outfall_id": "OUT1", "node_id": "OF1", "boundary_type": "free",
             "fixed_stage_m": np.nan, "geometry": where["OF1"]},
            {"outfall_id": "OUT2", "node_id": "OF2", "boundary_type": "normal",
             "fixed_stage_m": np.nan, "geometry": where["OF2"]},
            {"outfall_id": "OUT3", "node_id": "OF3", "boundary_type": "fixed",
             "fixed_stage_m": 7.7, "geometry": where["OF3"]},
        ],
        geometry="geometry",
        crs=crs,
    )
    excluded = gpd.GeoDataFrame(columns=["street_id", "reason", "geometry"], geometry="geometry", crs=crs)
    return UrbanDrainageNetwork(
        nodes=nodes, inlets=inlets, pipes=pipes, outfalls=outfalls,
        excluded_links=excluded, inlet_catchments=None, diagnostics={}, crs=crs,
    )


GRID_SHAPE = (16, 24)
GRID_DX_M = 5.0
GRID_X0_M = 500000.0
GRID_Y0_M = 4100000.0  # top edge


def grid_terrain_m() -> np.ndarray:
    """A street (gutter on row 7) with 2% cross slope, falling 0.5% to the east."""
    rows, cols = np.indices(GRID_SHAPE)
    return 10.0 + 0.02 * GRID_DX_M * np.abs(rows - 7) + 0.005 * GRID_DX_M * (GRID_SHAPE[1] - 1 - cols)


def grid_point(row: int, column: int):
    from shapely.geometry import Point

    return Point(GRID_X0_M + GRID_DX_M * (column + 0.5), GRID_Y0_M - GRID_DX_M * (row + 0.5))


def build_grid_case(repo: Path) -> Path:
    """Rasters and network of the MATLAB grid run (VAL-UD-001, part C)."""
    import geopandas as gpd
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import LineString

    from hydropol2d.urban_drainage import UrbanDrainageNetwork, write_drainage_geopackage

    folder = HERE / "Grid"
    folder.mkdir(exist_ok=True)
    terrain = grid_terrain_m()
    transform = from_origin(GRID_X0_M, GRID_Y0_M, GRID_DX_M, GRID_DX_M)
    profile = dict(driver="GTiff", height=GRID_SHAPE[0], width=GRID_SHAPE[1], count=1,
                   crs="EPSG:32633", transform=transform)
    for name, values, dtype in (("DEM.tif", terrain, "float64"),
                                ("LULC.tif", np.ones(GRID_SHAPE), "int16"),
                                ("SOIL.tif", np.ones(GRID_SHAPE), "int16")):
        with rasterio.open(folder / name, "w", dtype=dtype, **profile) as target:
            target.write(values.astype(dtype), 1)

    crs = "EPSG:32633"
    gutter = 7
    junction_cols = (3, 8, 13, 18)
    node_rows = []
    for index, col in enumerate(junction_cols, 1):
        rim = float(terrain[gutter, col])
        node_rows.append((f"J{index}", "junction", rim, rim - 1.6 - 0.05 * index, (gutter, col)))
    node_rows.append(("OF1", "outfall", float(terrain[gutter, 23]), float(terrain[gutter, 23]) - 2.0,
                      (gutter, 23)))
    hill = (2, 12)
    node_rows.append(("J5", "junction", float(terrain[hill]), float(terrain[hill]) - 1.2, hill))
    node_rows.append(("OF2", "outfall", float(terrain[1, 12]), float(terrain[hill]) - 1.5, (1, 12)))
    nodes = gpd.GeoDataFrame(
        [{"node_id": i, "node_type": k, "rim_elevation_m": rim, "invert_elevation_m": inv,
          "shaft_area_m2": 1.0, "initial_head_m": np.nan, "design_status": "valid",
          "geometry": grid_point(*cell)} for i, k, rim, inv, cell in node_rows],
        geometry="geometry", crs=crs)
    invert = {row[0]: row[3] for row in node_rows}
    cell_of = {row[0]: row[4] for row in node_rows}
    links = [("P1", "J1", "J2", 0.30), ("P2", "J2", "J3", 0.30), ("P3", "J3", "J4", 0.40),
             ("P4", "J4", "OF1", 0.45), ("P5", "J5", "OF2", 0.30)]
    pipe_records = []
    for pipe_id, start, end, diameter in links:
        (r0, c0), (r1, c1) = cell_of[start], cell_of[end]
        length = GRID_DX_M * float(np.hypot(r1 - r0, c1 - c0))
        up = invert[start]
        down = max(invert[end] + 0.02, up - 0.012 * length)
        pipe_records.append({"pipe_id": pipe_id, "upstream_node": start, "downstream_node": end,
                             "length_m": length, "diameter_m": diameter,
                             "slope": (up - down) / length, "upstream_invert_m": up,
                             "downstream_invert_m": down, "manning_n": 0.013,
                             "design_flow_m3_s": 0.0, "design_status": "valid",
                             "geometry": LineString([grid_point(*cell_of[start]), grid_point(*cell_of[end])])})
    pipes = gpd.GeoDataFrame(pipe_records, geometry="geometry", crs=crs)
    inlet_rows = [(f"I{i}", f"J{i}", (gutter, col), 0.10, True) for i, col in enumerate(junction_cols, 1)]
    inlet_rows.append(("I5", "J5", hill, 0.02, True))
    inlets = gpd.GeoDataFrame(
        [{"inlet_id": n, "identifier": i, "node_id": node, "row": cell[0], "column": cell[1],
          "placement_reason": "spacing", "surface_elevation_m": float(terrain[cell]),
          "crest_elevation_m": float(terrain[cell]), "opening_area_m2": area,
          "discharge_coefficient": 0.67, "reverse_flow_allowed": reverse,
          "geometry": grid_point(*cell)} for n, (i, node, cell, area, reverse) in enumerate(inlet_rows, 1)],
        geometry="geometry", crs=crs)
    outfalls = gpd.GeoDataFrame(
        [{"outfall_id": "OUT1", "node_id": "OF1", "boundary_type": "free", "fixed_stage_m": np.nan,
          "geometry": grid_point(*cell_of["OF1"])},
         {"outfall_id": "OUT2", "node_id": "OF2", "boundary_type": "fixed",
          "fixed_stage_m": float(terrain[hill]) + 0.03, "geometry": grid_point(*cell_of["OF2"])}],
        geometry="geometry", crs=crs)
    excluded = gpd.GeoDataFrame(columns=["street_id", "reason", "geometry"], geometry="geometry", crs=crs)
    network = UrbanDrainageNetwork(nodes=nodes, inlets=inlets, pipes=pipes, outfalls=outfalls,
                                   excluded_links=excluded, inlet_catchments=None, diagnostics={}, crs=crs)
    return write_drainage_geopackage(network, HERE / "network_grid.gpkg")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python-repo", required=True, type=Path)
    args = parser.parse_args()
    repo = args.python_repo.resolve()
    sys.path.insert(0, str(repo / "src"))
    commit = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    branch = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()

    from hydropol2d.core.backend import configure_backend

    configure_backend("cpu")
    import jax.numpy as jnp

    from hydropol2d.config.models import UrbanDrainageHydraulicsConfig
    from hydropol2d.core.engine import _compiled_urban_drainage_advance
    from hydropol2d.modules.urban_drainage import (
        initialize_urban_drainage_state,
        node_head_m,
        prepare_urban_drainage,
    )
    from hydropol2d.urban_drainage import read_drainage_geopackage, write_drainage_geopackage

    gpkg = HERE / "network_kernel.gpkg"
    write_drainage_geopackage(build_network(repo), gpkg)
    network = read_drainage_geopackage(gpkg)
    prepared = prepare_urban_drainage(
        network,
        terrain_elevation_m=TERRAIN_M,
        valid_mask=np.ones(TERRAIN_M.shape, dtype=bool),
        cell_area_m2=DX_M * DX_M,
        hydraulics=UrbanDrainageHydraulicsConfig(),
        maximum_velocity_m_s=10.0,
    )
    advance = _compiled_urban_drainage_advance(SimpleNamespace(urban_drainage=prepared))
    volume, discharge = initialize_urban_drainage_state(prepared)
    head = node_head_m(prepared, volume)
    n_inlets = len(prepared.inlet_ids)
    n_pipes = len(prepared.pipe_ids)
    n_nodes = len(prepared.node_ids)
    surface = np.zeros(TERRAIN_M.shape)
    surface[INITIAL_POND[0]] = INITIAL_POND[1]
    state = dict(
        surface=jnp.asarray(surface),
        volume=volume,
        head=head,
        q=discharge,
        inlet_q=jnp.zeros(n_inlets),
        capture=jnp.zeros(n_inlets),
        surcharge=jnp.zeros(n_inlets),
        max_head=head,
        duration=jnp.zeros(n_nodes),
        max_q=jnp.zeros(n_pipes),
        max_fill=jnp.zeros(n_pipes),
        reversals=jnp.zeros(n_pipes),
        peak=jnp.zeros(n_inlets),
    )
    cells = sorted({(int(r), int(c)) for r, c in zip(network.inlets["row"], network.inlets["column"])})
    header = ["step", "t_end_s", "source_00_m3", "source_12_m3", "substeps", "outflow_m3", "inflow_m3"]
    header += [f"h_{n}" for n in prepared.node_ids]
    header += [f"V_{n}" for n in prepared.node_ids]
    header += [f"q_{p}" for p in prepared.pipe_ids]
    header += [f"qi_{i}" for i in prepared.inlet_ids]
    header += [f"S_r{r}c{c}" for r, c in cells]
    rows = []
    for step in range(1, STEPS + 1):
        t_mid = (step - 0.5) * DT_S
        rates = source_rates_m3_s(t_mid)
        surface = np.asarray(state["surface"]).copy()
        added = [rate * DT_S for rate in rates]
        for (r, c), volume_m3 in zip(SOURCE_CELLS, added):
            surface[r, c] += volume_m3
        result = advance(
            jnp.asarray(surface), state["volume"], state["head"], state["q"], state["inlet_q"],
            state["capture"], state["surcharge"], state["max_head"], state["duration"],
            state["max_q"], state["max_fill"], state["reversals"], state["peak"],
            jnp.asarray(DT_S),
        )
        (remaining, count, surface_out, volume, head, q, inlet_q, capture, surcharge,
         max_head, duration, max_q, max_fill, reversals, peak, outflow, inflow,
         failed, _limit, _dt) = result
        if bool(failed) or float(remaining) > 1e-12:
            raise RuntimeError(f"Python drainage advance failed at step {step}.")
        state.update(surface=surface_out, volume=volume, head=head, q=q, inlet_q=inlet_q,
                     capture=capture, surcharge=surcharge, max_head=max_head, duration=duration,
                     max_q=max_q, max_fill=max_fill, reversals=reversals, peak=peak)
        s_out = np.asarray(surface_out)
        rows.append(
            [step, step * DT_S, added[0], added[1], int(count), float(outflow), float(inflow)]
            + list(np.asarray(head)) + list(np.asarray(volume)) + list(np.asarray(q))
            + list(np.asarray(inlet_q)) + [s_out[r, c] for r, c in cells]
        )
    # The grid run's network must also pass HydroPol2D-Python's own checks.
    grid_gpkg = build_grid_case(repo)
    from hydropol2d.urban_drainage import audit_drainage_network
    from hydropol2d.config.models import UrbanDrainageGenerationConfig

    grid_network = read_drainage_geopackage(grid_gpkg)
    audit = audit_drainage_network(
        grid_network, parameters=UrbanDrainageGenerationConfig(design_rainfall_intensity_mm_h=50.0)
    )
    if not audit["acceptable"]:
        raise RuntimeError(f"Grid network fails the Python audit: {audit['blocking_issues']}")
    prepare_urban_drainage(
        grid_network, terrain_elevation_m=grid_terrain_m(),
        valid_mask=np.ones(GRID_SHAPE, dtype=bool), cell_area_m2=GRID_DX_M**2,
        hydraulics=UrbanDrainageHydraulicsConfig(), maximum_velocity_m_s=10.0,
    )

    reference = HERE / "python_reference_kernel.csv"
    with reference.open("w", encoding="utf-8") as handle:
        handle.write(f"# HydroPol2D-Python {branch} {commit}\n")
        handle.write("# make_python_reference.py: _compiled_urban_drainage_advance, 400 x 5 s steps\n")
        handle.write(",".join(header) + "\n")
        for row in rows:
            handle.write(",".join(repr(float(v)) if not isinstance(v, int) else str(v) for v in row) + "\n")
    inputs = HERE / "python_reference_inputs.csv"
    with inputs.open("w", encoding="utf-8") as handle:
        handle.write(f"# HydroPol2D-Python {branch} {commit}\n")
        handle.write("# terrain_m of the 4 x 5 kernel grid (rows from the top), dx = 10 m\n")
        for row in TERRAIN_M:
            handle.write(",".join(repr(float(v)) for v in row) + "\n")
    print(f"Wrote {gpkg.name}, {reference.name}, {inputs.name} and the grid case from {commit[:10]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
