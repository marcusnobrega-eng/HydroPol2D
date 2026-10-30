"""Python reference for VAL-D8LI-001 (D8 local-inertial routing, MATLAB vs Python).

HydroPol2D-Python: /Users/mngomes/Documents/GitHub/HydroPol2D-Python-gui-polish,
branch release/studio-0.9.0rc2, commit 66cd12f (the D8 solver dates from 9b91aa1).

Run from this folder with the Python repository's environment:

    PYTHONPATH=<python repo>/src <python repo>/.venv/bin/python make_python_reference.py

The case is the 9 x 11, 20 m V-tilted catchment of the Python D8 tests
(tests/test_d8_local_inertial_vtilted.py): rivers above 0.0012 km2, width
5 A^0.2 and bank depth 0.15 A^0.2 (A in km2; shallower than the test's 0.5,
so the river leaves its banks), and a normal-depth outlet
(slope 0.004) at the bottom-centre cell. Rain is 47 mm/h for 10 min, then 0;
the run lasts 30 min. (47 mm/h, not 60: with whole-second steps 60 mm/h puts a
dry cell at exactly the 1 mm wet/dry depth, where rounding in the last bit
decides whether a link flows.)

Part 1 (kernel): the Python D8 kernels are driven step by step, as
run_local_inertial_window does, in hybrid mode (D4 floodplain + D8 river) and
full-domain mode. Land n = 0.05, channel n = 0.03 (in the prismatic tables);
the outlet uses the channel n (MATLAB's rule; see the README). Each step:
rain on every cell, [hybrid: D4 floodplain with the D8 stage], D8 links,
outlet. dt = min(2 s, link CFL with courant 0.5, time to the next record).
The step list is written so MATLAB can replay the same steps; the series
file holds the link CFL time step at each record without the 2 s cap.

Part 2 (engine): run_case on the same DEM with n = 0.04 everywhere, in both
modes, for the outlet hydrograph of the full Python engine.

Outputs (CSV, in this folder): reference_graph.csv, reference_tables_<mode>.csv,
reference_<mode>_steps.csv, reference_<mode>_links.csv,
reference_<mode>_cells.csv, reference_<mode>_series.csv,
reference_engine_<mode>.csv; and python_river_tables.nc (the prismatic river
tables of part 1, in Python's file layout).
"""

from __future__ import annotations

import csv
import shutil
import tempfile
from pathlib import Path
from types import SimpleNamespace

import jax

# Double precision before any HydroPol2D module builds an array (run_case does
# the same through configure_backend).
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
import numpy as np
import rasterio
from rasterio.transform import from_origin

from hydropol2d import preprocess_case, run_case
from hydropol2d.config.models import (
    BoundaryFace,
    CaseConfig,
    ForcingConfig,
    ForcingPoint,
    GridConfig,
    ModelFlags,
    NormalFlowBoundary,
    OutputConfig,
    RoutingConfig,
    StaticConfig,
    TerrainPreprocessingConfig,
    TimeConfig,
    VerticalHydrologyConfig,
)
from hydropol2d.routing.d8_river_local_inertial import (
    apply_d8_river_outlets_jax,
    build_d8_river_outlets,
    build_d8_river_topology,
    d4_faces_replaced_by_d8_river,
    d8_river_channel_depth_jax,
    d8_river_local_inertial_step_jax,
    d8_river_water_surface_jax,
    estimate_d8_river_timestep_jax,
    expand_d8_tables_to_full_domain,
    load_d8_river_subgrid_tables,
    prepared_d8_river_geometry,
    write_prismatic_d8_river_tables,
)
from hydropol2d.routing.local_inertial import local_inertial_step


HERE = Path(__file__).resolve().parent
NROWS, NCOLS, DX = 9, 11, 20.0
CENTER = NCOLS // 2
CELL_AREA = DX * DX
RAIN_MM_H = 47.0
RAIN_END_S = 600.0
END_S = 1800.0
RECORD_S = 120.0
MAX_DT_S = 2.0
COURANT = 0.5
MIN_DEPTH_M = 1.0e-3
OUTLET_SLOPE = 0.004
KERNEL_LAND_N = 0.05
KERNEL_CHANNEL_N = 0.03
ENGINE_N = 0.04


def dem() -> np.ndarray:
    rows, columns = np.meshgrid(np.arange(NROWS), np.arange(NCOLS), indexing="ij")
    return 0.004 * np.abs(columns - CENTER) * DX + 0.004 * (NROWS - 1 - rows) * DX


def write_raster(path: Path, values: np.ndarray) -> None:
    with rasterio.open(
        path, "w", driver="GTiff", height=NROWS, width=NCOLS, count=1, dtype="float64",
        crs="EPSG:3857", transform=from_origin(0.0, NROWS * DX, DX, DX), nodata=-9999.0,
    ) as dataset:
        dataset.write(values.astype(np.float64), 1)


def case(directory: Path, routing: str, manning_n: float, rainfall_series: bool) -> CaseConfig:
    directory.mkdir(parents=True, exist_ok=True)
    for name, values in {"dem.tif": dem(), "lulc.tif": np.ones((NROWS, NCOLS)),
                         "soil.tif": np.ones((NROWS, NCOLS))}.items():
        write_raster(directory / name, values)
    forcing = (
        ForcingConfig(rainfall_mm_h_series=(
            ForcingPoint(time_s=0.0, value=RAIN_MM_H),
            ForcingPoint(time_s=RAIN_END_S, value=0.0),
        ))
        if rainfall_series
        else ForcingConfig(rainfall_mm_h=RAIN_MM_H)
    )
    return CaseConfig(
        name=f"vtilted-{routing}",
        grid=GridConfig(nrows=NROWS, ncols=NCOLS, dx_m=DX),
        time=TimeConfig(end_s=END_S, dt_s=RECORD_S),
        flags=ModelFlags(canopy=False, infiltration=False, routing=routing,
                         structured_lookup_subgrid=True),
        forcing=forcing,
        vertical=VerticalHydrologyConfig(
            soil_depth_m=2.0, theta_initial=0.2, theta_r=0.05, theta_s=0.45,
            alpha_vg_m_inverse=2.5, n_vg=1.6, ksat_near_surface_mm_h=0.0),
        routing=RoutingConfig(
            default_manning_n=manning_n, courant=COURANT, min_dt_s=0.01, max_dt_s=MAX_DT_S,
            min_depth_m=MIN_DEPTH_M,
            structured_subgrid_table_path=Path("river_tables.nc"),
            hydraulic_table_coverage="rivers",
            normal_flow=(NormalFlowBoundary(
                faces=(BoundaryFace(row=NROWS - 1, column=CENTER, side="south"),),
                slope_m_m=OUTLET_SLOPE),)),
        terrain=TerrainPreprocessingConfig(
            minimum_stream_area_km2=0.0012, width_coefficient=5.0, width_exponent=0.2,
            depth_coefficient=0.15, depth_exponent=0.2),
        static=StaticConfig(dem_path=Path("dem.tif"), lulc_path=Path("lulc.tif"),
                            soil_path=Path("soil.tif")),
        outputs=OutputConfig(interval_s=RECORD_S, diagnostic_interval_s=RECORD_S,
                             write_geotiff=False, write_netcdf=False,
                             write_modeling_results=False),
        source_path=directory / "case.yaml",
    )


def write_csv(path: Path, header: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow([f"{value:.17g}" if isinstance(value, float) else value
                             for value in row])


def kernel_run(prepared, mode: str) -> None:
    full = mode == "full"
    terrain_products = prepared.terrain_products
    routing = prepared.routing
    config = prepared.config
    river_topology = build_d8_river_topology(
        terrain_products.stream_mask, terrain_products.d8.downstream_index, dx_m=DX, dy_m=DX)
    river_tables = load_d8_river_subgrid_tables(prepared, river_topology)
    topology = (
        build_d8_river_topology(routing.valid_mask, terrain_products.d8.downstream_index,
                                dx_m=DX, dy_m=DX)
        if full else river_topology)
    tables = (expand_d8_tables_to_full_domain(prepared, topology, river_topology, river_tables)
              if full else river_tables)
    outlets = build_d8_river_outlets(prepared)
    stream = np.asarray(terrain_products.stream_mask, dtype=bool)
    width, _ = prepared_d8_river_geometry(prepared)
    if full:
        width = np.where(stream, width, np.sqrt(CELL_AREA))
    # MATLAB's outlet rule: the channel n on river cells, land-cover n elsewhere.
    outlet_n = np.where(stream, KERNEL_CHANNEL_N, np.asarray(routing.manning_n))
    inactive_x, inactive_y = d4_faces_replaced_by_d8_river(topology)
    terrain = jnp.asarray(routing.terrain_elevation_m, dtype=jnp.float64)
    valid = np.asarray(routing.valid_mask, dtype=bool)

    # Tables as the solver uses them (full mode: the expanded link rows).
    levels_index = 60
    write_csv(HERE / f"reference_tables_{mode}.csv",
              ["link", "source_row", "source_col", "target_row", "target_col", "length_m",
               "diagonal", "face_datum_m", "area_at_level_60_m2", "conveyance_at_level_60",
               "conveyance_at_top"],
              [(k + 1, int(s // NCOLS) + 1, int(s % NCOLS) + 1, int(t // NCOLS) + 1,
                int(t % NCOLS) + 1, float(topology.edge_length_m[k]), int(topology.diagonal[k]),
                float(tables.faces.datum_m[k]), float(tables.faces.area_m2[k, levels_index]),
                float(tables.faces.conveyance[k, levels_index]),
                float(tables.faces.conveyance[k, int(tables.faces.count[k]) - 1]))
               for k, (s, t) in enumerate(zip(topology.source_index, topology.target_index))])

    volume = jnp.zeros((NROWS, NCOLS), dtype=jnp.float64)
    discharge = jnp.zeros(topology.edge_count, dtype=jnp.float64)
    qx = jnp.zeros((NROWS, NCOLS + 1), dtype=jnp.float64)
    qy = jnp.zeros((NROWS + 1, NCOLS), dtype=jnp.float64)
    time_s = 0.0
    next_record = RECORD_S
    cumulative_rain = 0.0
    cumulative_out = 0.0
    steps, link_rows, cell_rows, series_rows = [], [], [], []

    def record(t: float, dt_cfl: float) -> None:
        stage = np.asarray(d8_river_water_surface_jax(volume, terrain, CELL_AREA, tables))
        channel = np.asarray(d8_river_channel_depth_jax(volume, terrain, CELL_AREA, tables))
        vol = np.asarray(volume)
        surface = np.maximum(stage - np.asarray(terrain), 0.0)
        for k in range(topology.edge_count):
            link_rows.append((t, k + 1, float(discharge[k])))
        for r in range(NROWS):
            for c in range(NCOLS):
                cell_rows.append((t, r + 1, c + 1, float(vol[r, c]), float(stage[r, c]),
                                  float(surface[r, c]), float(channel[r, c])))
        series_rows.append((t, float(np.sum(vol)), cumulative_rain, cumulative_out,
                            float(last_outlet[0]), dt_cfl))

    last_outlet = [0.0]
    record(0.0, float(estimate_d8_river_timestep_jax(
        surface_volume_m3=volume, topology=topology, terrain_elevation_m=routing.terrain_elevation_m,
        tables=tables, cell_area_m2=CELL_AREA, courant=COURANT, max_dt_s=1.0e9)))
    step = 0
    while time_s < END_S - 1e-9:
        dt_cfl = float(estimate_d8_river_timestep_jax(
            surface_volume_m3=volume, topology=topology,
            terrain_elevation_m=routing.terrain_elevation_m, tables=tables,
            cell_area_m2=CELL_AREA, courant=COURANT, max_dt_s=MAX_DT_S))
        dt = min(MAX_DT_S, dt_cfl, next_record - time_s)
        rain_mm_h = RAIN_MM_H if time_s < RAIN_END_S - 1e-9 else 0.0
        rain_m3 = rain_mm_h / 1000.0 / 3600.0 * dt * CELL_AREA
        volume = jnp.where(valid, volume + rain_m3, 0.0)
        cumulative_rain += rain_m3 * float(np.sum(valid))
        if not full:
            water_surface = d8_river_water_surface_jax(volume, terrain, CELL_AREA, tables)
            floodplain = local_inertial_step(
                surface_volume_m3=volume, qx_previous_m3_s=qx, qy_previous_m3_s=qy,
                parameters=routing, config=config.routing, dx_m=DX, dy_m=DX, dt_s=dt,
                time_s=time_s, water_surface_m=water_surface, inactive_internal_x=inactive_x,
                inactive_internal_y=inactive_y, include_normal_flow=False)
            volume, qx, qy = floodplain.surface_volume_m3, floodplain.qx_m3_s, floodplain.qy_m3_s
        river = d8_river_local_inertial_step_jax(
            surface_volume_m3=volume, river_discharge_previous_m3_s=discharge,
            topology=topology, terrain_elevation_m=routing.terrain_elevation_m, tables=tables,
            cell_area_m2=CELL_AREA, dt_s=dt, min_depth_m=MIN_DEPTH_M)
        discharge = river.river_discharge_m3_s
        stage = d8_river_water_surface_jax(river.surface_volume_m3, terrain, CELL_AREA, tables)
        volume, outlet_q = apply_d8_river_outlets_jax(
            surface_volume_m3=river.surface_volume_m3, water_surface_m=stage,
            river_bed_elevation_m=tables.cells.datum_m, river_width_m=width,
            manning_n=outlet_n, outlets=outlets, dt_s=dt, min_depth_m=MIN_DEPTH_M)
        last_outlet[0] = float(outlet_q)
        cumulative_out += float(outlet_q) * dt
        step += 1
        steps.append((step, time_s, dt, rain_mm_h))
        time_s += dt
        if time_s >= next_record - 1e-9:
            time_s = next_record
            record(time_s, float(estimate_d8_river_timestep_jax(
                surface_volume_m3=volume, topology=topology,
                terrain_elevation_m=routing.terrain_elevation_m, tables=tables,
                cell_area_m2=CELL_AREA, courant=COURANT, max_dt_s=1.0e9)))
            next_record += RECORD_S

    write_csv(HERE / f"reference_{mode}_steps.csv", ["step", "t_start_s", "dt_s", "rain_mm_h"], steps)
    write_csv(HERE / f"reference_{mode}_links.csv", ["time_s", "link", "discharge_m3_s"], link_rows)
    write_csv(HERE / f"reference_{mode}_cells.csv",
              ["time_s", "row", "col", "volume_m3", "stage_m", "surface_depth_m",
               "channel_depth_m"], cell_rows)
    write_csv(HERE / f"reference_{mode}_series.csv",
              ["time_s", "storage_m3", "cumulative_rain_m3", "cumulative_outflow_m3",
               "outlet_m3_s", "cfl_dt_unclipped_s"], series_rows)
    print(f"kernel {mode}: {step} steps, links {topology.edge_count}, "
          f"outflow {cumulative_out:.6f} m3 of rain {cumulative_rain:.6f} m3")


def kernel_prepared(directory: Path, routing: str):
    prepared = preprocess_case(case(directory, routing, KERNEL_LAND_N, rainfall_series=False))
    # The prismatic tables with a channel n of its own (Python reads it from
    # config.subgrid.channel_manning_n when that is set).
    shim = SimpleNamespace(
        terrain_products=prepared.terrain_products, routing=prepared.routing,
        static_fields=prepared.static_fields,
        config=SimpleNamespace(grid=prepared.config.grid,
                               subgrid=SimpleNamespace(channel_manning_n=KERNEL_CHANNEL_N)))
    path = prepared.config.resolve_input_path(prepared.config.routing.structured_subgrid_table_path)
    write_prismatic_d8_river_tables(shim, path)
    return prepared, path


def main() -> None:
    work = Path(tempfile.mkdtemp(prefix="hp2d_d8_reference_"))
    try:
        prepared, table_path = kernel_prepared(work / "kernel_hybrid", "d8_local_inertial")
        shutil.copyfile(table_path, HERE / "python_river_tables.nc")
        products = prepared.terrain_products
        downstream = np.asarray(products.d8.downstream_index)
        width, depth = prepared_d8_river_geometry(prepared)
        rows = []
        for r in range(NROWS):
            for c in range(NCOLS):
                d = int(downstream[r, c])
                rows.append((r + 1, c + 1, float(dem()[r, c]),
                             d // NCOLS + 1 if d >= 0 else 0, d % NCOLS + 1 if d >= 0 else 0,
                             float(products.d8.contributing_cells[r, c]),
                             int(products.stream_mask[r, c]), float(width[r, c]),
                             float(depth[r, c])))
        write_csv(HERE / "reference_graph.csv",
                  ["row", "col", "elevation_m", "downstream_row", "downstream_col",
                   "contributing_cells", "river", "river_width_m", "river_depth_m"], rows)
        kernel_run(prepared, "hybrid")
        prepared_full, _ = kernel_prepared(work / "kernel_full", "full_d8_local_inertial")
        kernel_run(prepared_full, "full")

        for mode, routing in (("hybrid", "d8_local_inertial"), ("full", "full_d8_local_inertial")):
            engine_case = preprocess_case(case(work / f"engine_{mode}", routing, ENGINE_N,
                                               rainfall_series=True))
            write_prismatic_d8_river_tables(
                engine_case,
                engine_case.config.resolve_input_path(
                    engine_case.config.routing.structured_subgrid_table_path))
            result = run_case(engine_case)
            ds = result.dataset
            times = np.asarray(ds["time_s"].values, dtype=float)
            write_csv(HERE / f"reference_engine_{mode}.csv",
                      ["time_s", "outlet_m3_s", "cumulative_rainfall_m3",
                       "cumulative_outflow_m3", "surface_volume_m3", "residual_m3"],
                      [(float(t), float(q), float(p), float(o), float(s), float(e))
                       for t, q, p, o, s, e in zip(
                           times, ds["outlet_discharge_m3_s"].values,
                           ds["cumulative_rainfall_m3"].values,
                           ds["cumulative_boundary_outflow_m3"].values,
                           ds["total_surface_volume_m3"].values,
                           ds["water_balance_residual_m3"].values)])
            print(f"engine {mode}: peak {float(np.max(ds['outlet_discharge_m3_s'])):.6f} m3/s")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
