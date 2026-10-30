# D8 local-inertial routing vs HydroPol2D-Python (VAL-D8LI-001, VAL-D8LI-002)

Purpose: check the MATLAB port of HydroPol2D-Python's Wflow-style D8
local-inertial routing (`routing/d8_river_local_inertial.py`, Python commit
9b91aa1; reference made at 66cd12f). MATLAB selects it with
`flags.flag_d8_local_inertial = 1` (and `flag_inertial = 1`):

- `flag_d8_full_domain = 0`: a D8 river graph inside a D4 floodplain. Each
  step routes the D4 floodplain faces first (with the D8 stage, the faces of
  cardinal river links closed), then the D8 river links, then the outlets.
- `flag_d8_full_domain = 1`: every valid cell drains along its D8 link; there
  is no D4 step.

Each link carries `Q = (Q_prev - g A dt S) / (1 + g dt |Q_prev| A / K^2)`,
with flow area `A` and conveyance `K` from the link tables, no flow below
1 mm, a Froude clip `|Q| <= A sqrt(g h)` and a donor limiter. The step is
limited by the link CFL `alfa_min L / sqrt(g h)`.

## Case

The 9 x 11, 20 m V-tilted catchment of the Python D8 tests: rivers above
0.0012 km2 (3 cells), width `5 A^0.2`, bank depth `0.15 A^0.2` (A in km2, so
the river leaves its banks), prismatic river tables, a normal-depth outlet
(slope 0.004) at the bottom-centre cell. Rain is 47 mm/h for 10 min, then 0;
the run lasts 30 min. The river graph has 51 links (36 diagonal), the full
graph 98 (80 diagonal).

## VAL-D8LI-001 (`run_d8_local_inertial_check.m`, fast)

- The D8 receivers, contributing cells, river cells, widths and depths equal
  Python's terrain products (`reference_graph.csv`).
- The prismatic tables built in MATLAB equal the link rows Python uses
  (`reference_tables_*.csv`) and the tables MATLAB reads from Python's file
  (`python_river_tables.nc`, `hp2d_d8_read_tables.m`).
- `D8_Local_Inertial_Model` replays Python's steps (`reference_*_steps.csv`).
  Every 2 min the link discharges, cell volumes, surface and channel depths,
  outlet flow and link CFL step must be within 1e-9 of Python
  (`reference_*_links.csv`, `*_cells.csv`, `*_series.csv`); they agree to
  about 1e-14. Rain - outflow - storage closes to 1e-15.
- A third run adds two reservoir controls on a river cell
  (`hp2d_reservoir_releases`): one release leaves the model (destination
  (0, 0)) and is booked as outlet flow, the other goes to a downstream river
  cell. The water balance closes. Python refuses reservoirs with D8 routing,
  so this run has no Python reference.

Land n is 0.05 and the channel n 0.03 here, so the check also tests that the
outlet uses the channel n (see "Differences from Python").

The check fails when a term is dropped: no friction term (link discharge
error 54%), no Froude clip (3.7e-4), the D4 faces of cardinal river links
left open (8%), the land n at the outlet (15%), no floodplain step (75%).
The donor limiter does not bind in this case.

## VAL-D8LI-002 (`run_d8_full_model_check.m`, full profile, about 3 min)

`HydroPol2D_preprocessing` and `HydroPol2D_Main_While` run the same catchment
in bypass mode (`Config/input_data_bypass_script.m`), once per mode, with
n = 0.04 everywhere. The event ledger must close (to 1e-9 of the rain), the
D8 depth maps must be written, and the outlet hydrograph must follow Python's
`run_case` on the same DEM (`reference_engine_*.csv`): NSE >= 0.99, peak
within 5%, outflow volume within 2%. The two engines take different steps,
so they agree to the time-stepping error.

## Differences from Python

- Outlet roughness: Python's engine uses the land-cover n at the D8 outlets;
  MATLAB uses the channel n of the river tables on river cells (land-cover n
  elsewhere). The kernel reference passes the channel n to Python's outlet
  kernel.
- Contributing cells follow the D8 graph in topological order. Python's
  `_contributing_cell_count` processes equal elevations in row-major order,
  which misses the cells behind a flat link to the up-right neighbour. No
  flat occurs in this case.
- The link CFL step is snapped to `time_step_increments` and kept at or above
  `min_time_step`, like the other MATLAB rules; Python stops the run below
  `min_dt_s`.

## Regenerate the Python reference

```bash
cd Validation/D8LocalInertial/VTilted_PythonParity
PYTHONPATH=<HydroPol2D-Python>/src <HydroPol2D-Python>/.venv/bin/python make_python_reference.py
```

## Run

```matlab
run('Validation/D8LocalInertial/VTilted_PythonParity/run_d8_local_inertial_check.m')
run('Validation/D8LocalInertial/VTilted_PythonParity/run_d8_full_model_check.m')
```

Results: `Outputs/Validation/Pass_Fail.csv` and
`Outputs/Validation/FullModel/Pass_Fail.csv` (with the hydrographs).
