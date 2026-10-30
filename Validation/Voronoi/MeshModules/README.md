# Reservoirs, water quality and human risk on meshes (VAL-MESH-MODULES-001)

Purpose: the grid's reservoir controls, mass-based wash-off and human risk
index run on MATLAB mesh cells, as HydroPol2D-Python runs them
(`mesh_modules/reservoirs.py`, `quality.py`, `risk.py`).

A prepared mesh case switches each module on with a field of its `config`:

- `config.reservoirs`: the grid's `Reservoir_Data` table (`k1, h1, k2, k3, h2,
  k4` and the map coordinates `x_us, y_us, x_ds1, y_ds1, x_ds2, y_ds2`). Each
  point goes to the mesh cell that holds most of the raster cell under it
  (this needs the case's overlap file). Mesh cell ids can be given instead as
  `source_cell, ds1_cell, ds2_cell`. A destination of (0, 0), or cell 0,
  sends the release out of the model; NaN switches the outlet off. The release
  uses the grid's `hp2d_reservoir_releases`. Reservoirs are refused on meshes
  with river channels or sub-grid tables, as in Python.
- `config.water_quality`: `C_3` (kg/h per (m3/s)^C_4 per kg, MATLAB's unit),
  `C_4`, `min_Bt` (g/m2), and `initial_mass_kg` or `C_1, C_2, ADD`. Each edge a
  cell drains through gets W = C Q^C_4 B. Sub-steps sum to the step and a
  cell's outflow never exceeds its mass. The rating-curve wash-off is
  grid-only.
- `config.human_risk`: the grid's mode-1 parameters (`mu, Cd, weight_person,
  width1_person, width2_person, ro_water, gravity`). The cell speed comes from
  the polygon's edge discharges (Perot 2000). Mode 3 is grid-only, in both
  codes.

The modules run after surface routing, in the grid's order: reservoirs, then
wash-off, then risk. Reservoir water that leaves the model enters the mesh
mass ledger. `results.modules` holds the histories.

The check:

1. One module step on fixed volumes and edge discharges is compared with the
   HydroPol2D-Python kernels on the same numbers. MATLAB writes
   `kernel_cells.csv`, `kernel_edges.csv` and `kernel_outlets.csv`.
   `make_python_reference.py` turns them into `python_reference_*.csv` using
   HydroPol2D-Python at 66cd12f. Python reads `C_3` per second (audit
   section 4), so it gets `C_3/3600`.
2. A one-hour mesh run with rain, a critical-flow outlet and all three
   modules. The water ledger must close with the reservoir outflow, the
   pollutant mass must close with its export, and the risk stays in [0, 1].

## Run

```matlab
run('Validation/Voronoi/MeshModules/run_mesh_modules_check.m')
```

Acceptance:

- Kernel errors against Python are below 1e-9.
- The step and run balances close to 1e-9.
- The reservoir outflow and the pollutant export are positive.

Result: the largest error against Python is 7.1e-15. The run's water residual
is 4e-16 and its pollutant residual 4e-16. On the old code the mesh run ignored
the three fields: no reservoir released and no pollutant moved.

To remake the reference, run the check once so it writes the inputs, then:

```bash
PYTHONPATH=<HydroPol2D-Python>/src <venv>/bin/python Validation/Voronoi/MeshModules/make_python_reference.py
```
