# Urban drainage against the HydroPol2D-Python kernel (VAL-UD-001)

Purpose: check the storm-drain runtime ported from HydroPol2D-Python
(`modules/urban_drainage.py`, called from `core/engine.py` after the surface
routing) and its coupling to the MATLAB grid run.

The runtime (`HydroPol2D_Functions/hp2d_urban_drainage_*.m`):

- Inlets are two-way orifices, Q = Cd A sqrt(2 g dh).
  - Capture happens when the surface head exceeds max(sewer head, crest).
  - Surcharge goes back to the surface when the sewer head exceeds
    max(surface head, crest) and the inlet allows reverse flow.
  - Capture is capped at the cell's water, surcharge at the node's water.
- Pipes are local inertial with semi-implicit Manning friction. A and R are
  taken at the mean end depth.
- Nodes store water in the shaft and in half of each incident pipe. The head
  is found from the volume by bisection, or by safeguarded Newton steps from
  the previous head.
- Outfalls are free, fixed-stage or normal-depth.
  - A normal-depth outfall pipe carries Manning normal flow.
  - A free outfall takes no flow back.
- The network sub-cycles inside each surface step at
  C sqrt(L S / (g A)) per pipe, capped at max_dt_s.

One difference from Python is on by default. A free outfall sits at
min(critical, normal) depth of its pipe's flow, the SWMM rule.
HydroPol2D-Python puts it at the node invert (code review 2026-10-26, F2),
which throttles the pipe to half its upstream depth.
`hydraulics.free_outfall_depth = "invert"` gives the Python rule.

## Files

- `make_python_reference.py`: run with the HydroPol2D-Python environment.
  It writes:
  - `network_kernel.gpkg` (with Python's `write_drainage_geopackage`);
  - `python_reference_kernel.csv`, made by the engine's own sub-cycle
    `_compiled_urban_drainage_advance`;
  - `python_reference_inputs.csv` (the 4 x 5 terrain);
  - the grid-run inputs `network_grid.gpkg` and `Grid/*.tif`. This network
    passes Python's `audit_drainage_network` and `prepare_urban_drainage`
    on the same DEM.

  The Python commit is in the first line of each CSV.
- `run_urban_drainage_check.m`: parts 1 and 2 below.
- `run_urban_drainage_grid_run.m` with `Config/input_data_bypass_script.m`:
  part 3.

## Cases

1. **VAL-UD-001, kernel against Python.**
   - Setup:
     - The 4 x 5 grid of 10 m cells holds two components: J1-J2-J3 with a
       free and a normal-depth outfall, and J5-J4 on a fixed-stage outfall
       at 7.7 m. Six inlets; I1 and I6 share a cell, and I5 takes no
       reverse flow.
     - 400 steps of 5 s. Triangular sources go to two cells, and a 30 m3 pond
       starts on I5's cell.
   - The run has all the processes:
     - capture, and surcharge back to the surface at five inlets (down to
       -0.22 m3/s);
     - 193 pipe-flow reversals;
     - 140 m3 out of the outfalls and 47 m3 in through the fixed-stage
       outfall;
     - two sub-steps per surface step.
   - MATLAB runs it with the Python free-outfall rule. It must match Python
     at every step for node heads, node volumes, pipe and inlet discharges,
     inlet-cell volumes, outfall volumes and sub-step counts.
   - Every step must conserve water.
   - Both GeoPackage readers (`readgeotable` and Database Toolbox `sqlite`)
     must give the same tables.
   - Two altered runs must miss the reference: one with reverse flow off,
     and one with the default free-outfall rule.
2. **VAL-UD-001R, free-outfall rule.**
   - Setup: four 50 m pipes (0.5 m, slope 0.004) end at a free outfall and
     carry a steady 0.12 m3/s.
   - With the default rule, the nodes above the outfall reach Manning normal
     depth.
   - The Python rule leaves the last node 49% deeper.
   - The 400-step network run with the default rule also conserves water.
3. **VAL-UD-001F, grid run in bypass mode.**
   - Setup:
     - 16 x 24 cells of 5 m. Rain is 40 mm/h for 15 min and the run lasts
       30 min, local inertial, no infiltration.
     - Four gutter inlets are piped to a free outfall at the east edge.
     - A hillside inlet sits on a fixed-stage outfall 3 cm above the ground,
       which pushes water out onto the street.
   - The event ledger must close (`running_control.flag_system_mass_ledger`):
     - water in nodes and pipes counts as storage;
     - outfall outflow is booked as outlet flow at the active cell nearest
       the outfall;
     - fixed-stage inflow is boundary inflow.
   - Capture, surcharge, both outfall flows and the seven
     `Urban_Drainage_*.csv` outputs must be present.
   - The summary also gives the residual the ledger would show if any one
     booking were dropped.

## Run

```matlab
run('Validation/UrbanDrainage/NodeLinkNetwork_PythonKernel/run_urban_drainage_check.m')
run('Validation/UrbanDrainage/NodeLinkNetwork_PythonKernel/run_urban_drainage_grid_run.m')
```

Both write rows to `Outputs/Validation/Pass_Fail.csv`.

Acceptance:

| Case | Quantity | Limit |
|---|---|---|
| VAL-UD-001 | heads, node volumes | 1e-8 m, 1e-8 m3 |
| VAL-UD-001 | pipe and inlet discharge | 1e-9 m3/s |
| VAL-UD-001 | outfall volumes per step | 1e-9 m3 |
| VAL-UD-001 | sub-step counts | equal |
| VAL-UD-001 | step mass residual | 1e-9 m3 |
| VAL-UD-001R | depth above the outfall (default rule) | within 5% of normal depth |
| VAL-UD-001R | depth above the outfall (Python rule) | more than 20% above |
| VAL-UD-001F | ledger residual | below 1e-9 of rain plus outfall inflow |

To rebuild the reference:

```sh
PY=/Users/mngomes/Documents/GitHub/HydroPol2D-Python-gui-polish
$PY/.venv/bin/python make_python_reference.py --python-repo $PY
```
