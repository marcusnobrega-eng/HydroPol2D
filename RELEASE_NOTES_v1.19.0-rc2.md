# HydroPol2D MATLAB v1.19.0-rc2

## Purpose

This candidate brings the MATLAB model up to the HydroPol2D-Python `0.9.0-rc2` physics on the
regular grid and runs the shared optional modules on meshes. It stays coordinated with the
HydroBathyDEM `0.4.0rc1` mesh bundle.

## Main changes

- Urban drainage, detention ponds, green infrastructure and D8 local-inertial routing run on the
  regular grid with HydroPol2D-Python's equations. They are set up in the bypass script
  (`InputData_Bypass`); the case workbook has no inputs for them. A mesh run refuses drainage,
  ponds and green infrastructure. The old raster D8 (`flag_D8`) is retired.
- Reservoirs, water quality, human risk (mode 1) and snow run on Voronoi and quadtree meshes. A
  mesh run that asks for a module it would skip stops at the input check.
- Forcing times start an interval, as in HydroPol2D-Python. A rain, inflow or ET-map table written
  with interval ends must move each time one step earlier.
- Reservoir releases conserve mass, and a (0, 0) destination sends the water out of the model.
- Human risk mode 3: the lift lever arm uses the tangent of the slope angle, and a dry cell is
  stable.
- Snow: rain and snow are split before the canopy, and only rain is intercepted.
- Groundwater and ET: a gaining river takes only groundwater stored before the step's recharge;
  reference ET uses the climate stations' ground heat flux; groundwater time stays in double
  precision with single-precision heads; both soil schemes use one lateral update rule.
- The legacy water balance runs every accepted step again. The run summary reports the applied
  rain and the maximum depth of the whole run (before, the maximum came from the last chunk of
  saved maps only).
- Layered subsurface (opt-in, `flag_layered_subsurface = 1`): three soil layers over regolith
  slots down to the depth to bedrock, Feddes root-water stress in both soil schemes, an implicit
  water-table flux (`flag_implicit_water_table_flux`, default 1) and capillary refill up to field
  capacity, each with parity checks against HydroPol2D-Python. `DTB_path` is optional in layered
  mode.
- Deep-aquifer starting state from a starting baseflow or a mean recharge.
- Gauges record depth, water level and discharge across a line of cells.
- Nine optional soil-property maps override the soil-class values cell by cell.
- One case workbook, `HydroPol2D_Inputs.xlsx`, with settings, class-table and series tabs. Older
  `General_Data.xlsx` cases still run.
- Fix: the transpiration and evaporation map lists of `General_Data.xlsx` read their own times and
  paths. Runs with `flag_input_ETP_map = 1` made before 2026-10-20 used wrong ET inputs.
- The Stanford example uses a 30 m DEM from USGS 3DEP and ESA WorldCover 2020 land use. The README
  lists the examples you can run, their expected results and the data sources.

## Validation

- Each new module has a check against the HydroPol2D-Python kernel: VAL-UD-001 (urban drainage),
  VAL-POND-001 (detention ponds), VAL-GI-001 (green infrastructure) and the D8 local-inertial
  check.
- The fast release profile passes: 42 of 42 checks.
- Stanford 30 m, 100 mm/h for 1 h, 24 h simulated: 100.00 mm of rain, runoff coefficient 0.47,
  outlet peaks of 402.84 m³/s at 75 min and 473.69 m³/s at 3.5 h, maximum depth 5.87 m in the
  Searsville Lake basin, volume error -6.45e-12 of the inflow.

This is not a stable release. Drainage, ponds and green infrastructure on meshes, the
stage-hydrograph and face-flux inflow options, and open-water evaporation with ET maps are planned
for the next update. Current-candidate CUDA validation, Apple Metal validation and the Windows
packaged-Studio gate remain open, and `v1.17.0` stays the stable numerical reference.
