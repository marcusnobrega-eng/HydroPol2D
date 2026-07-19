# Canonical V-Tilted Routing Regression Benchmark

`run_vtilted_full_config_audit.m` executes the full HydroPol2D preprocessing
and routing workflow on the 20 m V-tilted catchment. It verifies the supplied
configuration before routing: two outlet cells, Manning values of `0.015` and
`0.150`, no hydrologic losses, and spatially invariant rainfall of `10.8 mm/h`
for 90 min. The simulation continues to 240 min.

The digitized hydrograph in `compare_flows.m` is associated with the supplied
local-inertial configuration. It is therefore used as a configuration-regression
benchmark, not as an analytical reference for all routing equations. The scorer
reports hydrograph shape, mass balance, outlet-volume difference, and final
surface storage separately. A mode satisfies the screening check only when all
three quantitative criteria pass.

The retired lookup-subgrid scripts remain in this folder for development
history. They are not active validation cases.

Run from MATLAB:

```matlab
cd Validation/Hydrodynamics/VTilted_Benchmark_Subgrid
run_vtilted_full_config_audit
```
