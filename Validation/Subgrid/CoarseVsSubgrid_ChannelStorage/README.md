# Retired Lookup-Subgrid Development Case

This folder is retained as internal development history for the lookup-table
subgrid pathway. It is not an active HydroPol2D validation case and must not
be used to support production or publication claims.

Historical diagnostics retained here:
- A synthetic rectangular incised channel has the exact analytical storage-depth and wetted-area curves.
- Shared-face lookup width, wetted perimeter, hydraulic radius, and conveyance match the rectangular-channel Manning reference for uniform roughness.
- Rainfall/infiltration-style areal depth changes update subgrid volume exactly, then invert back to representative depth.
- Local-inertial shared-face discharge matches the SFINCS-style one-step inertial equation used by the model, based on grid-average face depth and effective roughness.
- Full-momentum lookup storage reduces to the ordinary coarse-cell result for a flat table and conserves volume for an incised-channel diagnostic.
- A low-friction local-inertial diagnostic remains finite and volume conservative.
- Legacy `River_Width`/`River_Depth` plus `flag_overbanks=1` is not validated here and remains diagnostic/deprecated.

Required diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- `Figures/p1_subgrid_storage_curve.png`
- `Figures/p1_subgrid_conveyance.png`
- `Figures/p1_subgrid_full_momentum.png`

The historical unit diagnostics are preserved for traceability only. They do
not establish a valid coupled subgrid-routing implementation.
