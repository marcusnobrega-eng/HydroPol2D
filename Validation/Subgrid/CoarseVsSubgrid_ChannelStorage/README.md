# P1-SUBGRID-001: Lookup-table Subgrid Hydraulics

Purpose: validate the newer lookup-table subgrid hydraulic pathway for storage, inverse depth, shared-face hydraulic geometry/conveyance, areal source/sink conversion, and controlled coupling to local-inertial and full-momentum routing.

Expected behavior:
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

Acceptance thresholds: storage residual `< 1e-9 m3`, conveyance error `< 0.1%`, coupled mass error `< 0.1%`, and flat full-momentum equivalence NSE `> 0.999`.

Current-code result: passed. Volume RMSE is `3.28e-14 m3`, face geometry and conveyance errors are `0%`, source/sink volume residual is `2.26e-14 m3`, flat full-momentum depth RMSE is `1.96e-17 m`, incised-channel mass error is `3.16e-14%`, and the low-friction local-inertial closed-volume residual is `2.22e-16 m3`.
