# RES-001: Rating-curve Storage Depletion

Purpose: validate reservoir/control-structure storage and discharge behavior.

Expected behavior:
- Reservoir outflow follows the prescribed rating curve.
- Storage depletion equals integrated outflow after accounting for inflow and remaining storage.
- Downstream volume equals released volume within the water-balance tolerance.

Implemented diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- reservoir stage, storage, outflow, and downstream-volume figures

Current result: `pass`. Storage RMSE is `0.273 m3`, outflow RMSE is `4.96e-5 m3/s`, recovered rating coefficient error is `0%`, and the storage-release mass residual is `4.37e-13%`.

Acceptance threshold: storage RMSE `< 2 m3`, outflow RMSE `< 0.01 m3/s`, recovered `k` error `< 0.1%`, and mass residual `< 1e-8%`.
