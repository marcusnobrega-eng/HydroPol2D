# BC-INFLOW-001: Rectangular-domain Inflow Hydrograph

Purpose: validate imposed inflow hydrograph boundary conditions.

Expected behavior:
- The imposed hydrograph volume enters the model domain correctly.
- The routed outlet hydrograph conserves volume after storage accounting.
- Peak timing and peak magnitude errors are reported.
- Breach-style hydrographs, if used, are treated as inputs to this boundary-condition test.

Implemented diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- imposed and routed hydrograph figures

Current result: `pass`. The imposed volume is `9600 m3`, the exact imposed-volume error is `3.41e-13%`, the outlet-storage mass residual is `6.44e-13%`, and the downstream peak response is `7.377 m3/s`.

Acceptance threshold: imposed-volume error `< 0.01%`, water-balance error `< 0.1%`, and downstream response greater than `1e-4 m3/s`.
