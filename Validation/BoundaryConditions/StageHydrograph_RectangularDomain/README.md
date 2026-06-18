# BC-STAGE-001: Rectangular-domain Stage Hydrograph

Purpose: validate prescribed-stage boundary conditions.

Expected behavior:
- Boundary stage follows the prescribed hydrograph.
- Stage perturbation propagates through the domain consistently.
- Volume response is physically consistent with boundary stage and storage changes.

Required diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- boundary-stage and profile-propagation figures

Acceptance threshold: boundary stage at numerical tolerance and full-domain water-balance error `< 0.1%`.
