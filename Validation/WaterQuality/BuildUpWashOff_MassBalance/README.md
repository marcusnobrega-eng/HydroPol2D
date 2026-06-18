# WQ-001: Build-up/Wash-off Mass Balance

Purpose: validate water-quality build-up and wash-off under a prescribed runoff flux.

Expected behavior:
- Pollutant build-up follows the prescribed accumulation law.
- Wash-off follows the prescribed runoff-driven removal law.
- Exported, remaining, and residual pollutant mass close the analytical mass balance.
- Pollutant mass and concentration never become negative.

Implemented diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- pollutant mass and concentration time-series figures

Current result: `pass`. The test uses a TSS-like buildup scale of `56.04 kg/ha`, a `400 m2` cell, `30 mm/h` runoff, and a mass-based washoff coefficient selected for an initial concentration of `120 mg/L`. The single-cell mass RMSE is `4.998e-7 kg`; the two-cell mass RMSE is `6.272e-8 kg`; all mass residuals are zero and no negative masses occur.

Acceptance threshold: mass RMSE `< 1e-6 kg`, exported load error `< 0.01%`, mass residual `< 1e-6 kg` or `< 0.01%`, and no negative pollutant mass.
