# P1-HR-001: Human Risk and Instability

Purpose: validate HydroPol2D's deterministic human-risk and instability routines against independent algebraic force-balance references.

Expected behavior:
- `flag_human_instability = 1` returns the expected drag/friction instability index.
- `flag_human_instability = 3` returns the expected stable/sliding/toppling/drowning class.
- Exact geometry and drowning threshold values are classified consistently with the implemented open/inclusive boundary convention.
- No wet or dry test cell receives an undefined class.

Required diagnostics:
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Pass_Fail.csv`
- threshold grids, boundary checks, confusion table, and figures

Current result: report-ready. Both subcases pass with zero RMSE, 100% classification agreement, zero boundary errors, and zero unexplained `NaN`.

Notes:
- The validation exposed and fixed a division-by-zero edge case in the simple instability module: wet cells with no available friction now return full instability instead of `NaN`.
- The detailed classifier geometry split is now deterministic at the exact split depth by assigning the equality case to the lower-depth branch.
