# Standard Validation Diagnostic Schema

Every validation case should write diagnostics under:

`Outputs/Validation/`

Do not mark a case as `validated` or `report_ready` in `Validation_Cases.csv` until the following files exist and contain no unexplained `NaN`.

## Mass_Balance.csv

Required columns:

| column | description |
|---|---|
| `case_id` | Registry case identifier. |
| `control_volume` | Domain, cell, reservoir, canopy, snowpack, aquifer, or pollutant store. |
| `time_start_s` | Start time for the balance interval. |
| `time_end_s` | End time for the balance interval. |
| `initial_storage` | Initial water, snow, or pollutant mass/storage. |
| `inflow` | Total incoming flux over the interval. |
| `outflow` | Total outgoing flux over the interval. |
| `source` | Internal production or imposed source term. |
| `sink` | Internal loss term. |
| `final_storage` | Final storage at the end of the interval. |
| `residual` | `initial_storage + inflow + source - outflow - sink - final_storage`. |
| `relative_error_pct` | Residual normalized by the relevant input or throughput volume. |
| `units` | Units for the storage/flux columns. |
| `notes` | Explanation for any nonzero residual or `NaN`. |

## Metric_Summary.csv

Required columns:

| column | description |
|---|---|
| `case_id` | Registry case identifier. |
| `module` | Component being validated. |
| `benchmark_type` | Exact, analytical, independent_reference, current-code-rerun, or diagnostic. |
| `metric` | Metric name such as RMSE, MAE, max_error, relative_L2, NSE, peak_timing_error, or storage_residual. |
| `value` | Metric value. |
| `units` | Metric units. |
| `threshold` | Acceptance threshold from the registry or case README. |
| `pass` | `true`, `false`, or `pending`. |
| `notes` | Required when `pass=false`, `pass=pending`, or `value` is `NaN`. |

## Pass_Fail.csv

Required columns:

| column | description |
|---|---|
| `case_id` | Registry case identifier. |
| `module` | Component being validated. |
| `status` | `pass`, `fail`, `pending`, or `diagnostic_only`. |
| `primary_reason` | One-sentence explanation of the status. |
| `mass_balance_pass` | `true`, `false`, or `pending`. |
| `benchmark_pass` | `true`, `false`, or `pending`. |
| `nan_check_pass` | `true`, `false`, or `pending`. |
| `report_ready` | `true` only when the case can be cited in the LaTeX report as validation evidence. |
| `reviewer` | Person or script that approved the diagnostic. |
| `date_utc` | Date of the pass/fail statement. |

## Figure and table naming

Use stable names so the LaTeX report can reference generated artifacts without renaming:

- `Figures/<case_id>_timeseries.png`
- `Figures/<case_id>_profile.png`
- `Figures/<case_id>_map.png`
- `Tables/<case_id>_metrics.csv`
- `Tables/<case_id>_mass_balance.csv`

Cases without a spatial map may omit the map figure. Cases without a hydrograph/profile should replace the profile figure with the most relevant module diagnostic.
