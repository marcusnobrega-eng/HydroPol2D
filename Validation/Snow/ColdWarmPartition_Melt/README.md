# P1-SNOW-001: Configurable Snow Accumulation and Melt

This Phase 1 suite verifies the active snow pathway in HydroPol2D. Snow
parameters are defined per LULC class in `Input_Data_Sheets/LULC_parameters.xlsx`.
The test suite evaluates equation-level behavior and controlled snow-and-runoff
dynamics; it is not a field-calibration study.

## Test sequence

| Case | Test | Evidence |
|---|---|---|
| `P1-SNOW-001A` | Configurable cold, mixed, and warm partition | Exact prescribed partition and snow-storage balance |
| `P1-SNOW-001B` | Per-LULC parameter maps and initial snow states | No-raster, SWE-only, depth-only, both-raster, and inconsistent-raster paths |
| `P1-SNOW-001C` | Timestep refinement | 5, 15, and 60 min integrations of the same degree-day solution |
| `P1-SNOW-001D` | Controlled V-tilted snow/runoff coupling | Closed liquid-water and snow-storage ledger on a two-LULC grid |
| `P1-SNOW-001E` | Input-surface equivalence | Identical parameter arrays from Excel and bypass configuration |
| `P1-SNOW-001F` | Normal HydroPol2D V-tilted event | Regular preprocessing, meteorological forcing, snow routing, and event mass ledger |

## Run

Run the complete suite in MATLAB:

```matlab
run('Validation/Snow/ColdWarmPartition_Melt/run_snow_phase1_validation.m')
```

The driver writes `Metric_Summary.csv`, `Mass_Balance.csv`, and `Pass_Fail.csv`
under `Outputs/Validation/`. The full-model event writes its configuration,
timeseries, figures, and ledger below `FullModelRuns/P1-SNOW-001F/`.

## Inputs and interpretation

- Initial SWE and snow-depth rasters are optional. Where neither has a value,
  the cell starts snow-free with the LULC-specific initial density.
- Where only one initial-state raster is supplied, HydroPol2D derives the other
  variable from the LULC-specific density.
- Positive SWE and depth pairs must imply a physically valid density. Negative
  values and incompatible pairs stop preprocessing.
- Unmapped LULC or snow codes receive the area-weighted mean of mapped-class
  properties. `Input_Class_Code_Audit.csv` records this fallback.
- Snow routing needs air temperature, minimum temperature, and wind. The active
  implementation obtains those variables through the internal meteorological
  forcing pathway.

Acceptance requires finite states, nonnegative SWE/depth, bounded density, and
mass residuals below the stated test thresholds.
