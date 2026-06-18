# P1-SNOW-001: Cold/Warm Partition and Melt

Purpose: validate snow accumulation, rain/snow partitioning, melt, sublimation, and snowpack mass closure under prescribed forcing.

Phase 1 case-study domain: `HydroPol2D_Model/Validation/Phase1_VTilted_Catchment`. The snow routine is evaluated by v-tilted zone using representative left hillslope, channel-strip, and right-hillslope forcing sequences.

Expected behavior:
- Cold precipitation accumulates in snow storage.
- Transitional precipitation is partitioned using the implemented 4 to 7 C linear snow-fraction rule.
- Warm precipitation remains liquid.
- Melt is bounded by available snow water equivalent.
- Sublimation removes snow storage without generating runoff.

Required diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- snow storage, snow depth, melt, sublimation, rain/snow partition, and density time series

Acceptance threshold: storage residual `< 1e-6 m3` for unit-scale tests or `< 0.01%` of input volume.

## Phase 1 implementation

Generate the independent reference:

```bash
python3 HydroPol2D_Model/Validation/scripts/phase1_reference_solutions.py --case snow_degree_day --output HydroPol2D_Model/Validation/Reference_Outputs/Phase1
```

Run the actual HydroPol2D snow function in MATLAB:

```matlab
run('HydroPol2D_Model/Validation/Snow/ColdWarmPartition_Melt/run_snow_model.m')
```

Then compare:

```bash
python3 HydroPol2D_Model/Validation/Snow/ColdWarmPartition_Melt/compare_snow_model.py
```

Implementation note: `Snow_Model_Function.m` currently accepts a `T_thresh` argument, but the active precipitation partition uses hard-coded lower and upper transition temperatures of 4 C and 7 C. The Phase 1 reference mirrors the active implementation so the test is an exact formula/storage verification.
