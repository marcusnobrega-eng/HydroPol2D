# CANOPY-001: Single-cell Canopy Storage Balance

Purpose: validate canopy interception storage before coupled surface or subsurface routing.

Phase 1 case-study domain: `HydroPol2D_Model/Validation/Phase1_VTilted_Catchment`. The same canopy bucket reference is evaluated by v-tilted LAI zone, including the LAI=0 channel strip for bypass behavior.

Expected behavior:
- LAI=0 bypasses canopy interception.
- Nonzero LAI fills canopy storage to capacity.
- Throughfall starts only after canopy storage capacity is exceeded.
- Evaporation depletes canopy storage and cannot make storage negative.

Required diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- canopy storage and throughfall time-series figure

Acceptance threshold: storage residual `< 1e-6 m3` and bypass/throughfall-onset errors at numerical precision.

## Phase 1 implementation

Generate the independent reference:

```bash
python3 HydroPol2D_Model/Validation/scripts/phase1_reference_solutions.py --case canopy_bucket --output HydroPol2D_Model/Validation/Reference_Outputs/Phase1
```

Run the actual HydroPol2D interception function in MATLAB:

```matlab
run('HydroPol2D_Model/Validation/Canopy_Interception/SingleCell_Storage_Balance/run_canopy_interception_model.m')
```

Then compare:

```bash
python3 HydroPol2D_Model/Validation/Canopy_Interception/SingleCell_Storage_Balance/compare_canopy_interception.py
```

If MATLAB is unavailable, the comparison script can generate a Python mirror output for diagnostic scaffolding:

```bash
python3 HydroPol2D_Model/Validation/Canopy_Interception/SingleCell_Storage_Balance/compare_canopy_interception.py --generate-mirror-output
```
