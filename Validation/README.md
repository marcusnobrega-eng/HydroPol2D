# HydroPol2D Validation Campaign

This folder is the working validation package for HydroPol2D. The authoritative list of cases is `Validation_Cases.csv`; the companion report is `../../HydroPol2D-docs-main/latex/validation_report.tex`.

## Campaign phases

The campaign is split into two steps:

1. **Phase 1: formula and dynamics verification.** Use analytical solutions, independent reference calculations, manufactured/simple benchmark solutions, and synthetic inverse-modeling parameter-recovery cases. See `Phase1_Formula_Dynamics_Verification.md`.
2. **Phase 2: field and applied validation.** Use field observations, laboratory observations, split-sample calibration/validation, remote sensing, high-water marks, observed hydrographs, water-quality samples, and full applied cases. This phase will be designed after Phase 1 is stable.

Synthetic inverse-modeling cases in Phase 1 are parameter-recovery tests. They are useful for confidence and identifiability, but they are not a substitute for independent field validation.

Phase 1 is tracked separately in `Phase1_Cases.csv`. Reference outputs generated from analytical or bookkeeping truths are stored under `Reference_Outputs/Phase1/`.

The shared Phase 1 spatial case-study domain is `Phase1_VTilted_Catchment`. Use it as the HydroPol2D domain for analytical/reference equation tests whenever the benchmark does not require a separate geometry.

Generate or refresh Phase 1 reference outputs with:

```bash
python3 HydroPol2D_Model/Validation/scripts/phase1_reference_solutions.py --case all --output HydroPol2D_Model/Validation/Reference_Outputs/Phase1
```

## Validation order

1. Canopy interception
2. Snow
3. Infiltration
4. Evapotranspiration
5. Groundwater and recharge
6. Hydrodynamics
7. Reservoir and control structures
8. Inflow hydrograph boundary conditions
9. Stage hydrograph boundary conditions
10. Spatial rainfall rasters
11. Water quality
12. Human risk and instability

Routing alternatives are treated as a separate audit gate: kinematic and diffusive routing must be reviewed and corrected before any validation claim is made.

## Standard outputs

Each report-ready case should produce the same diagnostic bundle:

- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- `Outputs/Validation/Figures/`
- `Outputs/Validation/Tables/`

The required columns and naming conventions are defined in `Diagnostic_Schema.md`.

A case is report-ready only when it has no unexplained `NaN`, no undocumented residual, and a clear pass/fail statement tied to the registry threshold.

## Current status

The existing Ritter, non-breaking wave, and infiltration/plane cases provide the first reusable foundations. Existing groundwater Case 4, tilted-plane hydrodynamics, and routing alternatives remain diagnostic until their implementation and diagnostic issues are corrected and current-code reruns are archived.
