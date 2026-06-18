# V-Tilted Evapotranspiration Validation Suite

Purpose: validate HydroPol2D potential evapotranspiration and actual
evapotranspiration extraction on the Phase 1 V-tilted catchment context.

The suite uses a small V-tilted grid with left hillslope, channel strip, and
right hillslope cells. It validates the FAO-style potential ET formula and the
actual ET storage limiter without making a flux-tower or field-calibration
claim.

## Cases

- `P1-ET-VT-000`: potential/reference ET formula on the V-tilted elevation and albedo grid.
- `P1-ET-VT-001`: open-water evaporation capped by ponded surface storage.
- `P1-ET-VT-002`: soil ET demand is fully met when root-zone storage is abundant.
- `P1-ET-VT-003`: soil ET is storage-limited when accessible root-zone water is exhausted.
- `P1-ET-VT-004`: shallow-rooted vegetation can extract only the rooted fraction of the near-surface layer.
- `P1-ET-VT-005`: internal ETP mode uses `Hydro_States.ETP` for soil cells and `Hydro_States.Ep` for ponded cells.
- `P1-ET-VT-006`: dry, impervious, and invalid cells do not extract soil ET.

## Run

```matlab
run('HydroPol2D_Model/Validation/Evapotranspiration/Reference_ET_Extraction/run_vtilted_et_validation_suite.m')
```

## Outputs

- `Outputs/Validation/VTilted_ET_Diagnostics.csv`
- `Outputs/Validation/VTilted_ET_Pass_Fail.csv`
- `Outputs/Validation/CellResults/<case_id>_cells.csv`
- `Outputs/Validation/Figures/VTilted_ET_Overview.{png,pdf}`
- `Outputs/Validation/Figures/P1_ET_VT_*_et_maps.{png,pdf}`

Acceptance threshold: storage residual `< 1e-6 m3`; ET and storage errors at
machine precision for these deterministic Phase 1 cases.
