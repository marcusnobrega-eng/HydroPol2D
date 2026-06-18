# RAIN-MAP-001: Tiny Raster Rainfall Totals

Purpose: validate spatial rainfall raster ingestion and temporal accumulation with exact known cellwise totals.

Expected behavior:
- Each raster cell receives the prescribed rainfall depth.
- Temporal accumulation equals the exact sum over all rainfall rasters.
- Total rainfall volume equals cell area times accumulated depth.
- Runoff differences reflect spatial rainfall differences, not parser artifacts.

Implemented diagnostics:
- `Outputs/Validation/Mass_Balance.csv`
- `Outputs/Validation/Metric_Summary.csv`
- `Outputs/Validation/Pass_Fail.csv`
- rainfall-total map and cellwise error table

Current result: `pass`. The three-raster GeoTIFF stack has total rainfall volume `75.333333 m3`, cellwise accumulated-depth error `0 m`, step-volume error `0%`, and total-volume residual `0%`.

Acceptance threshold: cellwise rainfall parser error `< 1e-9 m` equivalent, step-volume error `< 0.1%`, and total-volume error `< 0.1%`.
