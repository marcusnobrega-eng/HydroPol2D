# Canonical V-Tilted Hydrograph and Subgrid Validation

This validation folder adds two V-tilted cases:

- `P1-HYDRO-LI-VTILT-BENCH-001`: a 20 m local-inertial run on the canonical Phase 1 V-tilted DEM compared with the benchmark hydrograph from `compare_flows.m`.
- `P1-SUBGRID-VTILT-001`: a fair local-inertial subgrid comparison using a 20 m ordinary-grid reference and a 60 m lookup-subgrid run.

The benchmark event uses `Reference/Rainfall_Intensity_Data.csv`, extracted from the supplied `Rainfall_Intensity_Data (5).xlsx` workbook. The series applies `10.8 mm/h` at `0, 15, 30, 45, 60, 75, and 90 min`, and the module-driver run uses the `180 min` bypass simulation window. Infiltration, ET, groundwater, and other losses are disabled by construction.

The current module-driver configuration also applies the bypass roughness pattern: side slopes use Manning `n = 0.015`, and the central strip uses `n = 0.15`. This is closer to the bypass files than the earlier uniform-roughness diagnostic, but it is still a lightweight module-driver run rather than a full HydroPol2D bypass-wrapper execution.

For the subgrid comparison, the 20 m DEM is cropped at the upstream end from 50 to 48 rows so the 3-by-3 aggregation to 60 m is exact while preserving the downstream outlet. The 60 m maximum water surface is projected back to the 20 m DEM with:

```text
h_projected_20m = max(eta_coarse_max - DEM_20m, 0)
```

Run from MATLAB:

```matlab
cd HydroPol2D_Model/Validation/Hydrodynamics/VTilted_Benchmark_Subgrid
run_vtilted_hydrograph_subgrid_validation
```
