# Voronoi CPU validation results

Validation began on `feature/voronoi-hybrid-routing` from HydroPol2D commit
`a50155e` and was consolidated through 2026-08-27 on
`integration/ecosystem-2026-09` with MATLAB R2025b and Python 3.13. Run
manifests record the exact component versions and input checksums used for each
result.

The structured-lookup roughness convention and explicit subgrid-flag
normalization were rerun on 2026-09-05 with MATLAB R2025b Update 3 and passed.

## Passed gates

- Legacy HydroPolMesh fixture suite: 10/10 tests passed (mesh QA, 20–1,000 m widths, urban masks,
  conservative overlap/GeoTIFF output, transitions, disconnected reaches,
  bends, bifurcations, and confluences).
- MATLAB Voronoi suite: all configuration, preflight, conservation,
  forcing/groundwater, boundary, lake-at-rest, channel-regime, hydrograph,
  transition, output, and fine-reference checks passed.
- Existing raster D4 fast release suite: 11/11 checks passed unchanged.
- MATLAB Code Analyzer: zero messages in the new Voronoi functions.

## Representative numerical results

| Check | Result |
|---|---:|
| Closed-domain relative mass residual | 3.75e-15 |
| Rainfall + groundwater relative residual | 6.50e-12 |
| Uniform-channel discharge error | 3.43e-15 |
| Uniform-channel depth error | 0 m |
| Fine resolved vs hybrid discharge error | 2.46% |
| Fine resolved vs hybrid volume error | 0.15% |
| Fine resolved vs hybrid bankfull-normalized level error | 3.55e-15 |
| Fine resolved vs hybrid inundated-area error | 0% |
| Coupled V-tilted variable-mesh cell reduction | 72.27% |
| Coupled surface-depth space-time relative L2 | 4.53% |
| Coupled soil-water space-time relative L2 | 0.26% |
| Coupled infiltration-volume error | 0.86% |
| Coupled actual-ET-volume error | 0.05% |
| Coupled recharge-volume error | 0.10% |
| Coupled capillary-volume error | 0.37% |
| Coupled river-exchange error | 0.16% |
| Coupled outlet-volume error | 1.46% |
| Coupled mass residual / rainfall | <1.5e-13 |
| Full output-slice mass residual / rainfall | 5.09e-14 |
| MATLAB-to-Python remapped-depth difference | 6.94e-18 m |
| MATLAB-to-Python remapped-volume difference | 0 m3 |
| Full output-slice GeoTIFF count | 171 |
| Full output-slice figures / videos | 9 / 20 |
| Pune standalone post-processing runtime | 73.75 s |
| Pune standalone output package | 82 files; 883.67 MB |

The coupled V-tilted case uses two fractional HRUs per polygon, internal
Penman-Monteith ET, canopy storage, layered vadose storage, infiltration,
recharge, capillary rise, 10-minute benchmark groundwater scheduling,
lateral Boussinesq flow, and riverbed exchange. Both the 1,464-cell fine mesh
and 406-cell variable mesh executed exactly 36 groundwater updates in six
hours. A separate gaining/losing channel test transferred +32.4 and -226.8
m3, respectively, with residuals below 5e-10 m3.

The fine reference represents a physical 30 m river with bank-aligned Voronoi
faces; the hybrid mesh retains the same 30 m width in its channel graph.

## Conventional output acceptance

The fully coupled V-tilted case was rerun through the canonical native archive
and the standalone MATLAB post-processor. It produced 171 GeoTIFFs, three PNG,
three PDF, three SVG, 20 MP4, five CSV, and one temporal manifest while retaining
the native `time_s` and raster `map_time_s` schedules. Python independently
opened and validated the MATLAB archive as
`hydropol2d-unstructured-output-1.0`, including 388 faces, 1,061 edges, one
gauge, 37 native output times, and seven raster output times.

The 1,032,414-cell Pune acceptance history was also reprocessed without rerunning
the solver. The standalone pass completed in 73.75 s and produced hourly depth,
velocity, WSE, hazard and instability stacks, final and temporal-maximum rasters,
PNG/PDF/SVG figures, six MP4 files, and diagnostic/water-balance tables. Its
legacy source archive contains only surface-routing states, so process-specific
infiltration, ET, groundwater, and channel products are correctly absent from
that particular package; those products were verified in the canonical coupled
V-tilted archive.

## Current solver and timestep checks

The following CPU checks were repeated on 2026-08-24 after the all-face and
channel-graph timestep updates. They verify the numerical implementation;
they do not constitute a production-scale performance qualification.

| Check | Result | Status |
|---|---:|---|
| Local-inertial 20 m V-tilted, 10.8 mm/h for 90 min | NSE 1.0000; peak error 0.0212%; outlet-volume error 0.00415%; mass error 0.00254% | pass |
| Kinematic explicit smoke, normal-flow boundary | 62 steps; adaptive timestep 15.56–300 s; closed mass ledger | pass |
| Diffusive explicit smoke, normal-flow boundary | 290 steps; adaptive timestep 1.08–300 s; closed mass ledger | pass |
| Full-momentum V-tilted CPU check | lake-at-rest error \(1.39\times10^{-17}\); rainfall mass error \(1.82\times10^{-14}\) | pass |
| Full-momentum adaptive versus fixed 1 s | final-depth difference 0.0248%; outlet-volume difference 0.1119% | pass |
| Neal graph-CFL self-contained fixture | first step 7.821 s versus 15.641 s single-link limit; minimum step 2.576 s; exact mass closure | pass |

The graph-CFL fixture has two 100 m links meeting at a node and a transition
to a resolved polygon. It therefore confirms that the active criterion sums
the incident channel and transition signal capacities rather than selecting
only the shortest individual link.

## Deliberately blocked

- GPU execution raises `HydroPol2D:VoronoiGPUNotValidated`.
- The device-only timestep loop and production GPU validation remain future
  gates. CPU Voronoi routing, the Python solver path, and HydroBathyDEM UGRID
  plus conservative-overlap exchange are implemented.
