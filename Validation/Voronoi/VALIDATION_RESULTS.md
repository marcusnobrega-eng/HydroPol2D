# Voronoi CPU validation results

Validated on 2026-08-23 with MATLAB R2025b and Python 3.13. The MATLAB work
is on `feature/voronoi-hybrid-routing`, based on HydroPol2D commit `a50155e`.

## Passed gates

- HydroPolMesh: 10/10 tests passed (mesh QA, 20–1,000 m widths, urban masks,
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

The coupled V-tilted case uses two fractional HRUs per polygon, internal
Penman-Monteith ET, canopy storage, layered vadose storage, infiltration,
recharge, capillary rise, 10-minute benchmark groundwater scheduling,
lateral Boussinesq flow, and riverbed exchange. Both the 1,464-cell fine mesh
and 406-cell variable mesh executed exactly 36 groundwater updates in six
hours. A separate gaining/losing channel test transferred +32.4 and -226.8
m3, respectively, with residuals below 5e-10 m3.

The fine reference represents a physical 30 m river with bank-aligned Voronoi
faces; the hybrid mesh retains the same 30 m width in its channel graph.

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
- India case generation and simulation are not permitted by this branch's
  preflight workflow yet.
- The 100,000/1,000,000-polygon Sherlock benchmarks, device-only timestep
  loop, Python HydroPol2D solver port, and HydroBathyDEM/GitHub integration
  remain the next rollout stage.
