# P1 Full-Momentum Hydrodynamics Validation

Purpose: validate the HydroPol2D full-momentum hydrodynamic solver before auditing local inertial, diffusive, kinematic, and CA/D8 alternatives.

Cases:
- `P1-HYDRO-FM-PLANE-001`: tilted-plane rainfall-runoff compared with the analytical kinematic-wave limit.
- `P1-HYDRO-FM-VTILT-001`: V-tilted catchment conservation and left/right symmetry check.
- `P1-HYDRO-FM-RITTER-001`: frictionless Ritter dam-break compared with the dry-bed analytical solution.
- `P1-BC-STAGE-FM-001`: prescribed stage-hydrograph boundary enforcement and exact volume bookkeeping.

Run:

```matlab
cd HydroPol2D_Model/Validation/Hydrodynamics/FullMomentum_Phase1
run_full_momentum_hydrodynamics_validation
```

Outputs are written to `Outputs/Validation`.
