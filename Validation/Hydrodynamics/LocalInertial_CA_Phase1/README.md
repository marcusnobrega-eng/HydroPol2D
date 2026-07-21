# P1 Local-Inertial and Cellular-Automata Hydrodynamics Validation

Purpose: repeat the Phase 1 hydrodynamics validation sequence for HydroPol2D local-inertial and cellular-automata routing after the full-momentum baseline.

Cases:
- tilted-plane rainfall-runoff against the analytical kinematic-wave limit;
- V-tilted conservation and left/right symmetry;
- Ritter dry-bed dam-break profile comparison;
- prescribed stage-hydrograph boundary bookkeeping.

Evidence level:
- Local inertial is evaluated as an approximation to shallow-water dynamics.
- Cellular automata is evaluated primarily for conservation, symmetry, and boundary bookkeeping; analytical wave/profile comparisons are reported as diagnostics because CA does not solve momentum.

Run:

```matlab
cd HydroPol2D_Model/Validation/Hydrodynamics/LocalInertial_CA_Phase1
run_local_inertial_ca_hydrodynamics_validation
```

Outputs are written to `Outputs/Validation`.
