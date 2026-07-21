# Analytical Hillslope Groundwater Benchmarks

This validation suite checks lateral groundwater dynamics against analytical hillslope solutions.

Cases:

- `VAL-GW-HILL-STEADY-001`: steady Dupuit hillslope with uniform recharge, no-flow divide, and fixed-head drain.
- `VAL-GW-HILL-TRANSIENT-001`: small-amplitude transient linearized Boussinesq decay mode with no-flow divide and fixed-head drain.

Run from MATLAB:

```matlab
run_hillslope_analytical_validation
```

Outputs are written to `Outputs/Validation`.
