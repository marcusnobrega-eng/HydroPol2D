# Boussinesq inner steps end at the requested time (VAL-GW-BOUSS-STEP-001)

Purpose: `Boussinesq_2D_explicit.m` advances exactly the step it is asked for
and reports all of its exfiltration.

The function splits a step where the Darcy-velocity Courant limit asks for
it. Its step counter was reset on every pass, so the last inner step was
never shortened. The loop ran past the requested step T, and the
exfiltration it returned covered only the last inner step (audit L14).
HydroPol2D-Python takes one predictor-corrector per requested step
(`modules/groundwater.py`). Now:

- the last inner step ends exactly at T;
- exfiltration and seepage are averaged over T;
- the mass error is summed over every inner step.

When the limit does not bind, this is Python's single step.

The case is a 3 x 8 strip of 10 m cells with a fast aquifer (K = 5e-3 m/s),
so the limit splits a two-day step into several inner steps. Recharge
pushes the table above the soil surface. The water balance over T must
close: recharge = storage change + exfiltration.

## Run

```matlab
run('Validation/Groundwater/BoussinesqInnerSteps/run_boussinesq_inner_step_check.m')
```

Acceptance: balance residual below 1e-9 of the recharge, positive
exfiltration. On the old code the residual was 17 %.
