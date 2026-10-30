# Lateral groundwater update interval (VAL-GW-LATERAL-DT-001)

Purpose: the lateral groundwater update runs when HydroPol2D-Python's rule
says it should, in both soil schemes.

Python runs the lateral update after max(min(target, stable, suggested),
min_dt) seconds:

- target is `groundwater_target_dt_s`, 86400 s;
- stable is the Boussinesq stability step;
- suggested is the adaptive solver's next step;
- min_dt is `groundwater_min_dt_s`, 60 s.

It also runs as soon as the head has moved by
`groundwater_max_head_change_m` (0.25 m) since the last update, and on the
last surface step (`modules/groundwater.py`, `core/layered_vertical.py`).

The audit found three MATLAB differences:

| Row | Scheme | MATLAB before this change |
|---|---|---|
| M9 | Legacy | min(stable, suggested): no one-day cap and no 60 s floor. With a long stable step, baseflow and exfiltration arrived in pulses. |
| M10 | Layered | Default interval 3600 s; the docs already said 86400 s. |
| M11 | Layered | No head-change trigger. |

`hp2d_groundwater_lateral_due` now holds the rule, and `Groundwater_Module.m`
and `Layered_Subsurface_Module.m` both call it. The legacy scheme keeps its
own tolerance-based head trigger as well (a known item, left as it was). The
layered default is now 86400 s. `groundwater_max_head_change_m` is a new
input (Flags sheet or bypass flags, default 0.25 m).

The check runs eight cases of the rule, checks the layered defaults, and
checks that both modules call the rule. The groundwater release checks
(`Validation/Groundwater/run_groundwater_release_checks.m`) and the
asynchronous groundwater and hillslope checks still pass.

## Run

```matlab
run('Validation/Groundwater/LateralUpdateInterval/run_lateral_update_interval_check.m')
```

Acceptance: all eight cases right, layered defaults 86400 s and 0.25 m, both
modules on the rule.
