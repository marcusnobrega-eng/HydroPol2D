# Legacy groundwater settings come from the case (VAL-GW-SETTINGS-001)

Purpose: the legacy groundwater scheme takes from the case the settings that
HydroPol2D-Python takes from its configuration (`core/engine.py`,
`modules/groundwater.py`).

| Row | Setting | MATLAB before | MATLAB now |
|---|---|---|---|
| L17 | Boussinesq Courant number | fixed at 0.25 | `groundwater_courant`, default 0.25, in both the stability estimate and the adaptive solver |
| L15 | Adaptive solver floor | min(dt, 1 s) | `groundwater_min_dt_s`, default 60 s |
| L15 | Stable step | no floor | at least 1 s, as in Python's `estimate_boussinesq_stable_dt_s` |
| L16 | Capillary-rise field-capacity head | fixed at 3.3 m | `capillary_field_capacity_head_m`, default 3.3 m |
| L16 | Capillary-rise Mualem l | fixed at 0.5 | the soil table's `l_vg` |

With the defaults nothing changes.

The check has three parts:

- A. `estimate_groundwater_stable_timestep` scales with the Courant number
  and is never below 1 s.
- B. `simulate_layered_capillary_rise` with l = −1 and a 1.0 m head matches
  `hp2d_vg_mualem_capillary_transfer` for one cell. That gives 0.223 mm,
  where the old fixed rule gave 0.128 mm.
- C. `Groundwater_Module` reads the settings.

## Run

```matlab
run('Validation/Groundwater/LegacySettingsFromCase/run_legacy_groundwater_settings_check.m')
```

The asynchronous groundwater cases (with capillary rise and the recharge
mound), the hillslope checks and the V-tilted legacy full run still pass.
