# Canopy evaporation potential; negative ET (VAL-ET-CANOPY-001)

Purpose: the canopy evaporates at the right potential, and a negative
Penman-Monteith value is no evaporation.

- **Canopy potential (audit M5).** HydroPol2D-Python evaporates canopy water
  at the reference ET (FAO-56), or at the prescribed evaporation +
  transpiration when ET comes from maps. MATLAB used the open-water rate
  `Ep`, which is set only in station mode. With ET maps the canopy therefore
  never emptied after its first fill. `interception_module.m` now uses the
  reference ET `ETP` in station mode and the current evaporation +
  transpiration maps in map mode. `Ep` stays the rate for ponded water.
- **Negative ET (audit L7).** Python clips negative reference ET and
  open-water evaporation at 0. MATLAB only replaced NaN, so on a cold, humid
  day a negative value added water to the canopy and to ponded cells.
  `update_spatial_BC.m` now clips `ETP` and `Ep` at 0.

The check runs these cases:

- A. Map mode with a full canopy (LAI 3, 0.6 mm), no rain, maps E = 1 and
  T = 3 mm/day, one hour. The canopy must lose 4/24 mm.
- B. Station mode with ETP 5 and Ep 7 mm/day. The canopy must lose 5/24 mm.
- C. `update_spatial_BC.m` on a station day whose Penman-Monteith values are
  negative. `mock/ETP_model.m` is a test double for the station model. ETP
  and Ep must come out clipped at 0.

## Run

```matlab
run('Validation/Evapotranspiration/CanopyPotential_NegativeET/run_canopy_potential_check.m')
```

Acceptance: canopy losses within 1e-12 mm, and no negative ETP or Ep.

On the old code:

| Case | Old result | Correct |
|---|---|---|
| Map mode loss | 0 mm | 0.1667 mm |
| Station mode loss | 0.2917 mm (from Ep) | 0.2083 mm |
| ETP and Ep | stayed at -0.4 and -0.6 mm/day | 0 |
