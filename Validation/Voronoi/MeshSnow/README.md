# Snow on mesh HRUs (VAL-MESH-SNOW-001)

Purpose: snow runs on MATLAB mesh cells, as HydroPol2D-Python runs
`snow_step` on its mesh HRU pieces (`mesh_modules/snow.py`).

A prepared mesh case switches snow on with `config.snow`. It takes the grid's
`Snow_Properties` parameters (`alpha, epsilon, C_e, DDF, T_snow_all,
T_rain_all, rho_snow_init, rho_max, k_t, k_swe, k_D`), each given as a
scalar, one value per cell, per HRU type (1 x h) or per HRU (n x h).
`initial_swe_mm` and `initial_depth_mm` are optional. The meteorology comes
from the case's `meteorology` forcing: `temperature_c`,
`minimum_temperature_c`, `wind_speed_m_s`, `latitude_deg` and
`day_of_year`.

In each step `hp2d_voronoi_snow_step` splits the gross precipitation of every
HRU with `Snow_Model_Function` before the canopy:

- The rain goes on to the canopy and the soil. Both mesh vertical steps now
  take one rain rate per HRU.
- The melt leaves the pack on the ground and goes straight to the surface.
- The pack is a store and sublimation a loss in the mesh water ledger.
  Sublimation is taken from the pack's own balance.

The check: a 2 x 2 mesh of 40 m cells, with two HRUs per cell that have
different degree-day factors and LAI. The ground is flat and the run lasts
12 h: 2 mm/h for 4 h at -5 C, then +6 C. It checks four things:

1. The ledger closes every step.
2. The recorded SWE equals `Snow_Model_Function` replayed offline on each
   HRU, with the run's own steps.
3. The canopy stays empty while it snows, and the pack melts away.
4. The same run on the layered subsurface closes its ledger and gives the
   same pack.

## Run

```matlab
run('Validation/Voronoi/MeshSnow/run_mesh_snow_check.m')
```

Acceptance: step residuals below 1e-9 of the rain, replay error below 1e-9 mm,
canopy storage 0 while snowing, final SWE 0, and the same pack on the layered
subsurface.

Result:

- Step residual 2.6e-13 (layered 1.8e-13).
- Replay error 0 mm.
- Canopy storage 0 while snowing; the peak SWE of 8 mm melts to 0.
- On the old code the mesh run ignored `config.snow`, so the snow fell as rain.

`run_layered_mesh_case` (Validation/Groundwater/Layered) gives the same
numbers before and after this change.
