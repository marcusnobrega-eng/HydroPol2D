# River availability, ground heat flux and single-precision heads

Three small checks of the grid groundwater and ET against HydroPol2D-Python.
Each one writes a row to `Outputs/Validation/Pass_Fail.csv`.

| Check | Audit row | MATLAB before | MATLAB now |
|---|---|---|---|
| VAL-GW-RIVAVAIL-001 | L18 | a gaining river could take this step's recharge | it takes at most the groundwater stored before this step's recharge, as in `core/engine.py` |
| VAL-ET-GHEAT-001 | L11 | `ETP_model` forced G = 0 | G comes from the stations (nearest station); with all stations at 0, nothing changes |
| VAL-GW-SINGLE-001 | none (found by the pond port) | `hp2d_boussinesq_adaptive` looped forever when the heads were single precision | time is kept in double precision, so the interval closes |

- **River:** one cell holds 2 mm of drainable groundwater, and 1.5 mm of it is
  this step's recharge. The riverbed is fast. The river may take 0.5 mm; the
  old rule took all 2 mm.
- **Heat flux:** one station has G = 1.5 MJ/m2/day. `ETP_model` must give the
  reference ET of `Evapotranspiration` with that G, and that ET must differ
  from the G = 0 result.
- **Single precision:** a 3 × 6 aquifer with single heads, over a 135.127 s
  interval. The call must return within 60 s and match the double-precision
  call to 1e-4 m.

## Run

```matlab
run('Validation/Groundwater/RiverAvailabilityHeatFlux/run_river_heat_flux_checks.m')
```
