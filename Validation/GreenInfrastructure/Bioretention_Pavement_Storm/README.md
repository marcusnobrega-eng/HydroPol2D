# VAL-GI-001: green infrastructure in a grid run, against HydroPol2D-Python

Driver: `run_green_infrastructure_check.m`. It writes `Outputs/Validation/Pass_Fail.csv`, `GI_Checks.csv` and `GI_Summary.csv`, and runs in about 4.5 min, so it is in the `full` release profile.

## Case

- **Domain.** A 10 × 7 tilted V of 10 m cells (x 0–70 m, y 0–100 m).
  - The valley is column 4, with slopes of 3 % across and 1 % down to the outlet (row 10).
  - Soil Ks is 10 mm/h, LAI is 1, and the initial soil moisture is 60 mm.
  - The rasters are written as double, so the model states are double where the model allows it.
- **Forcing.** 60 mm/h of rain for 45 min, then 55 min dry. Reference ET (about 4 mm/day) comes from station records.
- **Runs.** Two runs in bypass mode:
  - `legacy`: the three-layer soil.
  - `layered`: `flag_layered_subsurface`, with groundwater on.
- **Rain garden** (bioretention, 104 m²).
  - It is a polygon over three valley cells: 3 compartments with 2 lateral connections.
  - The crest is at 0.08 m, the media Ks is 80 mm/h and the inlet is 2 m wide.
  - Its underdrain, Q = 0.003 h^0.5, goes to the outlet cell.
- **Permeable pavement** (80 m²).
  - It covers two hillslope cells: 2 compartments with 1 connection.
  - Its underdrain, Q = 0.004 (h − 0.02)^0.5, leaves the model and is booked as outlet flow.

## What passes

Both runs, `legacy` and `layered`, must pass every item below. The values are from the committed reference.

| Check | Legacy | Layered | Tolerance |
|---|---|---|---|
| GI step books its transfers (host + soil + GI + out + ET = GI rain) | 9.8e-15 m³ | 8.9e-15 m³ | 1e-9 m³ |
| GI step on the model states (single precision in this build) | 2.3e-6 m³ | 2.2e-5 m³ | 1e-4 m³ |
| GI stores balance their fluxes | 1.3e-13 m³ | 1.1e-13 m³ | 1e-9 m³ |
| Event water balance (system ledger) | 7.6e-5 m³ | 6.1e-5 m³ | 1e-5 × rain (3.1e-3 m³) |
| Capture / rain-garden overflow | 83.8 / 74.4 m³ | 90.6 / 81.7 m³ | > 0 |
| Underdrain to the outlet cell / out of the model | 4.12 / 1.44 m³ | 3.96 / 1.34 m³ | > 0 |
| Native exfiltration into the soil | 2.48 m³ | 2.48 m³ | > 0 |
| Gross LAI restored after the run | 0 | 0 | 0 |
| Geometry vs Python `prepare_green_infrastructure` (areas, widths in m², m) | 9.7e-12 | 9.7e-12 | 1e-9 |
| Kernel replay vs Python kernel (states, cumulative fluxes, host volumes) | 1.4e-12 m³ | 1.4e-12 m³ | 1e-9 m³ |

Rain is 308 m³ in both runs. With GI off, the same case closes to -2.4e-4 m³ (legacy) and -1.7e-3 m³ (layered). The event residual comes from the existing model, not from the GI; the GI step itself conserves to round-off (first row).

**Controls.** These must fail, and they do:

- The ledger without the outlet booking: 1.44 m³.
- The ledger without the GI storage: 9.7 m³.
- A replay with the overflow coefficient 10 % low: off by 0.18 to 0.20 m³.

Mutating the code makes the check fail:

- Host forcing not scaled to the host share: the ledger is off by 8.39 m³.
- The outlet booking removed: the ledger is off by 1.44 m³.
- The kernel's ET dropped: the replay is off by 0.017 m³.

## Python reference

`Reference/make_python_reference.py` rebuilds the case in HydroPol2D-Python from `gi_case.json`. The commit is in `python_reference_info.csv`; it was 3cbd8dd on release/studio-0.9.0rc2. The script:

1. Runs `prepare_green_infrastructure` and writes `python_geometry.csv`.
2. Replays `green_infrastructure_step` on the GI inputs MATLAB recorded at every accepted step (`replay_inputs_<variant>.csv`). The inputs are:
   - the host-cell volumes before the GI step;
   - the gross rain and ET rates;
   - dt;
   - the fraction of the native exfiltration the MATLAB soil accepted.
3. Applies the engine coupling: what the soil does not take returns to the gravel.

The Python run replaces only `_limited_edge_transfer` with MATLAB's cap (see below). On this storm that cap never binds: the unpatched kernel gives the same states to 4e-22 m³ (`python_limiter_as_is.csv`).

**Limiter case** (`limiter_case.json`, `python_limiter_case.csv`). This is a 4 × 5 grid with the three default designs and 8 m³ in every host cell at each step, run for 20 steps of 30 s.

- Python caps each lateral edge by its donor alone. A compartment with several outflows goes negative in 18 of the 20 steps, and the clip at the next step creates 0.035 m³.
- MATLAB caps all outflows of a donor together, and gravel inflow by the free capacity. It conserves to 4e-14 m³ and matches Python run with that cap to round-off.

To regenerate: run the check with `HYDROPOL_GI_WRITE_REPLAY=1`, then run the Python script with the Python repo's venv.
