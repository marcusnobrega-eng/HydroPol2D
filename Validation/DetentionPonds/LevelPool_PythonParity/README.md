# VAL-POND-001: detention ponds against HydroPol2D-Python

Checks the MATLAB port of HydroPol2D-Python's level-pool detention ponds
(`routing/detention.py`, audit row H6) on the structured grid.

Run from MATLAB (the kernel part takes seconds; each full run 1-2 min on an idle laptop, longer under load):

```matlab
run('Validation/DetentionPonds/LevelPool_PythonParity/run_detention_ponds_check.m')
```

Writes `Outputs/Validation/Pass_Fail.csv` (`case_id,status,passed,report_ready,metric,notes`),
`Kernel_Parity.csv`, `FullModel_Summary.csv` and the pond time series of each full run.

## 1. Kernel parity (VAL-POND-001K)

`make_python_reference.py` builds a case on an 8 x 10 grid of 10 m cells with two ponds and
runs the public Python kernel functions for 600 steps of 10 s, in the engine's order:
`pond_forcing_step`, `add_groundwater_exfiltration_to_ponds`, `pond_exchange_step`.

- Pond A: DTM stage-storage over a 3 x 3 dug basin (polygon footprint), unlined (20 mm/h),
  with an orifice, a spillway and a pump (start 100.33 m, stop 100.18 m).
- Pond B: power law V = 150 (h - 100.3)^1.5, lined, with an orifice, a spillway and a
  disabled orifice.
- Forcing: 150 mm/h of rain from 5 to 25 min (varying across columns), 60 mm/day of
  open-water evaporation (varying across rows), 0.2 mm of groundwater exfiltration a step
  from 35 to 40 min, and 60 m3 poured on pond A's orifice outlet at 47.5 min (reverse flow).
- After the exchange every ordinary cell loses 1/60 of its volume a step: a stand-in for
  routing. Without it a destination cell fills until it is level with the pond, and there
  the explicit orifice law (Q ~ sqrt(dH), slope ~13 m2/s at dH = 17 um) overshoots every
  step and amplifies round-off by ~1.4 a step in both codes, so two correct codes drift
  apart.

The script writes `kernel_case.json` (the inputs, read by `detention_kernel_parity.m`, so
both codes see the same case), `python_reference_geometry.csv` and
`python_reference_steps.csv`. Regenerate with:

```bash
PYTHONPATH=<HydroPol2D-Python>/src <HydroPol2D-Python>/.venv/bin/python make_python_reference.py
```

The reference was made with HydroPol2D-Python commit 3cbd8ddb (release/studio-0.9.0rc2).

Pass: geometry (cells, faces, widths, sills, initial storage and stage) within 1e-10;
every step's storage, rain, evaporation, seepage, recharge, critical inflow, cumulative
device volume and total surface volume within 1e-9 m3, stage within 1e-10 m, device
discharge within 1e-10 m3/s; the same pump switches; the loop's water balance closes to
1e-9 m3; and the case exercises the pump stop, reverse flow and spillway flow.

## 2. Full model (VAL-POND-001G, VAL-POND-001N)

`run_detention_ponds_full_model.m` writes the rasters of a 12 x 16 plane of 10 m cells
(falling 0.01 to the south) and runs `HydroPol2D_preprocessing` and `HydroPol2D_Main_While`
in bypass mode (`Config/input_data_bypass_script.m`) for 60 min: rain 60 mm/h for 20 min
(N also has station ETP and infiltration).

- Pond A: DTM, dug 0.6 m (polygon footprint), unlined (10 mm/h); orifice and spillway onto a
  ditch that runs to the outlet, and a pump to (0, 0), which leaves the model.
- Pond B: power law, lined (cell-list footprint), one orifice.
- G: groundwater on (legacy module, local update, no lateral flow), infiltration off and the
  soil at theta_r, so the table moves only by pond seepage: the seepage recharges the aquifer
  under the pond. With soil drainage (infiltration on) the legacy event ledger of this domain is
  off by ~0.05 m3 with or without ponds (-0.047 m3 without ponds, -0.059 m3 with), so G leaves
  it out to close the ledger tightly.
- N: groundwater off, station ETP and infiltration: the seepage leaves the model and is
  booked as a loss.

Pass: the event mass ledger (`hp2d_update_system_mass_ledger`, with pond storage, pond
evaporation and pond seepage out) closes to 1e-9 of the rain in G (2e-13 m3 for the same
domain without ponds) and 1e-5 in N (where the domain without ponds is already off by 3e-6
of the rain, from the legacy soil bookkeeping); each pond's own balance
closes to 1e-9 m3; the pump's volume equals the volume booked as leaving the model; all
the seepage goes where the variant says; pond cells end with no surface water and no soil
water; every pond captured water, evaporated (when ETP is on) and every device moved water. The summary
also gives the ledger residual without the pond terms, to show they matter.

The rasters of the full runs are float64. With float32 heads, `hp2d_boussinesq_adaptive`
(lateral groundwater) can loop forever, with or without ponds: the suggested step it carries
over is `single`, so a step of single(interval) ends ~1e-6 s short of a double interval, and
every later step that small rounds away when added to the single `elapsed`, while the exit
tolerance is 10 eps of the double interval. Seen here at step 60 of the groundwater run
(interval 135.127 s). Not fixed in this check (it would change existing groundwater runs);
casting the interval and steps to double in that function removes it.

`run_detention_ponds_kernel_check.m` runs only part 1 (seconds; the release runner's fast
profile). `run_detention_ponds_check.m` runs both parts (several minutes; full profile).
