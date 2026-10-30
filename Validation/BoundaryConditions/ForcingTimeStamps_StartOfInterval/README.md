# Forcing time stamps (VAL-FORCING-STAMP-001)

Purpose: check how `update_spatial_BC.m` reads the times of a rate series.
Value i is the rate from t_i to t_i+1. The last value lasts one more
interval (the spacing of the last two times), and then the rate is 0. A
series with one time is a constant rate. Before the first time the rate is 0.
This is the convention of HydroPol2D-Python (`series_end_s` in
`forcing/gridded.py`).

The rule applies to lumped rain, rain gauges, rainfall maps, evaporation and
transpiration maps and inflow hydrographs. The stage hydrograph and the ETP
station records already held each value from its own time. Design storms
(alternated blocks, Huff) are restamped in `input_data_script.m` so each block
keeps its place.

The case has rain and inflow tables at 0, 10 and 20 min, read over model steps
of 7 min from 0 to 42 min. The check integrates the same rule exactly.

Before the fix MATLAB read value i as the rate for (t_i-1, t_i]. Row 1 of a
table that starts at 0 covered (-10, 0] and was never used: the event lost
1 of its 3.5 mm. Every rate also came one interval early.

## Run

```matlab
run('Validation/BoundaryConditions/ForcingTimeStamps_StartOfInterval/run_forcing_time_stamp_check.m')
```

Acceptance: every step's rain and inflow depth within 1e-12 mm of the
reference, and the total rain equal to the table's 3.5 mm.
