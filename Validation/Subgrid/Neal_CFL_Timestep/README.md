# Neal channel CFL time step (VAL-SUBGRID-NEAL-CFL-001)

Purpose: check that a Neal (2012) channel-subgrid run gets a CFL time step.
The Neal solver returns no face depth, so the velocity and wave rules in
`refreshing_timestep.m` saw still water and kept the step at
`max_time_step`. `hp2d_neal_stable_timestep.m` now limits the step:

    dt = alfa_min * dx / max(|v| + sqrt(g h))

- `h` is the depth above the channel bed.
- The cell discharge is the mean of its two faces in each direction, channel
  plus floodplain.
- `v` divides that discharge by the wetted section: channel width x depth in
  bank, cell width x depth over the banks.

This is the rule of `estimate_neal_subgrid_dt_s` in HydroPol2D-Python.

The case is a row of ten 50 m cells with a 5 m wide, 3 m deep channel. Cells
1-5 flow in bank (2 m, 20 m3/s); cells 6-10 flow over the banks (3.5 m,
20 m3/s in the channel and 15 m3/s on the floodplain). The check computes the
reference cell by cell, independently of the helper. With `alfa_min` 0.4 and
`max_time_step` 60 s the reference is 3.111 s. Before the fix the model
returned 60 s.

## Run

```matlab
run('Validation/Subgrid/Neal_CFL_Timestep/run_neal_cfl_timestep_check.m')
```

Acceptance: the model step is below `max_time_step` and within 0.001 s (the
`time_step_increments` floor) of the reference.
