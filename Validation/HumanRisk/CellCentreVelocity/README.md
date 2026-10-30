# Human risk uses the cell-centre speed of the step (VAL-HR-VEL-001)

Purpose: the speed in the human-risk index is the one HydroPol2D-Python uses
(audit M14).

Python evaluates both human-risk modes with the cell-centre speed, every
step (`core/engine.py`, `_cell_velocity_m_s`). That speed is the mean of the
net discharges through each pair of opposite faces, over the face width and
the cell depth, in cells deeper than 1 mm.

MATLAB used a different velocity for each mode, and both were refreshed only
when the time step was recomputed:

- mode 3: the largest face speed;
- mode 1: the cell velocity of `refreshing_timestep`, which also missed the
  west face (L4).

`hp2d_cell_centre_velocity` computes the cell-centre speed from the step's
face outflows. `Human_Instability_Module.m` now uses it for both modes
whenever the run has face discharges. A prescribed velocity field is still
accepted when there are none, so the equation check VAL-HR-001 still runs.

The case is a 1 x 3 strip of 10 m cells, 0.2 m deep. 0.2 m3/s enters the
middle cell through its west face and 0.4 m3/s leaves through its east face.
Its cell-centre speed is 0.15 m/s; the largest face speed is 0.2 m/s.

## Run

```matlab
run('Validation/HumanRisk/CellCentreVelocity/run_risk_cell_velocity_check.m')
```

Acceptance: the speeds are exact and the module's risk index matches them
(1e-12). VAL-HR-001A and VAL-HR-001B still pass.
