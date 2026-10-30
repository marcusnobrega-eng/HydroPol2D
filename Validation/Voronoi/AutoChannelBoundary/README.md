# Channel nodes take their boundary cell's outlet (VAL-MESH-CHANBC-001)

Purpose: a Neal channel that reaches the mesh boundary drains through it, as
in HydroPol2D-Python (audit M21).

Python applies the surface stage, normal-flow and critical-flow edges to the
Neal channel node of each boundary cell (`voronoi.py`,
`_channel_boundary_step`). A MATLAB mesh needed an explicit
`forcing.channel_boundary`; without one the node was closed. When the case
gives no channel boundary, `hp2d_voronoi_auto_channel_boundary` now builds
one from the surface boundary. Each node takes:

- the width of its widest link, at least √(plan area);
- the mean roughness of its links, 0.035 without links;
- a length of max(√(plan area), 1 m).

A stage given as a depth above the invert is measured from the channel bed.
MATLAB keeps its own channel stage formula, which carries the previous
discharge and friction; Python's form drops both.

The check has two parts:

- A. The helper on a 1 x 3 mesh with a channel node in every cell. A
  critical edge on cell 3 and a stage edge 0.4 m above the invert on cell 1
  give node 3 critical flow and node 1 a stage 0.4 m above its bed, with
  the right width and roughness.
- B. The same mesh with 0.5 m of water in the channel and a critical outlet
  on cell 3's east edge. The channel must drain through node 3 and the
  ledger must close. Without the outlet nothing leaves.

## Run

```matlab
run('Validation/Voronoi/AutoChannelBoundary/run_auto_channel_boundary_check.m')
```

Result: 29.96 of 30 m3 drains, residual 2.6e-16. The closed run keeps it all.
