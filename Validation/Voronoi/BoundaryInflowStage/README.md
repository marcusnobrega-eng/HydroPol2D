# Mesh inflow as cell volume; stage edges under every scheme (VAL-MESH-BC-002)

Purpose: mesh boundary conditions follow HydroPol2D-Python (`voronoi.py`).

| Audit row | HydroPol2D-Python | MATLAB before | MATLAB now |
|---|---|---|---|
| L22 | Adds an inflow hydrograph to its host cells as volume before the vertical step, so it can infiltrate in the step it arrives | Applied the inflow as an edge flux after routing (a ghost state in full momentum) | `hp2d_voronoi_split_inflow` puts each inflow edge's discharge × dt into the edge's cell before the vertical step. The edge is a wall for routing and the CFL. The volume is booked as boundary inflow. |
| L23 | Applies the local-inertial face formula at stage edges under every routing scheme | The diffusive solver used Manning at stage edges | The local-inertial formula under every scheme |
| L23 (12f494b) | A stage can be a depth above the edge cell's bed | — | `surface_boundary.stage_reference = "depth_above_local_invert"` gives that |

The check uses a 3 x 4 mesh of 10 m cells:

- A. A 0.05 m3/s inflow runs for 10 min onto a dry, infiltrating cell. That
  cell must infiltrate in the first step, and the ledger must close with the
  inflow volume.
- B. One 1 s diffusive step with a 0.30 m stage above a dry edge. The edge
  discharge must equal the local-inertial value g·h·dt·slope.
- C. The same stage given as a depth above the local invert must give the
  same discharge.

## Run

```matlab
run('Validation/Voronoi/BoundaryInflowStage/run_boundary_inflow_stage_check.m')
```

| Case | Old code | New code |
|---|---|---|
| A. First-step infiltration | none | 0.055 m3 |
| B. Stage edge discharge | 0.776 m2/s (Manning) | 0.088 m2/s (local inertial) |
| C. Invert-referenced stage | ignored | same discharge |
