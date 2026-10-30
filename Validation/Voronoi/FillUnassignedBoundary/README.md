# Unassigned mesh boundary edges get one condition (VAL-MESH-FILLBC-001)

Purpose: a mesh case can give every unlisted boundary edge the same
condition, as HydroPol2D-Python's `normal_flow_all_unassigned_boundary_edges`
does (audit M20).

MATLAB mesh cases had to list every outlet edge. A `surface_boundary` struct
can now carry `fill_unassigned_type` (for example `"normal_flow"`) and
`fill_unassigned_value` (its slope). `hp2d_voronoi_fill_boundary` then gives
that condition to every exterior edge the boundary does not list, except the
edges of inflow cells. Without these fields unlisted edges stay walls, as
before.

Python also builds outlet edges from the grid case's raster boundary faces.
That needs the grid routing configuration, which a MATLAB prepared mesh
case does not carry, so it is not ported. The prepared quadtree runner still
puts critical flow on every perimeter edge when the Builder wrote no
boundary plan.

The check has two parts:

- A. On a 3 x 4 mesh: with the fill, all 14 exterior edges get normal flow.
  An inflow cell's other exterior edges are left out. Without the fill
  nothing changes.
- B. Rain on a plane that falls to the west. With the fill, water leaves and
  the ledger closes; without it, the walls keep all of it.

## Run

```matlab
run('Validation/Voronoi/FillUnassignedBoundary/run_fill_unassigned_boundary_check.m')
```
