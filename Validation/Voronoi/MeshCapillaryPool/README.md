# HRUs of a mesh cell share its groundwater (VAL-MESH-CAP-001)

Purpose: capillary rise and open-water evaporation on meshes follow
HydroPol2D-Python (audit L27).

- **Capillary rise.** Python lets each HRU take the rise its deficit and
  conductivity allow. When a cell's HRUs together ask for more than the
  cell's groundwater holds, it scales them together (`voronoi.py`,
  `_capillary_rise`). MATLAB capped every HRU at the cell's available depth,
  so a dry HRU next to a wet one could take only part of a store the other
  did not need.
- **Open-water evaporation.** It now takes its volume over the cell's
  surface area, as in Python. Before, MATLAB used the cell area, which
  differs where a cell holds channel plan area.

The case is one 100 m2 cell with two HRUs of 50 %. HRU 1 is dry with a fast
soil; HRU 2 is at field capacity. The cell holds 2 mm of groundwater, less
than HRU 1 asks for. The rise must take the whole 2 mm (0.2 m3), all of it
into HRU 1.

## Run

```matlab
run('Validation/Voronoi/MeshCapillaryPool/run_mesh_capillary_pool_check.m')
```

Acceptance: 0.2 m3 taken (error below 1e-12), nothing into HRU 2. On the old
code HRU 1 got only 0.1 m3.
