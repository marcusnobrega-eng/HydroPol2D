# Rain maps on a quadtree mesh (VAL-MESH-RAINMAP-001)

Purpose: the prepared quadtree runner reads rain maps the way HydroPol2D-Python
does (audit M16).

Python samples each rain map on every HRU piece (one mesh-cell × raster-cell
overlap) and stops the rain one interval after the last map (`voronoi.py`).
`HydroPol2D_Run_Prepared_Quadtree` took the nearest map value at each mesh
cell's centre and held the last map to the end of the run.
`hp2d_quadtree_rainfall_rate` now:

- samples the map cell containing each land-use raster cell's centre
  (by rounding; `interp2` 'nearest' zeroed edge cells that sat a rounding
  error outside the map);
- averages over each mesh cell with the overlap weights;
- gives 0 after the last map.

`hp2d_overlap_raster_centres` gives the raster centres in the overlap's order.
Rain that differs between the HRUs of one mesh cell is not represented: the
mesh source is one rate per cell.

The same commit also changes how `HydroPol2D_Run_Prepared_Quadtree` reads the
case:

- M23: the routing scheme (`flags.routing`) and `flags.critical_flow` come
  from the case instead of being fixed at local inertial and on.
- L29: the wet/dry depth comes from `parameters.min_depth_m` when given
  (default 1 mm, as in Python).
- M18: the soil ET crop coefficient comes from the land-cover table's
  `crop_coefficient` when it has one (default 1).
- A land-cover row without a field keeps the fallback.

The check uses a 4 x 6 raster under two mesh cells, one raster cell outside
the domain, and two 30 min rain maps that vary by column.

## Run

```matlab
run('Validation/Voronoi/QuadtreeRainMaps/run_quadtree_rain_maps_check.m')
```

Acceptance: raster order right, both maps averaged exactly, rate 0 after
the last map.

A copy of a Builder package (Recife, 82,759 cells) was also run with two
uniform maps of 20 and 40 mm/h over 75 min. The new runner rains exactly
0.75 of the old volume: it stops at 60 min, where the old one held 40 mm/h
to the end.
