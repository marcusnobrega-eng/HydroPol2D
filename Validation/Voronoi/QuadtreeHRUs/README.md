# Quadtree cells keep their land-use and soil mix (VAL-MESH-HRU-001)

Purpose: a coarse quadtree cell that is part roof and part lawn must
infiltrate like that mix, not like its majority class.

HydroPol2D-Python gives every mesh-cell x raster-cell overlap its own HRU
(`voronoi.py`). `HydroPol2D_Run_Prepared_Quadtree` used to give each mesh
cell its majority land-use and soil class, with a 0/1 pervious fraction.
Roughness and the groundwater parameters also came from the majority class.

`hp2d_quadtree_hrus` now groups the overlaps of a cell by land-use x soil
pair. Each pair is one HRU with its share of the cell's active area. LAI,
initial moisture, depth to bedrock and the impervious share are area means
over the group. The class parameters are the same in every raster cell of a
pair, so this gives Python's result wherever the continuous rasters are
uniform within a pair. It needs far fewer HRU columns than one per overlap:
a 360 m cell on a 10 m raster has 1296 overlaps but only a few pairs.
Roughness, specific yield and groundwater conductivity are area means over
the cell, as in Python. The impervious share comes from
`ImperviousFraction.tif` when the Builder wrote one (percent when its largest
value exceeds 1, decided for the whole raster). Impervious classes count
as 1.

The check uses a 2 x 2 mesh of 40 m cells over an 8 x 8 raster of 10 m cells:

- A. The HRU shares and means equal the exact counts and means. A
  single-class cell has one HRU.
- B. A two-hour storm (40 mm/h, flat ground, infiltration and canopy) runs
  with one HRU per overlap (Python's construction), with the class HRUs, and
  with the old majority class.

## Run

```matlab
run('Validation/Voronoi/QuadtreeHRUs/run_quadtree_hru_check.m')
```

Acceptance: shares and means exact; the class HRUs within 1e-6 of the
per-overlap infiltration and within 1e-8 m of its depths; the majority run
more than 5 % off.

Result: share and mean errors 0. The class HRUs match the per-overlap run to
1.6e-16 in infiltration and 1.4e-17 m in depth. The majority class is 55 %
off.

The prepared quadtree runner was also run on a copy of a Builder package
(Recife, uniform 90 m mesh, 82,759 cells, 15 min, one synthetic rain map).
There every cell is one raster cell. The final depths equal the old runner's
exactly, and the mass residual is 3.4e-11.
