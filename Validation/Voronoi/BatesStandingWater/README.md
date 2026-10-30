# Mesh time step counts standing water (VAL-MESH-BATES-001)

Purpose: the mesh time step follows HydroPol2D-Python's rules (audit M22).

**Bates rule on square meshes.** On axis-aligned square meshes MATLAB runs
the structured Bates rule cell by cell (`bates_cartesian_timestep`). It kept
only cells moving at 0.01 m/s or more, so a deep lake at rest set no wave
limit and the step jumped to `max_dt_s`. Python's rule (`voronoi.py`) is

    alfa_min · width / (|q/h| + √(g h))

over every wet cell, where |q/h| is the cell's fastest interior face. The
grid has counted standing water this way since 20f2ded, and MATLAB now uses
the same rule on meshes.

**Unstructured rule with hydraulic tables.** The interior faces now take
their wave speed from the face's hydraulic depth A/P at the higher water
surface of the two cells, where either cell holds water, as in Python.
Before, MATLAB used the volume spread over the cell.

The check holds a 1 m lake at rest for 60 s on a 4 x 4 mesh of 10 m cells.
Every step must be 0.5 · 10 / √9.81 = 1.596 s, floored to the 1 ms
increment, and the lake must stay at rest.

## Run

```matlab
run('Validation/Voronoi/BatesStandingWater/run_mesh_bates_standing_water_check.m')
```

Acceptance: step error below 1e-12 s, depth error below 1e-12 m.

Result: the steps are 1.596 s and the depth error is 0. On the old code
every step was 30 s.

`run_voronoi_subgrid_validation` still passes with the same numbers. The
Recife Builder package, which uses hydraulic tables, also runs.
