# Mesh module refusal (VAL-MESH-REFUSE-001)

Purpose: a MATLAB mesh run must not silently leave out a module the user
switched on.

A mesh run (`flag_voronoi = 1` with a prepared `VoronoiCase`, or a Case
Builder quadtree package run by `HydroPol2D_Run_Prepared_Quadtree`) reads its
physics from the prepared case. The grid inputs are read only for the mesh
switch. Before this check, snow, reservoir controls, water quality, human
risk, dam break and the deep aquifer switched on in the grid inputs (or in
the Builder's `hydropol_config.json`) were dropped from the run without a
word. `hp2d_refuse_mesh_modules` now stops the run and names each module.
A module is let through only when the prepared case carries its data
(`config.reservoirs`, `config.water_quality`, `config.human_risk`,
`config.snow`). Human risk mode 3 is always refused on meshes. The deep
aquifer never runs on meshes, in MATLAB or in HydroPol2D-Python, which also
refuses it. Detention ponds, urban drainage and green infrastructure are in
the list for when the grid gets them.

The check covers the helper (each switch alone; all off; a carried and an
empty module), the switches `HydroPol2D_Voronoi_Selection` passes on in
bypass mode, and the quadtree runner refusing a Builder deep-aquifer switch
before it reads any raster.

## Run

```matlab
run('Validation/Voronoi/ModuleRefusal/run_mesh_module_refusal_check.m')
```

Acceptance: all 9 switches refused alone, the all-off case runs, a carried
module runs, the selection and the quadtree runner refuse.
On the old code no switch is refused (0 of 9).
