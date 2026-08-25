# Voronoi finite-volume validation

This directory validates the separate Voronoi finite-volume runner. It does
not alter or replace the raster D4 validation cases.

The suite currently covers:

- configuration rejection for retired raster D8 and conflicting solvers;
- closed-domain conservation and lake at rest;
- conservative raster forcing and lateral Boussinesq groundwater;
- exact hybrid channel/floodplain storage and runoff transfer;
- steady and unsteady channel flow, wetting/drying, backwater/reversal;
- both resolved/subgrid transition orientations and flow reversal;
- a fine resolved 30 m channel against a coarse Neal representation;
- prescribed surface/channel boundaries and native NetCDF output.
- fractional polygon HRUs with canopy interception, layered Darcy-vG
  infiltration, root-zone ET, vadose drainage, capillary rise, asynchronous
  lateral groundwater, seepage, and gaining/losing river exchange;
- direct Penman-Monteith forcing from polygon meteorological fields;
- fine-versus-variable V-tilted convergence for the fully coupled water
  balance.
- the user-facing `flag_voronoi` dispatch, prepared-case contract, native
  velocity/hydrology/groundwater NetCDF, conservative GeoTIFF remapping,
  and a compact six-panel state figure.

Generate the fixtures with the scripts in `HydroPolMesh/tests`, add
`HydroPol2D_Functions` and this directory to the MATLAB path, then run the
corresponding `run_voronoi_*_validation` functions. GPU execution deliberately
fails preflight until the CPU suite, performance benchmark, and device-only
timestep loop are frozen.

The coupled runner accepts `config.hydrology_enabled=true` and an
`n_cell`-by-`n_hru` `config.hydrology.hru_fraction` matrix. Soil, vegetation,
and hydraulic properties may be scalar, one value per HRU, one value per
cell, or full cell-by-HRU matrices. Fractions must sum to one in every
polygon. Meteorological forcing can be supplied as `forcing.meteorology`
with temperature, daily extrema, day of year, latitude, wind speed, and
optional humidity, Krs, albedo, and ground heat flux. Prescribed potential
ET remains available through `forcing.potential_et_m_s`.

Groundwater recharge and capillary exchange accumulate at the surface-flow
timestep. Lateral Boussinesq flow and river exchange run at
`groundwater_update_interval_s`, with conservative internal subcycling when
required by the groundwater stability limit. Default groundwater boundaries
are no-flow.

For the complete vertical slice, generate the V-tilted fixtures and run:

```matlab
run_voronoi_vtilted_end_to_end_validation(mesh_directory,output_directory)
```

This writes a reusable `vtilted-variable-case.mat`. In the ordinary launcher,
set `flag_voronoi=1` and set `Voronoi case file` (Excel) or
`InputData_Bypass.Voronoi.case_file` (script mode) to that MAT file. With
`flag_voronoi=0`, the existing raster D4 preprocessing and solver are used
unchanged. Unresolved rivers accept only `neal_subgrid` or `none`.

## Solver scope and stability limits

The UGRID runner is CPU-only. `local_inertial` supports both resolved 2D
flow and the generalized Neal channel graph. `kinematic`, explicit
`diffusive`, and `full_momentum` are currently resolved-2D solvers: the
preflight rejects an unresolved channel graph with any of those modes.

Each solver uses the active limiter below, in addition to the configured
minimum and maximum timestep and any forcing or groundwater update boundary:

- local inertial: \(\Delta t=C\min_i A_i/\sum_e L_e(|u_e|+\sqrt{gh_e})\);
- kinematic: \(\Delta t=C\min_i A_i/\sum_e L_e(5|q_e|/(3h_e))\), including
  normal- and critical-flow boundary outflow;
- explicit diffusive: \(\Delta t=C\min_i A_i/\sum_e K_e\), including stage,
  normal-flow, and critical-flow boundary conductance;
- full momentum: the all-face \(A/\sum L(|u|+c)\) CFL, capped at a Courant
  number of 0.45;
- Neal links and resolved/subgrid transitions: the minimum of the link
  travel limit \(L/(|u|+c)\) and the graph-node storage limit
  \(A_{node}/\sum W(|u|+c)\).

`run_voronoi_channel_cfl_validation` generates its own three-node Neal reach
and resolved transition. It verifies that the graph-node limiter is stricter
than an individual-link wave limit, transmits through both links and the
transition, and closes mass exactly. It is intentionally self-contained, so
it does not depend on an application mesh.
