# Voronoi finite-volume validation

This directory validates the separate Voronoi local-inertial and hybrid Neal
runner. It does not alter or replace the raster D4 validation cases.

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
