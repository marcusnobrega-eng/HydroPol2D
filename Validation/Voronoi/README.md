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

Generate the fixtures with the scripts in `HydroPolMesh/tests`, add
`HydroPol2D_Functions` and this directory to the MATLAB path, then run the
corresponding `run_voronoi_*_validation` functions. GPU execution deliberately
fails preflight until the CPU suite, performance benchmark, and device-only
timestep loop are frozen.
