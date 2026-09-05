# HydroPol2D MATLAB v1.18.0-rc1

## Purpose

This release candidate prepares the MATLAB model to consume the same HydroBathyDEM `0.3.0rc1`
mesh package as HydroPol2D-Python `0.8.0rc1`. The stable public MATLAB reference remains
`v1.17.0` until coordinated validation is complete.

## Main changes

- Strict readers and validators for mesh-contract `1.0` UGRID, overlap, and subgrid tables.
- Explicit zero-based file indexing with one-based conversion only inside MATLAB.
- Separate Neal channel, structured lookup subgrid, and Voronoi subgrid execution paths.
- Public configuration now uses `flag_neal_channel`, `flag_structured_lookup_subgrid`, and
  `flag_voronoi_subgrid`; legacy `flag_subgrid` and `flag_overbanks` remain migration inputs only.
- CPU Voronoi routing, hydrology, groundwater, boundary, output, and post-processing support.
- Courant limit of `0.2` by default and rejection above `0.3` for Voronoi runs.
- Correct output-time alignment and structured-lookup roughness convention checks.

## Release gates

MATLAB R2025b validation completed on the integration workstation: 23/23 raster checks passed,
and the 47-case registry completed with no errors and three documented warnings. Shared-bundle
Pune baseline and Voronoi-subgrid runs completed in MATLAB and Python at 2 h, 6 h, and 12 h.
The largest cumulative outlet-volume difference was `0.03779%` of total model input, below the
integrated-domain threshold of `0.1%`; both solvers independently closed mass to numerical
precision.

This candidate is not yet a stable release. Clean-clone verification, final Studio regression,
documentation review, and coordinated release-candidate packaging remain open.
