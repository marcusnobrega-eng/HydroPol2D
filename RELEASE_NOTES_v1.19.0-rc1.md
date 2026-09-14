# HydroPol2D MATLAB v1.19.0-rc1

## Purpose

This candidate consolidates the MATLAB runtime around the shared HydroBathyDEM `0.4.0rc1` mesh
bundle and the HydroPol2D-Python `0.9.0rc1` compatibility set.

## Main changes

- One canonical runtime path for the main loop, local-inertial solver, and inundation outputs.
- Historical `_old`, `_new`, `_Editing`, and testing variants moved outside the runtime path under
  a preserved rollback tag.
- Preallocated diagnostics and output histories with concise configurable progress reporting.
- Atomic structured checkpoints and final saves without duplicate workspace serialization.
- Reduced avoidable CPU/GPU transfers while preserving existing numerical behavior.
- Removed the arbitrary 10 m/s surface-velocity cap from regular-grid and Voronoi routing;
  stability remains controlled by the shared Courant timestep, wet/dry treatment, friction, and
  conservative draining limiters.
- Shared mesh acceptance for explicit zero-based file indexing, MATLAB one-based conversion,
  topology, initial storage, and mass closure.
- Native post-processing exports now preserve groundwater-only maps and use consistent gathered
  arrays for velocity, hazard, infiltration, and recharge products.
- The India high-resolution workflow now has portable preparation, preflight, calibration, Slurm,
  and output-summary utilities without hard-coded local runtime paths.

## Validation

- Controlled Voronoi routing closed mass with relative error `2.8623e-16` over 1,560 steps.
- Voronoi output produced six requested times and 35 variables.
- Irregular groundwater closed mass with relative error `3.793e-16`.
- Full-momentum validation retained all five regular-grid passes, and the Voronoi validation
  preserved a valid `12 m/s` state without clipping while closing rainfall-runoff mass to
  `1.21e-14`.
- Atomic checkpoint and preallocated output-history checks passed in MATLAB R2025b.
- The clean-wheel shared acceptance read 61 cells and 143 faces identically in MATLAB and Python
  and reported zero relative mass error.
- The portable India case configuration suite passes 16 checks covering environment overrides,
  preprocessing-only runs, output locations, and post-processing inputs.

This is not a stable release. Current-candidate CUDA validation, Apple Metal validation, and the
Windows packaged-Studio gate remain open. Voronoi execution remains CPU-only until a JAX/CUDA
implementation passes the coordinated backend parity suite.
