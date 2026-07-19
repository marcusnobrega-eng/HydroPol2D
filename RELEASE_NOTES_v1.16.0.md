# HydroPol2D v1.16.0

## Bundled Terrain Runtime

- Bundles the required TopoToolbox v2.4-based source subset under
  `third_party/topotoolbox_lite/`.
- Removes the user-facing external TopoToolbox setting from bypass
  configuration and spreadsheet inputs.
- Registers HydroPol2D functions and the bundled runtime through
  `hydropol2d_add_runtime_paths`.
- Keeps the existing `GRIDobj`, `FLOWobj`, and `STREAMobj` interfaces so
  terrain preprocessing results remain consistent with the prior workflow.
- Retains three documented, unmodified later upstream methods needed for
  current MATLAB smoothing compatibility; see
  `third_party/topotoolbox_lite/UPSTREAM.md`.

## Requirements and Licensing

- Mapping Toolbox and Image Processing Toolbox remain required.
- Optimization Toolbox is required only when DEM smoothing is enabled.
- Parallel Computing Toolbox remains optional for GPU execution.
- HydroPol2D and the bundled TopoToolbox source are distributed under GPLv3;
  see `LICENSE` and `third_party/NOTICE.md`.

## Upgrade Notes

Remove any previous `addpath` commands for an external TopoToolbox copy. Run
`HydroPol2D_V115.m` from the repository root, or use a supported case runner.
