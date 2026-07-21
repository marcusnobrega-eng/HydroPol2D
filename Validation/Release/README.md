# Terrain Runtime Verification

These utilities document the release checks for the bundled terrain runtime.
They are run from a clean MATLAB path so that `GRIDobj`, `FLOWobj`, and
`STREAMobj` must resolve inside `third_party/topotoolbox_lite/`.

## Utilities

- `trace_topotoolbox_lite_runtime.m` profiles the active terrain workflow and
  returns the bundled files reached at runtime.
- `review_topotoolbox_lite_dependencies.m` records MATLAB's broader static
  dependency view and checks it against `MANIFEST.txt`.
- `build_topotoolbox_lite_manifest.m` rebuilds the manifest from a local
  upstream snapshot and the flag-dependent terrain methods.
- `run_topotoolbox_lite_verification.m` compares a pre-release full snapshot
  with the curated runtime. The full snapshot is deliberately not part of a
  release checkout.

The runtime has a TopoToolbox v2.4 base. Three exact later upstream methods
are retained for behavior continuity and MATLAB `linprog` compatibility. The
method-level commits are recorded in `third_party/topotoolbox_lite/UPSTREAM.md`.
