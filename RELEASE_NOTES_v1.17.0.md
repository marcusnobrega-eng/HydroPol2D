# HydroPol2D v1.17.0

## Overview

HydroPol2D v1.17.0 consolidates the current MATLAB model, its input
workflows, and its reproducible validation suite in one self-contained
repository.

## Model and Inputs

- Retains the bundled GPLv3 TopoToolbox v2.4-derived runtime, so users do not
  need to configure an external TopoToolbox installation.
- Consolidates the current layered-soil, asynchronous groundwater, snow,
  evapotranspiration, routing, water-quality, calibration, and particle-filter
  workflows.
- Keeps spreadsheet and bypass-script inputs aligned with the active model
  configuration.

## Validation

- Standardizes validation paths, case identifiers, registry entries, and
  release commands under `Validation/`.
- Corrects cumulative-infiltration storage and map export, groundwater
  diagnostic drainage area, full-momentum hydrostatic reconstruction, and
  water-quality washoff time-step and outlet-mass handling.
- Provides one release command that runs 23 checks for interception, snow,
  infiltration, evapotranspiration, groundwater, hydrodynamics, routing,
  boundaries, reservoirs, rainfall, water quality, human risk, and the
  Neal (2012) channel extension.

Run the suite from the repository root:

```bash
python3 Validation/scripts/run_validation_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab
python3 Validation/scripts/audit_validation_registry.py
```

The local-inertial and diffusive Ritter comparisons are retained as
method-limited diagnostics because those approximations do not represent a
full-momentum dry-bed dam break. The Neal (2012) staged-overbank benchmark
documents its late-recession limitation.

## MATLAB Requirements

HydroPol2D requires Mapping Toolbox and Image Processing Toolbox. Optimization
Toolbox is required only when DEM smoothing is enabled.
