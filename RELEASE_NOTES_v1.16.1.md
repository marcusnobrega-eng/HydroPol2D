# HydroPol2D v1.16.1

## Purpose

This patch release consolidates the current HydroPol2D source code and the
Phase 1 formula-and-dynamics validation campaign into one runnable MATLAB
repository.

## Corrections

- Corrects cumulative-infiltration storage and map export to use metres
  consistently.
- Adds a drainage-area fallback to groundwater diagnostics for controlled
  domains that do not define a catchment-area field.
- Corrects the full-momentum hydrostatic reconstruction regression test and
  adds a lake-at-rest benchmark over stepped terrain.
- Corrects water-quality washoff time-step selection and outlet mass removal.
- Repairs the validation registry, runner, and analytical hillslope
  references so every public case points to an executable driver.

## Validation

The Phase 1 release suite was run from a clean MATLAB path with the bundled
TopoToolbox runtime only. It completed 23 checks covering interception, snow,
layered infiltration, evapotranspiration, groundwater and recharge,
hydrodynamics, routing, boundaries, reservoir dynamics, spatial rainfall,
water quality, human-risk classification, and the Neal (2012) channel
extension.

Run the same suite from the repository root:

```bash
python3 Validation/scripts/run_phase1_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab
```

Then check the registry:

```bash
python3 Validation/scripts/audit_validation_registry.py
```

The local-inertial and diffusive Ritter comparisons remain diagnostic, because
those approximations do not represent a full-momentum dry-bed dam break. The
staged Neal (2012) overbank case is reported with its late-recession
limitation. These results are retained explicitly and do not mask any
report-ready failure.

## Runtime

The repository remains self-contained with respect to TopoToolbox. It bundles
the required GPLv3 TopoToolbox v2.4-derived runtime in
`third_party/topotoolbox_lite/`; no external TopoToolbox path is configured by
users. Mapping Toolbox and Image Processing Toolbox remain required.
