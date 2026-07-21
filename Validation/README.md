# HydroPol2D Validation

This folder contains the HydroPol2D validation suite. The cases assess equations, numerical dynamics, mass conservation, and boundary handling against analytical solutions, independent reference calculations, benchmark hydrographs, or refined numerical references.

[Validation_Cases.csv](Validation_Cases.csv) is the executable registry. Each row identifies the model component, driver, reference source, acceptance criterion, and evidence status. The full report is maintained in the documentation repository.

## Run the validation suite

```bash
python3 Validation/scripts/run_validation_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab
```

The runner starts each test in a clean MATLAB process. It writes a compact summary to `Validation/Release/Outputs/Validation_Release_Summary.csv` and returns a nonzero status when a report-ready case fails.

Check the registry after a run:

```bash
python3 Validation/scripts/audit_validation_registry.py \
  Validation/Validation_Cases.csv
```

## Evidence status

- `report_ready`: meets its stated acceptance criterion.
- `diagnostic`: retained to document the behavior or limitation of an approximation; not used as acceptance evidence.
- `limited`: meets its primary metrics and has a stated limitation that accompanies the result.

## Coverage

The suite covers canopy interception, snow, layered infiltration, evapotranspiration, recharge and groundwater flow, full-momentum and local-inertial hydrodynamics, cellular automata, kinematic and diffusive routing, reservoir dynamics, inflow and stage boundaries, spatial rainfall, water quality, human-risk classification, and the Neal (2012) local-inertial channel extension.
