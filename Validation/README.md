# HydroPol2D Phase 1 Validation

This folder contains the controlled tests used to verify the current
HydroPol2D implementation. Phase 1 tests equations, numerical dynamics,
mass conservation, and boundary handling against analytical solutions,
independent reference calculations, or refined numerical references. It
does not constitute field calibration or field validation.

`Phase1_Cases.csv` is the release registry. Each row identifies the model
component, the driver, the truth source, the acceptance criterion, and the
evidence status. The detailed report is maintained separately in the
documentation repository.

## Run the release suite

From the repository root:

```bash
python3 Validation/scripts/run_phase1_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab
```

The runner starts every test in a clean MATLAB process. It writes a compact
summary to `Validation/Release/Outputs/Phase1_Release_Summary.csv` and
returns a nonzero exit status when a report-ready test fails. The V-tilted
infiltration formula suite and its full-domain reruns are executed separately:
the former supplies the strict criteria and the latter checks normal model
execution. Generated outputs are ignored by Git.

Run the registry check after the suite:

```bash
python3 Validation/scripts/audit_validation_registry.py Validation/Phase1_Cases.csv
```

## Evidence levels

- `report_ready`: passed the stated Phase 1 criterion.
- `diagnostic`: completed and retained to show the limits of an
  approximation; it is not a validation claim for that process.
- `limited`: passed its primary controlled metrics but has a stated
  limitation that must accompany any use of the result.

## Coverage

The release suite covers canopy interception, snow, layered infiltration,
evapotranspiration, recharge and groundwater flow, full-momentum and
local-inertial hydrodynamics, cellular automata, kinematic and diffusive
routing, reservoir dynamics, inflow and stage boundaries, spatial rainfall,
water quality, human-risk classification, and the Neal (2012) local-inertial
channel extension. The subgrid extension is tested only for its documented
channel geometries and local-inertial routing.
