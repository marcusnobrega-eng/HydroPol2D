# HydroPol2D Validation

The HydroPol2D validation suite evaluates the implemented equations, storage updates, boundary conditions, numerical dynamics, and parameter recovery under conditions with a documented reference.

The suite uses five evidence types:

1. Analytical or semi-analytical solutions.
2. Independent reference calculations.
3. Exact mass, water, and boundary-volume bookkeeping.
4. Documented benchmark hydrographs.
5. Refined numerical references for cases without an analytical solution.

The executable registry is [Validation_Cases.csv](Validation_Cases.csv). Each row identifies the module, driver, reference source, metrics, acceptance criteria, and evidence status. Reference outputs are stored under `Reference_Outputs/Validation/`.

## Shared V-tilted catchment

The V-tilted synthetic catchment provides a common spatial setting for hydrologic and groundwater validation. It includes terrain, soil, land cover, LAI, groundwater, channel, and boundary-condition inputs. A case still passes or fails against its analytical or independent reference, not against the catchment itself.

Ritter dam-break, non-breaking wave, and channel-subgrid benchmarks use their own geometries where required by the reference problem.

## Voronoi finite-volume mode

The Voronoi runner is a separate CPU-only mesh mode. Its local-inertial
solver supports the generalized Neal channel graph; kinematic, explicit
diffusive, and full-momentum modes currently require fully resolved 2D flow.
The specific tests, equations, acceptance checks, and measured results are
documented in [Voronoi/README.md](Voronoi/README.md) and
[Voronoi/VALIDATION_RESULTS.md](Voronoi/VALIDATION_RESULTS.md). GPU and
India-scale applications remain explicitly outside the validated scope.

## Acceptance requirements

A validation case must document its reference, units, diagnostics, metrics, acceptance threshold, and pass/fail result. It must close the relevant mass balance, contain no unexplained `NaN` values, and compare HydroPol2D outputs directly with the reference.

## Run the suite

```bash
python3 Validation/scripts/run_validation_release.py \
  --matlab /Applications/MATLAB_R2025b.app/bin/matlab

python3 Validation/scripts/audit_validation_registry.py \
  Validation/Validation_Cases.csv
```

The release runner starts each driver in a clean MATLAB process and records the case summary in `Validation/Release/Outputs/Validation_Release_Summary.csv`.
