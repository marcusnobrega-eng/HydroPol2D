# V-Tilted Infiltration Validation Suite

This suite validates meaningful infiltration regimes on the Phase 1 V-tilted
catchment context. The cases are not fallback tests. They are hydrologic
regime checks for rainfall supply, hydraulic intake capacity, finite vadose
storage, layer percolation, shallow groundwater coupling, and spatial soil
contrast.

The V-tilted catchment provides the spatial case-study framing. The first
diagnostic runner uses representative V-tilted soil columns so each regime has
a clear reference expectation before full HydroPol2D hydrograph runs are used.
The case registry uses physically interpretable sandy-loam, loam, and
clay-loam hydraulic properties. Conductivity multipliers are limited to
plausible structural contrasts such as compacted soil, structured topsoil, and
moderately slower transmission-zone flow.

## Cases

- `P1-INFIL-VT-000`: no-infiltration loam baseline.
- `P1-INFIL-VT-001`: supply-limited sandy-loam infiltration.
- `P1-INFIL-VT-002`: capacity-limited clay-loam infiltration.
- `P1-INFIL-VT-003`: shallow sandy-loam storage-excess transition.
- `P1-INFIL-VT-004`: moist sandy-loam layered percolation and recharge delay.
- `P1-INFIL-VT-005`: shallow-groundwater loam saturation-excess response.
- `P1-INFIL-VT-006`: sandy-loam conductivity contrast.

## Run

```matlab
run('HydroPol2D_Model/Validation/Infiltration/VTilted_Infiltration/run_vtilted_infiltration_validation_suite.m')
```

## Outputs

- `Outputs/Validation/VTilted_Infiltration_Diagnostics.csv`
- `Outputs/Validation/VTilted_Infiltration_Pass_Fail.csv`
- `Outputs/Validation/TimeSeries/<case_id>.csv`
- `Outputs/Validation/Figures/*hydrograph_infiltration.{png,pdf}`
- `Outputs/Validation/Figures/VTilted_Infiltration_*_Overview.{png,pdf}`
