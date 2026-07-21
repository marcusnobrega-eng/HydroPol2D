# Neal 2012 Simple Straight-Channel Examples

This folder contains simple, low-ambiguity benchmarks for the Neal et al. (2012) channel-subgrid formulation in HydroPol2D.

The goal is to isolate the effect of sub-resolution channel geometry without rainfall-runoff generation, infiltration, evapotranspiration, groundwater, or complex floodplain geometry.

## Model comparison

Each case compares:

1. A fine resolved one-row channel reference.
2. A coarse ordinary local-inertial model.
3. A coarse Neal-mode local-inertial model using:
   - `flag_subgrid = 1`
   - `flag_overbanks = 1`

## Geometry

- Fine channel width: `10 m`
- Coarse cell size: `30 m`
- Fine cell size: `10 m`
- Bank height: `1 m`
- Bed slope: `0.001 m/m`
- Manning `n`: `0.035`

## Cases

- `VAL-SUBGRID-NEAL-SIMPLE-001`
  - 5 coarse cells long
  - constant within-bank inflow

- `VAL-SUBGRID-NEAL-SIMPLE-002`
  - 30 coarse cells long
  - constant within-bank inflow

- `VAL-SUBGRID-NEAL-SIMPLE-003`
  - 10 coarse cells long
  - wet-start pulse that crosses bankfull

## Current interpretation

The first two examples do exactly what we wanted from a simple demonstration:

- the coarse Neal run is much closer to the fine reference than the coarse ordinary run;
- hydrograph timing, stage, and longitudinal profile all improve strongly;
- mass balance is essentially exact.

The third example remains diagnostic. It is useful because it shows that the current Neal-mode implementation still struggles once the pulse strongly crosses bankfull and drains back out.
