# Phase 1: Formula and Dynamics Verification

This phase intentionally stays simple. Its purpose is to validate the equations, storage updates, numerical dynamics, and parameter sensitivity before attempting complex field validation.

Strictly speaking, most Phase 1 cases are verification or benchmark tests rather than field validation. They still build confidence because they answer four foundational questions:

1. Does each formula reproduce a known analytical or independently computed reference?
2. Does each storage/dynamic module conserve mass?
3. Does each coupled process respond in the expected direction and magnitude?
4. Can simple inverse tests recover known synthetic parameters?

Phase 2 will handle field observations, split-sample calibration/validation, remote sensing, high-water marks, and full applied case studies.

The executable Phase 1 registry is `Phase1_Cases.csv`. Analytical/reference targets generated so far are stored in `Reference_Outputs/Phase1/`.

## Shared case-study domain

The shared Phase 1 case-study domain is the v-tilted synthetic catchment:

`HydroPol2D_Model/Validation/Phase1_VTilted_Catchment`

Use this domain for equation and dynamics tests whenever the analytical/reference solution can be embedded spatially. The analytical solution remains the truth source; the v-tilted catchment supplies a consistent HydroPol2D DEM, LAI, soil, land-cover, groundwater, and channel geometry context.

Exceptions:

- Ritter dam-break keeps its own idealized dam-break geometry.
- The non-breaking wave benchmark keeps its own wave-domain geometry.
- Any future benchmark with required laboratory/analytical geometry should remain separate and be cross-referenced from Phase 1.

## Evidence types allowed in Phase 1

| evidence type | meaning | how it should be used |
|---|---|---|
| Analytical solution | Closed-form or semi-analytical solution exists. | Best evidence for formulas/dynamics. |
| Independent reference calculation | Same formula computed outside HydroPol2D. | Good for algebraic modules such as ET, risk, rating curves, and water quality. |
| Method of manufactured solutions | Force a known solution into the governing equations. | Useful when no natural analytical solution exists. |
| Synthetic parameter recovery | Generate synthetic observations from known parameters, then recover them through inverse modeling. | Tests identifiability and implementation, not field truth. |
| Cross-model benchmark | Compare against a trusted simple solver for the same idealized problem. | Useful for coupled modules, but weaker than analytical solutions. |

## Recommended Phase 1 cases

| module | preferred Phase 1 truth source | recommended case |
|---|---|---|
| Canopy interception | Analytical/intermediate reference | Rutter/Gash-style storage capacity test: gross rainfall, canopy capacity, evaporation, throughfall, and drainage. |
| Snow | Analytical degree-day and storage balance | Prescribed temperature and precipitation sequence with known rain/snow partition, melt, sublimation, and SWE. |
| Infiltration | Analytical solution | Green-Ampt or Philip infiltration into a homogeneous soil column under constant rainfall. |
| Evapotranspiration | Independent reference calculation | FAO-56 reference ET plus constrained soil/open-water extraction tests. |
| Groundwater/recharge | Analytical or manufactured solution | 1D/2D Boussinesq-style recharge mound, linear reservoir recession, or synthetic water-table recovery. |
| Hydrodynamics | Analytical solution | Ritter dam-break, non-breaking wave, steady uniform flow, and gradually varied/backwater profile where applicable. |
| Routing alternatives | Analytical solution after audit | Kinematic-wave hydrograph and diffusive-wave attenuation on a plane with known slope/roughness. |
| Subgrid hydraulics | Independent reference/cross-model | Trapezoidal or rectangular channel rating/conveyance calculation compared with coarse and subgrid representations. |
| Reservoir/control structures | Analytical ODE | Storage depletion under rating curve `Q = a(h-h0)^b`, including exact or high-accuracy reference integration. |
| Inflow hydrograph BC | Exact bookkeeping/reference routing | Rectangular domain with prescribed inflow volume and simple expected travel/storage response. |
| Stage hydrograph BC | Analytical/simple reference | Prescribed boundary stage in a rectangular channel with expected profile or quasi-steady response. |
| Spatial rainfall | Exact raster bookkeeping | Tiny raster stack with known cellwise depth, duration, and domain volume. |
| Water quality | Analytical mass balance | Build-up/wash-off under fixed runoff flux with known pollutant mass history. |
| Human risk | Independent reference calculation | Deterministic depth/velocity grid checked against published or documented threshold equations. |

## Simple inverse-modeling cases

Synthetic inverse tests should be used sparingly and labeled as parameter-recovery tests.

Recommended inverse tests:

- Infiltration: recover saturated hydraulic conductivity or sorptivity from synthetic cumulative infiltration.
- Routing: recover Manning roughness from a synthetic outlet hydrograph after the routing implementation audit.
- Groundwater: recover recession constant or hydraulic diffusivity from synthetic water-table/baseflow recession.
- Reservoir: recover rating-curve coefficients from synthetic stage/outflow pairs.
- Water quality: recover wash-off coefficient from synthetic concentration/load data.

Rules for inverse tests:

- Generate synthetic observations from a known parameter set.
- Add one no-noise case and one controlled-noise case.
- Report recovered parameter, uncertainty, objective function, residual structure, and whether the true value lies inside the confidence/credible interval.
- Do not call parameter recovery field validation.

## Phase 1 acceptance standard

For a Phase 1 case to be accepted:

- the reference solution or synthetic truth must be documented;
- all units and conversions must be explicit;
- mass residual must pass the module-specific tolerance;
- no unexplained `NaN` may appear in summaries;
- plots must compare HydroPol2D against the reference;
- pass/fail must be tied to numerical thresholds.

## Reference-output generator

Several Phase 1 cases have simple analytical or bookkeeping references that can be generated before running HydroPol2D:

```bash
python3 HydroPol2D_Model/Validation/scripts/phase1_reference_solutions.py --case all --output HydroPol2D_Model/Validation/Reference_Outputs/Phase1
```

The generated CSVs are truth targets for comparison. They are not HydroPol2D outputs.

## Priority order

1. Hydrodynamics: Ritter, non-breaking wave, steady/backwater profile.
2. Infiltration: Green-Ampt or Philip column/plane case.
3. Reservoir and boundary conditions: exact volume/rating-curve tests.
4. ET, canopy, snow, water quality, human risk: formula/reference checks.
5. Groundwater and routing alternatives: analytical/manufactured tests after implementation audit.
6. Subgrid and spatial rainfall: bookkeeping/reference conveyance tests.

This order gives confidence in the numerical core before we spend effort on observational data, uncertain parameters, and complex applied cases.
