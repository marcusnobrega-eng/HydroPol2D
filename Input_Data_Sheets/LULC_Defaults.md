# LULC Default Parameter Basis

`LULC_parameters.xlsx` provides a consistent starting table for a new
HydroPol2D case. It does not provide universal physical constants. Terrain
resolution, vegetation structure, season, water depth, soil depth, and local
management all affect the effective parameters used by a distributed model.

## Surface Routing, Rooting, and ET

| LULC class | Manning n [s m^-1/3] | Root depth [m] | Kc [-] | Basis for the generic value |
| --- | ---: | ---: | ---: | --- |
| Tree cover | 0.100 | 1.50 | 1.05 | Forest roughness within standard 2D floodplain guidance; deep woody root zone. |
| Shrubland | 0.080 | 1.00 | 0.70 | Shrub/scrub roughness range; intermediate woody rooting. |
| Grassland | 0.060 | 0.60 | 0.85 | Conservative shallow-flow roughness; shallow herbaceous rooting. |
| Cropland | 0.050 | 1.00 | 1.00 | Upper end of a generic crop roughness range; mid-season crop ET starting value. |
| Built-up | 0.030 | 0.00 | 0.00 | Smooth sealed-surface value. Use a separate urban class with larger n where buildings or obstructions must be represented. |
| Bare / sparse vegetation | 0.030 | 0.10 | 0.30 | Barren land roughness range; minimal vegetation water uptake. |
| Snow and ice | 0.020 | 0.00 | 0.00 | Initial smooth snow/ice runoff surface; revise for rough ice, debris, or crevassed terrain. |
| Permanent water bodies | 0.035 | 0.00 | 0.00 | Within the standard open-water range. Ponded cells evaporate at `Ep`, not through Kc. |
| Herbaceous wetland | 0.070 | 0.40 | 1.00 | Emergent-wetland roughness range and shallow, water-table-limited rooting. |
| Mangroves | 0.150 | 1.00 | 1.05 | High-end woody-wetland roughness; revise using local vegetation density and tidal/channel geometry. |
| Moss and lichen | 0.080 | 0.05 | 0.40 | Conservative shallow-flow resistance and near-surface water uptake. |

The roughness values follow the class ranges in the USACE HEC-RAS 2D land-cover
guidance. That guidance also notes that its tabulated values are intended for
appreciable flow depths; shallow overland flow can require larger effective
values. Root-depth values are simple representative depths consistent with
global land-model parameterizations, rather than estimates of the deepest roots
at a site. `Kc` is applied only to internally computed reference ET and is an
initial, time-invariant value. Seasonal crop curves, remotely sensed ET, or
locally calibrated values are preferable when available.

`h0` and `d0` both default to zero. `h0` is used by the cellular-automata
routing option only; it is not a substitute for canopy storage or soil
infiltration. `d0` represents the initial surface-water state and should be set
from a restart or an initial-depth raster when the event begins on a wet
surface.

## Water Quality

The generic table sets `C1 = C2 = C3 = 0` and `C4 = 1`. This deliberately
represents no initial pollutant buildup and no washoff. There is no defensible
universal set of buildup/washoff coefficients by land-cover class alone. In the
current mass-based HydroPol2D formulation,

```
B0 = C1 [1 - exp(-C2 ADD)] A / 10000
W = C3 Q^C4 B
```

where `B0` is buildup mass, `ADD` is antecedent dry days, `A` is cell area in
m2, `Q` is runoff in m3 s^-1, and `B` is the available pollutant mass. `C1`
is in kg ha^-1, `C2` is in d^-1, `C4` is dimensionless, and the units of `C3`
depend on `C4`. Define these values for the pollutant, measurement units, and
case study, then calibrate them against concentration or load observations.

## Snow

The snow columns use broadly plausible initial values: snow albedo of 0.55-0.80,
emissivity of 0.98-0.99, a linear rain/snow transition of -1 to 2 degC,
fresh-snow density of 100 kg m^-3, and maximum seasonal-snow density of
450 kg m^-3 (550 kg m^-3 for the snow-and-ice class). Degree-day melt,
sublimation, and compaction coefficients are empirical rate parameters. They
are retained as explicit inputs because they should be fitted or evaluated
against local SWE, snow-depth, or runoff data; they should not be interpreted
as transferable LULC constants.

## Sources

- USACE HEC-RAS: [land cover and Manning n guidance](https://www.hec.usace.army.mil/confluence/rasdocs/r2dum/6.5/developing-a-terrain-model-and-geospatial-layers/creating-land-cover-mannings-n-values-and-impervious-layers).
- Zeng (2001): [Global Vegetation Root Distribution for Land Modeling](https://doi.org/10.1175/1525-7541(2001)002%3C0525:GVRDFL%3E2.0.CO;2).
- FAO-56: [crop evapotranspiration and crop-coefficient method](https://www.fao.org/4/x0490e/x0490e0a.htm).
- Hock (2003): [temperature-index melt modelling](https://doi.org/10.1016/S0022-1694(03)00257-9).
- Dai (2008): [rain-snow phase transition](https://doi.org/10.1029/2008GL033295).
- Vionnet et al. (2020): [seasonal snow-density evolution](https://doi.org/10.5194/tc-14-1829-2020).
- U.S. EPA SWMM: [Reference Manual, Volume III: Water Quality](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P100P2NY.TXT).
