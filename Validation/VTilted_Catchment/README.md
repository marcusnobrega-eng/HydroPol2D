# V-tilted Synthetic Catchment

This is the shared synthetic domain for HydroPol2D validation cases that require a spatial catchment setting. It provides a consistent terrain, soil, land-cover, LAI, groundwater, channel, and boundary-condition context. Each case is evaluated against its documented analytical or independent reference.

## Domain concept

The domain contains two vegetated hillslopes that drain laterally to a central channel strip. A longitudinal slope routes water toward the outlet. This setting is used for interception, snow, evapotranspiration, infiltration, groundwater, routing, boundary-condition, spatial-rainfall, water-quality, and risk tests.

Ritter dam-break and non-breaking-wave tests remain separate because their analytical solutions require dedicated geometries.

## Files

- `Config/Domain_Config.csv`: domain parameters.
- `Config/VTilted_Case_Map.csv`: use of the domain by each case.
- `generate_vtilted_domain.m`: MATLAB domain and raster generator.
- `generate_vtilted_domain.py`: Python domain and raster generator.
- `Static/`: generated static rasters.
- `Forcing/`: prescribed forcing inputs.
- `Outputs/Validation/`: shared diagnostics.

## Generate the static rasters

```bash
python3 Validation/VTilted_Catchment/generate_vtilted_domain.py
```

or in MATLAB:

```matlab
run('Validation/VTilted_Catchment/generate_vtilted_domain.m')
```
