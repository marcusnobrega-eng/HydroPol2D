# Phase 1 V-tilted Synthetic Catchment

This is the shared synthetic case-study domain for Phase 1 formula and dynamics verification.

The purpose is not to replace analytical solutions. The v-tilted catchment provides a consistent HydroPol2D spatial domain in which analytical/reference solutions can be embedded, sampled, and compared.

## Domain concept

The domain has three zones:

- left hillslope: vegetated pervious hillslope;
- channel strip: low-LAI channel/bare strip;
- right hillslope: denser vegetated pervious hillslope.

The terrain slopes laterally toward the central channel and longitudinally toward the outlet. This gives a simple spatial context for interception, snow, ET, infiltration, groundwater, routing, boundary-condition, rainfall-raster, water-quality, and risk checks.

## Phase 1 rule

Use this v-tilted catchment as the HydroPol2D case-study domain whenever a module can be tested spatially while still comparing to a documented analytical or independent reference solution.

Pure hydrodynamic benchmarks with their own required geometry, such as Ritter dam-break and non-breaking wave, remain separate benchmark domains.

## Files

- `Config/Domain_Config.csv`: canonical parameters.
- `Config/Phase1_VTilted_Case_Map.csv`: how each Phase 1 case uses the domain.
- `generate_phase1_vtilted_domain.m`: MATLAB raster/domain generator writing to this folder.
- `Static/`: generated static rasters.
- `Forcing/`: generated or hand-supplied Phase 1 forcing inputs.
- `Outputs/Validation/`: shared domain diagnostics.

## Generate static rasters

Run in Python from the repository root:

```bash
python3 HydroPol2D_Model/Validation/Phase1_VTilted_Catchment/generate_phase1_vtilted_domain.py
```

Or run in MATLAB:

```matlab
run('HydroPol2D_Model/Validation/Phase1_VTilted_Catchment/generate_phase1_vtilted_domain.m')
```

Both generators are local to this repository and do not write to external `/oak/...` paths.
