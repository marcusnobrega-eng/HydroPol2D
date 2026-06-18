# Generated Static Rasters

Run `../generate_phase1_vtilted_domain.py` with Python or `../generate_phase1_vtilted_domain.m` in MATLAB to generate the v-tilted Phase 1 static rasters here.

Expected outputs include:
- `DEM.tif`
- `LULC.tif`
- `SOIL.tif`
- `DTB.tif`
- `Albedo.tif`
- `LAI.tif`
- `Initial_SM.tif`
- `GW_depth.tif`
- `GW_table.tif`
- `Zone_ID.tif`
- `overview.png`

The Python generator uses `rasterio`; the MATLAB generator is provided for consistency with HydroPol2D workflows.
