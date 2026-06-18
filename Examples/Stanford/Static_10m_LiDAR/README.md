# Stanford 10 m Static Inputs

This folder contains the real 10 m Stanford terrain setup for current-model
reruns.

- `DEM.tif` is the USGS/3DEP 10 m DEM from
  `/Users/mngomes/Downloads/USGS_10m_DEM_clipped.tif`.
- `LULC.tif`, `SOIL.tif`, `DTB.tif`, `LAI.tif`, and `Albedo.tif` are the
  current Stanford ancillary rasters. They are coarser than the DEM but share
  the same domain bounds, and are aligned/resampled by HydroPol2D preprocessing
  for 10 m simulations.
- `Initial_Soil_Moisture.tif` is copied from the current Stanford static setup
  and is likewise aligned during preprocessing.

Do not use `Static_10m_Upsampled` for website or validation claims; that folder
is diagnostic only because its DEM was derived from the 30 m grid.
