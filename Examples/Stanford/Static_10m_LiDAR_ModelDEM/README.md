# Stanford 10 m LiDAR Static Setup

This folder contains the prepared Stanford 10 m terrain setup used for current-model QA after the raw clipped USGS DEM showed isolated single-cell terrain artifacts.

- `DEM.tif` is copied from `/Users/mngomes/Downloads/Static/DEM.tif`.
- The DEM is a real 10 m EPSG:3310 Stanford terrain raster from the downloaded Stanford setup.
- `LULC.tif`, `SOIL.tif`, `DTB.tif`, `LAI.tif`, `Albedo.tif`, and `Initial_Soil_Moisture.tif` are copied from the current Stanford example static rasters and are resampled/aligned by HydroPol2D preprocessing.

The raw clipped USGS/3DEP DEM is preserved separately in `Static_10m_LiDAR`. It is not used for website-ready claims until the isolated terrain artifacts are conditioned.
