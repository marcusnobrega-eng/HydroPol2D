# Stanford 10 m Static Inputs

This folder contains the real 10 m Stanford terrain setup for current-model
reruns (EPSG:3310, 10 m, nodata -9999, clipped to the 158.1 km2 catchment).

- `DEM.tif` is the example's earlier 10 m USGS 3DEP DEM, completed on 2026-11-03
  (`tools/stanford_dem/fill_10m_lidar_gap.py`):
  - its 5.9 km2 western headwater gap is filled from the USGS 3DEPElevation
    service at 10 m (downloaded 2026-11-03, bare earth, U.S. public domain),
    blended into the earlier file so the fill meets it without a step;
  - two 3DEP tile edges inside the filled area (10-15 m steps) are smoothed;
  - 26 single-cell pits and spikes (up to 82 m deep) take the median of their
    eight neighbours;
  - it is clipped to the 10 m catchment mask of `Static_10m_LiDAR_ModelDEM`.
- Outside the filled area the DEM is the earlier file. It differs from today's
  3DEP by 2.8 m RMS, most on the steep western slopes, so it is not the same
  download as the 30 m DEM in `Static`.
- `LULC.tif`, `SOIL.tif`, `DTB.tif`, `LAI.tif`, and `Albedo.tif` are the
  current Stanford ancillary rasters. They are coarser than the DEM but share
  the same domain bounds, and are aligned/resampled by HydroPol2D preprocessing
  for 10 m simulations.
- `Initial_Soil_Moisture.tif` is copied from the current Stanford static setup
  and is likewise aligned during preprocessing.

Do not use `Static_10m_Upsampled` for website or validation claims; that folder
is diagnostic only because its DEM was derived from the 30 m grid.
