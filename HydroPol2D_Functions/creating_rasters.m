% Generating Rasters

% V-Tilted
[dem] = V_Tilted_Plane_Watersheds(20,800,20,1000,0,0.05,0.02);

lulc = dem./dem;
lulc(:,41) = 2;

SaveAsciiRaster(dem,20,0,0,'DEM_V_Tilted.asc',-9999);
SaveAsciiRaster(lulc,20,0,0,'LULC_V_Tilted.asc',-9999);


% Using the Base Raster to create V-Tilted Rasters
