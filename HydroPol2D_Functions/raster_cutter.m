function [rainfall_raster] = raster_cutter(DEM_raster, rR, rr2, n)

    % === 1. Extract DEM CRS ===
    if isprop(DEM_raster.georef.SpatialRef, "ProjectedCRS")
        DEM_CRS = DEM_raster.georef.SpatialRef.ProjectedCRS;
    else
        error("DEM raster has no ProjectedCRS defined. You must define a CRS for alignment.");
    end

    % === 2. Extract Rainfall CRS ===
    if isa(rR, 'map.rasterref.GeographicCellsReference')
        isRainGeographic = true;
        Rain_CRS = geocrs(4326);  % WGS84 assumed
    elseif isprop(rR, "ProjectedCRS")
        isRainGeographic = false;
        Rain_CRS = rR.ProjectedCRS;
    else
        error("Rainfall raster reference has no CRS defined.");
    end

    crs_match = isequal(DEM_CRS, Rain_CRS);

    % === 3. Create rainfall coordinate grid ===
    xRain = rR.XWorldLimits(1) + (0:size(rr2, 2)-1) * rR.CellExtentInWorldX;
    yRain = rR.YWorldLimits(1) + (0:size(rr2, 1)-1) * rR.CellExtentInWorldY;
    [X_rain, Y_rain] = meshgrid(xRain, yRain);

    % === 4. Convert to lat/lon if needed ===
    if isRainGeographic
        lat = Y_rain;
        lon = X_rain;
    else
        [lat, lon] = projinv(Rain_CRS, X_rain, Y_rain);
    end

    % === 5. Generate DEM-aligned mesh ===
    xDEM = DEM_raster.georef.SpatialRef.XWorldLimits(1) + ...
           (0:DEM_raster.size(2)-1) * DEM_raster.georef.SpatialRef.CellExtentInWorldX;
    yDEM = DEM_raster.georef.SpatialRef.YWorldLimits(2) - ...
           (0:DEM_raster.size(1)-1) * DEM_raster.georef.SpatialRef.CellExtentInWorldY;
    [X_dem, Y_dem] = meshgrid(xDEM, yDEM);

    % === 6. Reproject and interpolate ===
    if ~crs_match
        disp('Reprojecting rainfall data to match DEM CRS...');
        [x_proj, y_proj] = projfwd(DEM_CRS, lat, lon);
        F = scatteredInterpolant(x_proj(:), y_proj(:), double(rr2(:)), 'linear', 'none');
        new_mesh = F(X_dem, Y_dem);
    else
        rr2_flipped = flipud(rr2);  % Match ascending yRain
        F = griddedInterpolant({yRain, xRain}, double(rr2_flipped), 'linear', 'none');
        new_mesh = F(Y_dem, X_dem);
    end

    % === 7. Build new raster reference aligned to DEM ===
    rasterSize = size(new_mesh);
    ref2 = maprefcells([min(xDEM), max(xDEM)], [min(yDEM), max(yDEM)], rasterSize);
    ref2.ProjectedCRS = DEM_CRS;
    ref2.ColumnsStartFrom = 'north';

    % === 8. Build rainfall raster object (clone from DEM) ===
    rr = DEM_raster;
    rr.Z = new_mesh;
    rr.name = 'rainfall_resampled';
    rr.size = size(new_mesh);
    rr.georef.SpatialRef = ref2;
    rr.georef.Height = rasterSize(1);
    rr.georef.Width = rasterSize(2);

    % Update RefMatrix explicitly to preserve alignment
    rr.refmat = [0 -DEM_raster.cellsize; DEM_raster.cellsize 0; ...
                 ref2.XWorldLimits(1) ref2.YWorldLimits(2)];
    rr.georef.RefMatrix = rr.refmat;

    % === 9. Crop to DEM spatial extent ===
    rainfall_raster = crop(rr, ...
        DEM_raster.georef.SpatialRef.XWorldLimits, ...
        DEM_raster.georef.SpatialRef.YWorldLimits);

    % === 10. Resample to match DEM grid (cell-by-cell) ===
    rainfall_raster = resample(rainfall_raster, DEM_raster);

    % === 11. Clip using DEM mask (remove NaNs in DEM) ===
    temp = DEM_raster;
    temp.Z = ~isnan(DEM_raster.Z);
    rainfall_raster = clip(rainfall_raster, temp);

    % === 12. Clean rainfall raster values ===
    rainfall_raster.Z(isnan(rainfall_raster.Z)) = 0;
    rainfall_raster.Z(rainfall_raster.Z > 300) = 0;
    rainfall_raster.Z(rainfall_raster.Z < 0) = 0;

end
