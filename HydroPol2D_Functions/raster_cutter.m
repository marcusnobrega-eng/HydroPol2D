function rainfall_raster = raster_cutter(DEM_raster, rR, rr2, n)
%RASTER_CUTTER  Resample/reproject an input raster (rainfall/ETP/etc) to DEM grid.
%
% Robust to:
%   - GeographicCellsReference (lon/lat)
%   - MapCellsReference (projected), even if ProjectedCRS is empty
%
% Output:
%   rainfall_raster : DEM-clone raster (same grid), with Z = resampled rr2
%
% NOTE: This function assumes DEM is projected (has ProjectedCRS).

% =========================================================================
% 1) DEM CRS + DEM grid
% =========================================================================
demRef = DEM_raster.georef.SpatialRef;

if ~isprop(demRef,"ProjectedCRS") || isempty(demRef.ProjectedCRS)
    error("DEM raster has no valid ProjectedCRS. DEM must be projected.");
end
DEM_CRS = demRef.ProjectedCRS;

% DEM grid (cell centers)
xDEM = demRef.XWorldLimits(1) + (0:DEM_raster.size(2)-1) * demRef.CellExtentInWorldX;
if isprop(demRef,"RowsStartFrom") && strcmpi(demRef.RowsStartFrom,"north")
    yDEM = demRef.YWorldLimits(2) - (0:DEM_raster.size(1)-1) * demRef.CellExtentInWorldY;
else
    yDEM = demRef.YWorldLimits(1) + (0:DEM_raster.size(1)-1) * demRef.CellExtentInWorldY;
end
[X_dem, Y_dem] = meshgrid(xDEM, yDEM);

% =========================================================================
% 2) Identify rainfall reference type correctly
% =========================================================================
isRainGeographic = isa(rR,'map.rasterref.GeographicCellsReference');
isRainProjected  = isa(rR,'map.rasterref.MapCellsReference') || isprop(rR,"ProjectedCRS");

Rain_CRS = [];
if isRainGeographic
    % If GeoTIFF doesn't specify GeographicCRS, assume WGS84
    Rain_CRS = geocrs(4326);
else
    % Projected rasterref
    if isprop(rR,"ProjectedCRS")
        Rain_CRS = rR.ProjectedCRS;
    end
end

% =========================================================================
% 3) Build rainfall coordinate grid (centers) in its own coordinate system
% =========================================================================
if isRainGeographic
    % Geographic: use lon/lat limits and cell extents
    lonVec = rR.LongitudeLimits(1) + (0:size(rr2,2)-1) * rR.CellExtentInLongitude;

    if isprop(rR,"RowsStartFrom") && strcmpi(rR.RowsStartFrom,"north")
        latVec = rR.LatitudeLimits(2) - (0:size(rr2,1)-1) * rR.CellExtentInLatitude;
    else
        latVec = rR.LatitudeLimits(1) + (0:size(rr2,1)-1) * rR.CellExtentInLatitude;
    end

    [LON, LAT] = meshgrid(lonVec, latVec);

    % Convert lon/lat -> DEM projected coords
    [x_proj, y_proj] = projfwd(DEM_CRS, LAT, LON);

else
    % Projected: use world limits/extents
    xVec = rR.XWorldLimits(1) + (0:size(rr2,2)-1) * rR.CellExtentInWorldX;

    if isprop(rR,"RowsStartFrom") && strcmpi(rR.RowsStartFrom,"north")
        yVec = rR.YWorldLimits(2) - (0:size(rr2,1)-1) * rR.CellExtentInWorldY;
    else
        yVec = rR.YWorldLimits(1) + (0:size(rr2,1)-1) * rR.CellExtentInWorldY;
    end

    [X_rain, Y_rain] = meshgrid(xVec, yVec);

    % If rainfall CRS is known and differs from DEM CRS, reproject
    if ~isempty(Rain_CRS) && ~isequal(Rain_CRS, DEM_CRS)
        [lat, lon] = projinv(Rain_CRS, X_rain, Y_rain);
        [x_proj, y_proj] = projfwd(DEM_CRS, lat, lon);

    else
        % CRS missing OR already matches DEM: assume X/Y are already in DEM coords
        % Validate overlap to avoid silent nonsense
        rainXlim = [min(X_rain(:)) max(X_rain(:))];
        rainYlim = [min(Y_rain(:)) max(Y_rain(:))];

        demXlim  = demRef.XWorldLimits;
        demYlim  = demRef.YWorldLimits;

        overlapX = ~(rainXlim(2) < demXlim(1) || rainXlim(1) > demXlim(2));
        overlapY = ~(rainYlim(2) < demYlim(1) || rainYlim(1) > demYlim(2));

        if ~(overlapX && overlapY)
            error(['Rainfall rasterref is projected but CRS is missing/unknown and extents do not overlap DEM. ' ...
                   'Fix GeoTIFF CRS metadata OR read as geographic lat/lon.']);
        end

        x_proj = X_rain;
        y_proj = Y_rain;
    end
end

% =========================================================================
% 4) Interpolate onto DEM grid
% =========================================================================
% IMPORTANT: rr2 orientation must match y direction used above.
% Our yVec/latVec were constructed to match RowsStartFrom.
% So we do NOT flip blindly; instead we interpolate as scattered points.
F = scatteredInterpolant(x_proj(:), y_proj(:), double(rr2(:)), 'linear', 'none');
new_mesh = F(X_dem, Y_dem);

% =========================================================================
% 5) Output raster (clone DEM grid exactly)
% =========================================================================
rainfall_raster = DEM_raster;
rainfall_raster.Z    = new_mesh;
rainfall_raster.name = 'rainfall_resampled';

% =========================================================================
% 6) Clip using DEM mask + clean
% =========================================================================
temp = DEM_raster;
temp.Z = ~isnan(DEM_raster.Z);
rainfall_raster = clip(rainfall_raster, temp);

rainfall_raster.Z(isnan(rainfall_raster.Z)) = 0;
rainfall_raster.Z(rainfall_raster.Z > 300)  = 0;
rainfall_raster.Z(rainfall_raster.Z < 0)    = 0;

end