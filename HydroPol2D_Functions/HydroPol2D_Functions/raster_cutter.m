function [rainfall_raster] = raster_cutter(DEM_raster,rR,rr2,n)
    if n==1
        % Reprojecting the latitude and longitude
        [x, y] = projfwd(projcrs(3857), rR.LatitudeLimits, rR.LongitudeLimits);
        % Average dimention of the pixel size acordding the ecuator
        % Accuracy could be affected according the latitude
        cs_deg = km2deg(DEM_raster.cellsize/1000);
        % resampling the rainfall mesh
        % Making base mesh with data
        [temp_x,temp_y] = ndgrid(1:rR.CellExtentInLatitude*1000000:abs(rR.LatitudeLimits(1)-rR.LatitudeLimits(2))*1000000, ...
            1:rR.CellExtentInLatitude*1000000:abs(rR.LongitudeLimits(1)-rR.LongitudeLimits(2))*1000000);
        F = griddedInterpolant(temp_x,temp_y,double(rr2),'nearest');
        % Making new mesh
        [nx,ny] = ndgrid(1:cs_deg*1000000:abs(rR.LatitudeLimits(1)-rR.LatitudeLimits(2))*1000000, ...
            1:cs_deg*1000000:abs(rR.LongitudeLimits(1)-rR.LongitudeLimits(2))*1000000);
        new_mesh = F(nx,ny);
    
        % New raster size, creating Gridobj
        rasterSize = [size(ny,1) size(nx,2)];
        % Creating the reference for the cropped area
        ref2 = maprefcells(x, y, rasterSize);
        ref2.ProjectedCRS = projcrs(3857);
        ref2.ColumnsStartFrom = 'north';
        % Making a copy of the DEM_raster to generate the Satellite raster
        rr = DEM_raster; rr.Z = new_mesh; rr.name = 'satellite_raster';
        rr.size = size(new_mesh); rr.georef.SpatialRef = ref2; rr.georef.Height = size(ny,1);
        rr.georef.Width = size(nx,2);
        rr.refmat = [0 -DEM_raster.cellsize; DEM_raster.cellsize 0; ref2.XWorldLimits(1) ref2.YWorldLimits(2)];
        rr.georef.RefMatrix = rr.refmat;
    
        % delete('spatial_rainfall_data\rainfall_3857.tif')
        rainfall_raster = crop(rr,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
        clear rr;
        % to match with the reference DEM array size
        temp = DEM_raster; temp.Z = ~isnan(DEM_raster.Z);
        rainfall_raster = resample(rainfall_raster,DEM_raster);
        rainfall_raster = clip(rainfall_raster,temp);
        
        %filtering NaNs and outliers
        rainfall_raster(isnan(rainfall_raster)) = 0;
        rainfall_raster(rainfall_raster>300) = 0;
    
    
    elseif n==0
        % resampling the rainfall mesh
        % Making base mesh with data
        % [temp_x,temp_y] = ndgrid(1:rR.CellExtentInWorldX:abs(rR.YWorldLimits(1)-rR.YWorldLimits(2)), ...
        %     1:rR.CellExtentInWorldX:abs(rR.XWorldLimits(1)-rR.XWorldLimits(2)));

        F = griddedInterpolant(rR.temp_x,rR.temp_y,double(rr2),'nearest');

        % Making new mesh
        % [nx,ny] = ndgrid(1:DEM_raster.cellsize:abs(rR.YWorldLimits(1)-rR.YWorldLimits(2)), ...
        %     1:DEM_raster.cellsize:abs(rR.XWorldLimits(1)-rR.XWorldLimits(2)));
        new_mesh = F(rR.nx,rR.ny);
        % New raster size, creating Gridobj
        rasterSize = [size(rR.ny,1) size(rR.nx,2)];
        % Creating the reference for the cropped area
        ref2 = maprefcells(rR.SpatialRef.XWorldLimits, rR.SpatialRef.YWorldLimits, rasterSize);
        ref2.ProjectedCRS = projcrs(3857);
        ref2.ColumnsStartFrom = 'north';
        % Making a copy of the DEM_raster to generate the Satellite raster
        rr = DEM_raster; rr.Z = new_mesh; rr.name = 'satellite_raster';
        rr.size = size(new_mesh); rr.georef.SpatialRef = ref2; rr.georef.Height = size(rR.ny,1);
        rr.georef.Width = size(rR.nx,2);
        rr.refmat = [0 -DEM_raster.cellsize; DEM_raster.cellsize 0; ref2.XWorldLimits(1) ref2.YWorldLimits(2)];
        rr.georef.RefMatrix = rr.refmat;
    
        % delete('spatial_rainfall_data\rainfall_3857.tif')
        rainfall_raster = crop(rr,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
        clear rr;
        % to match with the reference DEM array size
        temp = DEM_raster; temp.Z = ~isnan(DEM_raster.Z);
        rainfall_raster = resample(rainfall_raster,DEM_raster);
        rainfall_raster = clip(rainfall_raster,temp);

        %filtering NaNs and outliers
        rainfall_raster.Z(isnan(rainfall_raster.Z)) = 0;
        rainfall_raster.Z(rainfall_raster.Z>300) = 0;
    
    end
end
