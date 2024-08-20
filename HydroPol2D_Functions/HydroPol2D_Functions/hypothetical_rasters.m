%% Hypotetical scenarios maker
function [LULC_raster,DEM_raster,SOIL_raster,Initial_Water_Depth] = hypothetical_rasters(LULC_raster,DEM_raster,SOIL_raster,Resolution,Width,Length,Slope,Dam_Length,Dam_Height,Dam_Width,Dam_Bottom_Elevation)
    
    ny = Length/Resolution;
    nx = Width/Resolution;
    ny_dam = Dam_Length/Resolution;
    nx_dam = Dam_Width/Resolution;


    % Internal Points of the Dam
    idx = (round(nx/2) - round(nx_dam/2) : 1 : round(nx/2) + round(nx_dam/2));

    matrix = zeros(nx,ny);
    DEM = matrix;
    Water_Depth = matrix;

    % Elevations of the Dam
    for i = 1:ny
        if i <= ny_dam
            DEM(:,i) = Dam_Bottom_Elevation; % Inside the dam
            Water_Depth(idx,i) = Dam_Height;
        else
            DEM(:,i) = DEM(:,i-1) - Slope*Resolution; % Outside of the dam
        end
    end

    xlimits = [0 ny*Resolution];
    ylimits = [0 nx*Resolution];
    dx =  Resolution;
    dy = Resolution;
    refmat = [0  dy;  dx  0;  xlimits(1) - dx/2   ylimits(2) + dy/2];

    R = refmatToMapRasterReference(refmat, size(DEM));
            
%     % Soil and LULC Maps
    DEM_raster.Z = DEM;
    DEM_raster.georef.SpatialRef.CellExtentInWorldX = Resolution;
    DEM_raster.georef.SpatialRef.CellExtentInWorldY = Resolution;
    DEM_raster.cellsize = Resolution;
    DEM_raster.refmat = refmat;
    DEM_raster.size = size(DEM);
    DEM_raster.georef.SpatialRef = R;
    DEM_raster.georef.ProjectedCRS = projcrs(3857);
%     
    SOIL_raster = DEM_raster;
    LULC_raster = DEM_raster;
    Initial_Water_Depth = DEM_raster;
    SOIL_raster.Z = ones(nx,ny);
    LULC_raster.Z = ones(nx,ny);
    Initial_Water_Depth.Z = Water_Depth;

    % Export_Rasters
    GRIDobj2geotiff(DEM_raster,'DEM_Dam_Break.tif');
    GRIDobj2geotiff(SOIL_raster,'SOIL_Dam_Break.tif');
    GRIDobj2geotiff(LULC_raster,'LULC_Dam_Break.tif');
    GRIDobj2geotiff(Initial_Water_Depth,'Initial_Water_Depth_Dam_Break.tif');  
    

end