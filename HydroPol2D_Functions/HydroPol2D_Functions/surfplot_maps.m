% Map Scaled Plot
function [axis] = surfplot_maps(DEM_raster,matrix,coloramp,x_label,y_label,z_label,no_data_value,idx_nan,row,col,index,size_font)
    hold on;
    
  %% Creating the custom basemap
basemapName = "openstreetmap";
url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)
% Getting lat and lon from the study area

[lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
latlim = [lat(1) lat(2)];
lonlim = [lon(1) lon(2)];
% Retriving the basemap image
try 
    [A,RA,attribA] = readBasemapImage(basemapName,latlim,lonlim);
catch ME
    warning('You need matlab 2022a or higher to use basemaps in georeferenced plots.')
end

%% Creating a Shapefile from the DEM as reference

% Creating a binary mask
binaryMask = ~isnan(DEM_raster.Z);
boundaries = bwboundaries(binaryMask);
% Pre-allocate arrays to store combined X and Y coordinates
combinedX = [];
combinedY = [];

    % Combine all boundary coordinates into a single array
for k = 1:numel(boundaries)
    boundary = boundaries{k};
    X = boundary(:, 2);
    Y = boundary(:, 1);
    combinedX = [combinedX; X; NaN]; % Add NaN to separate polygons
    combinedY = [combinedY; Y; NaN]; % Add NaN to separate polygons
end
% Remove the trailing NaNs at the end (optional)
combinedX = combinedX(1:end-1);
combinedY = combinedY(1:end-1);
% making the geostruct to alocate the shapefile
% cheking if CRS of the project is on WGS84
web_mercator_crs = projcrs(3857);
no_plot=0;
if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= web_mercator_crs.Name;
    no_plot = 1;
else
   S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
   S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)]; % Calculate bounding box for each polygon
   S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
   S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
end

    %% Matrix
    strcat('ax',num2str(index)) = subplot(row,col,index);
    % set = strcat('ax',num2str(index));
    z = matrix; % Infiltration Map
    z(z<0)=nan; % Taking away negative values
    zmax = max(max(max(z))); zmin = min(min(min(z))); z(idx_nan) = no_data_value;
    % Using DEM as a reference
    F = DEM_raster; % Just to get a coordinate system
    F.Z = z;
    % Plotting Hillshade + Data
    if no_plot==0;
        try
            mapshow(A,RA,"AlphaData",0.45);hold on;
            mapshow(S_p,'FaceColor','n'); hold on;
        catch ME
        end
    end
    imagesc(hillshade(DEM_raster)); hold on;
    colormap(gca,'gray')
    surf(F);
    colormap(gca,coloramp)
    view(0,90);
    colorbar
    % imageschs(DEM_raster,F,'colormap',coloramp);

  
    if zmin == zmax
        zmax = zmin*2;
    end
    if zmin == 0 && zmin == zmax
        zmax = zmin + 10;
    end
    caxis([zmin zmax]);
    k = colorbar;
    k.TickDirection = 'out';
    ylabel(k,z_label,'Interpreter','Latex','FontSize',size_font)
    xlabel(x_label,'Interpreter','Latex','FontSize',size_font)
    ylabel(y_label,'Interpreter','Latex','FontSize',size_font) 
    zlabel (z_label,'Interpreter','Latex','FontSize',size_font)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    title(z_label,'Interpreter','Latex');
    grid off
    box on ;
    axis = gca;
end