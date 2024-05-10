% HydroPol2D Input Maps
% Developer: Marcus Nobrega
% Date 6/21/2023
% Goal: Plot Initial Maps

x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,2)]; y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,1)];
filename = 'Input_Maps';
set(gcf,'units','inches','position',[2,0,10,8])

%% Creating the custom basemap
web_mercator_crs = projcrs(3857);
no_plot=0;
try
if DEM_raster.georef.SpatialRef.ProjectedCRS.Name == web_mercator_crs.Name
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
    end
en
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
else
    no_plot = 1;


end
end

% -------- Manning  -------- %
t_title = 'Manning';
ax1 = subplot(3,2,1);
if no_plot==0;
    try
        ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax1 = mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
    end
end
axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
z = LULC_Properties.roughness; z(idx_nan) = nan;
idx = z < 0;
z(idx) = nan;
idx = isinf(z);
z(idx) = nan;
xmax = size(z,2);
xend = xmax;
ymax = size(z,1);
yend = ymax;
h_min = min(min(z));
F = z;
zmax = max(max(z(~isnan(z))));
if isempty(zmax) || isinf(zmax) || zmax == 0
    zmax = 0.1;
end
map = surf(x_grid,y_grid,F);
set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
title((t_title),'Interpreter','Latex','FontSize',12)
view(0,90)
if h_min == zmax
    zmax = 2*h_min;
end
caxis([h_min zmax]);
colormap(jet)
hold on
k = colorbar ;
ylabel(k,'$n$ ($\mathrm{s.m^{-1/3}}$)','Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax1, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

zlabel ('$n$ ($\mathrm{sm^{-1/3}}$)','Interpreter','Latex','FontSize',12)
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')
% ---------- h_0 --------------- %

ax2 = subplot(3,2,2);
if no_plot==0;
    try
        ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax1 = mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
    end
end
t_title = 'Initial Abstraction';
axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
z = LULC_Properties.h_0; z(idx_nan) = nan;
idx = z < 0;
z(idx) = nan;
idx = isinf(z);
z(idx) = nan;
xmax = size(z,2);
xend = xmax;
ymax = size(z,1);
yend = ymax;
h_min = min(min(z));
F = z;
zmax = max(max(z(~isnan(z))));
if isempty(zmax) || isinf(zmax) || zmax == 0
    zmax = 0.1;
end
map = surf(x_grid,y_grid,F);
set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
title((t_title),'Interpreter','Latex','FontSize',12)
view(0,90)
if h_min == zmax
    zmax = 2*h_min;
end
caxis([h_min zmax]);
colormap(jet)
hold on
k = colorbar ;
ylabel(k,'$h_0$ ($\mathrm{mm})$','Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax2, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

zlabel ('$h_0$ ($\mathrm{mm}$)','Interpreter','Latex','FontSize',12)
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')

% ----------  k_sat ------------- %
ax3 = subplot(3,2,3);
if no_plot==0;
    try
        ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax1 = mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
    end
end
t_title = 'Sat. Hyd. Conductivity';
axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
z = Soil_Properties.ksat; z(idx_nan) = nan;
idx = z < 0;
z(idx) = nan;
idx = isinf(z);
z(idx) = nan;
xmax = size(z,2);
xend = xmax;
ymax = size(z,1);
yend = ymax;
h_min = min(min(z));
F = z;
zmax = max(max(z(~isnan(z))));
if isempty(zmax) || isinf(zmax) || zmax == 0
    zmax = 0.1;
end
map = surf(x_grid,y_grid,F);
set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
title((t_title),'Interpreter','Latex','FontSize',12)
view(0,90)
if h_min == zmax
    zmax = 2*h_min;
end
caxis([h_min zmax]);
colormap(jet)
hold on
k = colorbar ;
ylabel(k,'$k_{sat}$ ($\mathrm{mm/h})$','Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax3, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

zlabel ('$k_{sat}$ ($\mathrm{mm/h}$)','Interpreter','Latex','FontSize',12)
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')

% ----------  dtheta ------------- %
ax4 = subplot(3,2,4);
if no_plot==0;
    try
        ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax1 = mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
    end
end
t_title = 'Moisture Deficit';
axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
z = (Soil_Properties.teta_sat - Soil_Properties.teta_i); z(idx_nan) = nan;
idx = z < 0;
z(idx) = nan;
idx = isinf(z);
z(idx) = nan;
xmax = size(z,2);
xend = xmax;
ymax = size(z,1);
yend = ymax;
h_min = min(min(z));
F = z;
zmax = max(max(z(~isnan(z))));
if isempty(zmax) || isinf(zmax) || zmax == 0
    zmax = 0.1;
end
map = surf(x_grid,y_grid,F);
set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
title((t_title),'Interpreter','Latex','FontSize',12)
view(0,90)
if h_min == zmax
    zmax = 2*h_min;
end
caxis([h_min zmax]);
colormap(jet)
hold on
k = colorbar ;
ylabel(k,'$\Delta \theta$ ($\mathrm{cm^3.cm^{-3}})$','Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax4, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

zlabel ('$\Delta \theta$ ($\mathrm{cm^3.cm^{-3}}$)','Interpreter','Latex','FontSize',12)
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')

% ----------  F_0 ------------- %
ax5 = subplot(3,2,5);
if no_plot==0;
    try
        ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax1 = mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
    end
end
t_title = 'Initial Soil Content';
z = Soil_Properties.I_0; z(idx_nan) = nan;
idx = z < 0;
z(idx) = nan;
idx = isinf(z);
z(idx) = nan;
xmax = size(z,2);
xend = xmax;
ymax = size(z,1);
yend = ymax;
h_min = min(min(z));
F = z;
zmax = max(max(z(~isnan(z))));
if isempty(zmax) || isinf(zmax) || zmax == 0
    zmax = 0.1;
end
map = surf(x_grid,y_grid,F);
set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
title((t_title),'Interpreter','Latex','FontSize',12)
view(0,90)
if h_min == zmax
    zmax = 2*h_min;
end
caxis([h_min zmax]);
colormap(jet)
hold on
k = colorbar ;
ylabel(k,'$F_0$ ($\mathrm{m}$)','Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax5, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

zlabel ('$I_0$ ($\mathrm{m}$)','Interpreter','Latex','FontSize',12)
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')

% ----------  D_0 ------------- %
ax6 = subplot(3,2,6);
if no_plot==0;
    try
        ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax1 = mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
    end
end
t_title = 'Initial Water Depth';
axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
z = depths.d_0/1000; z(idx_nan) = nan;
idx = z < 0;
z(idx) = nan;
idx = isinf(z);
z(idx) = nan;
xmax = size(z,2);
xend = xmax;
ymax = size(z,1);
yend = ymax;
h_min = min(min(z));
F = z;
zmax = max(max(z(~isnan(z))));
if isempty(zmax) || isinf(zmax) || zmax == 0
    zmax = 0.1;
end
map = surf(x_grid,y_grid,F);
set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
title((t_title),'Interpreter','Latex','FontSize',12)
view(0,90)
if h_min == zmax
    zmax = 2*h_min;
end
caxis([h_min zmax]);
colormap(jet)
hold on
k = colorbar ;
ylabel(k,'$d_{0}$ ($\mathrm{m}$)','Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax6, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

zlabel ('$d_{0}^0$ ($\mathrm{m}$)','Interpreter','Latex','FontSize',12)
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')

exportgraphics(gcf,fullfile(strcat('Modeling_Results','\','Input_Maps.TIF')),'ContentType','image','Colorspace','rgb','Resolution',600)


close all