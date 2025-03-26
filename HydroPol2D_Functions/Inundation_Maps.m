%%% Inundation Maps %%%
% Creates .GIF files from the saved data
% Last updated: 10/4/2021

close all
f = 1;
t_max = running_control.routing_time;
tfinal = t_max ; %t_max;
DEM_maps = gather(Elevation_Properties.elevation_cell);

%% Time Data Processing
% The following block processes time-related data for different components 
% of the hydrological-hydrodynamic model, such as spatial rainfall, 
% evapotranspiration (ETP), and time records. Time values are adjusted 
% if they are not already in datetime format.

% Check if elapsed time flag is not set
if flags.flag_elapsed_time ~= 1
    % If spatial rainfall is enabled, process the spatial rainfall duration
    if flags.flag_spatial_rainfall == 1
        % Convert rainfall spatial duration to datetime if not already in datetime format
        if ~isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration)
            % Convert from minutes to days and add start date
            Spatial_Rainfall_Parameters.rainfall_spatial_duration = ...
                Spatial_Rainfall_Parameters.rainfall_spatial_duration / 60 / 24 + date_begin;
        end
    end

    % If ETP (Evapotranspiration) is enabled, process the climatologic duration
    if flags.flag_ETP == 1
        % Convert climatologic spatial duration to datetime if not already in datetime format
        if ~isdatetime(ETP_Parameters.climatologic_spatial_duration)
            % Convert from minutes to days and add start date
            ETP_Parameters.climatologic_spatial_duration = ...
                ETP_Parameters.climatologic_spatial_duration / 60 / 24 + date_begin;
        end
    end

    % Convert time records to datetime format if not already
    if ~isdatetime(running_control.time_records)
        % Convert from minutes to days and add start date
        running_control.time_records = ...
            double(running_control.time_records / 60 / 24) + date_begin;
    end
end

%% Creating the Custom Basemap
% This section sets up a custom basemap using OpenStreetMap tiles and
% specifies the attribution for usage in the plot.

% Define basemap name and URL template for OpenStreetMap
basemapName = "openstreetmap"; % Custom name for the basemap
url = "c.tile.openstreetmap.org/${z}/${x}/${y}.png"; % URL template for fetching tiles
url2 = 'a'; % Placeholder, potentially unused
copyright = char(uint8(169)); % Copyright symbol
attribution = copyright + "OpenStreetMap contributors"; % Full attribution string
attribution_2 = "Stadia Maps @ OpenStreetMap contributors"; % Additional attribution info

% Add the custom basemap to the current figure
addCustomBasemap(basemapName, url, "Attribution", attribution); 
% The function 'addCustomBasemap' integrates the basemap into the map with the provided attribution

%% Getting Latitude and Longitude from the Study Area
% This block extracts the latitude and longitude limits from the georeferencing
% information of the DEM (Digital Elevation Model) raster. These limits define
% the geographical extent of the study area based on the projected coordinate system.

% Convert the projected coordinate system to geographic coordinates
% The 'projinv' function transforms the X and Y world limits of the DEM raster
% from the projected coordinate system into latitude and longitude (geographic coordinates)
[lat, lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS, ...
                     DEM_raster.georef.SpatialRef.XWorldLimits, ...
                     DEM_raster.georef.SpatialRef.YWorldLimits);

% Define latitude and longitude limits
latlim = [lat(1), lat(2)]; % Latitude range of the study area (min to max latitude)
lonlim = [lon(1), lon(2)]; % Longitude range of the study area (min to max longitude)

% These limits can now be used for plotting or other geospatial analysis


%% Retrieving the Basemap Image
% This block attempts to fetch the basemap image based on the latitude and longitude limits.

try
    % Read the basemap image for the specified geographic limits
    [A, RA, attribA] = readBasemapImage(basemapName, latlim, lonlim); 
    % A: The image array of the basemap
    % RA: The referencing matrix for the basemap
    % attribA: The attribution information for the basemap image
catch ME
    % If MATLAB version is below 2022a, warn the user about basemap support
    warning('You need MATLAB 2022a or higher to use basemaps in georeferenced plots.');
end

%% Creating a Shapefile from the DEM as Reference
% This block generates a shapefile representing the boundaries of valid 
% (non-NaN) data in the DEM raster. It first creates a binary mask, 
% then extracts the boundary coordinates, and finally creates a shapefile 
% in the appropriate coordinate reference system (CRS).

% Create a binary mask where non-NaN values in the DEM are set to 1
binaryMask = ~isnan(DEM_raster.Z); 
% Extract the boundaries of the valid (non-NaN) data
boundaries = bwboundaries(binaryMask); 

% Pre-allocate arrays to store combined X and Y coordinates
combinedX = [];
combinedY = [];

% Combine all boundary coordinates into a single array
% Each boundary (polygon) is stored in 'boundaries' as a cell array
for k = 1:numel(boundaries)
    boundary = boundaries{k};
    X = boundary(:, 2); % X-coordinates of the boundary
    Y = boundary(:, 1); % Y-coordinates of the boundary
    combinedX = [combinedX; X; NaN]; % Add NaN to separate different polygons
    combinedY = [combinedY; Y; NaN]; % Add NaN to separate different polygons
end

% Remove the trailing NaNs at the end (optional)
combinedX = combinedX(1:end-1);
combinedY = combinedY(1:end-1);

% Create a geospatial structure for storing the shapefile data
% Check if the CRS of the DEM raster is Web Mercator (EPSG:3857)
web_mercator_crs = projcrs(3857);
no_plot = 0; % Flag to check if plotting is required
if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= web_mercator_crs.Name
    % If CRS is not Web Mercator, skip plotting (no_plot = 1)
    no_plot = 1;
else
    % Prepare the structure for the shapefile data
    S_p = struct('Geometry', 'Polygon', ...
                 'BoundingBox', [], ...
                 'X', [], ...
                 'Y', [], ...
                 'fid', 1, 'DN', 0);  % Initialize the structure fields
             
    % Calculate the bounding box of the polygon (study area limits)
    S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), ...
                       DEM_raster.georef.SpatialRef.YWorldLimits(1,1); ...
                       DEM_raster.georef.SpatialRef.XWorldLimits(1,2), ...
                       DEM_raster.georef.SpatialRef.YWorldLimits(1,2)];

    % Adjust the X and Y coordinates based on the DEM's georeferencing
    % Convert the coordinates from raster grid to world coordinates
    S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + ...
             combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - ...
             DEM_raster.georef.SpatialRef.CellExtentInWorldX / 2)';
    S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - ...
             combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + ...
             DEM_raster.georef.SpatialRef.CellExtentInWorldX / 2)';
end

%% Set up HEC-RAS colors
hec_ras_colors = [52/255 85/255 132/255; 0 1 1; 0 128/255 1; 0 255/255 0; 1 1 0; 1 128/255 0; 1 0 0; 128/255 0 128/255];

%% Plot Elevation Model

% Set up figure for plotting
h = figure;
axis tight manual  % Ensures consistent size for getframe()
FileName_String = 'Elevation_Model.gif'; % Output file name for GIF
FileName = fullfile(folderName, FileName_String); % Full path to save the GIF

% Grid resolution and time steps
grid_resolution = Wshed_Properties.Resolution; % Grid resolution
tmax = 10;  % Number of frames for the animation (adjust as needed)

% Loop through each time step and generate frames
for t = 1:tmax
    % Extract DEM data and define grid dimensions
    z = DEM_maps; % DEM data
    xmax = length(z(1,:));
    ymax = length(z(:,1));
    
    % Define UTM coordinates
    x_grid = GIS_data.xulcorner + grid_resolution * (1:xmax);
    y_grid = GIS_data.yulcorner - grid_resolution * (1:ymax);
    
    % Set values <= 0 in the DEM to infinity (no data)
    z(z <= 0) = inf;
    
    % Define plot range based on DEM values
    h_min = min(z(~isinf(z))); % Minimum elevation value
    zmax = max(z(~isinf(z))); % Maximum elevation value
    zmax = round(zmax / 10, 0) * 10 * 1.05; % Adjust to nearest 10 and scale up

    % Extract the DEM data to be plotted
    F = DEM_maps(1:ymax, 1:xmax); 
    
    % Create the 3D surface plot
    kk = surf(x_grid, y_grid, F);
    set(kk, 'LineStyle', 'none'); % Remove grid lines
    set(gca, 'XTickLabel', x_grid, 'YTickLabel', y_grid);
    
    % Adjust axis limits based on DEM range
    axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid) h_min zmax]);
    
    % Set the viewing angle for rotation
    view_angle = (t - 1) * 360 / tmax; % Smooth rotation effect
    view(view_angle,(t-1)*90/tmax); % Adjust tilt

    % Color settings
    colorbar; 
    caxis([h_min, zmax]); % Set color scale for elevation
    colormap(Terrain_RAS); % Set terrain colormap
    
    % Plot aesthetics
    box on;
    title('Elevation', 'Interpreter', 'Latex', 'FontSize', 12);
    
    % Axis labels
    k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out'; % Colorbar for elevation
    ylabel(k, 'Elevation [m]', 'Interpreter', 'Latex', 'FontSize', 12);
    xlabel('Easting [m]', 'Interpreter', 'Latex', 'FontSize', 12);
    ylabel('Northing [m]', 'Interpreter', 'Latex', 'FontSize', 12);
    zlabel('Elevation [m]', 'Interpreter', 'Latex', 'FontSize', 12);
    
    % Set font properties
    set(gca, 'FontName', 'Garamond', 'FontSize', 12);
    
    % Increase border (axis) thickness
    ax = gca; 
    ax.LineWidth = 2; % Make axis lines thicker
    
    % Update the plot
    drawnow;

    % Capture the plot as an image frame for GIF
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256); % Convert to indexed image for GIF
    
    % Write to GIF file, append after the first frame
    if t == 1
        imwrite(imind, cm, FileName, 'gif', 'Loopcount', inf);
    else
        imwrite(imind, cm, FileName, 'gif', 'WriteMode', 'append');
    end
end

% Close the figure after GIF creation
clf; % Clear figure window

%% Plot Water Surface Elevation and Depths
% Adjusting the size

% Generate video showing water level profile over time
close all

Video_Name = 'WSE_Depths.mp4';

% Set up video
video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
% control the framerate
video.FrameRate = 10;
open(video);

% Set up HEC-RAS colors
hec_ras_colors = [52/255 85/255 132/255; 0 1 1; 0 128/255 1; 0 255/255 0; 1 1 0; 1 128/255 0; 1 0 0; 128/255 0 128/255];

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
FileName_String = 'WSE_and_Depths.gif';
FileName = fullfile(folderName,strcat('\',FileName_String));

set(gcf,'units','inches','position',[0,0,7,12])
set(gcf,'DefaultTextInterpreter','latex')


% createMatrixVideo(Maps.Hydro.d, running_control.time_records(1:5), 'test3', 'test2', 14, 'TEST', Terrain_RAS_ramp, 'Test', DEM_raster)

% This is absolutely memory expensive
% z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
% idx2 = Maps.Hydro.d < depths.depth_wse*1000;
% z1(z1<=0)=nan;
% z1(idx2) = nan;

% z1max = max(max(Max_depth.d/1000 + DEM_maps));
% z1min = min(min(Max_depth.d/1000 + DEM_maps));
% z1max = max(max(max(Maps.Hydro.d/1000 + DEM_maps)));
% z1min = min(min(min(Maps.Hydro.d/1000 + DEM_maps)));
z1max = max(max(Max_depth_d/1000 + DEM_maps)); % data comes from prost_processing
z1min = min(min(Max_depth_d/1000 + DEM_maps));

z2max = max(max(Max_depth_d/1000)); % data comes from prost_processing
z2min = min(min(Max_depth_d/1000));

% z2 = gather(Maps.Hydro.d/1000);
% z2(z2<=0)=nan;
% z2(idx2)=nan;

targetColor=[0.27 0.47 0.68];

store=1;
flag_loader=1;
for t = 1:f:length(running_control.time_records)
    clf
    ax1 = subplot(2,1,1);
    t_title = running_control.time_records(t);

    if t > saver_memory_maps*store
        store = store + 1;
        load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
        idx2 = Maps.Hydro.d < depths.depth_wse*1000;
        z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
        z2 = gather(Maps.Hydro.d)/1000;
        z1(z1<=0)= nan;
        z1(idx2) = nan;
        z2(z2<=0)= nan;
        z2(idx2) = nan;
    else
        if flag_loader == 1
          load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
          flag_loader=0;
          idx2 = Maps.Hydro.d < depths.depth_wse*1000;
          z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
          z2 = gather(Maps.Hydro.d)/1000;
          z1(z1<=0)= nan;
          z1(idx2) = nan;
          z2(z2<=0)= nan;
          z2(idx2) = nan;
        end
    end
    
    if isnan(z1max)
        zmax = 10;
    end
    if isnan(z1min)
        zmin = 0;
    end
    xmax = length(z(1,:));
    xend = xmax;
    ymax = length(z(:,1));
    yend = ymax;
    xbegin = 1;
    ybegin = 1;
    F = z1([ybegin:1:yend],[xbegin:1:xend],t - (store-1)*saver_memory_maps);
    F(idx2(:,:,t - (store-1)*saver_memory_maps)) = 0;
    F(F==0)=nan;
    % F = z1([ybegin:1:yend],[xbegin:1:xend],t);
    % F(idx2(:,:,t)) = 0;
    if no_plot == 0
        try
            mapshow(A,RA,"AlphaData",0.45);hold on;
        catch ME
            warning('You need matlab 2022 or higher to run basemaps.')
        end
    end
    surf(ax1,x_grid,y_grid,F);
    shading INTERP;

    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) z1min z1max])
    if flags.flag_elapsed_time == 1
        title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
    else
        title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
    end

    view(0,90);
    colorbar
    caxis([z1min z1max]);
    colormap(ax1,WSE_RAS)
    k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out'; k.TickDirection  = 'out';
    ylabel(k,'WSE [m]','Interpreter','Latex','FontSize',12)
    xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
    ylabel (' Northing [m] ','Interpreter','Latex','FontSize',12)
    zlabel ('WSE [m]','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    box on
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end

    ax2 = subplot(2,1,2);
    % zmax = max(max(max(z2(~isinf(z2)))));
    % zmin = min(min(min(z2)));
    % zmax = max(max(max(z2)));
    % zmin = min(min(min(z2)));
    if isnan(z2max)
        z2max = 10;
    end
    if isnan(z2min)
        z2min = 0;
    end
    F = z2([ybegin:1:yend],[xbegin:1:xend],t-(store-1)*saver_memory_maps);
    F(idx2(:,:,t-(store-1)*saver_memory_maps)) = 0;
    F(F==0)=nan;
    % F = z2([ybegin:1:yend],[xbegin:1:xend],t);
    % F(idx2(:,:,t)) = 0;
    % F(F==0)=nan;
    if no_plot==0
        try
            mapshow(A,RA,"AlphaData",0.25);hold on;
        catch ME
            warning('You need matlab 2022 or higher to run basemaps.')
        end
    end
    surf(ax2,x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) z2min z2max])

    shading INTERP;
    if flags.flag_elapsed_time == 1
        title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
    else
        title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
    end
    view(0,90);
    colorbar
    caxis([0 z2max]);
    colormap(ax2,Depth_Purple)
    k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';

    ylabel(k,'Depths [m]','Interpreter','Latex','FontSize',12)
    xlabel('Easting [m] ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
    zlabel ('WSE [m]','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end
    box on
    % Set background color and write to video
    frame = getframe(gcf);
    writeVideo(video,frame);
    hold off
    clf
end

% Close video writer    
close(video);
close all

%% Plot Water Surface Elevation and Depths
% Adjusting the size

% Generate video showing water level profile over time
close all

Video_Name = 'Depths.mp4';

% Set up video
video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
% control the framerate
video.FrameRate = 10;
open(video);

% Set up HEC-RAS colors
hec_ras_colors = [52/255 85/255 132/255; 0 1 1; 0 128/255 1; 0 255/255 0; 1 1 0; 1 128/255 0; 1 0 0; 128/255 0 128/255];

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
FileName_String = 'WSE_and_Depths.gif';
FileName = fullfile(folderName,strcat('\',FileName_String));

set(gcf,'units','inches','position',[0,0,7,12])
set(gcf,'DefaultTextInterpreter','latex')


% createMatrixVideo(Maps.Hydro.d, running_control.time_records(1:5), 'test3', 'test2', 14, 'TEST', Terrain_RAS_ramp, 'Test', DEM_raster)

% This is absolutely memory expensive
% z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
% idx2 = Maps.Hydro.d < depths.depth_wse*1000;
% z1(z1<=0)=nan;
% z1(idx2) = nan;

% z1max = max(max(Max_depth.d/1000 + DEM_maps));
% z1min = min(min(Max_depth.d/1000 + DEM_maps));
% z1max = max(max(max(Maps.Hydro.d/1000 + DEM_maps)));
% z1min = min(min(min(Maps.Hydro.d/1000 + DEM_maps)));
z1max = max(max(Max_depth_d/1000 + DEM_maps)); % data comes from prost_processing
z1min = min(min(Max_depth_d/1000 + DEM_maps));

z2max = max(max(Max_depth_d/1000)); % data comes from prost_processing
z2min = min(min(Max_depth_d/1000));

% z2 = gather(Maps.Hydro.d/1000);
% z2(z2<=0)=nan;
% z2(idx2)=nan;

targetColor=[0.27 0.47 0.68];

store=1;
flag_loader=1;
for t = 1:f:length(running_control.time_records)
    clf
    t_title = running_control.time_records(t);

    if t > saver_memory_maps*store
        store = store + 1;
        load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
        idx2 = Maps.Hydro.d < depths.depth_wse*1000;
        z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
        z2 = gather(Maps.Hydro.d)/1000;
        z1(z1<=0)= nan;
        z1(idx2) = nan;
        z2(z2<=0)= nan;
        z2(idx2) = nan;
    else
        if flag_loader == 1
          load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
          flag_loader=0;
          idx2 = Maps.Hydro.d < depths.depth_wse*1000;
          z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
          z2 = gather(Maps.Hydro.d)/1000;
          z1(z1<=0)= nan;
          z1(idx2) = nan;
          z2(z2<=0)= nan;
          z2(idx2) = nan;
        end
    end
    
    if isnan(z1max)
        zmax = 10;
    end
    if isnan(z1min)
        zmin = 0;
    end
    xmax = length(z(1,:));
    xend = xmax;
    ymax = length(z(:,1));
    yend = ymax;
    xbegin = 1;
    ybegin = 1;
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) z1min z1max])
    if flags.flag_elapsed_time == 1
        title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
    else
        title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
    end

    if isnan(z2max)
        z2max = 10;
    end
    if isnan(z2min)
        z2min = 0;
    end
    F = z2([ybegin:1:yend],[xbegin:1:xend],t-(store-1)*saver_memory_maps);
    F(idx2(:,:,t-(store-1)*saver_memory_maps)) = 0;
    F(F==0)=nan;
    % F = z2([ybegin:1:yend],[xbegin:1:xend],t);
    % F(idx2(:,:,t)) = 0;
    % F(F==0)=nan;
    if no_plot==0
        try
            mapshow(A,RA,"AlphaData",0.25);hold on;
        catch ME
            warning('You need matlab 2022 or higher to run basemaps.')
        end
    end
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) z2min z2max])

    shading INTERP;
    if flags.flag_elapsed_time == 1
        title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
    else
        title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
    end
    view(0,90);
    colorbar
    caxis([0 z2max]);
    colormap(Depth_Purple)
    k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';

    ylabel(k,'Depths [m]','Interpreter','Latex','FontSize',12)
    xlabel('Easting [m] ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
    zlabel ('WSE [m]','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end
    box on
    % Set background color and write to video
    frame = getframe(gcf);
    writeVideo(video,frame);
    hold off
    clf
end

% Close video writer    
close(video);
close all
%% WSE GIFS
% %% Plot Water Surface Elevation and Depths
% % Adjusting the size
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% FileName_String = 'WSE_and_Depths.gif';
% FileName = fullfile(folderName,strcat('\',FileName_String));
%
% set(gcf,'units','inches','position',[0,0,7,12])
% idx2 = Maps.Hydro.d < depths.depth_wse*1000;
% % This is absolutely memory expensive
% z1 = gather(Maps.Hydro.d)/1000 + DEM_maps;
% z1(z1<=0)=nan;
% z1(idx2) = nan;
%
% zmax = max(max(max(Maps.Hydro.d + DEM_maps)));
% zmin = min(min(min(Maps.Hydro.d + DEM_maps)));
%
% z2 = gather(Maps.Hydro.d/1000);
% z2(z2<=0)=nan;
% z2(idx2)=nan;
%
% targetColor=[0.27 0.47 0.68];
%
% for t = 1:f:length(running_control.time_records)
%     % clf
%     subplot(2,1,1)
%     t_title = running_control.time_records(t);
%
%     if isnan(zmax)
%         zmax = 10;
%     end
%     if isnan(zmin)
%         zmin = 0;
%     end
%     F = z1([ybegin:1:yend],[xbegin:1:xend],t);
%     F(idx2(:,:,t)) = 0;
%     if no_plot==0;
%         mapshow(A,RA,"AlphaData",0.45);hold on;
%         mapshow(S_p,'FaceColor','n'); hold on;
%     end
%     surf(x_grid,y_grid,F);
%     shading INTERP;
%
%     axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
%     if flags.flag_elapsed_time == 1
%         title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
%     else
%         title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
%     end
%
%     view(0,90);
%     colorbar
%     caxis([zmin zmax]);
%     colormap(depth_ramp)
%     k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
%     ylabel(k,'WSE [m]','Interpreter','Latex','FontSize',12)
%     xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
%     ylabel (' Northing [m] ','Interpreter','Latex','FontSize',12)
%     zlabel ('WSE [m]','Interpreter','Latex','FontSize',12)
%     set(gca,'FontName','Garamond')
%     ax = ancestor(gca, 'axes');
%     ax.XAxis.Exponent = 0;xtickformat('%.0f');
%     ax.YAxis.Exponent = 0;ytickformat('%.0f');
%
%     box on
%     drawnow
%
%     subplot(2,1,2)
%
%     % zmax = max(max(max(z2(~isinf(z2)))));
%     % zmin = min(min(min(z2)));
%     zmax = max(max(max(z2)));
%     zmin = min(min(min(z2)));
%     if isnan(zmax)
%         zmax = 10;
%     end
%     if isnan(zmin)
%         zmin = 0;
%     end
%     F = z2([ybegin:1:yend],[xbegin:1:xend],t);
%     F(idx2(:,:,t)) = 0;
%     F(F==0)=nan;
%     if no_plot==0;
%         mapshow(A,RA,"AlphaData",0.45);hold on;
%         mapshow(S_p,'FaceColor','n'); hold on;
%     end
%     surf(x_grid,y_grid,F);
%     axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
%
%     shading INTERP;
%     if flags.flag_elapsed_time == 1
%         title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
%     else
%         title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
%     end
%     view(0,90);
%     colorbar
%     caxis([0 zmax]);
%     colormap(depth_ramp)
%     k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
%
%     ylabel(k,'Depths [m]','Interpreter','Latex','FontSize',12)
%     xlabel('Easting [m] ','Interpreter','Latex','FontSize',12)
%     ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
%     zlabel ('WSE [m]','Interpreter','Latex','FontSize',12)
%     set(gca,'FontName','Garamond')
%     ax = ancestor(gca, 'axes');
%     ax.XAxis.Exponent = 0;xtickformat('%.0f');
%     ax.YAxis.Exponent = 0;ytickformat('%.0f');
%
%     box on
%     grid on
%     % Capture the plot as an image
%     frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if t == 1
%         imwrite(imind,cm,FileName,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,FileName,'gif','WriteMode','append');
%     end
%     clf
% end
% close all

%% Spatial Rainfall Isoietal
% time_step_rainfall = double(time_step_rainfall);
if flags.flag_spatial_rainfall == 1
    close all
    set(gcf,'units','inches','position',[2,2,6.5,4])
    if isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration(2))
        time_step_rainfall = hours(Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) -Spatial_Rainfall_Parameters.rainfall_spatial_duration(1)); % hours
    else
        time_step_rainfall = (Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) -Spatial_Rainfall_Parameters.rainfall_spatial_duration(1))/60; % hours
    end
    rain_total = rainfall_sum*time_step_rainfall; % data comes from post_processing
    time_total = days(date_end - date_begin);
    title_isoietal = strcat('Cumulative rainfall of',{' '},string(round(time_total,2)),' days ',' from ',{' '},cellstr(date_begin),' to ',{' '},cellstr(date_end));
    idx = isnan(Elevation_Properties.elevation_cell);
    rain_total(rain_total<=0) = nan;
    rain_total(idx) = nan;
    zmin = min(min(rain_total));
    zmax = max(max(rain_total));
    F = rain_total([ybegin:1:yend],[xbegin:1:xend]);
    if no_plot==0
        try
            mapshow(A,RA,"AlphaData",0.45);hold on;
            mapshow(S_p,'FaceColor','n'); hold on;
        catch ME

        end
    end
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
    shading INTERP;
    title(title_isoietal,'Interpreter','Latex','FontSize',12);
    colorbar
    caxis([zmin zmax]);
    colormap(Spectrum)
    k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
    ylabel(k,'Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
    xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)

    zlabel ('Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

    box on
    set(gca,'tickdir','out');
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    view(0,90)
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end
    exportgraphics(gcf,fullfile(folderName,'Isoietal_Map_Rainfall.png'),'ContentType','image','Colorspace','rgb','Resolution',1200)
    saveas(gcf,fullfile(folderName,'Isoeital_Rainfall.fig'))
    close all
end


%% Spatial ETR Isoietal
% time_step_rainfall = double(time_step_rainfall);
if flags.flag_ETP == 1
    close all
    set(gcf,'units','inches','position',[2,2,6.5,4])
    if isdatetime(ETP_Parameters.climatologic_spatial_duration(2))
        time_step_ETR = days(ETP_Parameters.climatologic_spatial_duration(2) -ETP_Parameters.climatologic_spatial_duration(1)); % hours
    else
        time_step_ETR = (ETP_Parameters.climatologic_spatial_duration(2) - ETP_Parameters.climatologic_spatial_duration(1))/60; % hours
    end
    ETR_Total = ETR_sum*time_step_ETR; % data comes from post_processing
    time_total = days(date_end - date_begin);
    title_isoietal = strcat('ETR of',{' '},string(round(time_total,2)),' days ',' from ',{' '},cellstr(date_begin),' to ',{' '},cellstr(date_end));
    idx = isnan(Elevation_Properties.elevation_cell);
    ETR_Total(ETR_Total<=0) = nan;
    ETR_Total(idx) = nan;
    zmin = min(min(ETR_Total));
    zmax = max(max(ETR_Total));
    F = ETR_Total([ybegin:1:yend],[xbegin:1:xend]);
    if no_plot==0
        try
            mapshow(A,RA,"AlphaData",0.45);hold on;
            mapshow(S_p,'FaceColor','n'); hold on;
        catch ME

        end
    end
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
    shading INTERP;
    title(title_isoietal,'Interpreter','Latex','FontSize',12);
    colorbar
    caxis([zmin zmax]);
    colormap(Spectrum)
    k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
    ylabel(k,'ETR Volume (mm)','Interpreter','Latex','FontSize',12)
    xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)

    zlabel ('ETR Volume (mm)','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

    box on
    set(gca,'tickdir','out');
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca,'FontName','Garamond','FontSize',12)
    % Increase border (axis) thickness
    ax = gca;          % Get current axis
    ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    view(0,90)
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end
    exportgraphics(gcf,fullfile(folderName,'Isoietal_Map_ETR.png'),'ContentType','image','Colorspace','rgb','Resolution',1200)
    saveas(gcf,fullfile(folderName,'Isoeital_ETR.fig'))
    close all
end

%% Dividing the Rainfall into n intervals
% zero_matrix = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2));
% if flags.flag_ETP == 1
%     store=1;
%     flag_loader=1;
%     ETR_sum = zeros(size(zero_matrix));
%     for i = 1:length(running_control.time_records)
%         if i > saver_memory_maps*store
%             store = store + 1;
%             load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
%         else
%             if flag_loader == 1
%                 load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
%                 flag_loader=0;
%             end
%         end
%         z = Maps.Hydro.ETR_save(:,:,i - ((store-1)*saver_memory_maps));
%         z(isnan(z)) = 0; % Attention here
%         ETR_sum = ETR_sum + z; % mm/h
%         ETP_Parameters.std_dev_ETR(i,1) = nanstd(z(:));
%     end
% end



% if flags.flag_spatial_rainfall == 1
%     n_plots = 8;
%     n_rainfall = length(running_control.time_records);
%     time_step_rainfall_aggregated = floor((n_rainfall/n_plots))*time_step_rainfall;
%     zmin = 0;
%     zmax = max(max(rainfall_sum*time_step_rainfall));
% 
%     for i = 1:n_plots
%         subplot(2,n_plots/2,i)
%         prior_range = floor(n_rainfall/n_plots*(i-1) + 1);
%         range_rainfall = floor(n_rainfall/n_plots*i);
%         rain_total = sum(Maps.Hydro.spatial_rainfall_maps(:,:,prior_range:range_rainfall)*time_step_rainfall,3);
%         time_total = days(date_end - date_begin);
%         zmin = max(min(min(rain_total)),0);
%         zmax = max(max(max(rain_total)),0);
%         time_step_rainfall = double(time_step_rainfall);
%         if zmin == zmax
%             zmax = zmin + 10;
%         end
%         if i == 1
%             t_prior = date_begin; % Time
%         else
%             t_prior = date_begin + floor(n_rainfall/n_plots*((i-1) + 1))*time_step_rainfall/24; % Time
%         end
%         t_range = date_begin + floor(n_rainfall/n_plots*((i)))*time_step_rainfall/24; % Time;
%         title_isoietal = strcat('Interval from',{' '},string(t_prior),' to ',{' '},string(t_range));
%         idx = isnan(Elevation_Properties.elevation_cell);
%         rain_total(rain_total<=0) = nan;
%         rain_total(idx) = nan;
%         F = rain_total([ybegin:1:yend],[xbegin:1:xend]);
%         if no_plot==0
%             try 
%                 mapshow(A,RA,"AlphaData",0.45);hold on;
%             catch ME
% 
%             end
%         end
%         surf(x_grid,y_grid,F);
%         axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
%         shading INTERP;
%         title(title_isoietal,'Interpreter','Latex','FontSize',12);
%         colorbar
%         caxis([zmin zmax]);
%         colormap(linspecer)
%         k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
%         ylabel(k,'Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
%         xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
%         ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
% 
%         zlabel ('Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
%         if no_plot == 0
%             mapshow(S_p,'FaceColor','n'); hold on;
%         end
%         set(gca,'FontName','Garamond')
% 
%         box on
%         set(gca,'tickdir','out');
%         set(gca, 'TickLength', [0.02 0.01]);
%         set(gca,'Tickdir','out')
%         set(gca,'FontName','Garamond');
%         ax = ancestor(gca, 'axes');
%         ax.XAxis.Exponent = 0;xtickformat('%.0f');
%         ax.YAxis.Exponent = 0;ytickformat('%.0f');
%         hold off
%         view(0,90)
%     end
% end
% close all
%% Spatial Rainfall Video
% Generate video showing water level profile over time

try 
    delete('Modeling_Results\Rainfall_Intensities.mp4')
end

%Calculating general Max and Min
if flags.flag_spatial_rainfall == 1
    store = 1;
    flag_loader=1;
    for i = 1:length(running_control.time_records)
        try
        if i > saver_memory_maps*store
            store = store + 1;
            load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');      
            Max_rains = max(max(Maps.Hydro.spatial_rainfall_maps,[],3),Max_rains);
            Min_rains = min(min(Maps.Hydro.spatial_rainfall_maps,[],3),Min_rains);
            if flags.flag_ETP
                Max_etp = max(max(Maps.Hydro.ETP_save,[],3),Max_etp);
                Min_etp = min(min(Maps.Hydro.ETP_save,[],3),Min_etp);
            end
        else
            if flag_loader == 1
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                flag_loader=0;
                Max_rains = max(Maps.Hydro.spatial_rainfall_maps,[],3);
                Min_rains = min(Maps.Hydro.spatial_rainfall_maps,[],3);
                if flags.flag_ETP
                    Max_etp = max(Maps.Hydro.ETP_save,[],3);
                    Min_etp = min(Maps.Hydro.ETP_save,[],3);
                end
            end
        end
        end
    end
end

if flags.flag_spatial_rainfall == 1
    close all
    Video_Name = 'Rainfall_Intensities.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=2;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[0,0,7,12])
    set(gcf,'DefaultTextInterpreter','latex')
    if flags.flag_spatial_rainfall == 1 && flags.flag_rainfall == 1 && flags.flag_input_rainfall_map ~= 1
        close all
        zmax = max(max(Max_rains));
        zmin = min(min(Min_rains));
        if zmin == zmax
            zmax = zmin + 10; % mm/h
        end
        store = 1;
        flag_loader=1;
        for t = 1:1:length(running_control.time_records)
            clf
            if t > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                    flag_loader=0;
                end
            end
            rain = Maps.Hydro.spatial_rainfall_maps(:,:,t - (store-1)*saver_memory_maps);
            idx = isnan(Elevation_Properties.elevation_cell);
            rain(rain<=0) = nan;
            rain(idx) = nan;
            spatial_rainfall = rain;
            if t > size(Spatial_Rainfall_Parameters.rainfall_spatial_duration,2)
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(end) + Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2);
                rainfall = zeros(size(Elevation_Properties.elevation_cell));
                z = rainfall;
                plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'r.', 'MarkerSize', 30)
            else
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(t);
                % Draw plot for d = x.^n
                z = rain; %
                z(idx_nan)=nan;
                if flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
                    rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(t,1:Spatial_Rainfall_Parameters.n_raingauges)'; % Values of rainfall at t for each rain gauge
                    % idx_rainfall = logical(isnan(rainfall) | rainfall == 0);
                    idx_rainfall = logical(isnan(rainfall));
                    Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1); % Coordinates (easting) of each rain gauge
                    Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2); % Coordinates (northing) of each rain gauge
                    Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = []; % Taking out nans
                    Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = []; % Taking out nans
                    rainfall(idx_rainfall) = []; % Taking out nans
                    plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'r.', 'MarkerSize', 30)
                end
            end
            hold on
            F = z([ybegin:1:yend],[xbegin:1:xend]);
            if no_plot==0
                try
                    mapshow(A,RA,"AlphaData",0.45);hold on;
                    mapshow(S_p,'FaceColor','n'); hold on;
                catch ME

                end
            end
            surf(x_grid,y_grid,F);
            axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])

            shading INTERP;

            %
            if flags.flag_elapsed_time == 1
                title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
            else
                title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
            end
            view(0,90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)
            k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
            ylabel(k,'Rainfall [mm/h]','Interpreter','Latex','FontSize',12)
            xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
            ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)

            zlabel ('Rainfall [mm/h]','Interpreter','Latex','FontSize',12)
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

            box on
            set(gca,'tickdir','out');
            set(gca, 'TickLength', [0.02 0.01]);
            set(gca,'Tickdir','out')
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
            ax = ancestor(gca, 'axes');
            ax.XAxis.Exponent = 0;xtickformat('%.0f');
            ax.YAxis.Exponent = 0;ytickformat('%.0f');
            hold off
            if no_plot == 0
                mapshow(S_p,'FaceColor','n'); hold on;
            end
            box on
            % Set background color and write to video
            frame = getframe(gcf);
            writeVideo(video,frame);
            hold off
            clf
        end
    end

    % Close video writer
    close(video);
    close all
end
%% Spatial Rainfall Video
if flags.flag_spatial_rainfall == 1 && flags.flag_rainfall == 1 && flags.flag_input_rainfall_map ~= 1

    close all
    Video_Name = 'Rainfall_Maps.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=2;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[0,0,6.5,4])
    set(gcf,'DefaultTextInterpreter','latex')


    zmax = max(max(max(Maps.Hydro.spatial_rainfall_maps)));
    zmin = min(min(min(Maps.Hydro.spatial_rainfall_maps)));
    if zmin == zmax
        zmax = zmin + 10; % mm/h
    end

    for t = 1:(size(Maps.Hydro.spatial_rainfall_maps,3)-1)
        clf
        rain = Maps.Hydro.spatial_rainfall_maps(:,:,t);
        idx = isnan(Elevation_Properties.elevation_cell);
        rain(rain<=0) = nan;
        rain(idx) = nan;
        spatial_rainfall = rain;
        if t > size(Spatial_Rainfall_Parameters.rainfall_spatial_duration,2)
            t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(end) + Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2);
            rainfall = zeros(size(Elevation_Properties.elevation_cell));
            z = rainfall;
            plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'r.', 'MarkerSize', 30)
        else
            t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(t);
            % Draw plot for d = x.^n
            z = rain; %
            z(idx_nan)=nan;
            if flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
                rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(t,1:Spatial_Rainfall_Parameters.n_raingauges)'; % Values of rainfall at t for each rain gauge
                % idx_rainfall = logical(isnan(rainfall) | rainfall == 0);
                idx_rainfall = logical(isnan(rainfall));
                Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1); % Coordinates (easting) of each rain gauge
                Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2); % Coordinates (northing) of each rain gauge
                Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = []; % Taking out nans
                Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = []; % Taking out nans
                rainfall(idx_rainfall) = []; % Taking out nans
                plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'r.', 'MarkerSize', 30)
            end
        end
        hold on
        F = z([ybegin:1:yend],[xbegin:1:xend]);
        if no_plot==0
            try
                mapshow(A,RA,"AlphaData",0.45);hold on;
                mapshow(S_p,'FaceColor','n'); hold on;
            catch ME

            end
        end
        surf(x_grid,y_grid,F);
        axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])

        shading INTERP;

        %
        if flags.flag_elapsed_time == 1
            title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        colorbar
        caxis([zmin zmax]);
        colormap(Spectrum)
        k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
        ylabel(k,'Rainfall [mm/h]','Interpreter','Latex','FontSize',12)
        xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)

        zlabel ('Rainfall [mm/h]','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond','FontSize',12)
        % Increase border (axis) thickness
        ax = gca;          % Get current axis
        ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

        box on
        set(gca,'tickdir','out');
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca,'FontName','Garamond','FontSize',12)
        % Increase border (axis) thickness
        ax = gca;          % Get current axis
        ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
        ax = ancestor(gca, 'axes');
        ax.XAxis.Exponent = 0;xtickformat('%.0f');
        ax.YAxis.Exponent = 0;ytickformat('%.0f');
        if no_plot == 0
            mapshow(S_p,'FaceColor','n'); hold on;
        end
        box on
        hold off
        % Set background color and write to video
        frame = getframe(gcf);
        writeVideo(video,frame);
        hold off
        clf
    end
end
% Close video writer
close(video);
close all
%% Flag_Spatial_Rainfall - Input Maps
if flags.flag_spatial_rainfall == 1 && flags.flag_rainfall == 1 && flags.flag_input_rainfall_map == 1
    % Generate video showing water level profile over time
    close all
    Video_Name = 'Rainfall_Maps.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=2;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[0,0,6.5,4])
    set(gcf,'DefaultTextInterpreter','latex')

    zmax = max(max(max(Maps.Hydro.spatial_rainfall_maps)));
    zmin = min(min(min(Maps.Hydro.spatial_rainfall_maps)));
    if zmin == zmax
        zmax = zmin + 10; % mm/h
    end

    for t = 1:(size(Maps.Hydro.spatial_rainfall_maps,3)-1)
        clf
        rain = Maps.Hydro.spatial_rainfall_maps(:,:,t);
        idx = isnan(Elevation_Properties.elevation_cell);
        rain(rain<0) = nan;
        rain(idx) = nan;
        spatial_rainfall = rain;
        if t > length(Spatial_Rainfall_Parameters.rainfall_spatial_duration)
            t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(end) + Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2);
            rainfall = zeros(size(Elevation_Properties.elevation_cell));
            z = rainfall;
        else
            t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(t);
            % Draw plot for d = x.^n
            z = rain; %
            z(idx_nan)=nan;
        end
        hold on
        F = z([ybegin:1:yend],[xbegin:1:xend]);
        if no_plot==0
            try
                mapshow(A,RA,"AlphaData",0.45);hold on;
                mapshow(S_p,'FaceColor','n'); hold on;
            catch ME

            end
        end
        surf(x_grid,y_grid,F);
        axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])

        shading INTERP;

        %
        if flags.flag_elapsed_time == 1
            title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        colorbar
        caxis([zmin zmax]);
        colormap(Spectrum)
        k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
        ylabel(k,'Rainfall [mm/h]','Interpreter','Latex','FontSize',12)
        xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
        zlabel ('Rainfall [mm/h]','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond','FontSize',12)
        % Increase border (axis) thickness
        ax = gca;          % Get current axis
        ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

        box on
        set(gca,'tickdir','out');
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca,'FontName','Garamond','FontSize',12)
        % Increase border (axis) thickness
        ax = gca;          % Get current axis
        ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
        ax = ancestor(gca, 'axes');
        ax.XAxis.Exponent = 0;xtickformat('%.0f');
        ax.YAxis.Exponent = 0;ytickformat('%.0f');
        if no_plot == 0
            mapshow(S_p,'FaceColor','n'); hold on;
        end
        box on
        % Set background color and write to video
        frame = getframe(gcf);
        writeVideo(video,frame);
        hold off
        clf
    end
end
% Close video writer
close(video);
close all


%% ETP  Video
try 
    delete('Modeling_Results\ETP_Intensities.mp4')
end

%Calculating general Max and Min
if flags.flag_ETP == 1
    store = 1;
    flag_loader=1;
    for i = 1:length(running_control.time_records)
        try
        if i > saver_memory_maps*store
            store = store + 1;
            load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');      
            if flags.flag_ETP
                Max_etp = max(max(Maps.Hydro.ETP_save,[],3),Max_etp);
                Min_etp = min(min(Maps.Hydro.ETP_save,[],3),Min_etp);
            end
        else
            if flag_loader == 1
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                flag_loader=0;
                if flags.flag_ETP
                    Max_etp = max(Maps.Hydro.ETP_save,[],3);
                    Min_etp = min(Maps.Hydro.ETP_save,[],3);
                end
            end
        end
        end
    end
end

if flags.flag_ETP == 1
    close all
    Video_Name = 'ETP_Intensities.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=2;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[0,0,7,12])
    set(gcf,'DefaultTextInterpreter','latex')
    if flags.flag_spatial_rainfall == 1 && flags.flag_rainfall == 1 && flags.flag_input_rainfall_map ~= 1
        close all
        zmax = max(max(Max_etp));
        zmin = min(min(Min_etp));
        if zmin == zmax
            zmax = zmin + 10; % mm/h
        end
        store = 1;
        flag_loader=1;
        for t = 1:1:length(running_control.time_records)
            clf
            if t > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                    flag_loader=0;
                end
            end
            ETP = Maps.Hydro.ETP_save(:,:,t - (store-1)*saver_memory_maps);
            idx = isnan(Elevation_Properties.elevation_cell);
            ETP(ETP<=0) = nan;
            ETP(idx) = nan;
            spatial_ETP = ETP;
            if t > size(Spatial_Rainfall_Parameters.rainfall_spatial_duration,2)
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(end) + Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2);
                spatial_ETP = zeros(size(Elevation_Properties.elevation_cell));
                z = spatial_ETP;
                plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'r.', 'MarkerSize', 30)
            else
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(t);
                % Draw plot for d = x.^n
                z = ETP; %
                z(idx_nan)=nan;
            end
            hold on
            F = z([ybegin:1:yend],[xbegin:1:xend]);
            if no_plot==0
                try
                    mapshow(A,RA,"AlphaData",0.45);hold on;
                    mapshow(S_p,'FaceColor','n'); hold on;
                catch ME

                end
            end
            surf(x_grid,y_grid,F);
            axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])

            shading INTERP;

            %
            if flags.flag_elapsed_time == 1
                title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
            else
                title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
            end
            view(0,90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)
            k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
            ylabel(k,'ETP [mm/day]','Interpreter','Latex','FontSize',12)
            xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
            ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)

            zlabel ('ETP [mm/day]','Interpreter','Latex','FontSize',12)
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

            box on
            set(gca,'tickdir','out');
            set(gca, 'TickLength', [0.02 0.01]);
            set(gca,'Tickdir','out')
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
            ax = ancestor(gca, 'axes');
            ax.XAxis.Exponent = 0;xtickformat('%.0f');
            ax.YAxis.Exponent = 0;ytickformat('%.0f');
            hold off
            if no_plot == 0
                mapshow(S_p,'FaceColor','n'); hold on;
            end
            box on
            % Set background color and write to video
            frame = getframe(gcf);
            writeVideo(video,frame);
            hold off
            clf
        end
    end

    % Close video writer
    close(video);
    close all
end

%% ETR  Video
try 
    delete('Modeling_Results\ETR_Intensities.mp4')
end

%Calculating general Max and Min
if flags.flag_ETP == 1
    store = 1;
    flag_loader=1;
    for i = 1:length(running_control.time_records)
        try
        if i > saver_memory_maps*store
            store = store + 1;
            load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');      
            if flags.flag_ETP
                Max_etr = max(max(Maps.Hydro.ETR_save,[],3),Max_etr);
                Min_etp = min(min(Maps.Hydro.ETR_save,[],3),Min_etp);
            end
        else
            if flag_loader == 1
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                flag_loader=0;
                if flags.flag_ETP
                    Max_etr = max(Maps.Hydro.ETR_save,[],3);
                    Min_etr = max(Maps.Hydro.ETR_save,[],3);
                end
            end
        end
        end
    end
end

if flags.flag_ETP == 1
    close all
    Video_Name = 'ETR_Intensities.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=2;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[0,0,7,12])
    set(gcf,'DefaultTextInterpreter','latex')
    if flags.flag_spatial_rainfall == 1 && flags.flag_rainfall == 1 && flags.flag_input_rainfall_map ~= 1
        close all
        zmax = max(max(Max_etr));
        zmin = min(min(Min_etr));
        if zmin == zmax
            zmax = zmin + 10; % mm/h
        end
        store = 1;
        flag_loader=1;
        for t = 1:1:length(running_control.time_records)
            clf
            if t > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                    flag_loader=0;
                end
            end
            ETR = Maps.Hydro.ETR_save(:,:,t - (store-1)*saver_memory_maps);
            idx = isnan(Elevation_Properties.elevation_cell);
            ETR(ETR<=0) = nan;
            ETR(idx) = nan;
            spatial_ETR = ETR;
            if t > size(Spatial_Rainfall_Parameters.rainfall_spatial_duration,2)
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(end) + Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2);
                spatial_ETR = zeros(size(Elevation_Properties.elevation_cell));
                z = spatial_ETR;
                plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'r.', 'MarkerSize', 30)
            else
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(t);
                % Draw plot for d = x.^n
                z = ETR; %
                z(idx_nan)=nan;
            end
            hold on
            F = z([ybegin:1:yend],[xbegin:1:xend]);
            if no_plot==0
                try
                    mapshow(A,RA,"AlphaData",0.45);hold on;
                    mapshow(S_p,'FaceColor','n'); hold on;
                catch ME

                end
            end
            surf(x_grid,y_grid,F);
            axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])

            shading INTERP;

            %
            if flags.flag_elapsed_time == 1
                title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
            else
                title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
            end
            view(0,90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)
            k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
            ylabel(k,'ETP [mm/day]','Interpreter','Latex','FontSize',12)
            xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
            ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)

            zlabel ('ETP [mm/day]','Interpreter','Latex','FontSize',12)
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 

            box on
            set(gca,'tickdir','out');
            set(gca, 'TickLength', [0.02 0.01]);
            set(gca,'Tickdir','out')
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
            ax = ancestor(gca, 'axes');
            ax.XAxis.Exponent = 0;xtickformat('%.0f');
            ax.YAxis.Exponent = 0;ytickformat('%.0f');
            hold off
            if no_plot == 0
                mapshow(S_p,'FaceColor','n'); hold on;
            end
            box on
            % Set background color and write to video
            frame = getframe(gcf);
            writeVideo(video,frame);
            hold off
            clf
        end
    end

    % Close video writer
    close(video);
    close all
end

%% GW Video
try 
    delete('Modeling_Results\GW_Depths.mp4')
end

% Calculating general Max and Min
if flags.flag_baseflow == 1
    store = 1;
    flag_loader = 1;
    for i = 1:length(running_control.time_records)
        try
            Max_gw = zeros(size(elevation));
            Min_gw = zeros(size(elevation));
            if i > saver_memory_maps * store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_', num2str(store)), 'Maps');      
                if flags.flag_ETP
                    Max_gw = max(max(Maps.Hydro.GWdepth_save, [], 3), Max_gw);
                    Min_gw = min(min(Maps.Hydro.GWdepth_save, [], 3), Min_gw);
                end
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_', num2str(store)), 'Maps');
                    flag_loader = 0;
                    if flags.baseflow
                        Max_gw = max(Maps.Hydro.GWdepth_save, [], 3);
                        Min_gw = max(Maps.Hydro.GWdepth_save, [], 3);
                    end
                end
            end
        end
    end
end

if flags.flag_baseflow == 1
    close all
    Video_Name = 'GW_Depths.mp4';
    
    % Set up video
    video = VideoWriter(fullfile(folderName, Video_Name), 'MPEG-4');
    video.FrameRate = 2;
    open(video);
    
    h = figure;
    set(gcf, 'units', 'inches', 'position', [0, 0, 10, 7]) % Larger figure size
    set(gcf, 'DefaultTextInterpreter', 'latex')
    set(gcf, 'PaperPositionMode', 'auto'); % Ensures high-quality rendering

    if flags.flag_baseflow == 1
        close all
        zmax = max(max(Max_gw));
        zmin = min(min(Min_gw));
        if zmin == zmax
            zmax = zmin + 10; % mm/h
        end
        store = 1;
        flag_loader = 1;
        for t = 1:length(running_control.time_records)
            clf
            if t > saver_memory_maps * store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_', num2str(store)), 'Maps');
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_', num2str(store)), 'Maps');
                    flag_loader = 0;
                end
            end
            GW = Maps.Hydro.GWdepth_save(:, :, t - (store - 1) * saver_memory_maps);
            idx = isnan(Elevation_Properties.elevation_cell);
            GW(GW < 0) = nan;
            GW(idx) = nan;
            spatial_GW = GW;
            
            if t > size(Spatial_Rainfall_Parameters.rainfall_spatial_duration, 2)
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(end) + Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2);
                spatial_GW = zeros(size(Elevation_Properties.elevation_cell));
                z = spatial_GW;
            else
                t_title = Spatial_Rainfall_Parameters.rainfall_spatial_duration(t);
                z = spatial_GW;
                z(idx_nan) = nan;
            end
            
            hold on
            F = z([ybegin:1:yend], [xbegin:1:xend]);
            if no_plot == 0
                try
                    mapshow(A, RA, "AlphaData", 0.45); hold on;
                    mapshow(S_p, 'FaceColor', 'n'); hold on;
                catch ME
                end
            end
            surf(x_grid, y_grid, F);
            axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])
            shading INTERP;
            
            if flags.flag_elapsed_time == 1
                title(sprintf('Time [h] = %4.2f', t_title / 60), 'Interpreter', 'Latex', 'FontSize', 14)
            else
                title(sprintf(string(t_title)), 'Interpreter', 'Latex', 'FontSize', 14);
            end
            view(0, 90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)
            k = colorbar; 
            k.FontName = 'Garamond'; 
            k.FontSize = 14; 
            k.TickDirection = 'out';
            ylabel(k, 'GW [m]', 'Interpreter', 'Latex', 'FontSize', 14)
            xlabel('Easting [m]', 'Interpreter', 'Latex', 'FontSize', 14)
            ylabel('Northing [m]', 'Interpreter', 'Latex', 'FontSize', 14)
            zlabel('GW [m]', 'Interpreter', 'Latex', 'FontSize', 14)
            set(gca, 'FontName', 'Garamond', 'FontSize', 14)

            % Increase border (axis) thickness
            ax = gca; 
            ax.LineWidth = 2;  
            box on
            set(gca, 'TickDir', 'out', 'TickLength', [0.02 0.01]);
            set(gca, 'FontName', 'Garamond', 'FontSize', 14)
            ax.XAxis.Exponent = 0; xtickformat('%.0f');
            ax.YAxis.Exponent = 0; ytickformat('%.0f');

            hold off
            if no_plot == 0
                mapshow(S_p, 'FaceColor', 'n'); hold on;
            end
            box on

            % Save high-resolution frame before writing to video
            frame_filename = 'temp_frame.png';
            print(gcf, '-dpng', '-r300', frame_filename); % 300 DPI output
            frame_img = imread(frame_filename);
            writeVideo(video, frame_img);
            delete(frame_filename); % Clean up temporary file

            hold off
            clf
        end
    end

    % Close video writer
    close(video);
    close all
end


%% Plot Risk Maps
% Adjusting the size
if flags.flag_human_instability > 0
    % Generate video showing water level profile over time
    close all
    Video_Name = 'Slide_Risk.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=1.5;
    video.Quality = 100;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','pixels','position',[0,0,1080,1080])
    set(gcf,'DefaultTextInterpreter','latex')
    set(gcf,'units','pixels','position',[0,0,1080,1080])

    % loop for video frames
    store=1;
    flag_loader=1;
    zzz = zeros(size(DEM_raster.Z,1),size(DEM_raster.Z,2),8);
    for i = 1:length(running_control.time_records)
        clf
        % subplot(2,1,2)
        t_title = running_control.time_records(i);
        if i > saver_memory_maps*store
            store = store + 1;
            load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
        else
            if flag_loader == 1
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                flag_loader=0;
            end
        end

        % subplot(2,1,1)
        subplot(2,1,1)
        set(gca,'Position',[0.06 0.85 0.9 0.1])
        if flags.flag_rainfall == 1
            if flags.flag_spatial_rainfall ~=1
                plot(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall), '-o','LineWidth',2);hold on
                ylabel(['Rainfall Intensity' newline '[mm/h]'],'interpreter','latex','FontSize',12);
                ylim([0 max(gather(Rainfall_Parameters.intensity_rainfall))*1.2])
                scatter(gather(Rainfall_Parameters.time_rainfall(i)),gather(Rainfall_Parameters.intensity_rainfall(i)),"filled",'MarkerFaceColor',[102 255 102]./255,'SizeData',100,'LineWidth',2.5,'MarkerEdgeColor',[0 0 0])
                set(gca,'ydir','reverse')
            else
                % bar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),gather(BC_States.average_spatial_rainfall),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
                % ylabel('Aerial Mean Rainfall Intensity [mm/h]','interpreter','latex');
                % hold on
                % er = errorbar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),BC_States.average_spatial_rainfall,Rainfall_Parameters.std_dev(1:(dim),1),Rainfall_Parameters.std_dev(1:(dim),1));
                % er.Color = [0 0 0];
                % er.LineStyle = 'none';
                % plot((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)))),gather(BC_States.average_spatial_rainfall),'LineWidth',1.5,'color','blue')
                % ylim([0 max(max(gather(BC_States.average_spatial_rainfall)))*6])
            end
        end

        if flags.flag_elapsed_time == 1
            title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end

        for j = 1:8
            zzz(:,:,j) = double(Maps.Hydro.(strcat('risk',Human_Instability_text.list{j}))(:,:,i - (store-1)*saver_memory_maps)>0)*find(Human_Instability.order==Human_Instability.order(j));
        end
        zzz_2 = max(zzz,[],3);
        zzz_2(zzz_2==0)=nan;

        F = zzz_2([ybegin:1:yend],[xbegin:1:xend]);
        % Plotting base map and study area boundaries
        subplot(2,1,2)
        set(gca,'Position',[0.035 0.05 0.9 0.75])
        if no_plot==0
            try
                mapshow(A,RA,"AlphaData",0.25);hold on;
                mapshow(S_p,'FaceColor','n'); hold on;
            catch ME

            end
        end

        cmap = [204 255 153; 102 255 102; 255 255 153; 255 255 0; 255 178 102; 255 128 0; 255 102 102; 255 51 51]./255;
        colormap(cmap);
        s= pcolor(x_grid,y_grid,F);
        s.EdgeColor = 'none';
        caxis([1 8]); % Set the colorbar scale to 1 to 8

        k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
        xlabel(' Easting [m] ','Interpreter','Latex','FontSize',14)
        ylabel ('Northing [m] ','Interpreter','Latex','FontSize',14)
        colorbar('Ticks',[1.4,1.6, 2.2,2.4, 3.1,3.3, 4.0,4.2, 4.9,5.1, 5.7,5.9, 6.6,6.8, 7.5,7.7], ...
            'TickLabels',Human_Instability_text.names,'FontSize',13, 'TickLength',0);
        set(gca,'FontName','Garamond','FontSize',12)
        % Increase border (axis) thickness
        ax = gca;          % Get current axis
        ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
        box on
        ax = ancestor(gca, 'axes');
        ax.XAxis.Exponent = 0;xtickformat('%.0f');
        ax.YAxis.Exponent = 0;ytickformat('%.0f');
        hold off
        % Set background color and write to video
        frame = getframe(gcf);

        % Save high-resolution frame before writing to video
        frame_filename = 'temp_frame.png';
        print(gcf, '-dpng', '-r300', frame_filename); % 300 DPI output
        frame_img = imread(frame_filename);
        writeVideo(video, frame_img);
        delete(frame_filename); % Clean up temporary file

        writeVideo(video,frame);
        box on
        hold off
        clf
    end
end
% Close video writer
close(video);
close all

%% Pollutant Concentration
if flags.flag_waterquality == 1
    % Adjusting the size
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    FileName_String = 'Pollutant_Concentration.gif';
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[3,3,6.5,5])

    % Replace all infs for nan
    z = Maps.WQ_States.Pol_Conc_Map; % Plotting concentration
    z(idx_nan) = nan;
    z(z<LULC_Properties.Pol_min) = nan;
    if ~isnan(max(max(max(z))))
        for t = 1:f:length(running_control.time_records)
            t_title = running_control.time_records(t);
            % Draw plot for d = x.^n
            zmax = max(max(z(:,:,t)));
            zmin = min(min(z(:,:,t)));
            zmax = max(zmax,0);
            zmin = max(zmin,0);
            if zmax == 0 || zmin == zmax
                zmin = min(min(min(z)));
                zmax = max(max(max(z)));
            end
            z(idx_nan) = nan;
            F = z([ybegin:1:yend],[xbegin:1:xend],t);
            if no_plot==0;
                try
                    mapshow(A,RA,"AlphaData",0.45);hold on;
                    mapshow(S_p,'FaceColor','n'); hold on;
                catch ME

                end
            end
            surf(x_grid,y_grid,F);
            axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax])

            shading INTERP;
            if flags.flag_elapsed_time == 1
                title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
            else
                title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
            end
            view(0,90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)

            k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
            ylabel(k,'Concentration (mg/L)','Interpreter','Latex','FontSize',12)
            xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
            ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
            zlabel ('Concentration (mg/L)','Interpreter','Latex','FontSize',12)
            set(gca,'FontName','Garamond','FontSize',12)
            % Increase border (axis) thickness
            ax = gca;          % Get current axis
            ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
            box on
            ax = ancestor(gca, 'axes');
            ax.XAxis.Exponent = 0;xtickformat('%.0f');
            ax.YAxis.Exponent = 0;ytickformat('%.0f');
            drawnow
            % Capture the plot as an image
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if t == 1
                imwrite(imind,cm,FileName,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,FileName,'gif','WriteMode','append');
            end
        end
    end
    
    clf
    close all

    % Pollutant Mass
    % Adjusting the size
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    FileName_String = 'Mass_of_pollutant.gif';
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[3,3,6.5,5])
    % Replace all infs for nan
    z = Maps.WQ_States.Pol_mass_map/Wshed_Properties.cell_area*1000; % Plotting pollutant mass in grams
    z = log(z);
    z(isinf(z)) = nan;
    zmax = max(max(max(z)));
    zmin = min(min(min(z)));
    for t = 1:f:length(running_control.time_records)
        t_title = running_control.time_records(t);
        % Draw plot for d = x.^n
        LULC_Properties.Pol_min = 0.01; % g/m2
        z(z<=LULC_Properties.Pol_min)=nan;
        F = z([ybegin:1:yend],[xbegin:1:xend],t);
        if no_plot==0;
            mapshow(A,RA,"AlphaData",0.45);hold on;
            mapshow(S_p,'FaceColor','n'); hold on;
        end
        surf(x_grid,y_grid,F);
        axis tight
        shading INTERP;
        if flags.flag_elapsed_time == 1
            title(sprintf('Time [h] = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        c = colorbar;
        caxis([zmin zmax]);
        colormap(Velocity_RAS)
        k = colorbar; k.FontName = 'Garamond'; k.FontSize = 12; k.TickDirection  = 'out';
        ylabel(k,'Log-scale Mass of pollutant ($\mathrm{g/m^2}$)','Interpreter','Latex','FontSize',12)
        xlabel(' Easting [m] ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing [m] ','Interpreter','Latex','FontSize',12)
        zlabel ('Mass of pollutant (kg/m^2)','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond','FontSize',12)
        % Increase border (axis) thickness
        ax = gca;          % Get current axis
        ax.LineWidth = 2;   % Set the line width to 2 (default is 0.5) 
        box on
        ax = ancestor(gca, 'axes');
        ax.XAxis.Exponent = 0;xtickformat('%.0f');
        ax.YAxis.Exponent = 0;ytickformat('%.0f');
        drawnow
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if t == 1
            imwrite(imind,cm,FileName,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,FileName,'gif','WriteMode','append');
        end
    end
end

close all