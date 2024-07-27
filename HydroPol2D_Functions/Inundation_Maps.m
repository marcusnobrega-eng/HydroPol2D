%%% Inundation Maps %%%
% Creates .GIF files from the saved data
% Last updated: 10/4/2021

close all
f = 1;
t_max = running_control.routing_time;
tfinal = t_max ; %t_max;
DEM_maps = gather(Elevation_Properties.elevation_cell);

%% Time Data
if flags.flag_elapsed_time ~=1.
    if flags.flag_spatial_rainfall == 1
        if ~isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration)
            Spatial_Rainfall_Parameters.rainfall_spatial_duration = Spatial_Rainfall_Parameters.rainfall_spatial_duration/60/24 + date_begin;
        end
    end
    if flags.flag_ETP == 1
        if ~isdatetime(ETP_Parameters.climatologic_spatial_duration)
            ETP_Parameters.climatologic_spatial_duration = ETP_Parameters.climatologic_spatial_duration/60/24 + date_begin;
        end
    end
    if ~isdatetime(running_control.time_records)
        running_control.time_records =  double(running_control.time_records/60/24) + date_begin;
    end
end

%% Creating the custom basemap
basemapName = "openstreetmap";
url = "c.tile.openstreetmap.org/${z}/${x}/${y}.png";
url2 = 'a';
copyright = char(uint8(169));
attribution = copyright + "OpenStreetMap contributors";
attribution_2 = "Stadia Maps @ OpenStreetMap contributors";
addCustomBasemap(basemapName,url,"Attribution",attribution)

% Getting lat and lon from the study area

[lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
latlim = [lat(1) lat(2)];
lonlim = [lon(1) lon(2)];
% Retriving the basemap image
try
    [A,RA,attribA] = readBasemapImage(basemapName,latlim,lonlim);
catch ME
    warning('You need matlab 2022a or higher to use basemaps in georeference plots.')
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

%% Set up HEC-RAS colors
hec_ras_colors = [52/255 85/255 132/255; 0 1 1; 0 128/255 1; 0 255/255 0; 1 1 0; 1 128/255 0; 1 0 0; 128/255 0 128/255];

%% Plot Elevation Model
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
FileName_String = 'Elevation_Model.gif';
FileName = fullfile(folderName,strcat('\',FileName_String));

a_grid = Wshed_Properties.Resolution;
b_grid = Wshed_Properties.Resolution;
tmax = 10;
for t = 1:tmax
    % Draw plot
    z = DEM_maps;
    xmax = length(z(1,:));
    xend = xmax;
    ymax = length(z(:,1));
    yend = ymax;
    xbegin = 1;
    ybegin = 1;
    max_h = max(max(max(z)));
    h_max = round(max_h/10,0)*10*1.05;
    % UTM Coordinates
    x_grid = GIS_data.xulcorner + a_grid*[xbegin:1:xend]; y_grid = GIS_data.yulcorner - a_grid*[ybegin:1:yend];
    %     y_grid = flip(y_grid); % Make sure we plot the graphs properly with y towards vertical top.
    z(z<=0)=inf;
    h_min = min(min(z));
    F = DEM_maps([ybegin:1:yend],[xbegin:1:xend]);
    zmax = max(max(max(z(~isinf(z)))));
    kk = surf(x_grid,y_grid,F);
    set(kk,'LineStyle','none');
    set(gca,'XTickLabel',x_grid)
    set(gca,'YTickLabel',y_grid)
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) h_min zmax])
    view(-(t-1)*360/tmax,(t-1)*90/tmax);
    colorbar
    caxis([h_min zmax]);
    colormap(Terrain_RAS)
    box on
    hold on
    title('Elevation','Interpreter','Latex','FontSize',12)
    k = colorbar ;
    ylabel(k,'Elevation (m)','Interpreter','Latex','FontSize',12)
    xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
    zlabel ('Elevation (m)','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond')
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
clf

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
    subplot(2,1,1)
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
    surf(x_grid,y_grid,F);
    shading INTERP;

    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) z1min z1max])
    if flags.flag_elapsed_time == 1
        title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
    else
        title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
    end

    view(0,90);
    colorbar
    caxis([z1min z1max]);
    colormap(WSE_RAS)
    k = colorbar;
    ylabel(k,'WSE (m)','Interpreter','Latex','FontSize',12)
    xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
    ylabel (' Northing (m) ','Interpreter','Latex','FontSize',12)
    zlabel ('WSE (m)','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond')
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    box on
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end

    subplot(2,1,2)
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
    surf(x_grid,y_grid,F);
    axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) z2min z2max])

    shading INTERP;
    if flags.flag_elapsed_time == 1
        title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
    else
        title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
    end
    view(0,90);
    colorbar
    caxis([0 z2max]);
    colormap(Depth_RAS)
    k = colorbar;

    ylabel(k,'Depths (m)','Interpreter','Latex','FontSize',12)
    xlabel('Easting (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
    zlabel ('WSE (m)','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond')
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
%         title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
%     else
%         title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
%     end
%
%     view(0,90);
%     colorbar
%     caxis([zmin zmax]);
%     colormap(depth_ramp)
%     k = colorbar;
%     ylabel(k,'WSE (m)','Interpreter','Latex','FontSize',12)
%     xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
%     ylabel (' Northing (m) ','Interpreter','Latex','FontSize',12)
%     zlabel ('WSE (m)','Interpreter','Latex','FontSize',12)
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
%         title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
%     else
%         title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
%     end
%     view(0,90);
%     colorbar
%     caxis([0 zmax]);
%     colormap(depth_ramp)
%     k = colorbar;
%
%     ylabel(k,'Depths (m)','Interpreter','Latex','FontSize',12)
%     xlabel('Easting (m) ','Interpreter','Latex','FontSize',12)
%     ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
%     zlabel ('WSE (m)','Interpreter','Latex','FontSize',12)
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
    rain_total = rainfall_sum*time_step_rainfall; % data comes form post_processing
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
    k = colorbar;
    ylabel(k,'Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
    xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)

    zlabel ('Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
    set(gca,'FontName','Garamond')

    box on
    set(gca,'tickdir','out');
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca,'FontName','Garamond');
    ax = ancestor(gca, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');
    view(0,90)
    if no_plot == 0
        mapshow(S_p,'FaceColor','n'); hold on;
    end
    exportgraphics(gcf,fullfile(folderName,'Isoietal_Map.png'),'ContentType','image','Colorspace','rgb','Resolution',1200)
    saveas(gcf,fullfile(folderName,'Isoeital.fig'))
    close all
end

%% Dividing the Rainfall into n intervals
% if flags.flag_spatial_rainfall == 1
%     n_plots = 8;
%     n_rainfall = size(Maps.Hydro.spatial_rainfall_maps*time_step_rainfall,3);
%     time_step_rainfall_aggregated = floor((n_rainfall/n_plots))*time_step_rainfall;
%     zmin = 0;
%     zmax = max(max(max((Maps.Hydro.spatial_rainfall_maps))))*time_step_rainfall_aggregated;
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
%         k = colorbar;
%         ylabel(k,'Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
%         xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
%         ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
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
close all
%% Spatial Rainfall Video
% Generate video showing water level profile over time

%Calculating general Max and Min
if flags.flag_spatial_rainfall == 1
    store = 1;
    flag_loader=1;
    for i = 1:length(running_control.time_records)
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
                title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
            else
                title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
            end
            view(0,90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)
            k = colorbar;
            ylabel(k,'Rainfall (mm/h)','Interpreter','Latex','FontSize',12)
            xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
            ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)

            zlabel ('Rainfall (mm/h)','Interpreter','Latex','FontSize',12)
            set(gca,'FontName','Garamond')

            box on
            set(gca,'tickdir','out');
            set(gca, 'TickLength', [0.02 0.01]);
            set(gca,'Tickdir','out')
            set(gca,'FontName','Garamond');
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
%% Spatial Rainfall GIFS Flag_Spatial_Rainfall
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
            title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        colorbar
        caxis([zmin zmax]);
        colormap(Spectrum)
        k = colorbar;
        ylabel(k,'Rainfall (mm/h)','Interpreter','Latex','FontSize',12)
        xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)

        zlabel ('Rainfall (mm/h)','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond')

        box on
        set(gca,'tickdir','out');
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca,'FontName','Garamond');
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
            title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        colorbar
        caxis([zmin zmax]);
        colormap(Spectrum)
        k = colorbar;
        ylabel(k,'Rainfall (mm/h)','Interpreter','Latex','FontSize',12)
        xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
        zlabel ('Rainfall (mm/h)','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond')

        box on
        set(gca,'tickdir','out');
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca,'FontName','Garamond');
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


%% Flag_ETP
if flags.flag_ETP == 1
    % Generate video showing water level profile over time
    close all
    Video_Name = 'ETP_Maps.mp4';
    % Set up video
    video = VideoWriter(fullfile(folderName,Video_Name),'MPEG-4');
    video.FrameRate=2;
    open(video);
    h = figure;
    FileName = fullfile(folderName,strcat('\',FileName_String));
    set(gcf,'units','inches','position',[0,0,6.5,4])
    set(gcf,'DefaultTextInterpreter','latex')

    zmax = max(max(max(Maps.Hydro.ETP_save)));
    zmin = min(min(min(Maps.Hydro.ETP_save)));
    if zmin == zmax
        zmax = zmin + 10; % mm/h
    end
    tmax_plot = min(days(date_end - date_begin),((size(Maps.Hydro.ETP_save,3))));
    for t = 1:tmax_plot
        clf
        ETP = Maps.Hydro.ETP_save(:,:,t);
        idx = isnan(Elevation_Properties.elevation_cell);
        ETP(ETP<=0) = nan;
        ETP(idx) = nan;
        spatial_ETP = ETP;
        t_title = ETP_Parameters.climatologic_spatial_duration(t);
        % Draw plot for d = x.^n
        z = ETP; %
        z(idx_nan)=nan;
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
        axis tight
        shading INTERP;

        %
        if flags.flag_elapsed_time == 1
            title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        colorbar
        caxis([zmin zmax]);
        colormap(Spectrum)
        k = colorbar;
        ylabel(k,'ETP (mm/day)','Interpreter','Latex','FontSize',12)
        xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
        zlabel ('ETP (mm/day)','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond')
        ax = ancestor(gca, 'axes');
        ax.XAxis.Exponent = 0;xtickformat('%.0f');
        ax.YAxis.Exponent = 0;ytickformat('%.0f');

        drawnow
        box on
        set(gca,'tickdir','out');
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca,'FontName','Garamond');
        if no_plot == 0
            mapshow(S_p,'FaceColor','n'); hold on;
        end
        box on
        hold off
        % Set background color and write to video
        frame = getframe(gcf);
        writeVideo(video,frame);
        clf
    end
end
% Close video writer
close(video);
close all

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
                ylabel(['Rainfall Intensity' newline '(mm/h)'],'interpreter','latex','FontSize',12);
                ylim([0 max(gather(Rainfall_Parameters.intensity_rainfall))*1.2])
                scatter(gather(Rainfall_Parameters.time_rainfall(i)),gather(Rainfall_Parameters.intensity_rainfall(i)),"filled",'MarkerFaceColor',[102 255 102]./255,'SizeData',100,'LineWidth',2.5,'MarkerEdgeColor',[0 0 0])
                set(gca,'ydir','reverse')
            else
                % bar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),gather(BC_States.average_spatial_rainfall),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
                % ylabel('Aerial Mean Rainfall Intensity (mm/h)','interpreter','latex');
                % hold on
                % er = errorbar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),BC_States.average_spatial_rainfall,Rainfall_Parameters.std_dev(1:(dim),1),Rainfall_Parameters.std_dev(1:(dim),1));
                % er.Color = [0 0 0];
                % er.LineStyle = 'none';
                % plot((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)))),gather(BC_States.average_spatial_rainfall),'LineWidth',1.5,'color','blue')
                % ylim([0 max(max(gather(BC_States.average_spatial_rainfall)))*6])
            end
        end

        if flags.flag_elapsed_time == 1
            title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
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

        k = colorbar;
        xlabel(' Easting (m) ','Interpreter','Latex','FontSize',14)
        ylabel ('Northing (m) ','Interpreter','Latex','FontSize',14)
        colorbar('Ticks',[1.4,1.6, 2.2,2.4, 3.1,3.3, 4.0,4.2, 4.9,5.1, 5.7,5.9, 6.6,6.8, 7.5,7.7], ...
            'TickLabels',Human_Instability_text.names,'FontSize',13, 'TickLength',0);
        set(gca,'FontName','Garamond')
        box on
        ax = ancestor(gca, 'axes');
        ax.XAxis.Exponent = 0;xtickformat('%.0f');
        ax.YAxis.Exponent = 0;ytickformat('%.0f');
        hold off
        % Set background color and write to video
        frame = getframe(gcf);
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
                title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
            else
                title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
            end
            view(0,90);
            colorbar
            caxis([zmin zmax]);
            colormap(Spectrum)

            k = colorbar;
            ylabel(k,'Concentration (mg/L)','Interpreter','Latex','FontSize',12)
            xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
            ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
            zlabel ('Concentration (mg/L)','Interpreter','Latex','FontSize',12)
            set(gca,'FontName','Garamond')
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
            title(sprintf('Time(h) = %4.2f',t_title/60),'Interpreter','Latex','FontSize',12)
        else
            title(sprintf(string(t_title)),'Interpreter','Latex','FontSize',12);
        end
        view(0,90);
        c = colorbar;
        caxis([zmin zmax]);
        colormap(Velocity_RAS)
        k = colorbar;
        ylabel(k,'Log-scale Mass of pollutant ($\mathrm{g/m^2}$)','Interpreter','Latex','FontSize',12)
        xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)
        zlabel ('Mass of pollutant (kg/m^2)','Interpreter','Latex','FontSize',12)
        set(gca,'FontName','Garamond')
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