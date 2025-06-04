%% Normalized Discharge gauges%%
normalized_table = table(gather(running_control.time_hydrograph),'VariableNames', {'Time'});
save('myWorkspace4graphs.mat')

zero_matrix = zeros(size(elevation,1),size(elevation,2));
if flags.flag_obs_gauges == 1 && flags.flag_rainfall == 1
    % Rainfall Std Deviation
    if flags.flag_spatial_rainfall == 1
        store=1;
        flag_loader=1;
        for i = 1:length(running_control.time_records)
            if i > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            else
                if flag_loader == 1
                  load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                  flag_loader=0;
                end
            end
            for j = 1:length(gauges.easting_obs_gauges)
                zgauges = Maps.Hydro.spatial_rainfall_maps(:,:,i - ((store-1)*saver_memory_maps)).*subcatchments{1,j};
                Rainfall_Parameters.std_dev_gauges{1,j}(i,1) = nanstd(zgauges(:));
            end
        end    
    end
    % Catchment area of each gauge and graphics
    for i = 1:length(gauges.easting_obs_gauges)
        gauges.catchment_area(i,1) = Wshed_Properties.fac_area(gauges.northing_obs_gauges(i,1),gauges.easting_obs_gauges(i,1)); % km2
    end

    %% Creating the folder
    mkdir(strcat(folderName,'\Normalized_Discharge_gauges'));
    % Specify the folder where the files live.
    myFolder_wd = strcat(pwd,'\',folderName,'\Normalized_Discharge_gauges'); % Current folder
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isfolder(myFolder_wd)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder_wd);
        uiwait(warndlg(errorMessage));
        return;
    end
    DEM_maps = gather(Elevation_Properties.elevation_cell);
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
    a_grid = Wshed_Properties.Resolution;
    b_grid = Wshed_Properties.Resolution;
    tmax = 10;
    x_grid = GIS_data.xulcorner + a_grid*[xbegin:1:xend]; y_grid = GIS_data.yulcorner - a_grid*[ybegin:1:yend];
    
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
if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= web_mercator_crs.Name
    no_plot = 1;
else
    S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
    S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)]; % Calculate bounding box for each polygon
    S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
    S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
end
%% Read the rainfall file
if flags.flag_spatial_rainfall == 1
    if isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration(2))
        time_step_rainfall = hours(Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) -Spatial_Rainfall_Parameters.rainfall_spatial_duration(1)); % hours
    else
        time_step_rainfall = (Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) -Spatial_Rainfall_Parameters.rainfall_spatial_duration(1))/60; % hours
    end
end

%% Plotting by gauges
labels_depth = gauges.labels_observed_string; labels_depth{gauges.num_obs_gauges+1} = 'Outlet';
labels_depth{gauges.num_obs_gauges + 2} = 'Rainfall Intensity';

for i = 1:length(gauges.easting_obs_gauges)
    % Initialize variables to track min and max values for scaling
    minRain = Inf;
    maxRain = -Inf;
    minDischarge = Inf;
    maxDischarge = -Inf;

    %Create a figure
    figure('Units', 'centimeters', 'Position', [2 2 29 12]); % ~18x9 cm is common for Nature half-width figures

    % Setup tiled layout for more control than subplot
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    %subplot(2, 1, 2);
    %set(gcf,'units','inches','position',[3,0,9,5]);

    %normalized discharge
    nexttile
    set(gca,'ycolor','black');
    specific_discharge = gauges.hydrograph_cell(:,i)/gauges.catchment_area(i,1) ; % m3/s/km2
    normalized_table.(strcat('Q_',labels_depth{i})) = specific_discharge; %save the values

    plot(gather(running_control.time_hydrograph),specific_discharge,'LineWidth',1.5,'linestyle',ls);
    hold on;

    xlabel('Elapsed Time (min)','interpreter','latex'); ylabel('Specific Discharge $(\mathrm{m^3/s/km^2})$','interpreter','latex'); set(gca,'FontSize',12);

    %Rainfall in a inverse y-axis
    yyaxis right; set(gca,'ydir','reverse','ycolor','black');
    if flags.flag_rainfall == 1
        if flags.flag_spatial_rainfall ~=1
            bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',[0 .55 .55],'EdgeColor',[0 .5 .5],'LineWidth',1.5)
            ylabel('Rainfall Intensity (mm/h)','interpreter','latex');
            ylim([0 max(gather(Rainfall_Parameters.intensity_rainfall))*6])

            Rain = gather(Rainfall_Parameters.intensity_rainfall);
            normalized_table.(strcat('P_',labels_depth{i})) = Rainfall_Parameters.intensity_rainfall;
        else
            bar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),gather(BC_States.average_spatial_rainfall_gauges{1,i}),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
            ylabel('Aerial Mean Rainfall Intensity (mm/h)','interpreter','latex'); set(gca,'FontSize',12);
            hold on
            try
                er = errorbar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),BC_States.average_spatial_rainfall_gauges{1,i},Rainfall_Parameters.std_dev_gauges{1,i}(1:(dim),1),Rainfall_Parameters.std_dev_gauges{1,i}(1:(dim),1));
                er.Color = [0 0 0];
                er.LineStyle = 'none';
            end
            plot((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)))),gather(BC_States.average_spatial_rainfall_gauges{1,i}),'LineWidth',1.5,'color','blue')
            ylim([0 max(max(gather(BC_States.average_spatial_rainfall_gauges{1,i})))*6])

            pivot_column = BC_States.average_spatial_rainfall_gauges{1,i};
            pivot_column(end+1) = 0;
            Rain = gather(BC_States.average_spatial_rainfall_gauges{1,i});

            normalized_table.(strcat('P_',labels_depth{i})) = pivot_column; %save the values
            normalized_table.(strcat('Std_dev_',labels_depth{i})) = Rainfall_Parameters.std_dev_gauges{1,i};%save the values
        end
    end

    % Update min and max values
    minRain = min(minRain, min(Rain));
    maxRain = max(maxRain, max(Rain));
    minDischarge = min(minDischarge, min(specific_discharge));
    maxDischarge = max(maxDischarge, max(specific_discharge));

    % Format X-axis

    xlim([min(gather(running_control.time_hydrograph)), max(gather(running_control.time_hydrograph))]);
    % Set limits and rescale discharge
    yyaxis left;
    ylim([0, 1.55*maxDischarge]); % Set limits for Discharge
    % Set limits and rescale Rain
    yyaxis right;
    scaledMinRain = 0;
    scaledMaxRain = (maxRain/(0.5*maxDischarge))*(1.55*maxDischarge);
    ylim([scaledMinRain, scaledMaxRain]); % Set scaled limits for Rain

    grid on;

    lgd_hydro = legend([strcat('SD on',{' '},labels_depth{i}),labels_depth(end),'Standar Dev'],'Interpreter','Latex','FontSize',8,'location','east');

    title(strcat('Specific Discharge on gauge',{' '},labels_depth{i}),'interpreter','latex','fontsize',12);
    box on;
    set(gca,'tickdir','out');
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca,'FontName','Garamond');

    %rainfall map
    if flags.flag_spatial_rainfall == 1
        nexttile
        hold on;
        rain_total = rainfall_sum*time_step_rainfall.*subcatchments{1,i}; % data comes form post_processing
        time_total = days(date_end - date_begin);
        title_isoietal = strcat('Cumulative rainfall of the event on the drainage area');
        idx = isnan(Elevation_Properties.elevation_cell);
        rain_total(rain_total<=0) = nan;
        rain_total(idx) = nan;
        zmin = min(min(rain_total));
        zmax = max(max(rain_total));
        F = rain_total([ybegin:1:yend],[xbegin:1:xend]);
        surf(x_grid,y_grid,F);hold on;
        axis([min(min(x_grid)) max(max(x_grid)) min(min(y_grid)) max(max(y_grid)) zmin zmax]); hold on;
        shading INTERP;
        title(title_isoietal,'Interpreter','Latex','FontSize',12);
        colorbar
        caxis([zmin zmax]);
        colormap(linspecer)
        k = colorbar;
        %F_subcatchment = F.*subcatchments{1,i};
        if no_plot==0
            try
                zLevel = zmin;  % or (zmin + zmax)/2 if you want to center it

                % For the raster map, use a surf workaround
                [xRaster, yRaster] = meshgrid(RA.XWorldLimits(1):RA.CellExtentInWorldX:RA.XWorldLimits(2)-RA.CellExtentInWorldX, ...
                    RA.YWorldLimits(2):-RA.CellExtentInWorldY:RA.YWorldLimits(1)+RA.CellExtentInWorldY);
                surf(xRaster, yRaster, zLevel * ones(size(A(:,:,1))), ...
                    'CData', A, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.45); hold on;
                
                % Polygon overlay (after raster)
                zLevelShape = zmin + 0.01;  % Slightly above raster surface
                X = S_p.X;
                Y = S_p.Y;

                % Handle NaN-separated parts
                nanIdx = isnan(X);
                partStart = [1, find(nanIdx) + 1];
                partEnd = [find(nanIdx) - 1, length(X)];

                for x1 = 1:length(partStart)
                    xi = X(partStart(x1):partEnd(x1));
                    yi = Y(partStart(x1):partEnd(x1));
                    if all(isnan(xi)) || all(isnan(yi))
                        continue;
                    end
                    fill3(xi, yi, repmat(zLevelShape, size(xi)), ...
                        'k', 'LineWidth', 1.2, 'FaceColor', 'none', 'EdgeAlpha',0.5); % Border only
                end
                
            catch ME

            end
        end
        
        ylabel(k,'Cumulative Rainfall Volume (mm)','Interpreter','Latex','FontSize',12)
        xlabel(' Easting (m) ','Interpreter','Latex','FontSize',12)
        ylabel ('Northing (m) ','Interpreter','Latex','FontSize',12)

        hold on;
        %plotting gauges location
        if flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
            rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(1,1:Spatial_Rainfall_Parameters.n_raingauges)'; % Values of rainfall at t for each rain gauge
            % idx_rainfall = logical(isnan(rainfall) | rainfall == 0);

            idx_rainfall = logical(isnan(rainfall));
            Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1); % Coordinates (easting) of each rain gauge
            Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2); % Coordinates (northing) of each rain gauge
            Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = []; % Taking out nans
            Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = []; % Taking out nans
            rainfall(idx_rainfall) = []; % Taking out nans
            plot3(Spatial_Rainfall_Parameters.x_coordinate, Spatial_Rainfall_Parameters.y_coordinate,zmax*ones(size(rainfall)), 'ko', 'MarkerSize', 8,'LineWidth',1.5, 'MarkerFaceColor','r' );
        end

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
    end
    % adjust the hydrograph legend
    lgd_hydro.Position = [0.2810, 0.65, 0.1205, 0.0958];
    %lgd_hydro.Position = [0.083, 0.55, 0.1205, 0.0958];
    try
        %exportgraphics(gcf,fullfile(myFolder_wd,'Specific_Discharge_Gauges_'+labels_depth{i}+'.pdf'),'ContentType','vector');
        exportgraphics(gcf,fullfile(myFolder_wd,'Specific_Discharge_Gauges_'+labels_depth{i}+'.png'),'ContentType','image','Colorspace','rgb','Resolution',1200);
    catch
        fprintf('Specific discharge gauges no exported, PDF export error')
    end
    saveas(gcf,fullfile(myFolder_wd,'Specific_Discharge_Gauges_'+labels_depth{i}+'.fig'))
    close all

end

%export the noramlized information    
FileName_String = 'normalized_table';    
FileName = fullfile(myFolder_wd,strcat('\',FileName_String,'.csv'));  
writetable(normalized_table,FileName);    
end    
   

