close all
%% Post-Processing - Graphs
simulation_time = toc;

%% Coloramp
[Spectrum,Depth_Purple,Terrain_RAS_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS,WSE_RAS] = coloramps();

%% Creating Modeling Results Folder
% Create the folder name
folderName = 'Modeling_Results';

try
    % If it doesn't exist, create the folder
    mkdir(folderName);
    disp('Folder "Modeling_Results" created successfully!');
catch ME
    disp('Impossible to create the folder for some reason');
end

%% Change directory if necessary
if flags.flag_rainfall_multiple_runs == 1
    input_rains = spreadsheetDatastore('rainfalls_cc.xlsx');
    mkdir(strcat('Modeling_Results','\',input_rains.VariableNames{Rainfall_Parameters.name_time}))
    folderName = strcat('Modeling_Results','\',input_rains.VariableNames{Rainfall_Parameters.name_time});
end


%% Time Data
if flags.flag_elapsed_time ~=1
    if flags.flag_spatial_rainfall == 1
        if ~isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration)
            Spatial_Rainfall_Parameters.rainfall_spatial_duration = double(Spatial_Rainfall_Parameters.rainfall_spatial_duration);
            Spatial_Rainfall_Parameters.rainfall_spatial_duration = Spatial_Rainfall_Parameters.rainfall_spatial_duration/60/24 + date_begin;
        end
    else
        if flags.flag_rainfall == 1
            try
                Rainfall_Parameters.time_rainfall =  double(Rainfall_Parameters.time_rainfall/60/24) + date_begin;
                Rainfall_Parameters.time_rainfall_saved = Rainfall_Parameters.time_rainfall;
                size_rain = time_rainfall >= Rainfall_Parameters.time_rainfall(1);
                Rainfall_Parameters.time_rainfall(size_rain~= 1) = [];
                if sum(size_rain) < length(Rainfall_Parameters.intensity_rainfall)
                    Rainfall_Parameters.intensity_rainfall(size_rain~= 1) = [];
                end
            catch
                Rainfall_Parameters.time_rainfall_saved = Rainfall_Parameters.time_rainfall;
                Rainfall_Parameters.time_rainfall = Rainfall_Parameters.time_rainfall_saved;
                size_rain = Rainfall_Parameters.time_rainfall >= Rainfall_Parameters.time_rainfall(1);
                Rainfall_Parameters.time_rainfall(size_rain~= 1) = [];
                if sum(size_rain) < length(Rainfall_Parameters.intensity_rainfall)
                    Rainfall_Parameters.intensity_rainfall(size_rain~= 1) = [];
                end
            end
        end
    end
    clear size_rain
    if flags.flag_ETP == 1
        if ~isdatetime(ETP_Parameters.climatologic_spatial_duration)
            ETP_Parameters.climatologic_spatial_duration = double(ETP_Parameters.climatologic_spatial_duration/60/24) + date_begin;
        end
    end
    try
        running_control.time_records = double(running_control.time_records/60/24) + date_begin;
    catch ME
        running_control.time_records = running_control.time_records;
    end

    try
        running_control.time_hydrograph =  double(running_control.time_hydrograph/60/24) + date_begin;
        running_control.time_hydrograph_save = running_control.time_hydrograph;
    catch
        running_control.time_hydrograph = running_control.time_hydrograph_save;
    end
    if flags.flag_inflow == 1 && flags.flag_elapsed_time ~= 1
        if ~isdatetime(Inflow_Parameters.time_inflow)
            Inflow_Parameters.time_inflow = double(Inflow_Parameters.time_inflow/60/24) + date_begin;
        elseif flags.flag_elapsed_time ==1
            Inflow_Parameters.time_inflow = Inflow_Parameters.time_inflow/60/24;
        end
    end
end


%% Outlet Hydrograph
if flags.flag_elapsed_time == 1
    line_plot(gather(running_control.time_hydrograph),'\mathrm{Elapsed~Time} ','min',gather(outlet_states.outlet_hydrograph),'\mathrm{Discharge} ','\mathrm{m^3 \cdot s^{-1}}',[],[],[],[],'Hydrograph',1,1);
else
    line_plot(gather(running_control.time_hydrograph),'\mathrm{Date} ','',gather(outlet_states.outlet_hydrograph),'\mathrm{Discharge} ','\mathrm{m^3 \cdot s^{-1}}',[],[],[],[],'Hydrograph',1,1);
end
yyaxis right; set(gca,'ydir','reverse','ycolor','black');
if flags.flag_rainfall == 1
    if flags.flag_spatial_rainfall ~=1
        Rainfall_Parameters.time_rainfall = Rainfall_Parameters.time_rainfall(1:length(Rainfall_Parameters.intensity_rainfall));
        bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',pallete.blue_colors(2,:),'EdgeColor',[0 .5 .5],'LineWidth',1.5)
        ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(gather(Rainfall_Parameters.intensity_rainfall)) + 200]);
    else
        if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
            dim = length(BC_States.average_spatial_rainfall);
            bar((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim))))',gather(BC_States.average_spatial_rainfall),'FaceColor',pallete.blue_colors(2,:),'EdgeColor',[0 .5 .5],'LineWidth',1.5);
            ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(BC_States.average_spatial_rainfall)*5]);
        else
            dim = length(BC_States.average_spatial_rainfall);
            bar((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim))))',gather(BC_States.average_spatial_rainfall),'FaceColor',pallete.blue_colors(2,:),'EdgeColor',[0 .5 .5],'LineWidth',1.5);
            ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(BC_States.average_spatial_rainfall)*5]);
        end
    end

end
% ylim([0 600]);
set(gcf,'units','inches','position',[3,3,6,4])
if flags.flag_inflow == 1
    hold on
    yyaxis left; set(gca,'ydir','normal','ycolor','black');
    if flags.flag_elapsed_time ~=1
        Inflow_Parameters.inflow_duration = minutes(Inflow_Parameters.time_inflow(end) - Inflow_Parameters.time_inflow(1));
    else
        Inflow_Parameters.inflow_duration = max(Inflow_Parameters.time_inflow);
    end
    if Inflow_Parameters.inflow_duration > running_control.routing_time
        if flags.flag_elapsed_time == 1
            tfinal_inflow = find(Inflow_Parameters.time_inflow >= running_control.routing_time,1,'first');
        else
            tfinal_inflow = find(Inflow_Parameters.time_inflow >= double(running_control.routing_time/60/24) + date_begin,1,'first');
        end
    else
        tfinal_inflow = length(Inflow_Parameters.time_inflow);
    end
    hold on
    plot(gather(Inflow_Parameters.time_inflow(1:tfinal_inflow,1)),gather(Inflow_Parameters.inflow_discharge(1:tfinal_inflow,:)),'LineWidth',1.5','color','black','LineStyle','--');
    ylim([0 max(max(max(gather(outlet_states.outlet_hydrograph))*1.2,1.5*max(max(gather(Inflow_Parameters.inflow_discharge)))))]);
    grid on
    if size(Inflow_Parameters.inflow_discharge,2) == 1
        legend('Outflow','Rainfall','Inflow','interpreter','latex');
    end
else
    legend('Outflow','Rainfall','interpreter','latex');
end
try
    exportgraphics(gcf,fullfile(folderName,'Hydrograph.pdf'),'ContentType','vector')
catch
    fprintf('No hydrograph exported, PDF export error')
end
saveas(gcf,fullfile(folderName,'Hydrograph.fig'))
close all

%% Stage and Hydrograph at the Outlet
if flags.flag_obs_gauges == 1
    set(gcf,'units','inches','position',[3,0,6.5,4])
    if gauges.num_obs_gauges == 1
        color_plot = linspecer(2);
    else
        color_plot = linspecer(gauges.num_obs_gauges);
    end
    font_size = 12;
    for i = 1:1
        fsize = 12;
        subplot(ceil(1),1,i)
        yyaxis left;
        set(gca,'ycolor',color_plot(1,:))
        plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph),'LineWidth',1.5','color',color_plot(1,:));
        xlabel('Elapsed Time [min]','Interpreter','latex');
        ylabel('Q $(\mathrm{m^3/s})$','interpreter','latex'); set(gca,'FontSize',12);
        ylabel('Q $(\mathrm{m^3/s})$','interpreter','latex');
        title('Outlet','FontName','Garamond')
        yyaxis right
        set(gca,'ycolor',color_plot(2,:))
        plot(gather(running_control.time_hydrograph),gather(outlet_states.depth_outlet/1000),'LineWidth',1.5','linestyle','--','color',color_plot(2,:));
        ylabel('h (m)','interpreter','latex')
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
        box on
    end
    try
        exportgraphics(gcf,fullfile(folderName,'Stage_Hydrograph_Outlet.pdf'),'ContentType','vector')
    catch
        fprintf('No Stage hydrographs exported, pdf export error')
    end
    saveas(gcf,fullfile(folderName,'Stage_Hydrograph_Outlet.fig'))
    close all
end


%% Stage Hydrograph %%
if flags.flag_obs_gauges == 1
    color_plots = linspecer(gauges.num_obs_gauges);
    set(gcf,'units','inches','position',[3,0,9,5])
    % Stage Hydrograph
    % plot(gather(running_control.time_hydrograph),gather(wse_outlet_2),'LineWidth',1.5','color','black'); xlabel('Time [min]','interpreter','latex'); ylabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex'); set(gca,'FontSize',12);
    for i = 1:gauges.num_obs_gauges
        if mod(i,3) == 0
            ls = '--';
        elseif mod(i,3) == 1
            ls = '-';
        else
            ls = ':';
        end
        set(gca,'ycolor','black')
        plot(gather(running_control.time_hydrograph),gather(gauges.depth_cell(:,i)),'LineWidth',1.5,'linestyle',ls,'color',color_plots(i,:));
        hold on
    end
    % Outlet
    plot(gather(running_control.time_hydrograph),gather(outlet_states.depth_outlet/1000),'LineWidth',1.5,'linestyle',ls,'marker','.','color','red');
    xlabel('Time ','interpreter','latex'); ylabel('Depth $(\mathrm{m})$','interpreter','latex'); set(gca,'FontSize',12);
    hold on
    yyaxis right; set(gca,'ydir','reverse','ycolor','black');
    if flags.flag_rainfall == 1
        if flags.flag_spatial_rainfall ~=1
            bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',[0 .55 .55],'EdgeColor',[0 .5 .5],'LineWidth',1.5)
            ylabel('Rainfall Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
        else
            ylabel('Mean Rainfall Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
            plot((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)))),gather(BC_States.average_spatial_rainfall),'LineWidth',1.5,'color','blue')
        end
    end
    labels_depth = gauges.labels_observed_string; labels_depth{gauges.num_obs_gauges+1} = 'Outlet';
    labels_depth{gauges.num_obs_gauges + 2} = 'Rainfall Intensity';
    legend(labels_depth,'FontName','Garamond','FontSize',8,'location','bestoutside')
    try
        ylim([0 max(max(gather(BC_States.average_spatial_rainfall)))*6]);
    catch
        ylim([0 10]);
    end
    title('Surface Runoff Depth','interpreter','latex','fontsize',12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca,'FontName','Garamond','FontSize',12)
    try
        exportgraphics(gcf,fullfile(folderName,'Stage_Hydrograph_Gauges.pdf'),'ContentType','vector')
    catch
        fprintf('No stage hydrograph gauges saved, PDF export error')
    end
    saveas(gcf,fullfile(folderName,'Stage_Hydrograph_Gauges.fig'))
    close all
end

%% Normalized Discharge %%
% Rainfall Std Deviation
zero_matrix = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2));
if flags.flag_spatial_rainfall == 1 && running_control.record_time_spatial_rainfall
    store=1;
    flag_loader=1;
    rainfall_sum = zeros(size(zero_matrix));
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
        z = Maps.Hydro.spatial_rainfall_maps(:,:,i - ((store-1)*saver_memory_maps));
        rainfall_sum = rainfall_sum + z; % mm/h
        Rainfall_Parameters.std_dev(i,1) = nanstd(z(:));
    end
    % for i = 1:size(Maps.Hydro.spatial_rainfall_maps,3)
    %     z = Maps.Hydro.spatial_rainfall_maps(:,:,i);
    %     Rainfall_Parameters.std_dev(i,1) = nanstd(z(:));
    % end
end

if flags.flag_obs_gauges == 1 && flags.flag_rainfall == 1
    % Catchment area of each gauge
    for i = 1:length(gauges.easting_obs_gauges)
        gauges.catchment_area(i,1) = Wshed_Properties.fac_area(gauges.northing_obs_gauges(i,1),gauges.easting_obs_gauges(i,1)); % km2
    end
    color_plots = linspecer(gauges.num_obs_gauges);
    set(gcf,'units','inches','position',[3,0,9,5])
    % Stage Hydrograph
    % plot(gather(running_control.time_hydrograph),gather(wse_outlet_2),'LineWidth',1.5','color','black'); xlabel('Time [min]','interpreter','latex'); ylabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex'); set(gca,'FontSize',12);
    for i = 1:gauges.num_obs_gauges
        if mod(i,3) == 0
            ls = '--';
        elseif mod(i,3) == 1
            ls = '-';
        else
            ls = ':';
        end
        set(gca,'ycolor','black')
        specific_discharge = gauges.hydrograph_cell(:,i)/gauges.catchment_area(i,1) ; % m3/s/km2
        plot(gather(running_control.time_hydrograph),specific_discharge,'LineWidth',1.5,'linestyle',ls,'color',color_plots(i,:));
        hold on
    end
    % Outlet
    plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph)/(Wshed_Properties.drainage_area/1000/1000),'LineWidth',1.5,'linestyle',ls,'marker','.','color','red');
    xlabel('Elapsed Time [min]','interpreter','latex'); ylabel('Specific Discharge $(\mathrm{m^3/s/km^2})$','interpreter','latex'); set(gca,'FontSize',12);
    hold on
    yyaxis right; set(gca,'ydir','reverse','ycolor','black');
    if flags.flag_rainfall == 1
        if flags.flag_spatial_rainfall ~=1
            bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',[0 .55 .55],'EdgeColor',[0 .5 .5],'LineWidth',1.5)
            ylabel('Rainfall Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
            ylim([0 max(gather(Rainfall_Parameters.intensity_rainfall))*6])
        else
            bar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),gather(BC_States.average_spatial_rainfall),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
            ylabel('Aerial Mean Rainfall Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
            hold on
            try
                er = errorbar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),BC_States.average_spatial_rainfall,Rainfall_Parameters.std_dev(1:(dim),1),Rainfall_Parameters.std_dev(1:(dim),1));
                er.Color = [0 0 0];
                er.LineStyle = 'none';
            end
            plot((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)))),gather(BC_States.average_spatial_rainfall),'LineWidth',1.5,'color','blue')
            ylim([0 max(max(gather(BC_States.average_spatial_rainfall)))*6])
        end
    end
    labels_depth = gauges.labels_observed_string; labels_depth{gauges.num_obs_gauges+1} = 'Outlet';
    labels_depth{gauges.num_obs_gauges + 2} = 'Rainfall Intensity';
    legend(labels_depth,'FontName','Garamond','FontSize',8,'location','bestoutside')

    title('Specific Discharge','interpreter','latex','fontsize',12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca,'FontName','Garamond','FontSize',12)
    try
        exportgraphics(gcf,fullfile(folderName,'Specific_Discharge_Gauges.pdf'),'ContentType','vector')
    catch
        fprintf('Specific discharge gauges no exported, PDF export error')
    end
    saveas(gcf,fullfile(folderName,'Specific_Discharge_Gauges.fig'))
    close all
end

%% Total ETR
zero_matrix = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2));
if flags.flag_ETP == 1
    store=1;
    flag_loader=1;
    ETR_sum = zeros(size(zero_matrix));
    for i = 1:length(running_control.time_records)
        try
            if i > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                    flag_loader=0;
                end
            end
            z = Maps.Hydro.spatial_rainfall_maps(:,:,i - ((store-1)*saver_memory_maps));
            z = Maps.Hydro.ETR_save(:,:,i - ((store-1)*saver_memory_maps));
            z(isnan(z)) = 0; % Attention here
            ETR_sum = ETR_sum + z; % mm/24h
            ETP_Parameters.std_dev_ETR(i,1) = nanstd(z(:));
        end
    end
end

%% Rating Curve - Outlet
close all
set(gcf,'units','inches','position',[3,3,6,4])
if size(Wshed_Properties.row_min,1) > 0
    % We have more than one outlet
    wse_outlet = Wshed_Properties.el_outlet(Wshed_Properties.row_min(1),Wshed_Properties.col_min(1)) + outlet_states.depth_outlet/1000;
    s = scatter(outlet_states.outlet_hydrograph,wse_outlet,'o','b'); xlabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex');
else
    % We have only one outlet
    wse_outlet = Wshed_Properties.el_outlet(Wshed_Properties.row_min,Wshed_Properties.col_min) + outlet_states.depth_outlet/1000; % wse
    s = scatter(gather(outlet_states.outlet_hydrograph),gather(wse_outlet),'o','b'); xlabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex');
end
ylabel('Water Surface Elevation $(\mathrm{m})$','interpreter','latex'); set(gca,'FontSize',12);
if max(wse_outlet) == Wshed_Properties.stage_min
    ylim([Wshed_Properties.stage_min Wshed_Properties.stage_min + 1])
else
    ylim([Wshed_Properties.stage_min max(wse_outlet)])
end
set(gca,'FontName','Garamond')
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
s.SizeData = 10;
box on
try
    exportgraphics(gcf,fullfile(folderName,'Rating_Curve_Outlet.pdf'),'ContentType','vector')
catch
    fprintf('No rating curge outlet export, PDF export error')
end
saveas(gcf,fullfile(folderName,'Rating_Curve_Outlet.fig'))
Rating_Curve_Data = table(gather(outlet_states.outlet_hydrograph),gather(wse_outlet),gather(outlet_states.depth_outlet)/1000,'VariableNames',{'Flow Discharge (m3/s)','WSE (m)','Depth (m)'});

FileName_String = 'Rating_Curve_Data';
FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));
writetable(Rating_Curve_Data,FileName);

% Hydrograph table
Outlet_Hydrograph_Data = table(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph),'VariableNames',{'Time ','Flow Discharge $(\mathrm{m^3/s})$'});

FileName_String = 'Outlet_Hydrograph_Data_Outlet';
FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));
writetable(Outlet_Hydrograph_Data,FileName);

% writetable(Outlet_Hydrograph_Data)
close all

%% Rating Curve - Specific Cell
if flags.flag_obs_gauges == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    color_plot = linspecer(gauges.num_obs_gauges);
    for i = 1:gauges.num_obs_gauges
        if gauges.num_obs_gauges > 5
            fsize = 8;
        else
            fsize = 12;
        end
        if gauges.num_obs_gauges > 3
            subplot(ceil(gauges.num_obs_gauges/3),3,i)
        elseif gauges.num_obs_gauges == 2
            subplot(ceil(gauges.num_obs_gauges/2),2,i)
        elseif gauges.num_obs_gauges == 1
            subplot(1,1,1)
        end
        s = scatter(gather(gauges.hydrograph_cell(:,i)),gather(gauges.wse_cell(:,i)),'o'); xlabel('Q $(\mathrm{m^3/s})$','interpreter','latex');
        s.LineWidth = 0.6;
        %     s.MarkerEdgeColor = 'b';
        s.MarkerFaceColor = color_plot(i,:);
        s.SizeData = 10;
        hold on
        ylabel('WSE $(\mathrm{m})$','interpreter','latex'); set(gca,'FontSize',fsize);
        title(gauges.labels_observed_string{i},'FontName','Garamond')
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca, 'FontName', 'Garamond', 'FontSize', fsize)
        box on
    end
    try
        exportgraphics(gcf,fullfile(folderName,'Rating_Curve_Gauges.pdf'),'ContentType','vector')
    catch
        fprintf('No rating curge gaugest export, PDF export error')
    end
    saveas(gcf,fullfile(folderName,'Rating_Curve_Gauges.fig'))
    Rating_Curve_Specifc_Cell = table(gather(running_control.time_hydrograph), gather(gauges.hydrograph_cell),gather(gauges.wse_cell),gather(gauges.depth_cell),'VariableNames',{'Time [min] or Date','Flow Discharge (m3/s)','WSE (m)','Water Depth (m)'});

    FileName_String = 'Rating_Curve_Gauges';
    FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));
    writetable(Rating_Curve_Specifc_Cell,FileName);


    % writetable(Rating_Curve_Specifc_Cell)
end
close all

%% Hydrographs - Specific Cell
if flags.flag_obs_gauges == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    if gauges.num_obs_gauges == 1
        color_plot = linspecer(10);
    else
        color_plot = linspecer(gauges.num_obs_gauges);
    end
    font_size = 12;
    for i = 1:gauges.num_obs_gauges
        if gauges.num_obs_gauges > 5
            fsize = 8;
        else
            fsize = 12;
        end
        if gauges.num_obs_gauges > 3
            subplot(ceil(gauges.num_obs_gauges/3),3,i)
        elseif gauges.num_obs_gauges == 2
            subplot(ceil(gauges.num_obs_gauges/2),2,i)
        elseif gauges.num_obs_gauges == 1
            subplot(1,1,1)
        end
        yyaxis left;
        set(gca,'ycolor',pallete.blue_colors(1,:))
        plot(gather(running_control.time_hydrograph),gather(gauges.hydrograph_cell(:,i)),'LineWidth',1.5','color',pallete.blue_colors(2,:));
        xlabel('Elapsed Time [min]','Interpreter','latex');
        ylabel('Q $(\mathrm{m^3/s})$','interpreter','latex'); set(gca,'FontSize',12);
        title(gauges.labels_observed_string{i},'FontName','Garamond')
        yyaxis right
        set(gca,'ycolor',color_plot(2,:))
        plot(gather(running_control.time_hydrograph),gather(gauges.depth_cell(:,i)),'LineWidth',1.5','linestyle','--','color',pallete.red_colors(3,:));
        ylabel('h (m)','interpreter','latex')
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
        box on
    end
    try
        exportgraphics(gcf,fullfile(folderName,'Hydrograph_Gauges.pdf'),'ContentType','vector')
    catch
        fprintf('No Hydrograph gauges export, PDF export error')
    end
    saveas(gcf,fullfile(folderName,'Hydrograph_Gauges.fig'))
end
close all

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
% Retriving the basemap imageshapefile
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
if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= web_mercator_crs.Name
    no_plot = 1;
else
    S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
    S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)]; % Calculate bounding box for each polygon
    S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
    S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
end

%% DEM with Streams and Obseved Points
FD = FLOWobj(DEM_raster);
area_km2 = GIS_data.min_area; % km2
area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
if no_plot==0
    try
        mapshow(A,RA,"AlphaData",0.45);hold on;
        mapshow(S_p,'FaceColor','n'); hold on;
    catch ME
        warning('You need matlab 2022 or higher to run basemaps')
    end
end

S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
ax1 = plot(S);
hold on
title(('Streams and Obseved Points'),'Interpreter','Latex','FontSize',12)
xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
ax = ancestor(ax1, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');
if flags.flag_obs_gauges == 1
    scatter(gauges.easting_obs_gauges_absolute,gauges.northing_obs_gauges_absolute)
    MS = STREAMobj2mapstruct(S);

    FileName_String = 'streamnetwork.shp';
    FileName = fullfile(folderName,strcat('\',FileName_String));
    shapewrite(MS,FileName);


    % Example x and y coordinates
    gauges.x_coord_gauges = gauges.easting_obs_gauges_absolute;
    gauges.y_coord_gauges = gauges.northing_obs_gauges_absolute;
end
if flags.flag_obs_gauges == 1
    labels_gauges = gauges.labels_observed_string;  % Labels for each point


    % Create a shapefile writer
    % shapefile = shapewrite('observed_gauges.shp');

    % Create a geoshape object with the points and labels
    points = mappoint(gauges.x_coord_gauges', gauges.y_coord_gauges', 'Label', labels_gauges');

    % Write the geoshape object to a shapefile
    FileName_String = 'observed_gauges.shp';
    FileName = fullfile(folderName,strcat('\',FileName_String));
    shapewrite(points, FileName); % This is in the same coordinate system of your DEM
end

close all

%% Water Quality Analysis
if flags.flag_waterquality == 1
    plot(gather(running_control.time_hydrograph),gather(WQ_States.outet_pollutograph),'LineWidth',1.5,'color','Red','Marker','o') ;
    xlabel('Time [min]','Interpreter','Latex');
    ylabel('Pol. Concentration (mg/L)','Interpreter','Latex')
    grid on
    hold on
    yyaxis right
    set(gca,'ycolor','black')
    ylabel('Load (kg/sec)','interpreter','latex')
    load_wq = 10^(-3)*gather(WQ_States.outet_pollutograph).*gather(outlet_states.outlet_hydrograph); % kg/sec
    plot(gather(running_control.time_hydrograph),load_wq,'LineWidth',1.5,'color','Blue','Marker','*')
    legend('Concentration','Load','interpreter','latex')
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
    box on
    exportgraphics(gcf,fullfile(folderName,'Pollutograph.pdf'),'ContentType','vector')
    saveas(gcf,fullfile(folderName,'Pollutograph.fig'))
    Outlet_Pollutograph_Data = table(gather(running_control.time_hydrograph),gather(WQ_States.outet_pollutograph),gather(load_wq),'VariableNames',{'Time [min]','Concentration (mg/L)','load (kg/sec)'});

    FileName_String = 'Outlet_Pollutograph_Data';
    FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));
    writetable(Outlet_Pollutograph_Data,FileName);

    %     writetable(Outlet_Pollutograph_Data)
    close all
    %% Histeresis Effect
    set(gcf,'units','inches','position',[3,3,6,4])
    plot(gather(running_control.time_hydrograph),outlet_states.outlet_hydrograph,'color','black','linewidth',1.5,'marker','*');
    xlabel('Time [min]','Interpreter','Latex');
    ylabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex')
    yyaxis right
    set(gca,'ycolor','black')
    plot(gather(running_control.time_hydrograph),gather(WQ_States.outet_pollutograph),'LineWidth',1.5,'color','Red','Marker','o') ;
    ylabel('Pol. Concentration (mg/L)','Interpreter','Latex')
    legend('Flow discharge','Pollutant Concentration','interpreter','latex')
    exportgraphics(gcf,fullfile(folderName,'Histeresis.pdf'),'ContentType','vector')
    saveas(gcf,fullfile(folderName,'Histeresis.fig'))
    close all
    % M(v) Curves
    m_M = WQ_States.mass_outlet_save./(max(max(WQ_States.mass_outlet)));
    v_V = WQ_States.vol_outlet_save./(max(max(WQ_States.vol_outlet)));
    plot(v_V,m_M,'LineWidth',1.5,'color','black','Marker','*');
    hold on
    plot(0:1,0:1,'LineWidth',1,'Color','black','LineStyle','--')
    ylabel('$m/m_{tot}$','Interpreter','Latex');
    xlabel('$V/V_{tot}$','Interpreter','Latex');
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
    box on
    exportgraphics(gcf,fullfile(folderName,'M(V)_Curve.pdf'),'ContentType','vector')
    saveas(gcf,fullfile(folderName,'M(V)_Curve.fig'))
    Outlet_M_V_Curve = table(gather(round(v_V,3)),gather(round(m_M,3)),'VariableNames',{'Normalized Volume','Normalized Polutant Mass'});

    FileName_String = 'Outlet_M_V_Curve';
    FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));
    writetable(Outlet_M_V_Curve,FileName);

    %     writetable(Outlet_M_V_Curve);

    close all

    %% EMC(Curve)
    plot(gather(running_control.time_hydrograph),WQ_States.EMC_outlet,'LineWidth',1.5,'color','black');
    ylabel('$\mathrm{EMC (mg/L)}$','Interpreter','Latex');
    xlabel('Elapsed Time [min]','Interpreter','Latex');
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
    box on
    exportgraphics(gcf,fullfile(folderName,'EMC Curve.pdf'),'ContentType','vector')
    saveas(gcf,fullfile(folderName,'EMC_Curve.fig'))
    Outlet_EMC_Curve = table(gather(running_control.time_hydrograph,3),gather(round(WQ_States.EMC_outlet,2)),'VariableNames',{'Normalized Volume','Normalized Polutant Mass'});

    FileName_String = 'Outlet_EMC_Curve';
    FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));
    writetable(Outlet_EMC_Curve,FileName);

    %     writetable(Outlet_EMC_Curve);
end
%% Exporting Rasters
if flags.flag_export_maps == 1
    % Delete previous rasters in the folder
    no_data_value = nan;

    % Create their own folder
    if flags.flag_waterquality==1
        mkdir(strcat(folderName,'\Water_Quality_Maps'));
        % Specify the folder where the files live.
        myFolder_wq = strcat(pwd,'\',folderName,'\Water_Quality_Maps\'); % Current folder
        % Get a list of all files in the folder with the desired file name pattern.
        filePattern = fullfile(myFolder_wq, '*.tif'); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        for k = 1 : length(theFiles)
            baseFileName = theFiles(k).name;
            fullFileName = fullfile(myFolder_wq, baseFileName);
            fprintf(1, 'Now deleting %s\n', fullFileName);
            delete(fullFileName);
        end
        filePattern = fullfile(myFolder_wq, '*.tif.aux'); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        for k = 1 : length(theFiles)
            baseFileName = theFiles(k).name;
            fullFileName = fullfile(myFolder_wq, baseFileName);
            fprintf(1, 'Now deleting %s\n', fullFileName);
            delete(fullFileName);
        end
    elseif flags.flag_human_instability>0
        % Specify the folder where the files live.
        myFolder_hr = strcat(pwd,'\',folderName,'\Human_Risk_Maps\'); % Current folder
        mkdir(strcat(folderName,'\Human_Risk_Maps'));
        % Get a list of all files in the folder with the desired file name pattern.
        filePattern = fullfile(myFolder_hr, '*.tif'); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        for k = 1 : length(theFiles)
            baseFileName = theFiles(k).name;
            fullFileName = fullfile(myFolder_hr, baseFileName);
            fprintf(1, 'Now deleting %s\n', fullFileName);
            delete(fullFileName);
        end
        filePattern = fullfile(myFolder_hr, '*.tif.aux'); % Change to whatever pattern you need.
        theFiles = dir(filePattern);
        for k = 1 : length(theFiles)
            baseFileName = theFiles(k).name;
            fullFileName = fullfile(myFolder_hr, baseFileName);
            fprintf(1, 'Now deleting %s\n', fullFileName);
            delete(fullFileName);
        end
    end

    mkdir(strcat(folderName,'\Water_Depths_Maps'));
    % Specify the folder where the files live.
    myFolder_wd = strcat(pwd,'\',folderName,'\Water_Depths_Maps\'); % Current folder
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isfolder(myFolder_wd)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder_wd);
        uiwait(warndlg(errorMessage));
        return;
    end
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolder_wd, '*.tif'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(myFolder_wd, baseFileName);
        fprintf(1, 'Now deleting %s\n', fullFileName);
        delete(fullFileName);
    end
    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolder_wd, '*.tif.aux'); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(myFolder_wd, baseFileName);
        fprintf(1, 'Now deleting %s\n', fullFileName);
        delete(fullFileName);
    end
    clc
    % Changing Nan Values
    % Maps.Hydro.d(isnan(Maps.Hydro.d)) = no_data_value;
    % Maps.Hydro.d(isinf(Maps.Hydro.d)) = no_data_value;
    % Maps.Hydro.d(idx_nan) = no_data_value;
    % Rasters - Depths, WSE, and Pol. Conc
    store=1;
    flag_loader=1;
    flag_loader_wq = 1;
    for i = 1:length(running_control.time_records)
        raster_exportion_percentage = i/length(running_control.time_records)*100
        %         flags.flag_elapsed_time = 1;

        if i > saver_memory_maps*store
            store = store + 1;
            load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            % Changing Nan Values
            Maps.Hydro.d(isnan(Maps.Hydro.d)) = no_data_value;
            Maps.Hydro.d(isinf(Maps.Hydro.d)) = no_data_value;
            Maps.Hydro.d(idx_nan) = no_data_value;
            Max_depth_d = max(max(Maps.Hydro.d,[],3),Max_depth_d);
        else
            if flag_loader == 1
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                flag_loader=0;
                % Changing Nan Values
                Maps.Hydro.d(isnan(Maps.Hydro.d)) = no_data_value;
                Maps.Hydro.d(isinf(Maps.Hydro.d)) = no_data_value;
                Maps.Hydro.d(idx_nan) = no_data_value;
                Max_depth_d = max(Maps.Hydro.d,[],3);
            end
        end

        idx_depth = Maps.Hydro.d(:,:,i - ((store-1)*saver_memory_maps)) < depths.depth_wse*1000;

        % idx_depth = Maps.Hydro.d(:,:,i) < depths.depth_wse*1000;
        if flags.flag_elapsed_time ~= 1
            time_map = datestr(running_control.time_records(i),'yyyy_mm_dd_hh_MM_ss');
        else
            time_map = running_control.time_records(i)/60; % hours
        end
        if flags.flag_wse == 0  % We will save only water surface depth
            if flags.flag_elapsed_time == 1
                FileName = strcat('Flood_Depths_t_', num2str(time_map),'_h.tif');
            else
                FileName = strcat('Flood_Depths_', string(time_map));
            end
            raster_exportion = Maps.Hydro.d(:,:,i - ((store-1)*saver_memory_maps))/1000;
            % raster_exportion = Maps.Hydro.d(:,:,i)/1000;
            raster_exportion(idx_nan) = no_data_value;
            raster_exportion(idx_depth) = no_data_value;
        else
            if flags.flag_elapsed_time == 1
                FileName =  strcat('Water_Surface_Elevation_t_', num2str(time_map),'_h.tif');
            else
                FileName = strcat('Water_Surface_Elevation_', string(time_map));
            end
            raster_exportion = Maps.Hydro.d(:,:,i - ((store-1)*saver_memory_maps))/1000 + double(idx_Elevation_Properties.elevation_cell).*Elevation_Properties.elevation_cell;
            % raster_exportion = Maps.Hydro.d(:,:,i)/1000 + double(idx_Elevation_Properties.elevation_cell).*Elevation_Properties.elevation_cell;
            raster_exportion(idx_nan) = no_data_value;
            raster_exportion(idx_depth) = no_data_value;
        end
        % Save Map
        FileName = fullfile(myFolder_wd,FileName);
        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
        %         SaveAsciiRaster(raster_exportion,Resolution,xllcorner,yllcorner,FileName,no_data_value)
        geotiffwrite(FileName,raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
        % Water Quality
        if flags.flag_waterquality == 1
            if flags.flag_elapsed_time == 1
                FileName =  strcat('Pollutant_Concentration', num2str(time_map),'min');
            else
                FileName = strcat('Pollutant_Concentration', string(time_map));
            end
            FileName = fullfile(myFolder_wq,FileName);

            raster_exportion = Maps.WQ_States.Pol_Conc_Map(:,:,i - (store-1)*saver_memory_maps);
            % raster_exportion = Maps.WQ_States.Pol_Conc_Map(:,:,i);
            idx_ = raster_exportion < LULC_Properties.Pol_min; % Finding values below Pol_min
            raster_exportion(idx_) = no_data_value;
            raster_exportion(isnan(raster_exportion)) = no_data_value;
            raster_exportion(isinf(raster_exportion)) = no_data_value;
            raster_exportion(raster_exportion < 0) = no_data_value;

            raster_to_export = DEM_raster; % Just to get the properties
            raster_to_export.Z = raster_exportion; % Putting the right values
            % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
            geotiffwrite(FileName,raster_to_export.Z,raster_to_export.georef.SpatialRef,...
                'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
        end
        if flags.flag_human_instability == 2
        elseif flags.flag_human_instability == 3
            list={'_cm','_tm','_am','_om','_cf','_tf','_af','_of'};
            if flags.flag_elapsed_time == 1
                FileName =  strcat('Human_instability_', num2str(time_map),'min');
            else
                FileName = strcat('Human_instability_', string(time_map));
            end
            FileName = fullfile(myFolder_hr,FileName);

            raster_exportion = zeros(size(DEM_raster.Z,1),size(DEM_raster.Z,2),8);
            for j = 1:8
                raster_exportion(:,:,j) = double(Maps.Hydro.(strcat('risk',list{j}))(:,:,i - (store-1)*saver_memory_maps)>0)*Human_Instability.order(j);
            end
            raster_exportion = max(raster_exportion,[],3);
            raster_exportion(isnan(raster_exportion)) = no_data_value;
            raster_exportion(isinf(raster_exportion)) = no_data_value;
            raster_exportion(raster_exportion < 0) = no_data_value;
            raster_exportion(raster_exportion==0) = no_data_value;
            geotiffwrite(FileName,raster_exportion,raster_to_export.georef.SpatialRef,...
                'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
        end
    end

    % Initial_Buildup Map
    if flags.flag_waterquality == 1
        FileName = 'Initial_Buildup_kg';
        FileName = fullfile(folderName,FileName);
        raster_exportion = Maps.WQ_States.initial_buildup_map;
        raster_exportion(isnan(raster_exportion)) = no_data_value;
        raster_exportion(isinf(raster_exportion)) = no_data_value;
        raster_exportion(raster_exportion < 0) = no_data_value;

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

        % Map of Total Washed Mass
        FileName = 'Total_Washed_Mass_Kg';
        FileName = fullfile(folderName,FileName);
        raster_exportion = WQ_States.Tot_Washed;
        raster_exportion(isnan(raster_exportion)) = no_data_value;
        raster_exportion(isinf(raster_exportion)) = no_data_value;
        raster_exportion(raster_exportion < 0) = no_data_value;

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
    end

    % Points of accumulation of Depths
    FileName = strcat('Accumulation_Areas_Larger_1m');
    FileName = fullfile(folderName,FileName);
    idx_depth = depths.d_t > 1*1000; % Larger than 1 m
    raster_exportion = no_data_value*ones(size(depths.d_t));
    raster_exportion(idx_depth) = 1;

    raster_to_export = DEM_raster; % Just to get the properties
    raster_to_export.Z = raster_exportion; % Putting the right values
    % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
    geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
        'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

    % Points of accumulation of pollutants
    if flags.flag_waterquality == 1
        FileName = strcat('Accumulation_Areas_10_g_m2');
        FileName = fullfile(folderName,FileName);
        pol_accumulation = 10; % g/m2
        zzz = WQ_States.B_t/Wshed_Properties.cell_area*1000; % g/m2
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;
        idx_bt = zzz < pol_accumulation; %
        zzz(idx_bt) = no_data_value;
        raster_exportion = zzz;

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

        % Final Polutant Mass
        FileName = strcat('Final_Pollutant_Mass_g_m2');
        FileName = fullfile(folderName,FileName);
        raster_exportion = no_data_value*ones(size(WQ_States.B_t));
        final_mass = Maps.WQ_States.Pol_mass_map(:,:,end)/Wshed_Properties.cell_area*1000; % g/m2
        idx_bt = final_mass > 0;
        raster_exportion(idx_bt) = final_mass(idx_bt);
        raster_exportion(isinf(raster_exportion)) = no_data_value;

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
    end
    % Maximum Velocity
    zzz = velocities.vmax_final; % m/s
    idx_wse = depths.dmax_final/1000 < depths.depth_wse; % Finding values below the threshold
    % Maximum_Depths
    FileName = strcat('Maximum_Velocity');
    FileName = fullfile(folderName,FileName);
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;
    zzz(idx_wse) = no_data_value;
    raster_exportion = zzz;
    raster_to_export = DEM_raster; % Just to get the properties
    raster_to_export.Z = raster_exportion; % Putting the right values
    % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
    geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
        'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

    % Risk Map
    if flags.flag_human_instability == 1
        store=1;
        flag_loader=1;
        for i = 1:size(length(running_control.time_records))
            if i > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                zzz = max(max(Maps.Hydro.risk,[],3),zzz); % m/s
            else
                if flag_loader == 1
                    load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                    flag_loader=0;
                    zzz = max(Maps.Hydro.risk,[],3); % m/s
                end
            end
        end
        idx_wse = depths.dmax_final/1000 < depths.depth_wse; % Finding values below the threshold
        % Maximum_Depths
        zzz = Human_Instability.max_risk; % Value of max risk        
        FileName = strcat('Maximum_Instability_Risk');
        FileName = fullfile(folderName,FileName);
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;
        zzz(idx_wse) = no_data_value;
        raster_exportion = zzz;
        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
    elseif flags.flag_human_instability == 2
    elseif flags.flag_human_instability == 3
        list={'_cm','_tm','_am','_om','_cf','_tf','_af','_of'};
        risk_summary = table('Size', [1 9], 'VariableNames',{'Risk','risk_cm', 'risk_tm', 'risk_am', 'risk_om', 'risk_cf', 'risk_tf', 'risk_af', 'risk_of'}, ...
            'VariableTypes', {'string','double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'});
        risk_summary.Risk(1) = 'Slide'; risk_summary.Risk(2) = 'Toppling'; risk_summary.Risk(3) = 'Drawing';
        for k = 1:8
            store=1;
            flag_loader=1;
            for i = 1:length(running_control.time_records)
                if i > saver_memory_maps*store
                    store = store + 1;
                    load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                else
                    if flag_loader == 1
                        load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                    end
                end
                zzz = double(Maps.Hydro.(strcat('risk',list{k}))(:,:,i - (store-1)*saver_memory_maps));
                if flag_loader ==1
                    zzz_2 = max(zzz,[],3);
                    flag_loader=0;
                else
                    zzz_2 = max(zzz,zzz_2);
                end
                risk_summary.(strcat('risk',list{k}))(1) = max(risk_summary.(strcat('risk',list{k}))(1),sum(sum(sum(zzz==1))));
                risk_summary.(strcat('risk',list{k}))(2) = max(risk_summary.(strcat('risk',list{k}))(2),sum(sum(sum(zzz==2))));
                risk_summary.(strcat('risk',list{k}))(3) = max(risk_summary.(strcat('risk',list{k}))(3),sum(sum(sum(zzz==3))));
            end
            zzz=zzz_2;

            % Marcus Edit
            zzz = Human_Instability.max_risk; % Value of max risk
            idx_wse = depths.dmax_final/1000 < depths.depth_wse; % Finding values below the threshold
            % Maximum_Depths
            FileName = strcat(strcat('Maximum_Instability_Risk',list{k}));
            FileName = fullfile(folderName,FileName);
            zzz(isinf(zzz)) = no_data_value;
            zzz(isnan(zzz)) = no_data_value;
            % zzz(idx_wse) = no_data_value;
            raster_exportion = zzz;
            raster_to_export = DEM_raster; % Just to get the properties
            raster_to_export.Z = raster_exportion; % Putting the right values
            geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
                'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
        end
        risk_summary(:,2:9)=risk_summary(:,2:9).*power(Wshed_Properties.Resolution,2)./1000000;
        writetable(risk_summary,strcat(folderName,"\",'Risk_summary.txt'),'Delimiter',',');
    end

    zzz = Soil_Properties.I_t; % m/s
    % Maximum_Depths
    FileName = strcat('Infiltrated_Depth');
    FileName = fullfile(folderName,FileName);
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;
    raster_exportion = zzz;
    raster_to_export = DEM_raster; % Just to get the properties
    raster_to_export.Z = raster_exportion; % Putting the right values
    % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
    geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
        'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)


    % Maximum
    zzz = depths.dmax_final/1000; % m
    idx_wse = zzz < depths.depth_wse; % Finding values below the threshold
    if flags.flag_wse == 0 % Save only max depth
        % Maximum_Depths
        FileName = strcat('Maximum_Depths');
        FileName = fullfile(folderName,FileName);
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;
        zzz(idx_wse) = no_data_value;
        raster_exportion = zzz;

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

    else
        raster_exportion = zzz + Elevation_Properties.elevation_cell; % Surface elevation
        raster_exportion(idx) = no_data_value;
        FileName = strcat('Max_Water_Surface_Elevation');
        FileName = fullfile(folderName,FileName);

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
    end

    % DEM
    if flags.flag_resample == 1
        FileName = strcat('DEM_resampled');
    else
        FileName = strcat('DEM_Treated');
    end
    FileName = fullfile(folderName,FileName);

    zzz = Elevation_Properties.elevation_cell;
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;
    raster_exportion = zzz;

    raster_to_export = DEM_raster; % Just to get the properties
    raster_to_export.Z = raster_exportion; % Putting the right values
    % Exporting the Map
    geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
        'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

    % LULC
    FileName = strcat('Land_Cover_Data');
    FileName = fullfile(folderName,FileName);
    zzz = LULC_Properties.LULC;
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;
    raster_exportion = zzz;

    raster_to_export = DEM_raster; % Just to get the properties
    raster_to_export.Z = raster_exportion; % Putting the right values
    % Exporting the Map
    geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
        'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)

    if flags.flag_waterquality == 1
        store=1;
        flag_loader=1;
        for i = 1:length(running_control.time_records)
            if i == length(running_control.time_records)
                % Water Quality - Mass of Pollutant
                time_map = running_control.time_records(i);
                FileName = 'Final_Mass_Of_Pollutant';
                FileName = fullfile(folderName,FileName);

                if i > saver_memory_maps*store
                    store = store + 1;
                    load(strcat('Temporary_Files\save_map_hydro',num2str(store)),'Maps');
                    Max_Pol_Conc_Map = max(max(Maps.WQ_States.Pol_Conc_Map,[],3),Max_Pol_Conc_Map);
                else
                    if flag_loader == 1
                        load(strcat('Temporary_Files\save_map_hydro',num2str(store)),'Maps');
                        flag_loader=0;
                        Max_Pol_Conc_Map = max(Maps.WQ_States.Pol_Conc_Map,[],3);
                    end
                end
                raster_exportion = Maps.WQ_States.Pol_mass_map(:,:,i - (store-1)*saver_memory_maps);
                % raster_exportion = Maps.WQ_States.Pol_mass_map(:,:,i);
                idx = raster_exportion < 0.01; % Finding values below Pol_min
                raster_exportion(idx) = no_data_value;
                raster_exportion = raster_exportion/Wshed_Properties.cell_area*1000; % g/m2

                raster_to_export = DEM_raster; % Just to get the properties
                raster_to_export.Z = raster_exportion; % Putting the right values
                % Exporting the Map
                geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
                    'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
            end
        end
        % Maximum
        FileName = strcat('Maximum_Pol_Conc','min');
        FileName = fullfile(folderName,FileName);

        % raster_exportion = Max_Pol_Conc_Map;
        raster_exportion = Max_Pol_Conc_Map;
        idx = raster_exportion < LULC_Properties.Pol_min; % Finding values below Pol_min
        raster_exportion(idx) = no_data_value;
        raster_exportion(isnan(raster_exportion)) = no_data_value;
        raster_exportion(isinf(raster_exportion)) = no_data_value;
        raster_exportion(raster_exportion < 0) = no_data_value;

        raster_to_export = DEM_raster; % Just to get the properties
        raster_to_export.Z = raster_exportion; % Putting the right values
        % GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
        geotiffwrite(strcat(cd,"\",FileName),raster_to_export.Z,raster_to_export.georef.SpatialRef,...
            'GeoKeyDirectoryTag',raster_to_export.georef.GeoKeyDirectoryTag)
        % Final Pollutant Mass
    end
end
%% Generate GIF Files of the Simulation
Inundation_Maps

%% Summary Table
%%% - Equivalent Width to Simulate in SWMM %%% - (LENHS, 2012)
% W = kc sqrt(A) / 1.12 * (1 - sqrt(1 - {1.128 / k_c}^2))
Wshed_Properties.width_SWMM = Wshed_Properties.compactness_coefficient*sqrt(Wshed_Properties.drainage_area)/1.12*(1 - sqrt(1 - (1.128/Wshed_Properties.compactness_coefficient)^2));
%%% - Runoff Coefficient -
Wshed_Properties.C_r = BC_States.outflow_volume/BC_States.inflow_volume;
%%% - Rainfall Volume
if flags.flag_spatial_rainfall ~=1 && flags.flag_rainfall == 1
    rainfall_vol = sum(sum(BC_States.delta_p));
elseif flags.flag_spatial_rainfall == 1 && flags.flag_rainfall == 1
    % tot_rain = rainfall_sum*running_control.record_time_spatial_rainfall/60/1000*Wshed_Properties.cell_area; % m3 for each cell
    tot_rain = sum(Maps.Hydro.spatial_rainfall_maps,3)*running_control.record_time_spatial_rainfall/60/1000*Wshed_Properties.cell_area; % m3 for each cell
    rainfall_vol = nansum(nansum(tot_rain))/Wshed_Properties.drainage_area*1000; % mm for the whole catchment
end

if flags.flag_waterquality == 0
    if flags.flag_rainfall == 0
        rainfall_vol = 0;
    end
    % Summary_Table = table(round(Wshed_Properties.drainage_area/1000/1000,3),rainfall_vol,round(Wshed_Properties.C_r,2),round(Wshed_Properties.impervious_area/1000/1000,3),round(Wshed_Properties.compactness_coefficient,3),round(Wshed_Properties.circularity_index,3),round(Wshed_Properties.width_SWMM,3),round(Wshed_Properties.form_factor,3),round(simulation_time/60,3),round(1/1000*max(max(Max_depth_d)),round(max(outlet_states.outlet_hydrograph),4),...
    Summary_Table = table(round(Wshed_Properties.drainage_area/1000/1000,3),rainfall_vol,round(Wshed_Properties.C_r,2),round(Wshed_Properties.impervious_area/1000/1000,3),round(Wshed_Properties.compactness_coefficient,3),round(Wshed_Properties.circularity_index,3),round(Wshed_Properties.width_SWMM,3),round(Wshed_Properties.form_factor,3),round(simulation_time/60,3),round(1/1000*max(max(max(Maps.Hydro.d))),2),round(max(outlet_states.outlet_hydrograph),4),...
        'VariableNames',...
        {'Drainage area (km2)','Rainfall Vol (mm)','Runoff Coefficient','Impervious Area (km2)','Compactness Coefficient','Circularity index','Equivalent Width (m)','Form Factor', ...
        'Simulation time (minutes)','Maximum Flood Depth (m)'...
        'Maximum Outflow (m^3/s)'});
else
    % Summary_Table = table(round(Wshed_Properties.drainage_area/1000/1000,3),rainfall_vol,round(Wshed_Properties.C_r,2),round(Wshed_Properties.impervious_area/1000/1000,3),round(Wshed_Properties.compactness_coefficient,3),round(Wshed_Properties.circularity_index,3),round(Wshed_Properties.width_SWMM,3),round(Wshed_Properties.form_factor,3),round(simulation_time/60,3),round(1/1000*max(max(Max_depth_d)),round(max(outlet_states.outlet_hydrograph),4),...
    % round(max(max(Max_Pol_Conc_Map))),1000*(1/Wshed_Properties.cell_area)*round(max(max(WQ_States.B_t(~isinf(WQ_States.B_t)))),4),round(initial_mass/1000,4),round(sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t))))/1000,4),1-round(sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t))))/1000,4)/round(initial_mass/1000,4),round(WQ_States.EMC_outlet(end,1),2),'VariableNames',...
    Summary_Table = table(round(Wshed_Properties.drainage_area/1000/1000,3),rainfall_vol,round(Wshed_Properties.C_r,2),round(Wshed_Properties.impervious_area/1000/1000,3),round(Wshed_Properties.compactness_coefficient,3),round(Wshed_Properties.circularity_index,3),round(Wshed_Properties.width_SWMM,3),round(Wshed_Properties.form_factor,3),round(simulation_time/60,3),round(1/1000*max(max(max(Maps.Hydro.d))),2),round(max(outlet_states.outlet_hydrograph),4),...
        round(max(max(max(Maps.WQ_States.Pol_Conc_Map))),3),1000*(1/Wshed_Properties.cell_area)*round(max(max(WQ_States.B_t(~isinf(WQ_States.B_t)))),4),round(initial_mass/1000,4),round(sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t))))/1000,4),1-round(sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t))))/1000,4)/round(initial_mass/1000,4),round(WQ_States.EMC_outlet(end,1),2),'VariableNames',...
        {'Drainage area (km2)','Rainfall Vol (mm)','Runoff Coefficient','Impervious Area (km2)','Compactness Coefficient','Circularity index','Equivalent Width (m)','Form Factor','Simulation time (minutes)','Maximum Flood Depth (m)','Maximum Outflow (m^3/s)',...
        'Maximum Concentration (mg/L)','Maximum Stored Mass of Pollutant  (g/m2)','Initial Pollutant Mass  (ton)','Final Pollutant Mass  (ton)','Wash-off Ratio','EMC (mg/L)'});
end
% writetable(Summary_Table)
FileName_String = 'Summary_Table';
FileName = fullfile(folderName,strcat('\',FileName_String,'.csv'));

writetable(Summary_Table,FileName);


% Show Summary Results
fprintf('Drainage area (km2) = %f\n',round(Wshed_Properties.drainage_area/1000/1000,3))
fprintf('Rainfall Vol (mm) = %f\n',rainfall_vol)
fprintf('Runoff Coefficient = %f\n',round(Wshed_Properties.C_r,2))
fprintf('Impervious Area (km2) = %f\n',round(Wshed_Properties.impervious_area/1000/1000,3))
fprintf('Compactness Coefficient = %f\n',round(Wshed_Properties.compactness_coefficient,3))
fprintf('Equivalent Width (m) = %f\n',round(Wshed_Properties.width_SWMM,3));
fprintf('Form Factor = %f\n',round(Wshed_Properties.form_factor,3))
fprintf('Circularity Index = %f\n',round(Wshed_Properties.circularity_index,3));
fprintf('Simulation time (minutes) = %f\n', round(simulation_time/60,3));
% Summary
% fprintf('Maximum Flood Depth (m) = %f\n', round(1/1000*max(max(Max_depth_d))));
fprintf('Maximum Flood Depth (m) = %f\n', round(1/1000*max(max(max(Maps.Hydro.d))),2));
fprintf('Maximum Outflow (m^3/s) = %f\n', round(max(outlet_states.outlet_hydrograph),4));
if flags.flag_waterquality == 1
    % fprintf('Maximum Concentration (mg/L) = %f\n', round(max(max(Max_Pol_Conc_Map))));
    fprintf('Maximum Concentration (mg/L) = %f\n', round(max(max(max(Maps.WQ_States.Pol_Conc_Map))),3));
    fprintf('Maximum Stored Mass of Pollutant  (g/m2) = %f\n', 1000*(1/Wshed_Properties.cell_area)*round(max(max(WQ_States.B_t(~isinf(WQ_States.B_t)))),4));
    fprintf('Initial Pollutant Mass  (ton) = %f\n', round(initial_mass/1000,4));
    fprintf('Final Pollutant Mass  (ton) = %f\n', round(sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t))))/1000,4));
    fprintf('Wash-off Ratio = %f\n', 1-round(sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t))))/1000,4)/round(initial_mass/1000,4))
    fprintf('EMC (mg/L) = %f\n',round(WQ_States.EMC_outlet(end,1),2));
end


%% Maps and Hydrographs
close all
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];
size_font = 10;
if flags.flag_waterquality == 1
    size_font = 12;
    [axis1] = surfplot_maps(DEM_raster,depths.dmax_final/1000,Spectrum,'Easting (m)','Northing (m)','Depth (m)',no_data_value,idx_nan,3,3,1,size_font);
    min_washed = 1e-4;
    % Tot_Washed
    % Min_Washed
    z = WQ_States.Tot_Washed;
    z(z<=min_washed)=nan;
    [axis2] = surfplot_maps(DEM_raster,z,Spectrum,'Easting (m)','Northing (m)','Total Washed Mass (kg)',no_data_value,idx_nan,3,3,2,size_font);
    % Infiltration
    [axis3] = surfplot_maps(DEM_raster,Soil_Properties.I_t,WSE_RAS,'Easting (m)','Northing (m)','Infiltration (mm)',no_data_value,idx_nan,3,3,3,size_font);
    % Velocity
    [axis4] = surfplot_maps(DEM_raster,velocities.vmax_final,Velocity_RAS,'Easting (m)','Northing (m)','Max. Velocity (m/s)',no_data_value,idx_nan,3,3,4,size_font);
    % B_f
    [axis5] = surfplot_maps(DEM_raster,Maps.WQ_States.Pol_mass_map(:,:,end)/Wshed_Properties.cell_area*1000,Velocity_RAS,'Easting (m)','Northing (m)','Pollutant Mass at the end ($\mathrm{g/m^2}$)',no_data_value,idx_nan,3,3,5,size_font);
    % Concentration
    [axis6] = surfplot_maps(DEM_raster,Maps.WQ_States.Pol_mass_map(:,:,end)/Wshed_Properties.cell_area*1000,Spectrum,'Easting (m)','Northing (m)','Pollutant Mass at the end ($\mathrm{mg/L}$)',no_data_value,idx_nan,3,3,6,size_font);

    subplot(3,3,7)
    % Hydrographs
    plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph),'LineWidth',1.5','color','black'); xlabel('Time [min]','interpreter','latex'); ylabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex'); set(gca,'FontSize',size_font);
    hold on
    yyaxis right; set(gca,'ydir','reverse','ycolor','black');
    bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
    ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(gather(Rainfall_Parameters.intensity_rainfall)) + 200]);
    if flags.flag_inflow == 1
        hold on
        yyaxis left; set(gca,'ydir','normal','ycolor','black');
        plot(gather(Inflow_Parameters.time_inflow(1:tfinal_inflow,1)),gather(Inflow_Parameters.inflow_discharge(1:tfinal_inflow,:)),'LineWidth',1.5','color','black','LineStyle','--'); xlabel('Time [min]','interpreter','latex');
        grid on
        legend('Outflow','Inflow','Rainfall','interpreter','latex');
    else
        legend('Outflow','Rainfall','interpreter','latex');
    end
    subplot(3,3,8)
    % Load and Pollutograph
    plot(gather(running_control.time_hydrograph),gather(WQ_States.outet_pollutograph),'LineWidth',1.5,'color',[255,140,0]/255) ;
    xlabel('Time [min]','Interpreter','Latex');
    ylabel('Pol. Concentration (mg/L)','Interpreter','Latex')
    grid on
    hold on
    yyaxis right
    set(gca,'ycolor','black')
    ylabel('Load (kg/sec)','interpreter','latex')
    load_wq = 10^(-3)*gather(WQ_States.outet_pollutograph).*gather(outlet_states.outlet_hydrograph); % kg/sec
    plot(gather(running_control.time_hydrograph),load_wq,'LineWidth',1.5,'color',[34,139,34]/255)
    legend('Concentration','Load','interpreter','latex')

    subplot(3,3,9)
    % First Flush
    m_M = WQ_States.mass_outlet_save./(max(max(WQ_States.mass_outlet)));
    v_V = WQ_States.vol_outlet_save./(max(max(WQ_States.vol_outlet)));
    plot(v_V,m_M,'LineWidth',1.5,'color','black','Marker','*');
    hold on
    plot(0:1,0:1,'LineWidth',1,'Color','black','LineStyle','--')
    ylabel('$m/m_{\mathrm{tot}}$','Interpreter','Latex');
    xlabel('$V/V_{\mathrm{tot}}$','Interpreter','Latex');
    exportgraphics(gcf,fullfile(folderName,'Mapas.png'),'ContentType','image','Colorspace','rgb','Resolution',300)
    % exportgraphics(gcf,fullfile(folderName,'Mapas.pdf'),'ContentType','vector')
    % saveas(gcf,fullfile(folderName,'Mapas.fig'))
else %%%%%%%%%%%%%%%%%%%%
    % No water quality
    size_font = 12;
    no_data_value = nan;
    [axis1] = surfplot_maps(DEM_raster,depths.dmax_final/1000,Spectrum,'Easting (m)','Northing (m)','Maximum Depth (m)',no_data_value,idx_nan,2,2,1,size_font);

    %%% Infiltration

    size_font = 12;
    [axis2] = surfplot_maps(DEM_raster,Soil_Properties.I_t,Spectrum,'Easting (m)','Northing (m)','Infiltration (mm)',no_data_value,idx_nan,2,2,2,size_font);


    %%%% Velocity
    size_font = 12;
    [axis3] = surfplot_maps(DEM_raster,velocities.vmax_final,Velocity_RAS,'Easting (m)','Northing (m)','Max. Velocity (m/s)',no_data_value,idx_nan,2,2,3,size_font);

    ax4 = subplot(2,2,4);
    % Hidrographs
    plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph),'LineWidth',1.5','color','black','marker','*'); xlabel('Time [min]','interpreter','latex'); ylabel('Flow Discharge $(\mathrm{m^3/s}$)','interpreter','latex'); set(gca,'FontSize',size_font);
    hold on
    yyaxis right; set(gca,'ydir','reverse','ycolor','black');
    if flags.flag_spatial_rainfall ~= 1 && flags.flag_rainfall == 1
        bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
        ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(gather(Rainfall_Parameters.intensity_rainfall)) + 200]);
    elseif flags.flag_rainfall == 1
        bar(gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1:(length(BC_States.average_spatial_rainfall)))),gather(BC_States.average_spatial_rainfall),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
        ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(gather(BC_States.average_spatial_rainfall)) + 200]);
    end
    if flags.flag_inflow == 1
        hold on
        yyaxis left; set(gca,'ydir','normal','ycolor','black')
        tfinal_inflow = length(Inflow_Parameters.time_inflow);
        plot(gather(Inflow_Parameters.time_inflow(1:tfinal_inflow,1)),gather(Inflow_Parameters.inflow_discharge(1:tfinal_inflow,:)),'LineWidth',1.5','color','black','LineStyle','--'); xlabel('Time [min]','interpreter','latex');
        ylim([0 max(max(gather(outlet_states.outlet_hydrograph))*1.2,max(max(gather(Inflow_Parameters.inflow_discharge))))]);
        grid on
        legend('Outflow','Inflow','Rainfall');
    else
        legend('Outflow','Rainfall');
    end
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    title('Hydrograph','Interpreter','Latex');
    set(gca,'Tickdir','out')

    % exportgraphics(gcf,fullfile(folderName,'Mapas.png'),'ContentType','image','Colorspace','rgb','Resolution',800)
    % exportgraphics(gcf,fullfile(folderName,'Mapas.pdf'),'ContentType','vector')
    % saveas(gcf,fullfile(folderName,'Mapas.fig'))
end

%% Input Data Maps
% DEM, LULC, SOIL, n, h0, ksat,
close all
% --- DEM --- %
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];

% MAIN INPUT MAPS
size_font = 12;
[axis1] = surfplot_maps(DEM_raster,Elevation_Properties.elevation_cell,Terrain_RAS_ramp,'Easting (m)','Northing (m)','Elevation (m)',no_data_value,idx_nan,1,3,1,size_font);

% --- LULC --- %
size_font = 12;
[axis2] = surfplot_maps(DEM_raster,LULC_Properties.LULC,linspecer(LULC_Properties.n_lulc),'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,1,3,2,size_font);

% --- SOIL --- %
size_font = 12;
[axis3] = surfplot_maps(DEM_raster,Soil_Properties.soil_matrix,linspecer(Soil_Properties.n_soil),'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,1,3,3,size_font);

exportgraphics(gcf,fullfile(strcat(cd,"\",folderName),'Input_Data_Maps.pdf'),'ContentType','vector')
saveas(gcf,fullfile(folderName,'Input_Data_Maps.fig'))
close all

fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];


% --- Manning --- %
if flags.flag_waterquality ~= 1
    [axis1] = surfplot_maps(DEM_raster,LULC_Properties.roughness,Terrain_RAS_ramp,'Easting (m)','Northing (m)','$\mathrm{n~(sm^{-1/3})}    $',no_data_value,idx_nan,4,2,1,size_font);

    % --- h0 --- %
    [axis2] = surfplot_maps(DEM_raster,LULC_Properties.h_0,linspecer,'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,4,2,2,size_font);
else
    [axis1] = surfplot_maps(DEM_raster,LULC_Properties.roughness,Terrain_RAS_ramp,'Easting (m)','Northing (m)','$\mathrm{n~(sm^{-1/3})}$',no_data_value,idx_nan,4,2,3,size_font);

    % --- h0 --- %
    [axis2] = surfplot_maps(DEM_raster,LULC_Properties.h_0,linspecer,'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,4,2,4,size_font);

    % --- C_1 --- %
    [axis3] = surfplot_maps(DEM_raster,LULC_Properties.C_1,linspecer,'Easting (m)','Northing (m)','$C_1~(\mathrm{kg.ha^{-1}}$)',no_data_value,idx_nan,4,2,5,size_font);

    % --- C_2 --- %
    [axis4] = surfplot_maps(DEM_raster,LULC_Properties.C_2,linspecer,'Easting (m)','Northing (m)','$C_2~(\mathrm{day^{-1}}$)',no_data_value,idx_nan,4,2,6,size_font);

    % --- C_3 --- %
    [axis5] = surfplot_maps(DEM_raster,LULC_Properties.C_3,linspecer,'Easting (m)','Northing (m)','$C_3$~(-)',no_data_value,idx_nan,4,2,7,size_font);

    % --- C_4 --- %
    [axis6] = surfplot_maps(DEM_raster,LULC_Properties.C_4,linspecer,'Easting (m)','Northing (m)','$C_4~(-)',no_data_value,idx_nan,4,2,8,size_font);

end
% exportgraphics(gcf,fullfile(folderName,'LULC_Parameters.png'),'ContentType','image','Colorspace','rgb','Resolution',600);
% saveas(gcf,fullfile(folderName,'LULC_Parameters.fig'))

close all
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];
% Ksat
[axis1] = surfplot_maps(DEM_raster,Soil_Properties.ksat,Spectrum,'Easting (m)','Northing (m)','$\mathrm{k_{sat}~(mm.h^{-1})}    $',no_data_value,idx_nan,1,3,1,size_font);

% Dtheta
dheta = Soil_Properties.teta_sat - Soil_Properties.teta_i;
[axis2] = surfplot_maps(DEM_raster,dheta,Spectrum,'Easting (m)','Northing (m)','$\mathrm{\Delta \theta~(-)}    $',no_data_value,idx_nan,1,3,2,size_font);

% Psi
[axis3] = surfplot_maps(DEM_raster,Soil_Properties.psi,Spectrum,'Easting (m)','Northing (m)','$\mathrm{\psi~(mm)}    $',no_data_value,idx_nan,1,3,3,size_font);

exportgraphics(gcf,fullfile(folderName,'SOIL_Parameters.png'),'ContentType','image','Colorspace','rgb','Resolution',150);
% saveas(gcf,fullfile(folderName,'SOIL_Parameters.fig'))
close all

%% DEM with Rivers
close all
tiledlayout(2,1)
set(gcf,'units','inches','position',[2,0,8,6])
ax1=nexttile;

% subplot(2,1,1);
if no_plot==0
    mapshow(A,RA,"AlphaData",0.45);hold on;
    mapshow(S_p,'FaceColor','n'); hold on;
end
imagesc(hillshade(DEM_raster)); hold on;
colormap('gray')
surf(DEM_raster);
colormap(ax1,'landcolor')
% imageschs(DEM_raster,[],'colormap','landcolor','ticklabels','nice'); % DEM
h = colorbar;
h.Label.String = 'Elevation (m)';
h.TickLabelInterpreter = 'latex';
xlabel('Easting (m)','interpreter','latex')
ylabel('Northing (m)','interpreter','latex')
ax = ancestor(ax1, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');

ax2 = nexttile;
% subplot(2,1,2)
FD = FLOWobj(DEM_raster);
As  = flowacc(FD);

if no_plot==0
    mapshow(A,RA,"AlphaData",0.45);hold on;
    mapshow(S_p,'FaceColor','n'); hold on;
end
imagesc(hillshade(DEM_raster)); hold on;
colormap(ax2,'gray')
F = dilate(sqrt(As),ones(5));
F.Z(F.Z==1)=NaN;
surf(F);
colormap(ax2,flowcolor);
h = colorbar;
h.Label.String = '(sqrt(# of pixels)';
h.TickLabelInterpreter = 'latex';
% imageschs(DEM_raster,dilate(sqrt(As),ones(5)),'colormap',flowcolor,...
%     'colorbarylabel','Flow accumulation (sqrt(# of pixels))',...
%     'ticklabel','nice');
xlabel('Easting (m)','interpreter','latex')
ylabel('Northing (m)','interpreter','latex')
ax = ancestor(ax2, 'axes');
ax.XAxis.Exponent = 0;xtickformat('%.0f');
ax.YAxis.Exponent = 0;ytickformat('%.0f');
% exportgraphics(gcf,fullfile(folderName,'Flow_Accumulation.png'),'ContentType','image','Colorspace','rgb','Resolution',600);
close all
%
% dev=Rainfall_Parameters.std_dev;
% rainfall_time=Spatial_Rainfall_Parameters.rainfall_spatial_duration;
% rainfall=BC_States.average_spatial_rainfall;
% streamflow_time=running_control.time_hydrograph;
% streamflow=outlet_states.outlet_hydrograph;
% save('outlet_data.mat',"dev","rainfall_time","rainfall","streamflow_time","streamflow")

clearvars a_grid area_cells area_km2 b_grid baseFileName C cm color_plot color_plots depth_accumulation Depth_RAS elevation f FileName FileName_String filePattern FolderName font_size frame fsize fullFileName h h_max h_min i idx2 idx3 idx_depth idx_i_a idx_wse im imind labels_depth labels_gauges ls max_depth max_h max_inf max_v MS myFolder no_data_value nx_max ny_max Out_Conc points raster_exportion raster_exportion_percentage s size_font t t_max t_previous t_save t_store t_title theFiles topo_path x_grid xbrgin xend xmax y_grid ybegin yend ymax z z1 z2 zmax zmin

disp('Thank you for using HydroPol2D. Results are exported in Modeling Results folder.')

%% Deleting temporary files

% files_to_delete = dir('Temporary_Files');
% for k = 1 : length(files_to_delete)
%     baseFileName = files_to_delete(k).name;
%     fullFileName = fullfile('Temporary_Files/', baseFileName);
%     fprintf(1, 'Now deleting %s\n', fullFileName);
%     delete(fullFileName);
% end




