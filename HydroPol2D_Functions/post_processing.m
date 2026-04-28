%% ========================================================================
% HydroPol2D — post_processing (WHAT THIS SCRIPT *ACTUALLY* EXPORTS)
%
% This post_processing script produces:
%   • Figures (PDF + optional .fig)
%   • Tables (CSV)
%   • Shapefiles (streams + observed gauges)
%   • GeoTIFF rasters (time-varying + static)
%   • (Optionally) GIF/MP4 animations via Inundation_Maps (depends on that function)
%
% ------------------------------------------------------------------------
% 0) EXPORT ROOT (comes from main)
% ------------------------------------------------------------------------
% REQUIRED workspace variable:
%   ExportRootDir   (absolute path chosen by user in the v155 main file)
%
% This script writes everything under:
%   resultsDir = fullfile(ExportRootDir, "Modeling_Results");
%
% It also expects temporary map chunks at:
%   tempDir   = fullfile(ExportRootDir, "Temporary_Files");
%   (save_map_hydro_1.mat, save_map_hydro_2.mat, ...)
%
% ------------------------------------------------------------------------
% 1) FOLDER STRUCTURE CREATED BY THIS SCRIPT
% ------------------------------------------------------------------------
% ExportRootDir
%   ├── Modeling_Results
%   │    ├── Figures_PDF
%   │    ├── Figures_FIG
%   │    ├── Tables_CSV
%   │    ├── Rasters_Water_Depths         % time-varying depth OR WSE (depending on flags.flag_wse)
%   │    ├── Rasters_WSE                  % created, but not used by current code (reserved)
%   │    ├── Rasters_Static               % static rasters (max, totals, DEM, LULC, etc.)
%   │    ├── Rasters_WQ                   % time-varying water quality rasters (if flag_waterquality==1)
%   │    ├── Rasters_Human_Risk           % time-varying human risk rasters (if flag_human_instability>0)
%   │    ├── GIFs_MP4                     % created; used only if Inundation_Maps writes here
%   │    └── Shapefiles
%   └── Temporary_Files                   % NOT created here; must exist from main / model loop
%
% NOTE: This script does NOT (currently) copy General_Data_used.xlsx nor write run_log.txt.
%       If you want those, do it in MAIN (or add it here explicitly).
%
% ------------------------------------------------------------------------
% 2) FIGURES EXPORTED (PDF + FIG)
% ------------------------------------------------------------------------
% Modeling_Results/Figures_PDF
%   • Hydrograph.pdf
%   • Stage_Hydrograph_Outlet.pdf                 (only if flag_obs_gauges==1)
%   • Stage_Hydrograph_Gauges.pdf                 (only if flag_obs_gauges==1)
%   • Specific_Discharge_Gauges.pdf               (only if flag_obs_gauges==1 AND flag_rainfall==1)
%   • Rating_Curve_Outlet.pdf
%   • Rating_Curve_Gauges.pdf                     (only if flag_obs_gauges==1)
%   • Hydrograph_Gauges.pdf                       (only if flag_obs_gauges==1)
%   • DEM_Streams_Observed_Gauges.pdf             (streams always; gauges only if flag_obs_gauges==1)
%
% If Water Quality is enabled (flag_waterquality==1), also:
%   • Pollutograph.pdf
%   • Hysteresis.pdf
%   • M(V)_Curve.pdf
%   • EMC_Curve.pdf
%
% Modeling_Results/Figures_FIG
%   • Same names as above, with .fig extension (where implemented)
%
% ------------------------------------------------------------------------
% 3) TABLES EXPORTED (CSV)
% ------------------------------------------------------------------------
% Modeling_Results/Tables_CSV
%   • Rating_Curve_Data.csv
%   • Outlet_Hydrograph_Data_Outlet.csv
%   • Rating_Curve_Gauges.csv                      (only if flag_obs_gauges==1)
%   • Summary_Table.csv
%
% If Water Quality is enabled (flag_waterquality==1), also:
%   • Outlet_Pollutograph_Data.csv
%   • Outlet_M_V_Curve.csv
%   • Outlet_EMC_Curve.csv
%
% If flag_human_instability == 3:
%   • Risk_summary.csv
%
% ------------------------------------------------------------------------
% 4) SHAPEFILES EXPORTED
% ------------------------------------------------------------------------
% Modeling_Results/Shapefiles
%   • streamnetwork.shp/.shx/.dbf/.prj            (written if flag_obs_gauges==1 in current code path)
%   • observed_gauges.shp/.shx/.dbf/.prj          (only if flag_obs_gauges==1)
%
% NOTE: streamnetwork is currently written inside the flag_obs_gauges==1 block.
%       If you want streams ALWAYS, move shapewrite(MS,...) outside that if.
%
% ------------------------------------------------------------------------
% 5) RASTERS EXPORTED (GeoTIFF)
% ------------------------------------------------------------------------
% A) Time-varying rasters (only if flags.flag_export_maps==1)
%   Modeling_Results/Rasters_Water_Depths
%     • Flood_Depths_YYYY_MM_DD_hh_MM_ss.tif                (if flags.flag_wse==0 and elapsed_time==0)
%     • Flood_Depths_t_<hours>_h.tif                        (if flags.flag_wse==0 and elapsed_time==1)
%     • Water_Surface_Elevation_YYYY_MM_DD_hh_MM_ss.tif     (if flags.flag_wse==1 and elapsed_time==0)
%     • Water_Surface_Elevation_t_<hours>_h.tif             (if flags.flag_wse==1 and elapsed_time==1)
%
%   Modeling_Results/Rasters_WQ   (only if flag_waterquality==1)
%     • Pollutant_Concentration_<timestamp>.tif
%
%   Modeling_Results/Rasters_Human_Risk (only if flag_human_instability==1)
%     • Human_instability_<timestamp>.tif
%
% B) Static rasters (only if flags.flag_export_maps==1)
%   Modeling_Results/Rasters_Static
%     • Accumulation_Areas_Larger_1m.tif
%     • Maximum_Velocity.tif
%     • Infiltrated_Depth.tif
%     • Maximum_Depths.tif                          (if flags.flag_wse==0)
%     • Max_Water_Surface_Elevation.tif              (if flags.flag_wse==1)
%     • DEM_resampled.tif OR DEM_Treated.tif         (depends on flags.flag_resample)
%     • Land_Cover_Data.tif
%
%   If flags.flag_infiltration==1:
%     • Cumulative_Infiltration.tif
%
%   If flags.flag_snow_modeling==1:
%     • Maximum_Snowpack.tif
%
%   If flags.flag_groundwater_modeling==1:
%     • Max_GW_depth.tif
%
%   If flag_waterquality==1:
%     • Initial_Buildup_kg.tif
%     • Total_Washed_Mass_Kg.tif
%     • Accumulation_Areas_10_g_m2.tif
%     • Final_Pollutant_Mass_g_m2.tif
%     • Final_Mass_Of_Pollutant.tif
%     • Maximum_Pol_Conc_min.tif
%
%   If flag_human_instability==1:
%     • Maximum_Instability_Risk.tif
%
%   If flag_human_instability==3:
%     • Maximum_Instability_Risk_<suffix>.tif        (suffix list: _cm,_tm,_am,_om,_cf,_tf,_af,_of)
%
% ------------------------------------------------------------------------
% 6) ANIMATIONS
% ------------------------------------------------------------------------
% This script calls:
%   Inundation_Maps
% ========================================================================

close all
%% Post-Processing - Graphs
simulation_time = toc;

%% Coloramp
[Spectrum,Depth_Purple,Terrain_RAS_ramp,blue_ramp,blues_2,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS,WSE_RAS] = coloramps();

%% Fig Sizes
set(groot,'DefaultFigureWindowStyle','normal');

show_figures = true;   % while tuning
paper_bg = 'white';

figsize.hydrograph            = [6 4];
figsize.stage_outlet          = [6 4];
figsize.stage_gauges          = [9 5];
figsize.specific_discharge    = [9 5];
figsize.rating_outlet         = [5 3];
figsize.rating_gauges         = [9 5];
figsize.hydrograph_gauges     = [9 5];
figsize.pollutograph          = [8 5];
figsize.hysteresis            = [6 4];
figsize.mvcurve               = [6 4];
figsize.emc                   = [6 4];
%% ========= EXPORT ROOT (MUST come from v155 main file) =========
% In v155 you must define and carry this variable into workspace:
%   ExportRootDir  (absolute path chosen by user)
%
% Example in v155:
%   ExportRootDir = fullfile(ProjectDir, "My_Run_2026_03_05");
%
% Here in post-processing we ONLY use it.

ExportRootDir = Paths.Root;

if ~exist('ExportRootDir','var') || isempty(ExportRootDir)
    error(['Post-processing: ExportRootDir is missing. ' ...
        'Define ExportRootDir in the v155 main file and pass it here.']);
end

resultsDir = Paths.Results;

if ~isfolder(resultsDir), mkdir(resultsDir); end

% folderName = base folder where your current code is writing files
folderName = Paths;

%% ========= Subfolders inside resultsDir =========
Dirs.FigPDF   = Paths.FigPDF;
Dirs.FigFIG   = Paths.FigFIG;
Dirs.Tables   = Paths.Tables;

Dirs.RastersWD = Paths.RastersWD;
Dirs.RastersWSE = Paths.RastersWSE;
Dirs.RastersStatic = Paths.RastersStatic;

Dirs.WQMaps   = Paths.WQMaps;
Dirs.HRMaps   = Paths.HRMaps;

Dirs.Anim     = Paths.Anim;

Dirs.Shapes = Paths.Shapes;
if ~isfolder(Dirs.Shapes), mkdir(Dirs.Shapes); end

mapFileOf = @(s) fullfile(tempDir, sprintf('save_map_hydro_%d.mat', s));

% Create all folders
fn = fieldnames(Dirs);
for k = 1:numel(fn)
    if ~isfolder(Dirs.(fn{k})), mkdir(Dirs.(fn{k})); end
end

% Temporary files MUST also be tied to ExportRootDir (NOT pwd)
tempDir = Paths.Temp;
if ~isfolder(tempDir)
    warning("Temporary_Files folder not found at: %s", tempDir);
end
addpath(tempDir);

disp("Post-processing export root:");
disp(folderName);

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
close all force
fig = createStyledFigure(show_figures, paper_bg, figsize.hydrograph);

if flags.flag_elapsed_time == 1
    line_plot(gather(running_control.time_hydrograph),'\mathrm{Elapsed~Time} ','min', ...
        gather(outlet_states.outlet_hydrograph),'\mathrm{Discharge} ','\mathrm{m^3 \cdot s^{-1}}', ...
        [],[],[],[],'Hydrograph',1,1);
else
    line_plot(gather(running_control.time_hydrograph),'\mathrm{Date} ','', ...
        gather(outlet_states.outlet_hydrograph),'\mathrm{Discharge} ','\mathrm{m^3 \cdot s^{-1}}', ...
        [],[],[],[],'Hydrograph',1,1);
end

set(fig,'CurrentAxes',gca);
yyaxis right; set(gca,'ydir','reverse','ycolor','black');
drawnow;
set(fig,'Units','inches','Position',[1 1 figsize.hydrograph]);
drawnow;


if flags.flag_satellite_rainfall == 1
    Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg = Spatial_Rainfall_Parameters.rainfall_spatial_duration;
end

if flags.flag_spatial_rainfall == 1
    if isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1)) && ...
            ~isdatetime(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg)

        % Use the first rainfall timestamp as the reference datetime
        baseDate = Spatial_Rainfall_Parameters.rainfall_spatial_duration(1);

        % Force numeric type to double, then convert minutes -> duration
        dt_minutes = double(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg);

        % Build datetime array correctly
        Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg = baseDate + minutes(dt_minutes);
    end
end

if flags.flag_rainfall == 1
    if flags.flag_spatial_rainfall ~=1
        Rainfall_Parameters.time_rainfall = Rainfall_Parameters.time_rainfall(1:length(Rainfall_Parameters.intensity_rainfall));
        bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',pallete.blue_colors(2,:),'EdgeColor',[0 .5 .5],'LineWidth',1.5)
        ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex'); ylim([0,max(gather(Rainfall_Parameters.intensity_rainfall)) + 200]);
    else
        if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
            dim = length(BC_States.average_spatial_rainfall);
            bar((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(1,1:(dim))))',gather(BC_States.average_spatial_rainfall),'FaceColor',pallete.blue_colors(2,:),'EdgeColor',[0 .5 .5],'LineWidth',1.5);
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
    exportgraphics(gcf, fullfile(Dirs.FigPDF,'Hydrograph.pdf'), 'ContentType','vector');
    saveas(gcf, fullfile(Dirs.FigFIG,'Hydrograph.fig'));
catch
    fprintf('No hydrograph exported, PDF export error')
end
saveas(gcf, fullfile(Dirs.FigFIG,'Hydrograph.fig'));
close all

%% Stage and Hydrograph at the Outlet
if flags.flag_obs_gauges == 1
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.stage_outlet);
    ax = axes('Parent', fig);

    if gauges.num_obs_gauges == 1
        color_plot = linspecer(2);
    else
        color_plot = linspecer(gauges.num_obs_gauges);
    end

    font_size = 12;

    % Left axis: Discharge
    yyaxis(ax, 'left');
    set(ax, 'YColor', color_plot(1,:));
    plot(ax, gather(running_control.time_hydrograph), gather(outlet_states.outlet_hydrograph), ...
        'LineWidth', 1.5, 'Color', color_plot(1,:));
    xlabel(ax, 'Elapsed Time', 'Interpreter', 'latex');
    ylabel(ax, '$Q~(\mathrm{m^3/s})$', 'Interpreter', 'latex');

    % Right axis: Stage
    yyaxis(ax, 'right');
    set(ax, 'YColor', color_plot(2,:));
    plot(ax, gather(running_control.time_hydrograph), gather(outlet_states.depth_outlet / 1000), ...
        'LineWidth', 1.5, 'LineStyle', '--', 'Color', color_plot(2,:));
    ylabel(ax, '$h~(\mathrm{m})$', 'Interpreter', 'latex');

    % Styling
    set(ax, 'FontSize', font_size, ...
        'FontName', 'Helvetica', ...
        'TickDir', 'out', ...
        'LineWidth', 2, ...
        'TickLength', [0.02 0.01]);
    box(ax, 'on');
    grid(ax, 'on');
    title(ax, 'Outlet', 'Interpreter', 'latex', 'FontSize', font_size);

    drawnow;

    % Export to vector PDF
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Stage_Hydrograph_Outlet.pdf'), 'ContentType', 'vector');
    catch
        fprintf('No Stage hydrographs exported, PDF export error\n');
    end

    % Save .fig file
    saveas(fig, fullfile(Dirs.FigFIG, 'Stage_Hydrograph_Outlet.fig'));

    % close(fig);   % leave commented while tuning figure sizes
end

%% Stage Hydrograph %%
if flags.flag_obs_gauges == 1
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.stage_gauges);
    ax = axes('Parent', fig);

    color_plots = linspecer(gauges.num_obs_gauges);
    font_size = 14;
    hold(ax, 'on');

    % Plot gauge depths
    for i = 1:gauges.num_obs_gauges
        if mod(i,3) == 0
            ls = '--';
        elseif mod(i,3) == 1
            ls = '-';
        else
            ls = ':';
        end

        plot(ax, gather(running_control.time_hydrograph), gather(gauges.depth_cell(:,i)), ...
            'LineWidth', 1.5, 'LineStyle', ls, 'Color', color_plots(i,:));
    end

    % Plot outlet depth
    plot(ax, gather(running_control.time_hydrograph), gather(outlet_states.depth_outlet / 1000), ...
        'LineWidth', 1.5, 'LineStyle', '-', 'Marker', '.', 'Color', 'red');

    xlabel(ax, 'Time [min]', 'Interpreter', 'latex');
    ylabel(ax, 'Depth $[\mathrm{m}]$', 'Interpreter', 'latex');

    % Formatting: left y-axis
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', font_size, ...
        'TickDir', 'out', ...
        'LineWidth', 2, ...
        'TickLength', [0.02 0.01]);
    box(ax, 'on');
    grid(ax, 'on');

    % Right y-axis: rainfall
    yyaxis(ax, 'right');
    set(ax, 'YDir', 'reverse', 'YColor', 'black');

    if flags.flag_rainfall == 1
        if flags.flag_spatial_rainfall ~= 1
            bar(ax, gather(Rainfall_Parameters.time_rainfall), gather(Rainfall_Parameters.intensity_rainfall), ...
                'FaceColor', [0 0.55 0.55], 'EdgeColor', [0 0.5 0.5], 'LineWidth', 1.5);
            ylabel(ax, 'Rainfall Intensity $[\mathrm{mm \cdot h^{-1}}]$', 'Interpreter', 'latex');
        else
            dim = length(BC_States.average_spatial_rainfall);
            plot(ax, gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1, 1:dim)), ...
                gather(BC_States.average_spatial_rainfall), ...
                'LineWidth', 1.5, 'Color', 'blue');
            ylabel(ax, 'Mean Rainfall Intensity $[\mathrm{mm \cdot h^{-1}}]$', 'Interpreter', 'latex');
        end
    end

    % Legend
    if isempty(gauges.labels_observed_string)
        gauges.labels_observed_string = string(1:1:length(gauges.x_coord_gauges));
    end
    labels_depth = gauges.labels_observed_string;
    labels_depth{gauges.num_obs_gauges + 1} = 'Outlet';
    labels_depth{gauges.num_obs_gauges + 2} = 'Rainfall Intensity';
    legend(ax, labels_depth, 'FontName', 'Helvetica', 'FontSize', 10, 'Location', 'bestoutside');

    % Axis limit for rainfall
    try
        ylim(ax, [0, max(max(gather(BC_States.average_spatial_rainfall))) * 6]);
    catch
        ylim(ax, [0, 10]);
    end

    % Title
    title(ax, 'Surface Runoff Depth', 'Interpreter', 'latex', 'FontSize', font_size);

    % Export
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Stage_Hydrograph_Gauges.pdf'), 'ContentType', 'vector');
    catch
        fprintf('No stage hydrograph gauges saved, PDF export error\n');
    end

    saveas(fig, fullfile(Dirs.FigFIG, 'Stage_Hydrograph_Gauges.fig'));

    % close(fig);   % leave commented while tuning figure sizes
end

%% Normalized Discharge %%
% Rainfall Std Deviation
zero_matrix = zeros(size(Elevation_Properties.elevation_cell,1), size(Elevation_Properties.elevation_cell,2));

%-----------------------------------------
% Load maps in chunks (simple + robust)
%-----------------------------------------
store = 1;
flag_loader = 1;
i = 1;   % initialize first chunk load

if i > saver_memory_maps * store
    store = store + 1;
    flag_loader = 1; % force load for the new store
end

if flag_loader == 1
    mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d.mat', store));
    if ~isfile(mapFile)
        error("Missing map file: %s", mapFile);
    end
    load(mapFile, 'Maps');
    flag_loader = 0;

    % Clean once per loaded chunk
    no_data_value = nan;
    Maps.Hydro.d(isnan(Maps.Hydro.d)) = no_data_value;
    Maps.Hydro.d(isinf(Maps.Hydro.d)) = no_data_value;
    Maps.Hydro.d(idx_nan)             = no_data_value;

    % Update max depth from this chunk
    if exist('Max_depth_d','var') && ~isempty(Max_depth_d)
        Max_depth_d = max(Max_depth_d, max(Maps.Hydro.d, [], 3));
    else
        Max_depth_d = max(Maps.Hydro.d, [], 3);
    end
end

if flags.flag_spatial_rainfall == 1 && running_control.record_time_spatial_rainfall
    store = 1;
    flag_loader = 1;
    rainfall_sum = zeros(size(zero_matrix));
    for i = 1:length(running_control.time_records)
        if i > saver_memory_maps * store
            store = store + 1;
            load(fullfile(tempDir, ['save_map_hydro_' num2str(store)]), 'Maps');
        else
            if flag_loader == 1
                load(fullfile(tempDir, ['save_map_hydro_' num2str(store)]), 'Maps');
                flag_loader = 0;
            end
        end
        z = Maps.Hydro.spatial_rainfall_maps(:,:,i - ((store-1)*saver_memory_maps));
        rainfall_sum = rainfall_sum + z;
        Rainfall_Parameters.std_dev(i,1) = nanstd(z(:));
    end
end

if flags.flag_obs_gauges == 1 && flags.flag_rainfall == 1
    % Compute catchment area per gauge
    for i = 1:length(gauges.easting_obs_gauges)
        gauges.catchment_area(i,1) = Wshed_Properties.fac_area( ...
            gauges.northing_obs_gauges(i,1), gauges.easting_obs_gauges(i,1)); % km²
    end

    % Set up figure
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.specific_discharge);
    ax = axes('Parent', fig);
    color_plots = linspecer(gauges.num_obs_gauges);
    hold(ax, 'on');

    font_size = 12;

    % Plot specific discharge for each gauge
    for i = 1:gauges.num_obs_gauges
        if mod(i,3) == 0
            ls = '--';
        elseif mod(i,3) == 1
            ls = '-';
        else
            ls = ':';
        end

        specific_discharge = gauges.hydrograph_cell(:,i) ./ gauges.catchment_area(i,1); % m³/s/km²
        plot(ax, gather(running_control.time_hydrograph), specific_discharge, ...
            'LineWidth', 1.5, 'LineStyle', ls, 'Color', color_plots(i,:));
    end

    % Plot outlet
    outlet_sd = gather(outlet_states.outlet_hydrograph) / (Wshed_Properties.drainage_area / 1e6); % km²
    plot(ax, gather(running_control.time_hydrograph), outlet_sd, ...
        'LineWidth', 1.5, 'LineStyle', '-', 'Marker', '.', 'Color', 'red');

    xlabel(ax, 'Elapsed Time [min]', 'Interpreter', 'latex');
    ylabel(ax, 'Specific Discharge $[\mathrm{m^3/s/km^2}]$', 'Interpreter', 'latex');

    % Style left axis
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', font_size, ...
        'TickDir', 'out', ...
        'LineWidth', 2, ...
        'TickLength', [0.02 0.01]);
    box(ax, 'on');
    grid(ax, 'on');

    % Right axis: Rainfall
    yyaxis(ax, 'right');
    set(ax, 'YDir', 'reverse', 'YColor', 'black');

    if flags.flag_rainfall == 1
        if flags.flag_spatial_rainfall ~= 1
            bar(ax, gather(Rainfall_Parameters.time_rainfall), gather(Rainfall_Parameters.intensity_rainfall), ...
                'FaceColor', [0 .55 .55], 'EdgeColor', [0 .5 .5], 'LineWidth', 1.5);
            ylabel(ax, 'Rainfall Intensity $[\mathrm{mm \cdot h^{-1}}]$', 'Interpreter', 'latex');
            ylim(ax, [0, max(gather(Rainfall_Parameters.intensity_rainfall)) * 6]);
        else
            dim = length(BC_States.average_spatial_rainfall);
            bar(ax, Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:dim), ...
                gather(BC_States.average_spatial_rainfall), ...
                'FaceColor', [0 .5 .5], 'EdgeColor', [0 .55 .55], 'LineWidth', 1.5);
            ylabel(ax, 'Areal Mean Rainfall Intensity $[\mathrm{mm \cdot h^{-1}}]$', 'Interpreter', 'latex');
            hold(ax, 'on');
            try
                er = errorbar(ax, Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:dim), ...
                    BC_States.average_spatial_rainfall, ...
                    Rainfall_Parameters.std_dev(1:dim,1), ...
                    Rainfall_Parameters.std_dev(1:dim,1));
                er.Color = [0 0 0];
                er.LineStyle = 'none';
            catch
            end
            plot(ax, gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:dim)), ...
                gather(BC_States.average_spatial_rainfall), ...
                'LineWidth', 1.5, 'Color', 'blue');
            ylim(ax, [0, max(max(gather(BC_States.average_spatial_rainfall))) * 6]);
        end
    end

    % Legend
    labels_depth = gauges.labels_observed_string;
    labels_depth{gauges.num_obs_gauges + 1} = 'Outlet';
    labels_depth{gauges.num_obs_gauges + 2} = 'Rainfall Intensity';
    legend(ax, labels_depth, 'FontName', 'Helvetica', 'FontSize', 10, 'Location', 'bestoutside');

    % Title
    title(ax, 'Specific Discharge', 'Interpreter', 'latex', 'FontSize', font_size);

    % Export
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Specific_Discharge_Gauges.pdf'), 'ContentType', 'vector');
    catch
        fprintf('Specific discharge gauges not exported, PDF export error\n');
    end
    saveas(fig, fullfile(Dirs.FigFIG, 'Specific_Discharge_Gauges.fig'));

    % close(fig);   % leave commented while tuning figure sizes
end

%% Total ETR and ETP
% This assumes maps are extracted at consistent spatial resolution.
% ETP is typically daily; use rainfall time step to integrate over time.

zero_matrix = zeros(size(Elevation_Properties.elevation_cell, 1), ...
    size(Elevation_Properties.elevation_cell, 2));

if flags.flag_ETP == 1
    store = 1;
    flag_loader = 1;
    ETR_sum = zeros(size(zero_matrix));
    ETP_sum = zeros(size(zero_matrix));

    for i = 1:length(running_control.time_records)
        try
            % Load new batch of saved maps when needed
            if i > saver_memory_maps * store
                store = store + 1;
                load(fullfile(tempDir, ['save_map_hydro_' num2str(store)]), 'Maps');
            elseif flag_loader == 1
                load(fullfile(tempDir, ['save_map_hydro_' num2str(store)]), 'Maps');
                flag_loader = 0;
            end

            % Get index for this frame
            idx = i - ((store - 1) * saver_memory_maps);

            % ---- ETR ----
            z_etr = Maps.Hydro.ETR_save(:, :, idx);
            z_etr(isnan(z_etr)) = 0;
            ETR_sum = ETR_sum + z_etr;
            ETP_Parameters.std_dev_ETR(i, 1) = nanstd(z_etr(:));

            % ---- ETP ----
            z_etp = Maps.Hydro.ETP_save(:, :, idx);
            z_etp(isnan(z_etp)) = 0;
            ETP_sum = ETP_sum + z_etp;
            ETP_Parameters.std_dev_ETP(i, 1) = nanstd(z_etp(:));

        catch err
            warning('Error processing ETP/ETR map at step %d: %s', i, err.message);
        end
    end
end

%% Rating Curve – Outlet
close all force
fig = createStyledFigure(show_figures, paper_bg, figsize.rating_outlet);
ax = axes('Parent', fig);

% Determine WSE at outlet
if size(Wshed_Properties.row_min, 1) > 0
    % Multiple outlets
    wse_outlet = Wshed_Properties.el_outlet(Wshed_Properties.row_min(1), Wshed_Properties.col_min(1)) + outlet_states.depth_outlet / 1000;
    s = scatter(ax, outlet_states.outlet_hydrograph, wse_outlet, 'o');
else
    % Single outlet
    wse_outlet = Wshed_Properties.el_outlet(Wshed_Properties.row_min, Wshed_Properties.col_min) + outlet_states.depth_outlet / 1000;
    s = scatter(ax, gather(outlet_states.outlet_hydrograph), gather(wse_outlet), 'o');
end

% Axis labels
xlabel(ax, 'Flow Discharge $[\mathrm{m^3/s}]$', 'Interpreter', 'latex');
ylabel(ax, 'Water Surface Elevation $[\mathrm{m}]$', 'Interpreter', 'latex');

% Formatting
set(ax, 'FontSize', 14, ...
    'FontName', 'Helvetica', ...
    'TickDir', 'out', ...
    'LineWidth', 2);
box(ax, 'on');
grid(ax, 'on');

% Adjust Y-limits
if max(wse_outlet) == Wshed_Properties.stage_min
    ylim(ax, [Wshed_Properties.stage_min, Wshed_Properties.stage_min + 1]);
else
    ylim(ax, [Wshed_Properties.stage_min, max(wse_outlet)]);
end

% Style scatter points
s.LineWidth = 0.75;
s.MarkerEdgeColor = [0 0.4 0.7];
s.MarkerFaceColor = [0 0.6 0.6];
s.SizeData = 20;

% Export figure
try
    exportgraphics(fig, fullfile(Dirs.FigPDF, 'Rating_Curve_Outlet.pdf'), 'ContentType', 'vector');
catch
    fprintf('No rating curve outlet exported – PDF export error\n');
end
saveas(fig, fullfile(Dirs.FigFIG, 'Rating_Curve_Outlet.fig'));

% close(fig);   % leave commented while tuning figure sizes

% Export rating curve data
Rating_Curve_Data = table( ...
    gather(outlet_states.outlet_hydrograph), ...
    gather(wse_outlet), ...
    gather(outlet_states.depth_outlet) / 1000, ...
    'VariableNames', {'Flow Discharge (m3/s)', 'WSE (m)', 'Depth (m)'});
writetable(Rating_Curve_Data, fullfile(Dirs.Tables, 'Rating_Curve_Data.csv'));

% Export hydrograph table
Outlet_Hydrograph_Data = table( ...
    gather(running_control.time_hydrograph), ...
    gather(outlet_states.outlet_hydrograph), ...
    'VariableNames', {'Time (min)', 'Flow Discharge (m3/s)'});
writetable(Outlet_Hydrograph_Data, fullfile(Dirs.Tables, 'Outlet_Hydrograph_Data_Outlet.csv'));

%% Rating Curve – Specific Cell
if flags.flag_obs_gauges == 1
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.rating_gauges);
    color_plot = linspecer(gauges.num_obs_gauges);

    for i = 1:gauges.num_obs_gauges
        % Font size based on number of gauges
        if gauges.num_obs_gauges > 5
            fsize = 8;
        else
            fsize = 14;
        end

        % Subplot layout
        if gauges.num_obs_gauges > 3
            temp = double(ceil(gauges.num_obs_gauges / 3));
            subplot(temp, 3, double(i), 'Parent', fig);
        elseif gauges.num_obs_gauges == 2
            subplot(1, 2, double(i), 'Parent', fig);
        else
            subplot(1, 1, 1, 'Parent', fig);
        end

        % Scatter plot of Q vs WSE
        s = scatter(gather(gauges.hydrograph_cell(:, i)), gather(gauges.wse_cell(:, i)), 'o');
        s.MarkerFaceColor = color_plot(i,:);
        s.MarkerEdgeColor = [0.3 0.3 0.3];
        s.SizeData = 20;
        s.LineWidth = 0.75;

        xlabel('$Q~(\mathrm{m^3/s})$', 'Interpreter', 'latex');
        ylabel('$\mathrm{WSE~(m)}$', 'Interpreter', 'latex');

        % Styling
        title(gauges.labels_observed_string{i}, 'Interpreter', 'latex', 'FontSize', fsize);
        set(gca, 'FontSize', fsize, ...
            'FontName', 'Helvetica', ...
            'TickDir', 'out', ...
            'TickLength', [0.02 0.01], ...
            'LineWidth', 2);
        box on;
        grid on;
    end

    % Export figure
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Rating_Curve_Gauges.pdf'), 'ContentType', 'vector');
    catch
        fprintf('No rating curve gauges exported – PDF export error\n');
    end

    saveas(fig, fullfile(Dirs.FigFIG, 'Rating_Curve_Gauges.fig'));

    % close(fig);   % leave commented while tuning figure sizes

    % Export data
    Rating_Curve_Specific_Cell = table( ...
        gather(running_control.time_hydrograph), ...
        gather(gauges.hydrograph_cell), ...
        gather(gauges.wse_cell), ...
        gather(gauges.depth_cell), ...
        'VariableNames', {'Time [min] or Date', 'Flow Discharge (m3/s)', 'WSE (m)', 'Water Depth (m)'});

    writetable(Rating_Curve_Specific_Cell, fullfile(Dirs.Tables, 'Rating_Curve_Gauges.csv'));
end

%% Hydrographs – Specific Cell
if flags.flag_obs_gauges == 1
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.hydrograph_gauges);

    % Color handling
    if gauges.num_obs_gauges == 1
        color_plot = linspecer(10);
    else
        color_plot = linspecer(gauges.num_obs_gauges);
    end

    for i = 1:gauges.num_obs_gauges
        % Font size by layout density
        if gauges.num_obs_gauges > 5
            fsize = 8;
        else
            fsize = 14;
        end

        % Subplot layout
        if gauges.num_obs_gauges > 3
            temp = double(ceil(gauges.num_obs_gauges / 3));
            subplot(temp, 3, double(i), 'Parent', fig);
        elseif gauges.num_obs_gauges == 2
            subplot(1, 2, double(i), 'Parent', fig);
        else
            subplot(1, 1, 1, 'Parent', fig);
        end

        % Left axis: Discharge
        yyaxis left;
        set(gca, 'YColor', pallete.blue_colors(1,:));
        plot(gather(running_control.time_hydrograph), gather(gauges.hydrograph_cell(:,i)), ...
            'LineWidth', 1.5, 'Color', pallete.blue_colors(2,:));
        ylabel('$Q~(\mathrm{m^3/s})$', 'Interpreter', 'latex');

        % Right axis: Depth
        yyaxis right;
        set(gca, 'YColor', pallete.red_colors(3,:));
        plot(gather(running_control.time_hydrograph), gather(gauges.depth_cell(:,i)), ...
            'LineWidth', 1.5, 'LineStyle', '--', 'Color', pallete.red_colors(3,:));
        ylabel('$h~(\mathrm{m})$', 'Interpreter', 'latex');

        % Common formatting
        xlabel('Elapsed Time [min]', 'Interpreter', 'latex');
        title(gauges.labels_observed_string{i}, 'Interpreter', 'latex', 'FontSize', fsize);
        set(gca, 'FontSize', fsize, ...
            'FontName', 'Helvetica', ...
            'TickDir', 'out', ...
            'LineWidth', 2, ...
            'TickLength', [0.02 0.01]);
        box on;
        grid on;
    end

    % Export figure
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Hydrograph_Gauges.pdf'), 'ContentType', 'vector');
    catch
        fprintf('No Hydrograph gauges exported – PDF export error\n');
    end
    saveas(fig, fullfile(Dirs.FigFIG, 'Hydrograph_Gauges.fig'));

    % close(fig);   % leave commented while tuning figure sizes
end

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

%% DEM with Streams and Observed Points

% Flow direction
FD = FLOWobj(DEM_raster);
area_km2 = GIS_data.min_area;
area_cells = area_km2 / ((DEM_raster.cellsize / 1000)^2); % pixels

% Optional map display
if no_plot == 0
    try
        mapshow(A, RA, 'AlphaData', 0.45); hold on;
        mapshow(S_p, 'FaceColor', 'none'); hold on;
    catch
        warning('You need MATLAB 2022 or higher to run basemaps');
    end
end

% Create stream network
S = STREAMobj(FD, 'minarea', area_cells);
ax1 = plot(S); hold on;

% Labels and formatting
title('Streams and Observed Points', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x~(\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y~(\mathrm{m})$', 'Interpreter', 'latex', 'FontSize', 14);

% Format axis if not a map axis
ax = ancestor(ax1, 'axes');
if isprop(ax, 'XAxis') % avoid error in map axes
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', 14, ...
        'TickDir', 'out', ...
        'LineWidth', 2);
    box on; grid on;
end

% Gauges
if flags.flag_obs_gauges == 1
    scatter(gauges.easting_obs_gauges_absolute, gauges.northing_obs_gauges_absolute, ...
        40, 'filled', 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerEdgeColor', 'k');

    % Export streams
    MS = STREAMobj2mapstruct(S);
    shapewrite(MS, fullfile(Dirs.Shapes, 'streamnetwork.shp')); % ok to keep at root
    % OR keep them in a subfolder:
    % Dirs.Shapes = fullfile(folderName,"Shapefiles"); if ~isfolder(Dirs.Shapes), mkdir(Dirs.Shapes); end
    % shapewrite(MS, fullfile(Dirs.Shapes,'streamnetwork.shp'));

    % Gauges to shapefile
    gauges.x_coord_gauges = gauges.easting_obs_gauges_absolute;
    gauges.y_coord_gauges = gauges.northing_obs_gauges_absolute;
    labels_gauges = gauges.labels_observed_string;

    points = mappoint(gauges.x_coord_gauges', gauges.y_coord_gauges', 'Label', labels_gauges');
    shapewrite(points, fullfile(Dirs.Shapes, 'observed_gauges.shp'));
end

% ==== Save Figure ====
try
    % Save as FIG (MATLAB format)
    saveas(gcf, fullfile(Dirs.FigFIG, 'DEM_Streams_Observed_Gauges.fig'));

    % Export as high-quality vector PDF
    exportgraphics(gcf, fullfile(Dirs.FigPDF, 'DEM_Streams_Observed_Gauges.pdf'), 'ContentType', 'vector');

catch err
    warning('Could not export DEM with streams and gauges')
end


close all;

%% Water Quality Analysis
if flags.flag_waterquality == 1
    %% Pollutograph – Concentration and Load at Outlet
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.pollutograph);
    ax = axes('Parent', fig);

    font_size = 14;

    % Left Y-axis: Pollutant concentration
    yyaxis(ax, 'left');
    plot(ax, gather(running_control.time_hydrograph), gather(WQ_States.outet_pollutograph), ...
        'LineWidth', 1.5, 'Color', 'r', 'Marker', 'o');
    xlabel(ax, 'Time [min]', 'Interpreter', 'latex');
    ylabel(ax, 'Concentration $[\mathrm{mg/L}]$', 'Interpreter', 'latex');
    set(ax, 'YColor', 'r');

    % Right Y-axis: Pollutant load
    yyaxis(ax, 'right');
    load_wq = 1e-3 * gather(WQ_States.outet_pollutograph) .* gather(outlet_states.outlet_hydrograph); % kg/s
    plot(ax, gather(running_control.time_hydrograph), load_wq, ...
        'LineWidth', 1.5, 'Color', 'b', 'Marker', '*');
    ylabel(ax, 'Load $[\mathrm{kg/s}]$', 'Interpreter', 'latex');
    set(ax, 'YColor', 'b');

    % Styling
    legend(ax, {'Concentration', 'Load'}, 'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'best');
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', font_size, ...
        'TickDir', 'out', ...
        'LineWidth', 2, ...
        'TickLength', [0.02 0.01]);
    grid(ax, 'on');
    box(ax, 'on');

    % Export figure
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Pollutograph.pdf'), 'ContentType', 'vector');
    catch
        warning('Pollutograph PDF export failed.');
    end
    saveas(fig, fullfile(Dirs.FigFIG, 'Pollutograph.fig'));

    % Export data
    Outlet_Pollutograph_Data = table( ...
        gather(running_control.time_hydrograph), ...
        gather(WQ_States.outet_pollutograph), ...
        gather(load_wq), ...
        'VariableNames', {'Time [min]', 'Concentration (mg/L)', 'Load (kg/s)'});

    writetable(Outlet_Pollutograph_Data, fullfile(Dirs.Tables, 'Outlet_Pollutograph_Data.csv'));

    %% Hysteresis Effect
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.hysteresis);
    ax = axes('Parent', fig);
    font_size = 14;

    % Left Y-axis: Flow discharge
    yyaxis(ax, 'left');
    plot(ax, gather(running_control.time_hydrograph), gather(outlet_states.outlet_hydrograph), ...
        'Color', 'black', 'LineWidth', 1.5, 'Marker', '*');
    xlabel(ax, 'Time [min]', 'Interpreter', 'latex');
    ylabel(ax, 'Flow Discharge $[\mathrm{m^3/s}]$', 'Interpreter', 'latex');

    % Right Y-axis: Pollutant concentration
    yyaxis(ax, 'right');
    set(ax, 'YColor', 'r');
    plot(ax, gather(running_control.time_hydrograph), gather(WQ_States.outet_pollutograph), ...
        'LineWidth', 1.5, 'Color', 'red', 'Marker', 'o');
    ylabel(ax, 'Concentration $[\mathrm{mg/L}]$', 'Interpreter', 'latex');

    % Formatting
    legend(ax, {'Flow Discharge', 'Pollutant Concentration'}, 'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'best');
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', font_size, ...
        'TickDir', 'out', ...
        'LineWidth', 2, ...
        'TickLength', [0.02 0.01]);
    box(ax, 'on');
    grid(ax, 'on');

    % Export hysteresis plot
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'Hysteresis.pdf'), 'ContentType', 'vector');
    catch
        warning('Failed to export Hysteresis PDF');
    end
    saveas(fig, fullfile(Dirs.FigFIG, 'Hysteresis.fig'));

    %% M(V) Curve
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.mvcurve);
    ax = axes('Parent', fig);

    % Normalized data
    m_M = WQ_States.mass_outlet_save ./ max(WQ_States.mass_outlet(:)); % m/m_tot
    v_V = WQ_States.vol_outlet_save ./ max(WQ_States.vol_outlet(:));   % V/V_tot

    % Plot M(V) curve
    plot(ax, gather(v_V), gather(m_M), ...
        'LineWidth', 1.5, 'Color', 'black', 'Marker', '*');
    hold(ax, 'on');

    % 1:1 reference line
    plot(ax, [0 1], [0 1], 'LineWidth', 1, 'Color', 'black', 'LineStyle', '--');

    xlabel(ax, '$V/V_{tot}$', 'Interpreter', 'latex');
    ylabel(ax, '$m/m_{tot}$', 'Interpreter', 'latex');

    % Formatting
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', font_size, ...
        'TickDir', 'out', ...
        'LineWidth', 2, ...
        'TickLength', [0.02 0.01]);
    box(ax, 'on');
    grid(ax, 'on');

    % Export M(V) curve
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'M(V)_Curve.pdf'), 'ContentType', 'vector');
    catch
        warning('Failed to export M(V) curve PDF');
    end
    saveas(fig, fullfile(Dirs.FigFIG, 'M(V)_Curve.fig'));

    % Export CSV data
    Outlet_M_V_Curve = table(round(gather(v_V), 3), round(gather(m_M), 3), ...
        'VariableNames', {'Normalized Volume', 'Normalized Pollutant Mass'});

    writetable(Outlet_M_V_Curve, fullfile(Dirs.Tables, 'Outlet_M_V_Curve.csv'));

    %% EMC (Event Mean Concentration) Curve
    close all force
    fig = createStyledFigure(show_figures, paper_bg, figsize.emc);
    ax = axes('Parent', fig);

    font_size = 14;

    % Plot EMC curve
    plot(ax, gather(running_control.time_hydrograph), gather(WQ_States.EMC_outlet), ...
        'LineWidth', 1.5, 'Color', 'black');

    xlabel(ax, 'Elapsed Time [min]', 'Interpreter', 'latex');
    ylabel(ax, 'EMC $[\mathrm{mg/L}]$', 'Interpreter', 'latex');

    % Styling
    set(ax, 'FontName', 'Helvetica', ...
        'FontSize', font_size, ...
        'TickDir', 'out', ...
        'TickLength', [0.02 0.01], ...
        'LineWidth', 2);
    box(ax, 'on');
    grid(ax, 'on');

    % Export figure
    try
        exportgraphics(fig, fullfile(Dirs.FigPDF, 'EMC_Curve.pdf'), 'ContentType', 'vector');
    catch
        warning('Failed to export EMC curve PDF');
    end
    saveas(fig, fullfile(Dirs.FigFIG, 'EMC_Curve.fig'));

    % Export EMC data table
    Outlet_EMC_Curve = table( ...
        gather(running_control.time_hydrograph), ...
        round(gather(WQ_States.EMC_outlet), 2), ...
        'VariableNames', {'Elapsed Time [min]', 'EMC (mg/L)'});

    writetable(Outlet_EMC_Curve, fullfile(Dirs.Tables, 'Outlet_EMC_Curve.csv'));
end

%% Exporting Rasters (cross-platform) + ONE NetCDF PER TIME-SERIES VARIABLE
if flags.flag_export_maps == 1

    no_data_value = nan;

    % =========================================================
    % NetCDF export controls
    % =========================================================
    % NetCDF will be exported ONLY for time-varying rasters,
    % and each variable will be stored in ONE single 3D file:
    %   (x, y, time)
    export_netcdf_timeseries = false;

    % Compression level for NetCDF4 files (0 = none, 9 = max)
    nc_deflate_level = 4;

    %-----------------------------
    % Helper: ensure folder exists
    %-----------------------------
    ensureFolder = @(p) ( ~isfolder(p) && mkdir(p) );

    %=========================================
    % 0) Make sure main export folders exist
    %=========================================
    ensureFolder(Dirs.RastersStatic);
    ensureFolder(Dirs.Tables);

    %========================
    % 1) Create/Clean WQ/HR folders
    %========================
    myFolder_wq = "";
    myFolder_hr = "";

    if flags.flag_waterquality == 1
        myFolder_wq = Dirs.WQMaps;
        myFolder_hr = Dirs.HRMaps;

        ensureFolder(myFolder_wq);
        ensureFolder(myFolder_hr);

        % Remove old TIFFs
        deleteMatchingFiles(myFolder_wq, {'*.tif', '*.tif.aux'});
        deleteMatchingFiles(myFolder_hr, {'*.tif', '*.tif.aux'});

        % Remove old aggregated NetCDF files
        deleteMatchingFiles(myFolder_wq, {'*_timeseries.nc'});
        deleteMatchingFiles(myFolder_hr, {'*_timeseries.nc'});

    elseif flags.flag_human_instability > 0
        myFolder_hr = Dirs.HRMaps;
        ensureFolder(myFolder_hr);

        deleteMatchingFiles(myFolder_hr, {'*.tif', '*.tif.aux'});
        deleteMatchingFiles(myFolder_hr, {'*_timeseries.nc'});
    end

    %========================
    % 2) Create/Clean Water Depths folder
    %========================
    myFolder_wd = Dirs.RastersWD;
    ensureFolder(myFolder_wd);

    if ~isfolder(myFolder_wd)
        errorMessage = sprintf("Error: The following folder does not exist:\n%s", myFolder_wd);
        uiwait(warndlg(errorMessage));
        return;
    end

    deleteMatchingFiles(myFolder_wd, {'*.tif', '*.tif.aux'});
    deleteMatchingFiles(myFolder_wd, {'*_timeseries.nc'});

    %=========================================================
    % 3) PREPARE TIME VECTOR FOR AGGREGATED NETCDF FILES
    %=========================================================
    nTimes = length(running_control.time_records);

    tvec = running_control.time_records(:);

    if flags.flag_elapsed_time == 1
        % --------------------------------------------------
        % ELAPSED TIME (assumed numeric or duration)
        % --------------------------------------------------
        if isduration(tvec)
            time_values_nc = hours(tvec);   % ✅ FIX
        else
            time_values_nc = double(tvec) / 60; % if stored as minutes
        end

        time_units_nc  = 'hours since simulation start';
        time_calendar_nc = 'none';

    else
        % --------------------------------------------------
        % ABSOLUTE TIME (datetime case)
        % --------------------------------------------------
        if isdatetime(tvec)
            t0_nc = tvec(1);

            % Convert to days since t0
            time_values_nc = days(tvec - t0_nc);   % ✅ FIX

            time_units_nc  = ['days since ' datestr(t0_nc, 'yyyy-mm-dd HH:MM:ss')];
            time_calendar_nc = 'proleptic_gregorian';

        elseif isduration(tvec)
            % Rare case: duration but not elapsed flag
            time_values_nc = days(tvec);  % fallback
            time_units_nc  = 'days since simulation start';
            time_calendar_nc = 'none';

        else
            % Numeric (old behavior)
            t0_nc = tvec(1);
            time_values_nc = double(tvec - t0_nc);

            time_units_nc  = 'days since simulation start';
            time_calendar_nc = 'none';
        end
    end

    %=========================================================
    % 4) INITIALIZE AGGREGATED NETCDF FILES (once)
    %=========================================================
    ncFile_hydro = '';
    ncVar_hydro  = '';
    ncFile_wq    = '';
    ncVar_wq     = '';
    ncFile_hr    = '';
    ncVar_hr     = '';

    ncHydroInitialized = false;
    ncWQInitialized    = false;
    ncHRInitialized    = false;

    if export_netcdf_timeseries
        if flags.flag_wse == 0
            ncFile_hydro = fullfile(myFolder_wd, 'Flood_Depths_timeseries.nc');
            ncVar_hydro  = 'flood_depth';
            ncLong_hydro = 'Flood depth';
            ncUnits_hydro = 'm';
        else
            ncFile_hydro = fullfile(myFolder_wd, 'Water_Surface_Elevation_timeseries.nc');
            ncVar_hydro  = 'water_surface_elevation';
            ncLong_hydro = 'Water surface elevation';
            ncUnits_hydro = 'm';
        end

        if flags.flag_waterquality == 1
            ncFile_wq = fullfile(myFolder_wq, 'Pollutant_Concentration_timeseries.nc');
            ncVar_wq  = 'pollutant_concentration';
            ncLong_wq = 'Pollutant concentration';
            ncUnits_wq = 'kg m-3';
        end

        if flags.flag_human_instability == 1
            ncFile_hr = fullfile(myFolder_hr, 'Human_Instability_timeseries.nc');
            ncVar_hr  = 'human_instability_risk';
            ncLong_hr = 'Human instability risk';
            ncUnits_hr = 'class';
        end
    end

    %========================
    % 5) Export time-varying rasters
    %========================
    store = 1;
    flag_loader = 1;
    mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));

    for i = 1:length(running_control.time_records)

        raster_exportion_percentage = i / length(running_control.time_records) * 100 

        %-----------------------------------------
        % Load maps in chunks (platform independent)
        %-----------------------------------------
        if i > saver_memory_maps * store
            store = store + 1;

            mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
            load(mapFile, 'Maps');

            % Changing NaN / Inf values
            Maps.Hydro.d(isnan(Maps.Hydro.d)) = no_data_value;
            Maps.Hydro.d(isinf(Maps.Hydro.d)) = no_data_value;
            Maps.Hydro.d(idx_nan) = no_data_value;

            Max_depth_d = max(max(Maps.Hydro.d, [], 3), Max_depth_d);

        else
            if flag_loader == 1
                load(mapFile, 'Maps');
                flag_loader = 0;

                % Changing NaN / Inf values
                Maps.Hydro.d(isnan(Maps.Hydro.d)) = no_data_value;
                Maps.Hydro.d(isinf(Maps.Hydro.d)) = no_data_value;
                Maps.Hydro.d(idx_nan) = no_data_value;

                Max_depth_d = max(Maps.Hydro.d, [], 3);
            end
        end

        local_i = i - ((store-1) * saver_memory_maps);

        % Threshold mask
        idx_depth = Maps.Hydro.d(:,:,local_i) < 0 * 1000; % currently deactivated

        if flags.flag_elapsed_time ~= 1
            time_map = datestr(running_control.time_records(i), 'yyyy_mm_dd_hh_MM_ss');
            time_stamp_str = datestr(running_control.time_records(i), 'yyyy-mm-dd HH:MM:ss');
        else
            time_map = running_control.time_records(i) / 60; % hours
            time_stamp_str = sprintf('Elapsed time = %.6f h', time_map);
        end

        %-----------------------------------------
        % HYDRO MAP: build filename + raster values
        %-----------------------------------------
        if flags.flag_wse == 0  % save water depth
            if flags.flag_elapsed_time == 1
                baseName = sprintf('Flood_Depths_t_%s_h', num2str(time_map));
            else
                baseName = sprintf('Flood_Depths_%s', string(time_map));
            end

            raster_exportion = Maps.Hydro.d(:,:,local_i) / 1000; % m
            raster_exportion(idx_nan) = no_data_value;
            raster_exportion(idx_depth) = no_data_value;

        else % save WSE
            if flags.flag_elapsed_time == 1
                baseName = sprintf('Water_Surface_Elevation_t_%s_h', num2str(time_map));
            else
                baseName = sprintf('Water_Surface_Elevation_%s', string(time_map));
            end

            raster_exportion = Maps.Hydro.d(:,:,local_i)/1000 + ...
                double(idx_Elevation_Properties.elevation_cell) .* Elevation_Properties.elevation_cell;

            raster_exportion(idx_nan) = no_data_value;
            raster_exportion(idx_depth) = no_data_value;
        end

        outTif = fullfile(myFolder_wd, baseName + ".tif");

        %-----------------------------------------
        % Export HYDRO (subgrid vs regular)
        %-----------------------------------------
        if flags.flag_subgrid ~= 1
            raster_to_export = DEM_raster;
            raster_to_export.Z = raster_exportion;
        else
            raster_to_export = DEM_raster_high_resolution;
            wse = raster_exportion + Subgrid_Properties.invert_el;
            % high_res_flood_map = ProjectFloodMap(DEM_raster_high_resolution, DEM_raster, wse);
            high_res_flood_map = ProjectFloodMap(DEM_raster_high_resolution, DEM_raster, wse);
            raster_to_export.Z = high_res_flood_map;
        end

        % GeoTIFF
        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

        % Aggregated NetCDF for HYDRO
        if export_netcdf_timeseries
            if ~ncHydroInitialized
                initializeTimeSeriesNetCDF( ...
                    ncFile_hydro, ...
                    raster_to_export.Z, ...
                    raster_to_export.georef.SpatialRef, ...
                    ncVar_hydro, ...
                    ncLong_hydro, ...
                    ncUnits_hydro, ...
                    no_data_value, ...
                    time_values_nc, ...
                    time_units_nc, ...
                    time_calendar_nc, ...
                    nc_deflate_level, ...
                    raster_to_export.georef.GeoKeyDirectoryTag);
                ncHydroInitialized = true;
            end

            appendTimeSliceToNetCDF(ncFile_hydro, ncVar_hydro, raster_to_export.Z, i);
        end

        %========================
        % Water Quality export (time series)
        %========================
        if flags.flag_waterquality == 1
            if flags.flag_elapsed_time == 1
                baseNameWQ = sprintf('Pollutant_Concentration_%smin', num2str(time_map));
            else
                baseNameWQ = sprintf('Pollutant_Concentration_%s', string(time_map));
            end

            outTifWQ = fullfile(myFolder_wq, baseNameWQ + ".tif");

            raster_exportion = Maps.WQ_States.Pol_Conc_Map(:,:,local_i);
            idx_ = raster_exportion < LULC_Properties.Pol_min;
            raster_exportion(idx_) = no_data_value;
            raster_exportion(isnan(raster_exportion)) = no_data_value;
            raster_exportion(isinf(raster_exportion)) = no_data_value;
            raster_exportion(raster_exportion < 0) = no_data_value;

            raster_to_export = DEM_raster;
            raster_to_export.Z = raster_exportion;

            geotiffwrite(outTifWQ, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
                'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

            if export_netcdf_timeseries
                if ~ncWQInitialized
                    initializeTimeSeriesNetCDF( ...
                        ncFile_wq, ...
                        raster_to_export.Z, ...
                        raster_to_export.georef.SpatialRef, ...
                        ncVar_wq, ...
                        ncLong_wq, ...
                        ncUnits_wq, ...
                        no_data_value, ...
                        time_values_nc, ...
                        time_units_nc, ...
                        time_calendar_nc, ...
                        nc_deflate_level, ...
                        raster_to_export.georef.GeoKeyDirectoryTag);
                    ncWQInitialized = true;
                end

                appendTimeSliceToNetCDF(ncFile_wq, ncVar_wq, raster_to_export.Z, i);
            end
        end

        %========================
        % Human instability export (time series)
        %========================
        if flags.flag_human_instability == 1
            if flags.flag_elapsed_time == 1
                baseNameHR = sprintf('Human_instability_%smin', num2str(time_map));
            else
                baseNameHR = sprintf('Human_instability_%s', string(time_map));
            end

            outTifHR = fullfile(myFolder_hr, baseNameHR + ".tif");

            raster_exportion = double(Maps.Hydro.risk(:,:,local_i));
            raster_exportion(isnan(raster_exportion)) = no_data_value;
            raster_exportion(isinf(raster_exportion)) = no_data_value;
            raster_exportion(raster_exportion < 0) = no_data_value;
            raster_exportion(raster_exportion == 0) = no_data_value;

            raster_to_export = DEM_raster;
            raster_to_export.Z = raster_exportion;

            geotiffwrite(outTifHR, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
                'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

            if export_netcdf_timeseries
                if ~ncHRInitialized
                    initializeTimeSeriesNetCDF( ...
                        ncFile_hr, ...
                        raster_to_export.Z, ...
                        raster_to_export.georef.SpatialRef, ...
                        ncVar_hr, ...
                        ncLong_hr, ...
                        ncUnits_hr, ...
                        no_data_value, ...
                        time_values_nc, ...
                        time_units_nc, ...
                        time_calendar_nc, ...
                        nc_deflate_level, ...
                        raster_to_export.georef.GeoKeyDirectoryTag);
                    ncHRInitialized = true;
                end

                appendTimeSliceToNetCDF(ncFile_hr, ncVar_hr, raster_to_export.Z, i);
            end
        end
    end

    %=========================================================
    % 6) STATIC EXPORTS -> TIFF ONLY (NO NETCDF)
    %=========================================================

    %--------------------------------------
    % Initial Buildup Map
    %--------------------------------------
    if flags.flag_waterquality == 1
        outTif = fullfile(Dirs.RastersStatic, "Initial_Buildup_kg.tif");

        raster_exportion = Maps.WQ_States.initial_buildup_map;
        raster_exportion(isnan(raster_exportion)) = no_data_value;
        raster_exportion(isinf(raster_exportion)) = no_data_value;
        raster_exportion(raster_exportion < 0) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = raster_exportion;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

        % Total Washed Mass
        outTif = fullfile(Dirs.RastersStatic, "Total_Washed_Mass_Kg.tif");

        raster_exportion = WQ_States.Tot_Washed;
        raster_exportion(isnan(raster_exportion)) = no_data_value;
        raster_exportion(isinf(raster_exportion)) = no_data_value;
        raster_exportion(raster_exportion < 0) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = raster_exportion;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end

    %--------------------------------------
    % Points of accumulation of Depths
    %--------------------------------------
    outTif = fullfile(Dirs.RastersStatic, "Accumulation_Areas_Larger_1m.tif");

    idx_depth = depths.d_t > 1 * 1000;
    raster_exportion = no_data_value * ones(size(depths.d_t));
    raster_exportion(idx_depth) = 1;

    raster_to_export = DEM_raster;
    raster_to_export.Z = raster_exportion;

    geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
        'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

    %--------------------------------------
    % Points of accumulation of pollutants
    %--------------------------------------
    if flags.flag_waterquality == 1
        outTif = fullfile(Dirs.RastersStatic, "Accumulation_Areas_10_g_m2.tif");

        pol_accumulation = 10; % g/m2
        zzz = WQ_States.B_t / Wshed_Properties.cell_area * 1000; % g/m2
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;
        idx_bt = zzz < pol_accumulation;
        zzz(idx_bt) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = zzz;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

        % Final Pollutant Mass
        outTif = fullfile(Dirs.RastersStatic, "Final_Pollutant_Mass_g_m2.tif");

        raster_exportion = no_data_value * ones(size(WQ_States.B_t));
        final_mass = Maps.WQ_States.Pol_mass_map(:,:,end) / Wshed_Properties.cell_area * 1000; % g/m2
        idx_bt = final_mass > 0;
        raster_exportion(idx_bt) = final_mass(idx_bt);
        raster_exportion(isinf(raster_exportion)) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = raster_exportion;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end

    %--------------------------------------
    % Maximum Velocity
    %--------------------------------------
    outTif = fullfile(Dirs.RastersStatic, "Maximum_Velocity.tif");

    zzz = velocities.vmax_final; % m/s
    idx_wse = depths.dmax_final/1000 < depths.depth_wse;
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;
    zzz(idx_wse) = no_data_value;

    raster_to_export = DEM_raster;
    raster_to_export.Z = zzz;

    geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
        'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

    %--------------------------------------
    % Risk Map
    %--------------------------------------
    if flags.flag_human_instability == 1
        store = 1;
        flag_loader = 1;

        zzz = -inf(size(DEM_raster.Z,1), size(DEM_raster.Z,2));

        for i = 1:length(running_control.time_records)
            if i > saver_memory_maps * store
                store = store + 1;
                mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
                load(mapFile, 'Maps');
            else
                if flag_loader == 1
                    mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
                    load(mapFile, 'Maps');
                    flag_loader = 0;
                end
            end

            zzz = max(zzz, max(Maps.Hydro.risk, [], 3));
        end

        zzz = Human_Instability.max_risk;

        idx_wse = depths.dmax_final/1000 < depths.depth_wse;
        outTif = fullfile(Dirs.RastersStatic, "Maximum_Instability_Risk.tif");

        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;
        zzz(idx_wse) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = zzz;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

    elseif flags.flag_human_instability == 2
        % Leave as-is / your implementation

    elseif flags.flag_human_instability == 3
        list = {'_cm','_tm','_am','_om','_cf','_tf','_af','_of'};
        risk_summary = table('Size', [3 9], ...
            'VariableNames', {'Risk','risk_cm', 'risk_tm', 'risk_am', 'risk_om', 'risk_cf', 'risk_tf', 'risk_af', 'risk_of'}, ...
            'VariableTypes', {'string','double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'});

        risk_summary.Risk(1) = "Slide";
        risk_summary.Risk(2) = "Toppling";
        risk_summary.Risk(3) = "Drowning";

        for k = 1:8
            store = 1;
            flag_loader = 1;

            for i = 1:length(running_control.time_records)
                if i > saver_memory_maps * store
                    store = store + 1;
                    mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
                    load(mapFile, 'Maps');
                else
                    if flag_loader == 1
                        mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
                        load(mapFile, 'Maps');
                        flag_loader = 0;
                    end
                end

                local_i = i - ((store-1) * saver_memory_maps);
                zzz = double(Maps.Hydro.(strcat('risk', list{k}))(:,:,local_i));

                risk_summary.(strcat('risk',list{k}))(1) = max(risk_summary.(strcat('risk',list{k}))(1), sum(zzz(:) == 1));
                risk_summary.(strcat('risk',list{k}))(2) = max(risk_summary.(strcat('risk',list{k}))(2), sum(zzz(:) == 2));
                risk_summary.(strcat('risk',list{k}))(3) = max(risk_summary.(strcat('risk',list{k}))(3), sum(zzz(:) == 3));
            end

            zzz = Human_Instability.max_risk;
            outTif = fullfile(Dirs.RastersStatic, "Maximum_Instability_Risk" + string(list{k}) + ".tif");

            zzz(isinf(zzz)) = no_data_value;
            zzz(isnan(zzz)) = no_data_value;

            raster_to_export = DEM_raster;
            raster_to_export.Z = zzz;

            geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
                'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
        end

        risk_summary(:,2:9) = risk_summary(:,2:9) .* (Wshed_Properties.Resolution^2) ./ 1e6;
        writetable(risk_summary, fullfile(Dirs.Tables, 'Risk_summary.csv'));
    end

    %--------------------------------------
    % Infiltrated Depth
    %--------------------------------------
    outTif = fullfile(Dirs.RastersStatic, "Infiltrated_Depth.tif");

    zzz = Soil_Properties.I_t;
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;

    raster_to_export = DEM_raster;
    raster_to_export.Z = zzz;

    geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
        'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

    %--------------------------------------
    % Maximum Depths / Max WSE
    %--------------------------------------
    if flags.flag_subgrid ~= 1
        raster_to_export = DEM_raster;
        zzz = depths.dmax_final / 1000;
        idx_wse = zzz < depths.depth_wse;
    else
        raster_to_export = DEM_raster_high_resolution;
        zzz = depths.dmax_final / 1000;
        wse = zzz + Subgrid_Properties.invert_el;
        % high_res_flood_map = ProjectFloodMap(DEM_raster_high_resolution, DEM_raster, wse);
        high_res_flood_map = ProjectFloodMap(DEM_raster_high_resolution, DEM_raster, wse);
        raster_to_export.Z = high_res_flood_map;
        zzz = raster_to_export.Z;
        idx_wse = zzz < depths.depth_wse;
    end

    if flags.flag_wse == 0
        outTif = fullfile(Dirs.RastersStatic, "Maximum_Depths.tif");

        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;
        zzz(idx_wse) = no_data_value;

        raster_to_export.Z = zzz;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    else
        outTif = fullfile(Dirs.RastersStatic, "Max_Water_Surface_Elevation.tif");

        raster_exportion = zzz + Elevation_Properties.elevation_cell;
        raster_exportion(idx_wse) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = raster_exportion;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end

    %--------------------------------------
    % Maximum Snowpack
    %--------------------------------------
    if flags.flag_snow_modeling == 1
        outTif = fullfile(Dirs.RastersStatic, "Maximum_Snowpack.tif");

        zzz = max_Hsnow / 1000;
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = zzz;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end

    %--------------------------------------
    % Total Infiltration
    %--------------------------------------
    if flags.flag_infiltration == 1
        outTif = fullfile(Dirs.RastersStatic, "Cumulative_Infiltration.tif");

        zzz = cumulative_infiltration / 1000;
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = zzz;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end

    %--------------------------------------
    % Maximum GW Depth
    %--------------------------------------
    if flags.flag_groundwater_modeling == 1
        outTif = fullfile(Dirs.RastersStatic, "Max_GW_depth.tif");

        zzz = max_GW_depth;
        zzz(isinf(zzz)) = no_data_value;
        zzz(isnan(zzz)) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = zzz;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end

    %--------------------------------------
    % DEM
    %--------------------------------------
    if flags.flag_resample == 1
        outTif = fullfile(Dirs.RastersStatic, "DEM_resampled.tif");
    else
        outTif = fullfile(Dirs.RastersStatic, "DEM_Treated.tif");
    end

    zzz = Elevation_Properties.elevation_cell;
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;

    raster_to_export = DEM_raster;
    raster_to_export.Z = zzz;

    geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
        'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

    %--------------------------------------
    % LULC
    %--------------------------------------
    outTif = fullfile(Dirs.RastersStatic, "Land_Cover_Data.tif");

    zzz = LULC_Properties.LULC;
    zzz(isinf(zzz)) = no_data_value;
    zzz(isnan(zzz)) = no_data_value;

    raster_to_export = DEM_raster;
    raster_to_export.Z = zzz;

    geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
        'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);

    %--------------------------------------
    % Final WQ maps
    %--------------------------------------
    if flags.flag_waterquality == 1
        store = 1;
        flag_loader = 1;

        for i = 1:length(running_control.time_records)
            if i == length(running_control.time_records)

                outTif = fullfile(Dirs.RastersStatic, "Final_Mass_Of_Pollutant.tif");

                if i > saver_memory_maps * store
                    store = store + 1;
                    mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
                    load(mapFile, 'Maps');
                    Max_Pol_Conc_Map = max(max(Maps.WQ_States.Pol_Conc_Map, [], 3), Max_Pol_Conc_Map);
                else
                    if flag_loader == 1
                        mapFile = fullfile(tempDir, sprintf('save_map_hydro_%d', store));
                        load(mapFile, 'Maps');
                        flag_loader = 0;
                        Max_Pol_Conc_Map = max(Maps.WQ_States.Pol_Conc_Map, [], 3);
                    end
                end

                local_i = i - ((store-1) * saver_memory_maps);
                raster_exportion = Maps.WQ_States.Pol_mass_map(:,:,local_i);

                idx = raster_exportion < 0.01;
                raster_exportion(idx) = no_data_value;

                raster_exportion = raster_exportion / Wshed_Properties.cell_area * 1000; % g/m2

                raster_to_export = DEM_raster;
                raster_to_export.Z = raster_exportion;

                geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
                    'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
            end
        end

        outTif = fullfile(Dirs.RastersStatic, "Maximum_Pol_Conc_min.tif");

        raster_exportion = Max_Pol_Conc_Map;

        idx = raster_exportion < LULC_Properties.Pol_min;
        raster_exportion(idx) = no_data_value;
        raster_exportion(isnan(raster_exportion)) = no_data_value;
        raster_exportion(isinf(raster_exportion)) = no_data_value;
        raster_exportion(raster_exportion < 0) = no_data_value;

        raster_to_export = DEM_raster;
        raster_to_export.Z = raster_exportion;

        geotiffwrite(outTif, raster_to_export.Z, raster_to_export.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', raster_to_export.georef.GeoKeyDirectoryTag);
    end
end


%% Generate GIF Files of the Simulation
Inundation_Maps

%% Summary Table
%%% - Equivalent Width to Simulate in SWMM %%% - (LENHS, 2012)
% W = kc sqrt(A) / 1.12 * (1 - sqrt(1 - {1.128 / k_c}^2))
Wshed_Properties.width_SWMM = Wshed_Properties.compactness_coefficient*sqrt(Wshed_Properties.drainage_area)/1.12*(1 - sqrt(1 - (1.128/Wshed_Properties.compactness_coefficient)^2));
%%% - Runoff Coefficient -

try
if flags.flag_inflow == 1
    Wshed_Properties.C_r = BC_States.outflow_volume/BC_States.inflow_volume;
elseif flags.flag_spatial_rainfall == 1
    Wshed_Properties.C_r = BC_States.outflow_volume/(sum(BC_States.average_spatial_rainfall / 1000)*time_step_rainfall * Wshed_Properties.drainage_area);
else
    Wshed_Properties.C_r = BC_States.outflow_volume/(sum(BC_States.delta_p) / 1000 * Wshed_Properties.drainage_area);
end
catch
    Wshed_Properties.C_r = nan;
end
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
FileName = fullfile(Dirs.Tables, [FileName_String '.csv']);

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


% %% Maps and Hydrographs
% close all
%
% % ---- Use a FIXED pixel size (faster + consistent) instead of fullscreen ----
% fig = figure('Color','w','Units','pixels','Position',[50 50 1600 1000], ...
%     'MenuBar','none','ToolBar','none');
% set(fig,'DefaultTextInterpreter','latex');
%
% size_font = 10;
%
% if flags.flag_waterquality == 1
%     size_font = 12;
%
%     [axis1] = surfplot_maps(DEM_raster,depths.dmax_final/1000,Spectrum, ...
%         'Easting (m)','Northing (m)','Depth (m)',no_data_value,idx_nan,3,3,1,size_font);
%
%     min_washed = 1e-4;
%
%     z = WQ_States.Tot_Washed;
%     z(z<=min_washed)=nan;
%     [axis2] = surfplot_maps(DEM_raster,z,Spectrum, ...
%         'Easting (m)','Northing (m)','Total Washed Mass (kg)',no_data_value,idx_nan,3,3,2,size_font);
%
%     [axis3] = surfplot_maps(DEM_raster,Soil_Properties.I_t,WSE_RAS, ...
%         'Easting (m)','Northing (m)','Infiltration (mm)',no_data_value,idx_nan,3,3,3,size_font);
%
%     [axis4] = surfplot_maps(DEM_raster,velocities.vmax_final,Velocity_RAS, ...
%         'Easting (m)','Northing (m)','Max. Velocity (m/s)',no_data_value,idx_nan,3,3,4,size_font);
%
%     [axis5] = surfplot_maps(DEM_raster,Maps.WQ_States.Pol_mass_map(:,:,end)/Wshed_Properties.cell_area*1000,Velocity_RAS, ...
%         'Easting (m)','Northing (m)','Pollutant Mass at the end ($\mathrm{g/m^2}$)',no_data_value,idx_nan,3,3,5,size_font);
%
%     [axis6] = surfplot_maps(DEM_raster,Maps.WQ_States.Pol_mass_map(:,:,end)/Wshed_Properties.cell_area*1000,Spectrum, ...
%         'Easting (m)','Northing (m)','Pollutant Mass at the end ($\mathrm{mg/L}$)',no_data_value,idx_nan,3,3,6,size_font);
%
%     subplot(3,3,7)
%     plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph), ...
%         'LineWidth',1.5,'color','black');
%     xlabel('Time [min]','interpreter','latex');
%     ylabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex');
%     set(gca,'FontSize',size_font);
%     hold on
%     yyaxis right; set(gca,'ydir','reverse','ycolor','black');
%     bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall), ...
%         'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
%     ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
%     ylim([0,max(gather(Rainfall_Parameters.intensity_rainfall)) + 200]);
%
%     if flags.flag_inflow == 1
%         hold on
%         yyaxis left; set(gca,'ydir','normal','ycolor','black');
%         plot(gather(Inflow_Parameters.time_inflow(1:tfinal_inflow,1)), ...
%              gather(Inflow_Parameters.inflow_discharge(1:tfinal_inflow,:)), ...
%              'LineWidth',1.5,'color','black','LineStyle','--');
%         xlabel('Time [min]','interpreter','latex');
%         grid on
%         legend('Outflow','Inflow','Rainfall','interpreter','latex');
%     else
%         legend('Outflow','Rainfall','interpreter','latex');
%     end
%
%     subplot(3,3,8)
%     plot(gather(running_control.time_hydrograph),gather(WQ_States.outet_pollutograph), ...
%         'LineWidth',1.5,'color',[255,140,0]/255);
%     xlabel('Time [min]','Interpreter','Latex');
%     ylabel('Pol. Concentration (mg/L)','Interpreter','Latex')
%     grid on
%     hold on
%     yyaxis right
%     set(gca,'ycolor','black')
%     ylabel('Load (kg/sec)','interpreter','latex')
%     load_wq = 1e-3*gather(WQ_States.outet_pollutograph).*gather(outlet_states.outlet_hydrograph); % kg/sec
%     plot(gather(running_control.time_hydrograph),load_wq,'LineWidth',1.5,'color',[34,139,34]/255)
%     legend('Concentration','Load','interpreter','latex')
%
%     subplot(3,3,9)
%     m_M = WQ_States.mass_outlet_save./(max(max(WQ_States.mass_outlet)));
%     v_V = WQ_States.vol_outlet_save./(max(max(WQ_States.vol_outlet)));
%     plot(v_V,m_M,'LineWidth',1.5,'color','black','Marker','*');
%     hold on
%     plot(0:1,0:1,'LineWidth',1,'Color','black','LineStyle','--')
%     ylabel('$m/m_{\mathrm{tot}}$','Interpreter','Latex');
%     xlabel('$V/V_{\mathrm{tot}}$','Interpreter','Latex');
%
%     % ---- FAST adaptive export (caps output megapixels) ----
%     exportFast(gcf, fullfile(folderName,'Maps.png'), 300, 8);  % 8 MP cap
%
% else
%     % No water quality
%     size_font = 12;
%     no_data_value = nan;
%
%     [axis1] = surfplot_maps(DEM_raster,depths.dmax_final/1000,Spectrum, ...
%         'Easting (m)','Northing (m)','Maximum Depth (m)',no_data_value,idx_nan,2,2,1,size_font);
%
%     [axis2] = surfplot_maps(DEM_raster,Soil_Properties.I_t,Spectrum, ...
%         'Easting (m)','Northing (m)','Infiltration (mm)',no_data_value,idx_nan,2,2,2,size_font);
%
%     [axis3] = surfplot_maps(DEM_raster,velocities.vmax_final,Velocity_RAS, ...
%         'Easting (m)','Northing (m)','Max. Velocity (m/s)',no_data_value,idx_nan,2,2,3,size_font);
%
%     ax4 = subplot(2,2,4);
%     plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph), ...
%         'LineWidth',1.5,'color','black','marker','*');
%     xlabel('Time [min]','interpreter','latex');
%     ylabel('Flow Discharge $(\mathrm{m^3/s})$','interpreter','latex');
%     set(gca,'FontSize',size_font);
%     hold on
%     yyaxis right; set(gca,'ydir','reverse','ycolor','black');
%
%     if flags.flag_spatial_rainfall ~= 1 && flags.flag_rainfall == 1
%         bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall), ...
%             'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
%         ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
%         ylim([0,max(gather(Rainfall_Parameters.intensity_rainfall)) + 200]);
%     elseif flags.flag_rainfall == 1
%         bar(gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1:(length(BC_States.average_spatial_rainfall)))), ...
%             gather(BC_States.average_spatial_rainfall), ...
%             'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
%         ylabel('Intensity [$\mathrm{mm \cdot h^{-1}}$]','interpreter','latex');
%         ylim([0,max(gather(BC_States.average_spatial_rainfall)) + 200]);
%     end
%
%     if flags.flag_inflow == 1
%         hold on
%         yyaxis left; set(gca,'ydir','normal','ycolor','black')
%         tfinal_inflow = length(Inflow_Parameters.time_inflow);
%         plot(gather(Inflow_Parameters.time_inflow(1:tfinal_inflow,1)), ...
%              gather(Inflow_Parameters.inflow_discharge(1:tfinal_inflow,:)), ...
%              'LineWidth',1.5,'color','black','LineStyle','--');
%         xlabel('Time [min]','interpreter','latex');
%
%         % Robust ylim (avoid invalid 1x2)
%         y1 = max(gather(outlet_states.outlet_hydrograph))*1.2;
%         y2 = max(gather(Inflow_Parameters.inflow_discharge(:)));
%         yTop = max([y1 y2]);
%         if ~isfinite(yTop) || yTop <= 0, yTop = 1; end
%         ylim([0 yTop]);
%
%         grid on
%         legend('Outflow','Inflow','Rainfall');
%     else
%         legend('Outflow','Rainfall');
%     end
%
%     set(gca, 'FontName', 'Garamond', 'FontSize', 12)
%     set(gca, 'TickLength', [0.02 0.01]);
%     title('Hydrograph','Interpreter','Latex');
%     set(gca,'Tickdir','out')
%
%     % ---- FAST export for the no-WQ summary too (optional but recommended) ----
%     exportFast(gcf, fullfile(folderName,'Maps.png'), 250, 8);
% end
%
%
% %% Input Data Maps
% % DEM, LULC, SOIL, n, h0, ksat,
% close all
%
% % ---- Fixed pixel size instead of fullscreen ----
% fig = figure('Color','w','Units','pixels','Position',[50 50 1600 700], ...
%     'MenuBar','none','ToolBar','none');
% set(fig,'DefaultTextInterpreter','latex');
%
% % MAIN INPUT MAPS
% size_font = 12;
% [axis1] = surfplot_maps(DEM_raster,Elevation_Properties.elevation_cell,Terrain_RAS_ramp, ...
%     'Easting (m)','Northing (m)','Elevation (m)',no_data_value,idx_nan,1,3,1,size_font);
%
% [axis2] = surfplot_maps(DEM_raster,LULC_Properties.LULC,linspecer(LULC_Properties.n_lulc), ...
%     'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,1,3,2,size_font);
%
% [axis3] = surfplot_maps(DEM_raster,Soil_Properties.soil_matrix,linspecer(Soil_Properties.n_soil), ...
%     'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,1,3,3,size_font);
%
% % Vector PDF is fine here (usually), keep it:
% exportgraphics(gcf, fullfile(folderName, 'Input_Data_Maps.pdf'), 'ContentType', 'vector');
% saveas(gcf,fullfile(folderName,'Input_Data_Maps.fig'))
% close all
%
%
% fig = figure('Color','w','Units','pixels','Position',[50 50 1600 900], ...
%     'MenuBar','none','ToolBar','none');
% set(fig,'DefaultTextInterpreter','latex');
%
% % --- Manning / h0 / WQ params ----
% if flags.flag_waterquality ~= 1
%     [axis1] = surfplot_maps(DEM_raster,LULC_Properties.roughness,Terrain_RAS_ramp, ...
%         'Easting (m)','Northing (m)','$\mathrm{n~(sm^{-1/3})}    $',no_data_value,idx_nan,4,2,1,size_font);
%
%     [axis2] = surfplot_maps(DEM_raster,LULC_Properties.h_0,linspecer, ...
%         'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,4,2,2,size_font);
% else
%     [axis1] = surfplot_maps(DEM_raster,LULC_Properties.roughness,Terrain_RAS_ramp, ...
%         'Easting (m)','Northing (m)','$\mathrm{n~(sm^{-1/3})}$',no_data_value,idx_nan,4,2,3,size_font);
%
%     [axis2] = surfplot_maps(DEM_raster,LULC_Properties.h_0,linspecer, ...
%         'Easting (m)','Northing (m)','Classification',no_data_value,idx_nan,4,2,4,size_font);
%
%     [axis3] = surfplot_maps(DEM_raster,LULC_Properties.C_1,linspecer, ...
%         'Easting (m)','Northing (m)','$C_1~(\mathrm{kg.ha^{-1}}$)',no_data_value,idx_nan,4,2,5,size_font);
%
%     [axis4] = surfplot_maps(DEM_raster,LULC_Properties.C_2,linspecer, ...
%         'Easting (m)','Northing (m)','$C_2~(\mathrm{day^{-1}}$)',no_data_value,idx_nan,4,2,6,size_font);
%
%     [axis5] = surfplot_maps(DEM_raster,LULC_Properties.C_3,linspecer, ...
%         'Easting (m)','Northing (m)','$C_3$~(-)',no_data_value,idx_nan,4,2,7,size_font);
%
%     [axis6] = surfplot_maps(DEM_raster,LULC_Properties.C_4,linspecer, ...
%         'Easting (m)','Northing (m)','$C_4~(-)',no_data_value,idx_nan,4,2,8,size_font);
% end
%
% % If you want an image export here, use fast cap (optional):
% % exportFast(gcf, fullfile(folderName,'LULC_Parameters.png'), 250, 8);
% close all
%
%
% fig = figure('Color','w','Units','pixels','Position',[50 50 1600 700], ...
%     'MenuBar','none','ToolBar','none');
% set(fig,'DefaultTextInterpreter','latex');
%
% % Ksat
% [axis1] = surfplot_maps(DEM_raster,Soil_Properties.ksat,Spectrum, ...
%     'Easting (m)','Northing (m)','$\mathrm{k_{sat}~(mm.h^{-1})}    $',no_data_value,idx_nan,1,3,1,size_font);
%
% % Dtheta
% dheta = Soil_Properties.teta_sat - Soil_Properties.teta_i;
% [axis2] = surfplot_maps(DEM_raster,dheta,Spectrum, ...
%     'Easting (m)','Northing (m)','$\mathrm{\Delta \theta~(-)}    $',no_data_value,idx_nan,1,3,2,size_font);
%
% % Psi
% [axis3] = surfplot_maps(DEM_raster,Soil_Properties.psi,Spectrum, ...
%     'Easting (m)','Northing (m)','$\mathrm{\psi~(mm)}    $',no_data_value,idx_nan,1,3,3,size_font);
%
% % Keep this one lighter; you already used 150:
% exportFast(gcf, fullfile(folderName,'SOIL_Parameters.png'), 150, 6);
% close all
%
%
% %% DEM with Rivers
% close all
% tiledlayout(2,1)
% set(gcf,'units','inches','position',[2,0,8,6])
%
% ax1 = nexttile;
%
% if no_plot==0
%     mapshow(A,RA,"AlphaData",0.45);hold on;
%     mapshow(S_p,'FaceColor','n'); hold on;
% end
%
% imagesc(hillshade(DEM_raster)); hold on;
% colormap(ax1,'gray')
% surf(DEM_raster);
% colormap(ax1,'landcolor')
% h = colorbar;
% h.Label.String = 'Elevation (m)';
% h.TickLabelInterpreter = 'latex';
% xlabel('Easting (m)','interpreter','latex')
% ylabel('Northing (m)','interpreter','latex')
% ax = ancestor(ax1, 'axes');
% ax.XAxis.Exponent = 0; xtickformat('%.0f');
% ax.YAxis.Exponent = 0; ytickformat('%.0f');
%
% ax2 = nexttile;
% FD = FLOWobj(DEM_raster);
% As = flowacc(FD);
%
% if no_plot==0
%     mapshow(A,RA,"AlphaData",0.45);hold on;
%     mapshow(S_p,'FaceColor','n'); hold on;
% end
%
% imagesc(hillshade(DEM_raster)); hold on;
% colormap(ax2,'gray')
% F = dilate(sqrt(As),ones(5));
% F.Z(F.Z==1)=NaN;
% surf(F);
% colormap(ax2,flowcolor);
% h = colorbar;
% h.Label.String = '(sqrt(# of pixels)';
% h.TickLabelInterpreter = 'latex';
% xlabel('Easting (m)','interpreter','latex')
% ylabel('Northing (m)','interpreter','latex')
% ax = ancestor(ax2, 'axes');
% ax.XAxis.Exponent = 0; xtickformat('%.0f');
% ax.YAxis.Exponent = 0; ytickformat('%.0f');
% close all
%
%
% clearvars a_grid area_cells area_km2 b_grid baseFileName C cm color_plot color_plots depth_accumulation Depth_RAS elevation f FileName FileName_String filePattern FolderName font_size frame fsize fullFileName h h_max h_min i idx2 idx3 idx_depth idx_i_a idx_wse im imind labels_depth labels_gauges ls max_depth max_h max_inf max_v MS myFolder no_data_value nx_max ny_max Out_Conc points raster_exportion raster_exportion_percentage s size_font t t_max t_previous t_save t_store t_title theFiles topo_path x_grid xbrgin xend xmax y_grid ybegin yend ymax z z1 z2 zmax zmin
%
% disp('Thank you for using HydroPol2D. Results are exported in Modeling Results folder.')

%% Deleting temporary files

% files_to_delete = dir('Temporary_Files');
% for k = 1 : length(files_to_delete)
%     baseFileName = files_to_delete(k).name;
%     fullFileName = fullfile(tempDir, baseFileName);;
%     fprintf(1, 'Now deleting %s\n', fullFileName);
%     delete(fullFileName);
% end

%% ===== Local helper function(s) =====

function exportFast(fig, outPath, baseDPI, maxMegapixels)
% exportFast  Adaptive exportgraphics that caps output megapixels.
%
% fig: figure handle
% outPath: full output path (png/jpg/etc)
% baseDPI: typical 300 (or 150/200 for faster)
% maxMegapixels: cap output size (e.g., 6–10)

if nargin < 3 || isempty(baseDPI), baseDPI = 300; end
if nargin < 4 || isempty(maxMegapixels), maxMegapixels = 8; end

% Figure size in pixels (screen)
oldUnits = fig.Units;
fig.Units = 'pixels';
pos = fig.Position;  % [x y w h]
fig.Units = oldUnits;

w = max(1, pos(3));
h = max(1, pos(4));

% Approximate scaling: export DPI relative to ~96dpi screen
screenDPI = 96;
scale = baseDPI / screenDPI;

outW = w * scale;
outH = h * scale;
mp = (outW * outH) / 1e6;

if mp > maxMegapixels
    scale2 = sqrt(maxMegapixels / mp);
    dpi = max(72, round(baseDPI * scale2));
else
    dpi = baseDPI;
end

exportgraphics(fig, outPath, ...
    'ContentType','image', ...
    'Colorspace','rgb', ...
    'Resolution', dpi);
end


function deleteMatchingFiles(folderPath, patterns)
% Delete files matching wildcard patterns in a folder.

if ~isfolder(folderPath)
    return;
end

if ischar(patterns) || isstring(patterns)
    patterns = cellstr(patterns);
end

for p = 1:numel(patterns)
    files = dir(fullfile(folderPath, patterns{p}));
    for k = 1:numel(files)
        thisFile = fullfile(folderPath, files(k).name);
        if exist(thisFile, 'file')
            try
                delete(thisFile);
            catch
                warning('Could not delete file: %s', thisFile);
            end
        end
    end
end
end

function writeRasterNetCDF(ncFile, Z, R, varName, longName, units, fillValue, timeStampStr, deflateLevel, geoKeyTag)
% =========================================================
% Write one 2D raster to NetCDF
%
% INPUTS
%   ncFile        -> output NetCDF path
%   Z             -> 2D raster matrix
%   R             -> raster reference object
%   varName       -> NetCDF variable name
%   longName      -> descriptive variable name
%   units         -> variable units
%   fillValue     -> missing value
%   timeStampStr  -> descriptive timestamp string
%   deflateLevel  -> NetCDF compression level (0-9)
%   geoKeyTag     -> GeoTIFF CRS metadata (stored as text attribute if possible)
%
% NOTES
%   - This function is self-contained, so it does not depend on any other
%     helper function.
%   - The raster is written as [x,y] dimensions using Z' so that the first
%     dimension matches x and the second matches y.
% =========================================================

if exist(ncFile, 'file')
    delete(ncFile);
end

Z = double(Z);
[nRows, nCols] = size(Z);

% -----------------------------------------------------
% Detect whether raster reference is geographic/projected
% -----------------------------------------------------
isGeographic = isa(R, 'map.rasterref.GeographicCellsReference') || ...
    isa(R, 'map.rasterref.GeographicPostingsReference');

if isGeographic
    xName   = 'lon';
    yName   = 'lat';
    crsName = 'geographic';
    crsUnits = 'degrees';
else
    xName   = 'x';
    yName   = 'y';
    crsName = 'projected';
    crsUnits = 'm';
end

% -----------------------------------------------------
% Build coordinate vectors from raster reference
% -----------------------------------------------------
% Column centers
colCenters = 1:nCols;

% Row centers
rowCenters = 1:nRows;

% X coordinates from first row
[xCoords, ~] = intrinsicToWorld(R, colCenters, ones(1, nCols));

% Y coordinates from first column
[~, yCoords] = intrinsicToWorld(R, ones(1, nRows), rowCenters);

xCoords = double(xCoords(:));
yCoords = double(yCoords(:));

% -----------------------------------------------------
% Create dimensions and coordinate variables
% -----------------------------------------------------
nccreate(ncFile, xName, ...
    'Dimensions', {xName, numel(xCoords)}, ...
    'Datatype', 'double');

nccreate(ncFile, yName, ...
    'Dimensions', {yName, numel(yCoords)}, ...
    'Datatype', 'double');

% -----------------------------------------------------
% Create data variable
% -----------------------------------------------------
nccreate(ncFile, varName, ...
    'Dimensions', {xName, numel(xCoords), yName, numel(yCoords)}, ...
    'Datatype', 'double', ...
    'FillValue', fillValue, ...
    'Format', 'netcdf4_classic', ...
    'DeflateLevel', deflateLevel, ...
    'Shuffle', true);

% -----------------------------------------------------
% Write coordinates and raster
% -----------------------------------------------------
ncwrite(ncFile, xName, xCoords);
ncwrite(ncFile, yName, yCoords);
ncwrite(ncFile, varName, Z.');

% -----------------------------------------------------
% Coordinate attributes
% -----------------------------------------------------
if isGeographic
    ncwriteatt(ncFile, xName, 'standard_name', 'longitude');
    ncwriteatt(ncFile, xName, 'long_name', 'longitude');
    ncwriteatt(ncFile, xName, 'units', 'degrees_east');

    ncwriteatt(ncFile, yName, 'standard_name', 'latitude');
    ncwriteatt(ncFile, yName, 'long_name', 'latitude');
    ncwriteatt(ncFile, yName, 'units', 'degrees_north');
else
    ncwriteatt(ncFile, xName, 'standard_name', 'projection_x_coordinate');
    ncwriteatt(ncFile, xName, 'long_name', 'x coordinate of projection');
    ncwriteatt(ncFile, xName, 'units', crsUnits);

    ncwriteatt(ncFile, yName, 'standard_name', 'projection_y_coordinate');
    ncwriteatt(ncFile, yName, 'long_name', 'y coordinate of projection');
    ncwriteatt(ncFile, yName, 'units', crsUnits);
end

% -----------------------------------------------------
% Data variable attributes
% -----------------------------------------------------
ncwriteatt(ncFile, varName, 'long_name', longName);
ncwriteatt(ncFile, varName, 'units', units);
ncwriteatt(ncFile, varName, 'coordinates', sprintf('%s %s', xName, yName));

% -----------------------------------------------------
% Global attributes
% -----------------------------------------------------
ncwriteatt(ncFile, '/', 'title', longName);
ncwriteatt(ncFile, '/', 'source', 'HydroPol2D export');
ncwriteatt(ncFile, '/', 'history', ['Created on ' datestr(now, 'yyyy-mm-dd HH:MM:ss')]);
ncwriteatt(ncFile, '/', 'timestamp', timeStampStr);
ncwriteatt(ncFile, '/', 'crs_name', crsName);
ncwriteatt(ncFile, '/', 'crs_units', crsUnits);

% Try to preserve GeoTIFF CRS info as text
try
    ncwriteatt(ncFile, '/', 'GeoKeyDirectoryTag', evalc('disp(geoKeyTag)'));
catch
    % do nothing
end
end

function initializeTimeSeriesNetCDF(ncFile, Z, R, varName, longName, units, ...
    fillValue, timeValues, timeUnits, timeCalendar, deflateLevel, geoKeyTag)
% initializeTimeSeriesNetCDF
% ------------------------------------------------------------
% Creates one aggregated NetCDF file for a time-varying raster:
% dimensions = (x, y, time)
%
% INPUTS
%   ncFile        : full output path to .nc file
%   Z             : one sample 2D raster [nRows x nCols]
%   R             : raster reference object
%   varName       : NetCDF variable name
%   longName      : descriptive variable name
%   units         : variable units
%   fillValue     : missing value
%   timeValues    : vector of time values already converted to numeric
%   timeUnits     : e.g. 'hours since simulation start'
%   timeCalendar  : e.g. 'none' or 'proleptic_gregorian'
%   deflateLevel  : 0-9 compression
%   geoKeyTag     : GeoTIFF CRS metadata to store as text attribute
%
% NOTES
%   - Data variable is stored as:
%         varName(x, y, time)
%   - Each raster slice is written as Z.' with start = [1 1 i]
%   - x corresponds to columns
%   - y corresponds to rows
%   - CRS is taken directly from the DEM raster reference
% ------------------------------------------------------------

    if exist(ncFile, 'file')
        delete(ncFile);
    end

    Z = double(Z);
    [nRows, nCols] = size(Z);
    nTimes = numel(timeValues);
    timeValues = double(timeValues(:));

    % -----------------------------------------------------
    % Detect raster reference type
    % -----------------------------------------------------
    isGeographic = isa(R, 'map.rasterref.GeographicCellsReference') || ...
                   isa(R, 'map.rasterref.GeographicPostingsReference');

    % -----------------------------------------------------
    % Variable names and default metadata
    % -----------------------------------------------------
    if isGeographic
        xName    = 'lon';
        yName    = 'lat';
        xStdName = 'longitude';
        yStdName = 'latitude';
        xLong    = 'longitude';
        yLong    = 'latitude';
        xUnits   = 'degrees_east';
        yUnits   = 'degrees_north';
        crsName  = 'WGS 84';
        epsgCode = 'EPSG:4326';
        gridMappingName = 'latitude_longitude';
        wktString = '';
    else
        xName    = 'x';
        yName    = 'y';
        xStdName = 'projection_x_coordinate';
        yStdName = 'projection_y_coordinate';
        xLong    = 'x coordinate of projection';
        yLong    = 'y coordinate of projection';
        xUnits   = 'm';
        yUnits   = 'm';
        crsName  = 'projected';
        epsgCode = '';
        gridMappingName = 'transverse_mercator'; % sensible default for UTM
        wktString = '';
    end

    % -----------------------------------------------------
    % Build coordinate vectors from raster reference
    % -----------------------------------------------------
    colCenters = 1:nCols;
    rowCenters = 1:nRows;

    [xCoords, ~] = intrinsicToWorld(R, colCenters, ones(1, nCols));
    [~, yCoords] = intrinsicToWorld(R, ones(1, nRows), rowCenters);

    xCoords = double(xCoords(:));
    yCoords = double(yCoords(:));

    % -----------------------------------------------------
    % Store y increasing for GIS/NetCDF compatibility
    % If original raster has descending y, append function
    % must flip rows before writing.
    % -----------------------------------------------------
    flipY = false;
    if numel(yCoords) >= 2 && yCoords(2) < yCoords(1)
        yCoords = flipud(yCoords);
        flipY = true;
    end

    % -----------------------------------------------------
    % Try to infer CRS metadata from MATLAB CRS object
    % -----------------------------------------------------
    if isGeographic
        try
            if isprop(R, 'GeographicCRS') && ~isempty(R.GeographicCRS)
                g = R.GeographicCRS;
                if isprop(g, 'Name') && ~isempty(g.Name)
                    crsName = char(g.Name);
                end
            end
        catch
        end
    else
        try
            if isprop(R, 'ProjectedCRS') && ~isempty(R.ProjectedCRS)
                p = R.ProjectedCRS;

                if isprop(p, 'Name') && ~isempty(p.Name)
                    crsName = char(p.Name);
                end

                % Try direct authority
                try
                    if isprop(p, 'Authority') && isprop(p, 'AuthorityCode') && ...
                            ~isempty(p.Authority) && ~isempty(p.AuthorityCode)
                        epsgCode = sprintf('%s:%s', ...
                            char(string(p.Authority)), char(string(p.AuthorityCode)));
                    end
                catch
                end

                % Try direct EPSG code
                if isempty(epsgCode)
                    try
                        if isprop(p, 'EPSGCode') && ~isempty(p.EPSGCode)
                            epsgCode = sprintf('EPSG:%d', p.EPSGCode);
                        end
                    catch
                    end
                end

                % Infer from name + geographic CRS if needed
                if isempty(epsgCode)
                    try
                        epsgCode = inferEPSGFromCRS(p);
                    catch
                    end
                end

                % WKT
                try
                    wktString = p.WKT;
                catch
                    try
                        wktString = evalc('disp(p)');
                    catch
                        wktString = '';
                    end
                end

                % Projection method -> CF grid_mapping_name
                try
                    if isprop(p, 'ProjectionMethod') && ~isempty(p.ProjectionMethod)
                        pm = lower(char(string(p.ProjectionMethod)));

                        if contains(pm, 'transverse mercator')
                            gridMappingName = 'transverse_mercator';
                        elseif contains(pm, 'lambert conformal conic')
                            gridMappingName = 'lambert_conformal_conic';
                        elseif contains(pm, 'mercator')
                            gridMappingName = 'mercator';
                        elseif contains(pm, 'albers')
                            gridMappingName = 'albers_conical_equal_area';
                        elseif contains(pm, 'polar stereographic')
                            gridMappingName = 'polar_stereographic';
                        elseif contains(pm, 'stereographic')
                            gridMappingName = 'stereographic';
                        end
                    end
                catch
                end
            end
        catch
        end
    end

    % -----------------------------------------------------
    % Create coordinate variables
    % -----------------------------------------------------
    nccreate(ncFile, xName, ...
        'Dimensions', {xName, numel(xCoords)}, ...
        'Datatype', 'double', ...
        'Format', 'netcdf4');

    nccreate(ncFile, yName, ...
        'Dimensions', {yName, numel(yCoords)}, ...
        'Datatype', 'double', ...
        'Format', 'netcdf4');

    nccreate(ncFile, 'time', ...
        'Dimensions', {'time', nTimes}, ...
        'Datatype', 'double', ...
        'Format', 'netcdf4');

    % CRS variable
    nccreate(ncFile, 'crs', ...
        'Datatype', 'int32', ...
        'Format', 'netcdf4');

    % Data variable
    nccreate(ncFile, varName, ...
        'Dimensions', {xName, numel(xCoords), yName, numel(yCoords), 'time', nTimes}, ...
        'Datatype', 'single', ...
        'FillValue', single(fillValue), ...
        'Format', 'netcdf4', ...
        'DeflateLevel', deflateLevel, ...
        'Shuffle', true);

    % -----------------------------------------------------
    % Write coordinate vectors
    % -----------------------------------------------------
    ncwrite(ncFile, xName, xCoords);
    ncwrite(ncFile, yName, yCoords);
    ncwrite(ncFile, 'time', timeValues);
    ncwrite(ncFile, 'crs', int32(0));

    % -----------------------------------------------------
    % Coordinate attributes
    % -----------------------------------------------------
    ncwriteatt(ncFile, xName, 'standard_name', xStdName);
    ncwriteatt(ncFile, xName, 'long_name', xLong);
    ncwriteatt(ncFile, xName, 'units', xUnits);
    ncwriteatt(ncFile, xName, 'axis', 'X');

    ncwriteatt(ncFile, yName, 'standard_name', yStdName);
    ncwriteatt(ncFile, yName, 'long_name', yLong);
    ncwriteatt(ncFile, yName, 'units', yUnits);
    ncwriteatt(ncFile, yName, 'axis', 'Y');

    ncwriteatt(ncFile, 'time', 'standard_name', 'time');
    ncwriteatt(ncFile, 'time', 'long_name', 'time');
    ncwriteatt(ncFile, 'time', 'units', timeUnits);
    ncwriteatt(ncFile, 'time', 'calendar', timeCalendar);
    ncwriteatt(ncFile, 'time', 'axis', 'T');

    % -----------------------------------------------------
    % CRS attributes
    % -----------------------------------------------------
    ncwriteatt(ncFile, 'crs', 'grid_mapping_name', gridMappingName);

    if ~isempty(crsName)
        ncwriteatt(ncFile, 'crs', 'long_name', crsName);
    end

    if ~isempty(epsgCode)
        ncwriteatt(ncFile, 'crs', 'epsg_code', epsgCode);
    end

    if ~isempty(wktString)
        try
            ncwriteatt(ncFile, 'crs', 'spatial_ref', wktString);
        catch
        end
        try
            ncwriteatt(ncFile, 'crs', 'crs_wkt', wktString);
        catch
        end
    end

    if ~isGeographic
        % Projection parameters if available
        try
            p = R.ProjectedCRS;

            if isprop(p, 'ProjectionParameters') && ~isempty(p.ProjectionParameters)
                pp = p.ProjectionParameters;

                tryWriteProjectionParam(ncFile, 'crs', pp, 'LatitudeOfNaturalOrigin',      'latitude_of_projection_origin');
                tryWriteProjectionParam(ncFile, 'crs', pp, 'LongitudeOfNaturalOrigin',     'longitude_of_central_meridian');
                tryWriteProjectionParam(ncFile, 'crs', pp, 'ScaleFactorAtNaturalOrigin',   'scale_factor_at_central_meridian');
                tryWriteProjectionParam(ncFile, 'crs', pp, 'FalseEasting',                 'false_easting');
                tryWriteProjectionParam(ncFile, 'crs', pp, 'FalseNorthing',                'false_northing');
                tryWriteProjectionParam(ncFile, 'crs', pp, 'StandardParallel',             'standard_parallel');
            end
        catch
        end

        % Geographic ellipsoid info
        try
            g = R.ProjectedCRS.GeographicCRS;
            if isprop(g, 'Spheroid') && ~isempty(g.Spheroid)
                sph = g.Spheroid;
                if isprop(sph, 'SemimajorAxis') && ~isempty(sph.SemimajorAxis)
                    ncwriteatt(ncFile, 'crs', 'semi_major_axis', double(sph.SemimajorAxis));
                end
                if isprop(sph, 'InverseFlattening') && ~isempty(sph.InverseFlattening)
                    ncwriteatt(ncFile, 'crs', 'inverse_flattening', double(sph.InverseFlattening));
                end
            end
            if isprop(g, 'PrimeMeridian') && ~isempty(g.PrimeMeridian)
                ncwriteatt(ncFile, 'crs', 'longitude_of_prime_meridian', double(g.PrimeMeridian));
            end
        catch
        end
    else
        % Basic WGS84 geographic info
        ncwriteatt(ncFile, 'crs', 'semi_major_axis', 6378137.0);
        ncwriteatt(ncFile, 'crs', 'inverse_flattening', 298.257223563);
        ncwriteatt(ncFile, 'crs', 'longitude_of_prime_meridian', 0.0);
    end

    % -----------------------------------------------------
    % Data variable attributes
    % -----------------------------------------------------
    ncwriteatt(ncFile, varName, 'long_name', longName);
    ncwriteatt(ncFile, varName, 'units', units);
    ncwriteatt(ncFile, varName, 'coordinates', sprintf('%s %s time', xName, yName));
    ncwriteatt(ncFile, varName, 'grid_mapping', 'crs');

    % -----------------------------------------------------
    % Global attributes
    % -----------------------------------------------------
    ncwriteatt(ncFile, '/', 'title', [longName ' time series']);
    ncwriteatt(ncFile, '/', 'source', 'HydroPol2D export');
    ncwriteatt(ncFile, '/', 'history', ['Created on ' datestr(now, 'yyyy-mm-dd HH:MM:ss')]);
    ncwriteatt(ncFile, '/', 'Conventions', 'CF-1.8');
    ncwriteatt(ncFile, '/', 'crs_name', crsName);

    if ~isempty(epsgCode)
        ncwriteatt(ncFile, '/', 'coordinate_reference_system', epsgCode);
    end

    if flipY
        ncwriteatt(ncFile, '/', 'note_on_y_axis', ...
            'Y coordinates were stored in increasing order; raster rows must be flipped before writing time slices.');
    end

    try
        ncwriteatt(ncFile, '/', 'GeoKeyDirectoryTag', evalc('disp(geoKeyTag)'));
    catch
        % do nothing
    end
end


function epsgCode = inferEPSGFromCRS(projectedCRS)
% inferEPSGFromCRS
% Try to infer EPSG code from MATLAB projcrs/geocrs information.

    epsgCode = '';

    % 1) Direct authority code
    try
        if isprop(projectedCRS, 'Authority') && isprop(projectedCRS, 'AuthorityCode') && ...
                ~isempty(projectedCRS.Authority) && ~isempty(projectedCRS.AuthorityCode)
            epsgCode = sprintf('%s:%s', ...
                char(string(projectedCRS.Authority)), ...
                char(string(projectedCRS.AuthorityCode)));
            return
        end
    catch
    end

    % 2) Direct EPSGCode property
    try
        if isprop(projectedCRS, 'EPSGCode') && ~isempty(projectedCRS.EPSGCode)
            epsgCode = sprintf('EPSG:%d', projectedCRS.EPSGCode);
            return
        end
    catch
    end

    % 3) Infer from name + geographic CRS + projection method
    try
        pName = lower(string(projectedCRS.Name));
    catch
        pName = "";
    end

    try
        gName = lower(string(projectedCRS.GeographicCRS.Name));
    catch
        gName = "";
    end

    try
        pMethod = lower(string(projectedCRS.ProjectionMethod));
    catch
        pMethod = "";
    end

    % Example: "WGS 84 / UTM zone 43N"
    try
        tok = regexp(char(pName), 'utm zone\s+(\d{1,2})([ns])', 'tokens', 'once');
        if ~isempty(tok) && contains(gName, 'wgs 84') && contains(pMethod, 'transverse mercator')
            zoneNum = str2double(tok{1});
            hemi = lower(tok{2});

            if hemi == 'n'
                epsgCode = sprintf('EPSG:%d', 32600 + zoneNum);
            else
                epsgCode = sprintf('EPSG:%d', 32700 + zoneNum);
            end
            return
        end
    catch
    end
end


function tryWriteProjectionParam(ncFile, crsVarName, projParams, matlabField, cfField)
% Safely write a projection parameter if it exists.

    try
        if isprop(projParams, matlabField)
            val = projParams.(matlabField);
            if ~isempty(val) && isnumeric(val)
                ncwriteatt(ncFile, crsVarName, cfField, double(val));
            end
        end
    catch
    end
end


function appendTimeSliceToNetCDF(ncFile, varName, Z, timeIndex)
% appendTimeSliceToNetCDF
% Writes one raster slice into:
%   varName(x, y, time)
%
% IMPORTANT:
%   Y coordinates were stored in increasing order,
%   so raster rows are flipped before transpose.

    Z = single(Z);
    Z = flipud(Z);
    ncwrite(ncFile, varName, Z.', [1 1 timeIndex]);
end

function fig = createStyledFigure(show_figures, paper_bg, figsize_inches)

if show_figures
    vis = 'on';
else
    vis = 'off';
end

fig = figure( ...
    'Color', paper_bg, ...
    'Visible', vis, ...
    'WindowStyle', 'normal', ...
    'Units', 'inches', ...
    'Position', [1 1 figsize_inches(1) figsize_inches(2)], ...
    'Resize', 'off', ...
    'NumberTitle', 'off', ...
    'IntegerHandle', 'off', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 figsize_inches(1) figsize_inches(2)], ...
    'PaperSize', [figsize_inches(1) figsize_inches(2)] );

drawnow;
set(fig,'Units','inches');
set(fig,'Position',[1 1 figsize_inches(1) figsize_inches(2)]);
drawnow;
end
