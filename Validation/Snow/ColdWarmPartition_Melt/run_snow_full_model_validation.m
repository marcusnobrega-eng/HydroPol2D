% P1-SNOW-001F: full HydroPol2D V-tilted snow-and-runoff integration.
%
% This script intentionally uses the normal preprocessing and main routing
% loop. It is kept separate from run_snow_model.m because preprocessing
% clears the script workspace as part of the regular HydroPol2D workflow.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
base_static_dir = fullfile(repo_root, 'Validation', 'Phase1_VTilted_Catchment', 'Static');
full_root = fullfile(case_dir, 'FullModelRuns', 'P1-SNOW-001F');
forcing_dir = fullfile(full_root, 'Forcing');
output_root = fullfile(full_root, 'Outputs');
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
config_dir = fullfile(case_dir, 'Config');

if ~exist(forcing_dir, 'dir'); mkdir(forcing_dir); end
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end
if exist(output_root, 'dir'); rmdir(output_root, 's'); end

addpath(functions_dir, '-begin');
hydropol2d_add_runtime_paths(hydropol2d_find_root(case_dir));
addpath(config_dir, '-begin');

etp_path = fullfile(forcing_dir, 'ETP_input_data.xlsx');
write_station_forcing(etp_path);

Paths = make_output_paths(output_root);
InputPaths = make_input_paths(base_static_dir, etp_path);
input_data_bypass_script_path = fullfile(config_dir, 'input_data_bypass_script.m');
use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
clean_output_folder = true;
run_postprocessing = false;
enable_logging = false;
model_folder = '';
resultsDir = Paths.Results;
export_root_dir = Paths.Root;

HydroPol2D_preprocessing;
% This validation requests the non-intrusive event ledger. It is off by
% default for regular HydroPol2D simulations.
running_control.flag_system_mass_ledger = true;
HydroPol2D_Main_While;

valid = ~idx_nan;
area_m2 = nansum(C_a(valid), 'all');
if ~exist('system_mass_ledger', 'var')
    error('P1-SNOW-001F did not return the requested event mass ledger.');
end
rain_volume_m3 = system_mass_ledger.cumulative_precipitation_m3;
snow_volume_m3 = nansum(Snow_Properties.SWE_t(valid) .* C_a(valid) / 1000, 'all');
surface_volume_m3 = nansum(depths.d_t(valid) .* C_a(valid) / 1000, 'all');
soil_volume_m3 = nansum(Soil_Properties.I_t(valid) .* C_a(valid) / 1000, 'all');
outlet_volume_m3 = system_mass_ledger.cumulative_outlet_m3;
initial_storage_m3 = system_mass_ledger.initial_storage_m3;
final_storage_m3 = system_mass_ledger.final_storage_m3;
event_net_flux_m3 = system_mass_ledger.cumulative_precipitation_m3 + ...
    system_mass_ledger.cumulative_boundary_inflow_m3 + ...
    system_mass_ledger.cumulative_prescribed_recharge_m3 - ...
    system_mass_ledger.cumulative_canopy_evaporation_m3 - ...
    system_mass_ledger.cumulative_etr_m3 - ...
    system_mass_ledger.cumulative_open_water_evaporation_m3 - ...
    system_mass_ledger.cumulative_snow_sublimation_m3 - ...
    system_mass_ledger.cumulative_outlet_m3;
system_mass_residual_m3 = (final_storage_m3 - initial_storage_m3) - event_net_flux_m3;
system_mass_residual_pct = 100 * abs(system_mass_residual_m3) / max(rain_volume_m3, eps);
snow_module_residual_m3 = abs(mass_balance_history.cum_errors_m3(end, 2));
finite_states = all(isfinite([Snow_Properties.SWE_t(valid); Snow_Properties.H_snow_t(valid); ...
    Snow_Properties.rho_snow(valid); depths.d_t(valid)]));

Summary = table("P1-SNOW-001F", "Full V-tilted snow and runoff", rain_volume_m3, ...
    snow_volume_m3, surface_volume_m3, soil_volume_m3, outlet_volume_m3, ...
    initial_storage_m3, final_storage_m3, event_net_flux_m3, system_mass_residual_m3, ...
    system_mass_residual_pct, snow_module_residual_m3, finite_states, ...
    'VariableNames', {'case_id','case_name','rain_volume_m3','final_swe_volume_m3', ...
    'final_surface_volume_m3','final_soil_storage_m3','outlet_volume_m3', ...
    'initial_storage_m3','final_storage_m3','event_net_flux_m3', ...
    'system_mass_residual_m3','system_mass_residual_pct','snow_module_residual_m3', ...
    'finite_states'});
Summary.passed = Summary.system_mass_residual_pct < 0.1 && ...
    Summary.snow_module_residual_m3 < 1e-6 && Summary.finite_states;
Summary.status = strings(height(Summary),1);
Summary.status(Summary.passed) = "pass";
Summary.status(~Summary.passed) = "fail";
writetable(Summary, fullfile(summary_dir, 'VTilted_Snow_FullModel_Summary.csv'));

temp_ts_path = fullfile(Paths.Temp, 'current_time_series_outputs.mat');
Temp = load(temp_ts_path, 'TempTS');
Hydrograph = table(Temp.TempTS.time.time_hydrograph_min(:), ...
    Temp.TempTS.outlet.hydrograph_m3_s(:), ...
    'VariableNames', {'time_min','outlet_discharge_m3s'});
writetable(Hydrograph, fullfile(summary_dir, 'VTilted_Snow_FullModel_Hydrograph.csv'));
Ledger = struct2table(system_mass_ledger);
writetable(Ledger, fullfile(summary_dir, 'VTilted_Snow_FullModel_MassLedger.csv'));
write_full_model_timeseries(Maps, running_control, summary_dir);
append_standard_outputs(Summary, summary_dir);
disp(Summary);

function InputPaths = make_input_paths(static_dir, etp_path)
InputPaths = struct();
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_dir, 'DTB.tif');
InputPaths.GW_table_path = fullfile(static_dir, 'GW_table.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_dir, 'Initial_SM.tif');
InputPaths.LAI_path = fullfile(static_dir, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_dir, 'Albedo.tif');
InputPaths.Initial_SWE_path = '';
InputPaths.Initial_Snow_Depth_path = '';
InputPaths.ETP_input_spreadsheet = etp_path;
InputPaths.Subgrid_DEM_path = '';
InputPaths.RiverWidths_path = '';
InputPaths.RiverDepths_path = '';
InputPaths.Warmup_Depth_path = '';
InputPaths.Initial_Buildup_path = '';
InputPaths.B1_path = '';
InputPaths.B2_path = '';
InputPaths.W1_path = '';
InputPaths.W2_path = '';
InputPaths.Rainfall_Rasters_Folder = '';
InputPaths.Transpiration_Rasters_Folder = '';
InputPaths.Evaporation_Rasters_Folder = '';
InputPaths.Inflow_Hydrograph_CSV = '';
InputPaths.Stage_Hydrograph_CSV = '';
InputPaths.Observed_Gauges_CSV = '';
InputPaths.Rainfall_Timeseries_File = '';
InputPaths.Outlet_Cells_CSV = '';
end

function Paths = make_output_paths(root_dir)
Paths = struct();
Paths.Root = root_dir;
Paths.Results = fullfile(root_dir, 'Modeling_Results');
Paths.Temp = fullfile(root_dir, 'Temporary_Files');
Paths.Logs = fullfile(root_dir, 'Logs');
Paths.FigPDF = fullfile(Paths.Results, 'Figures_PDF');
Paths.FigFIG = fullfile(Paths.Results, 'Figures_FIG');
Paths.Tables = fullfile(Paths.Results, 'Tables_CSV');
Paths.RastersWD = fullfile(Paths.Results, 'Rasters_Water_Depths');
Paths.RastersWSE = fullfile(Paths.Results, 'Rasters_WSE');
Paths.RastersStatic = fullfile(Paths.Results, 'Rasters_Static');
Paths.RastersVelocity = fullfile(Paths.Results, 'Rasters_Velocity');
Paths.RastersHazard = fullfile(Paths.Results, 'Rasters_Hazard');
Paths.WQMaps = fullfile(Paths.Results, 'Rasters_WQ');
Paths.HRMaps = fullfile(Paths.Results, 'Rasters_Human_Risk');
Paths.Anim = fullfile(Paths.Results, 'GIFs_MP4');
Paths.Shapes = fullfile(Paths.Results, 'Shapefiles');
fields = fieldnames(Paths);
for i = 1:numel(fields)
    if ~exist(Paths.(fields{i}), 'dir')
        mkdir(Paths.(fields{i}));
    end
end
end

function write_station_forcing(path)
times = datetime(2025, 1, 30, 0, 0, 0) + minutes([0; 60; 120; 135]);
data = [-4 -8 -6 2 70 0; -2 -6 -4 2 70 0; 5 1 3 2 70 0; 5 1 3 2 70 0];
cells = cell(7, 8);
cells(1,:) = {'Index','Date','Tmax_C','Tmin_C','Tavg_C','U2_m_s','RH_pct','G_MJ_m2_day'};
cells(2,:) = {NaN,NaT,NaN,NaN,NaN,0,NaN,0};
cells(3,:) = {NaN,NaT,NaN,NaN,NaN,NaN,NaN,NaN};
for i = 1:numel(times)
    cells{3+i,1} = i;
    cells{3+i,2} = times(i);
    cells(3+i,3:8) = num2cell(data(i,:));
end
writecell(cells, path);
end

function write_full_model_timeseries(Maps, running_control, summary_dir)
if ~isfield(Maps.Hydro, 'Snowpack') || ~isfield(Maps.Hydro, 'd')
    error('P1-SNOW-001F did not retain snowpack and surface-water map series.');
end
n = min([size(Maps.Hydro.Snowpack, 3), size(Maps.Hydro.d, 3), ...
    numel(running_control.time_records)]);
t = running_control.time_records(1:n)';
mean_snow_depth_mm = squeeze(mean(Maps.Hydro.Snowpack(:,:,1:n), [1 2], 'omitnan'));
mean_surface_water_mm = squeeze(mean(Maps.Hydro.d(:,:,1:n), [1 2], 'omitnan'));
T = table(t, mean_snow_depth_mm(:), mean_surface_water_mm(:), ...
    'VariableNames', {'time_min','mean_snow_depth_mm','mean_surface_water_mm'});
writetable(T, fullfile(summary_dir, 'VTilted_Snow_FullModel_TimeSeries.csv'));

blue = [11 129 162] / 255;
red = [226 87 89] / 255;
fig = figure('Color', 'w', 'Position', [100 100 820 390], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');
plot(ax, T.time_min / 60, T.mean_snow_depth_mm, 'Color', blue, 'LineWidth', 2.0, ...
    'DisplayName', 'Snow depth');
plot(ax, T.time_min / 60, T.mean_surface_water_mm, 'Color', red, 'LineWidth', 2.0, ...
    'DisplayName', 'Surface water depth');
xlabel(ax, 'Time [h]');
ylabel(ax, 'Domain-mean depth [mm]');
legend(ax, 'Location', 'northwest', 'Box', 'off');
grid(ax, 'on'); box(ax, 'on');
set(ax, 'FontName', 'Helvetica', 'FontSize', 11, 'LineWidth', 1.6);
exportgraphics(fig, fullfile(summary_dir, 'snow_full_model_vtilted.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(summary_dir, 'snow_full_model_vtilted.pdf'), 'ContentType', 'vector');
close(fig);
end

function append_standard_outputs(Summary, summary_dir)
metric_path = fullfile(summary_dir, 'Metric_Summary.csv');
pass_path = fullfile(summary_dir, 'Pass_Fail.csv');
metric_row = table("P1-SNOW-001F", "Full V-tilted snow and runoff", ...
    Summary.system_mass_residual_pct, Summary.system_mass_residual_m3, NaN, 0, ...
    Summary.passed, "Normal HydroPol2D preprocessing and routing; max_abs_error is the system mass residual [%].", ...
    'VariableNames', {'case_id','case_name','max_abs_error','mass_residual_m3', ...
    'timestep_spread_mm','fallback_cells','passed','notes'});
pass_row = table("P1-SNOW-001F", Summary.status, Summary.passed, Summary.passed, ...
    "Full V-tilted snow-and-runoff integration with internal meteorological forcing.", ...
    'VariableNames', {'case_id','status','passed','report_ready','notes'});

if isfile(metric_path)
    T = readtable(metric_path, 'TextType', 'string');
    T(T.case_id == metric_row.case_id, :) = [];
    writetable([T; metric_row], metric_path);
else
    writetable(metric_row, metric_path);
end
if isfile(pass_path)
    T = readtable(pass_path, 'TextType', 'string');
    T(T.case_id == pass_row.case_id, :) = [];
    writetable([T; pass_row], pass_path);
else
    writetable(pass_row, pass_path);
end
end
