% VAL-UD-001F: urban drainage inside a full HydroPol2D grid run (bypass mode).
%
% Normal preprocessing and main loop on Grid/DEM.tif with network_grid.gpkg
% (the GeoPackage HydroPol2D-Python accepts on the same DEM). Rain 40 mm/h
% for 15 min, 30 min run, local inertial, no infiltration. The event ledger
% must close with the drainage storage counted and the outfall water booked
% as outlet flow (free outfall) or boundary inflow (fixed-stage outfall).
%
% A driver of its own, because preprocessing clears the workspace. It
% writes its row to Outputs/Validation/Pass_Fail.csv and keeps the rows of
% run_urban_drainage_check.m.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
addpath(fullfile(repo_root, 'HydroPol2D_Functions'), '-begin');
hydropol2d_add_runtime_paths(hydropol2d_find_root(case_dir));
config_dir = fullfile(case_dir, 'Config');
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
output_root = fullfile(case_dir, 'Outputs', 'GridRun');
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end
if exist(output_root, 'dir'); rmdir(output_root, 's'); end

Paths = struct('Root', output_root);
Paths.Results = fullfile(output_root, 'Modeling_Results');
Paths.Temp = fullfile(output_root, 'Temporary_Files');
Paths.Logs = fullfile(output_root, 'Logs');
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
path_fields = fieldnames(Paths);
for i = 1:numel(path_fields)
    if ~exist(Paths.(path_fields{i}), 'dir'); mkdir(Paths.(path_fields{i})); end
end

grid_dir = fullfile(case_dir, 'Grid');
InputPaths = struct();
InputPaths.DEM_path = fullfile(grid_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(grid_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(grid_dir, 'SOIL.tif');
empty_paths = {'DTB_path', 'GW_table_path', 'Initial_Soil_Moisture_path', 'LAI_path', ...
    'Albedo_path', 'Initial_SWE_path', 'Initial_Snow_Depth_path', 'ETP_input_spreadsheet', ...
    'Subgrid_DEM_path', 'RiverWidths_path', 'RiverDepths_path', 'Warmup_Depth_path', ...
    'Initial_Buildup_path', 'B1_path', 'B2_path', 'W1_path', 'W2_path', ...
    'Rainfall_Rasters_Folder', 'Transpiration_Rasters_Folder', 'Evaporation_Rasters_Folder', ...
    'Inflow_Hydrograph_CSV', 'Stage_Hydrograph_CSV', 'Observed_Gauges_CSV', ...
    'Rainfall_Timeseries_File', 'Outlet_Cells_CSV'};
for i = 1:numel(empty_paths)
    InputPaths.(empty_paths{i}) = '';
end
input_data_bypass_script_path = fullfile(config_dir, 'input_data_bypass_script.m');
use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
clean_output_folder = true;
run_postprocessing = false;
enable_logging = false;
model_folder = '';
resultsDir = Paths.Results;
export_root_dir = Paths.Root;
drainage_flag_for_run = 1;
setenv('HYDROPOL2D_SKIP_INITIAL_PLOTS', '1');   % no input-map figures

HydroPol2D_preprocessing;
running_control.flag_system_mass_ledger = true;
running_control.report_every_percent = 25;
HydroPol2D_Main_While;
setenv('HYDROPOL2D_SKIP_INITIAL_PLOTS', '');

% ---------------------------------------------------------------------
% Checks
% ---------------------------------------------------------------------
L = system_mass_ledger;
S = Drainage_States.State;
rain_m3 = L.cumulative_precipitation_m3;
net_flux_m3 = L.cumulative_precipitation_m3 + L.cumulative_boundary_inflow_m3 + ...
    L.cumulative_prescribed_recharge_m3 - L.cumulative_canopy_evaporation_m3 - ...
    L.cumulative_etr_m3 - L.cumulative_open_water_evaporation_m3 - ...
    L.cumulative_snow_sublimation_m3 - L.cumulative_outlet_m3;
residual_m3 = (L.final_storage_m3 - L.initial_storage_m3) - net_flux_m3;
drainage_storage_m3 = sum(S.node_volume_m3);
capture_m3 = sum(S.inlet_cumulative_capture_m3);
surcharge_m3 = sum(S.inlet_cumulative_surcharge_m3);
outfall_out_m3 = S.cumulative_outfall_outflow_m3;
outfall_in_m3 = S.cumulative_outfall_inflow_m3;
% Ledger outcomes if one booking were dropped (what the check would see).
residual_without_storage_m3 = residual_m3 - (drainage_storage_m3 - S.initial_storage_m3);
residual_without_outfall_booking_m3 = residual_m3 - outfall_out_m3;
residual_without_inflow_booking_m3 = residual_m3 + outfall_in_m3;
scale_m3 = rain_m3 + outfall_in_m3;
tolerance_m3 = 1e-9 * scale_m3;
csv_files = dir(fullfile(Paths.Tables, 'Urban_Drainage_*.csv'));
records = numel(Drainage_States.History.time_min);
passed = abs(residual_m3) < tolerance_m3 && ...
    abs(L.cumulative_drainage_outfall_outflow_m3 - outfall_out_m3) < 1e-9 * scale_m3 && ...
    abs(L.cumulative_drainage_outfall_inflow_m3 - outfall_in_m3) < 1e-9 * scale_m3 && ...
    capture_m3 > 0.1 * rain_m3 && outfall_out_m3 > 0.1 * rain_m3 && ...
    outfall_in_m3 > 1 && surcharge_m3 > 1 && drainage_storage_m3 > 0 && ...
    L.cumulative_outlet_m3 >= outfall_out_m3 && numel(csv_files) == 7 && ...
    records >= floor(running_control.routing_time / running_control.record_time_hydrographs) && ...
    min(abs([residual_without_storage_m3, residual_without_outfall_booking_m3, ...
    residual_without_inflow_booking_m3])) > 1e3 * tolerance_m3;

Summary = table("VAL-UD-001F", rain_m3, capture_m3, surcharge_m3, outfall_out_m3, outfall_in_m3, ...
    drainage_storage_m3, L.cumulative_outlet_m3, L.cumulative_boundary_inflow_m3, residual_m3, ...
    residual_without_storage_m3, residual_without_outfall_booking_m3, ...
    residual_without_inflow_booking_m3, S.total_substeps, S.surface_steps, S.max_substeps, ...
    records, numel(csv_files), passed, 'VariableNames', {'case_id', 'rain_m3', 'capture_m3', ...
    'surcharge_m3', 'outfall_outflow_m3', 'outfall_inflow_m3', 'final_drainage_storage_m3', ...
    'ledger_outlet_m3', 'ledger_boundary_inflow_m3', 'ledger_residual_m3', ...
    'residual_without_drainage_storage_m3', 'residual_without_outfall_booking_m3', ...
    'residual_without_inflow_booking_m3', 'drainage_substeps', 'surface_steps', ...
    'max_substeps_per_step', 'drainage_records', 'drainage_csv_files', 'passed'});
writetable(Summary, fullfile(summary_dir, 'UrbanDrainage_GridRun_Summary.csv'));
writetable(struct2table(L), fullfile(summary_dir, 'UrbanDrainage_GridRun_MassLedger.csv'));
disp(Summary);

status = "fail";
if passed; status = "pass"; end
row = table("VAL-UD-001F", status, double(passed), residual_m3, outfall_out_m3, outfall_in_m3, ...
    "Full grid run in bypass mode: event ledger closes with drainage storage and outfall flows.", ...
    'VariableNames', {'case_id', 'status', 'report_ready', 'max_error', 'outfall_outflow_m3', ...
    'outfall_inflow_m3', 'notes'});
pass_path = fullfile(summary_dir, 'Pass_Fail.csv');
if isfile(pass_path)
    T = readtable(pass_path, 'TextType', 'string', 'Delimiter', ',');
    T(T.case_id == row.case_id, :) = [];
    T = [T; row];
else
    T = row;
end
writetable(T, pass_path);
disp(T);
