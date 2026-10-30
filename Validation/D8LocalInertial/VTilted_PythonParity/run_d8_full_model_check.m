% VAL-D8LI-002: full HydroPol2D runs with D8 local-inertial routing.
%
% HydroPol2D_preprocessing + HydroPol2D_Main_While run the V-tilted
% catchment of VAL-D8LI-001 in bypass mode (Config/input_data_bypass_script.m),
% once as a D8 river graph in a D4 floodplain and once on the full-domain D8
% graph: 47 mm/h for 10 min, 30 min in all, n = 0.04. The event ledger
% (rain, outlet, storage) must close, the D8 maps must be written, and the
% outlet hydrograph must follow the Python engine (run_case on the same DEM,
% make_python_reference.py part 2): NSE >= 0.99, peak within 5 %, outflow
% volume within 2 %. The two codes take different time steps, so they agree
% to the time-stepping error, not to round-off.
%
% Results: Outputs/Validation/FullModel/Pass_Fail.csv and the hydrographs.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(case_dir, '..', '..', '..', 'HydroPol2D_Functions'), '-begin');
hydropol2d_add_runtime_paths(hydropol2d_find_root(case_dir));
summary_dir = fullfile(case_dir, 'Outputs', 'Validation', 'FullModel');
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end

RunSummary = table();
Cases = {'hybrid', 'full'};
for icase = 1:numel(Cases)
    % Start each run from a clean workspace (preprocessing reuses flags).
    clearvars -except RunSummary Cases icase case_dir summary_dir
    ValidationCase = struct('mode', Cases{icase});
    full_root = fullfile(case_dir, 'Outputs', 'FullModelRuns', Cases{icase});
    static_dir = fullfile(case_dir, 'Outputs', 'Static');
    output_root = fullfile(full_root, 'Run');
    config_dir = fullfile(case_dir, 'Config');
    if exist(output_root, 'dir'); rmdir(output_root, 's'); end
    d8_write_static_rasters(static_dir);
    Paths = d8_output_paths(output_root);
    InputPaths = d8_input_paths(static_dir);
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
    running_control.flag_system_mass_ledger = true;
    HydroPol2D_Main_While;

    SummaryRow = d8_model_row(case_dir, summary_dir, ValidationCase.mode, system_mass_ledger, ...
        fullfile(Paths.Temp, 'current_time_series_outputs.mat'), D8_States, Maps, depths, idx_nan);
    RunSummary = [RunSummary; SummaryRow]; %#ok<AGROW>
end

writetable(RunSummary, fullfile(summary_dir, 'Pass_Fail.csv'));
disp(RunSummary(:, {'case_id', 'status', 'report_ready', 'metric', 'value', 'threshold'}));
failed = unique(RunSummary.case_id(RunSummary.status ~= "PASS"));
if ~isempty(failed)
    error('VAL-D8LI-002 failed: %s.', strjoin(failed, ', '));
end
fprintf('VAL-D8LI-002 passed.\n');

%% ------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function Row = d8_model_row(case_dir, summary_dir, mode, ledger, series_path, D8_States, Maps, depths, idx_nan)
rain = ledger.cumulative_precipitation_m3;
net = rain + ledger.cumulative_boundary_inflow_m3 - ledger.cumulative_outlet_m3;
residual = abs((ledger.final_storage_m3 - ledger.initial_storage_m3) - net) / rain;
Temp = load(series_path, 'TempTS');
t_min = double(Temp.TempTS.time.time_hydrograph_min(:));
q = double(Temp.TempTS.outlet.hydrograph_m3_s(:));
keep = [true; diff(t_min) > 0];
E = readtable(fullfile(case_dir, "reference_engine_" + mode + ".csv"));
q_matlab = interp1(t_min(keep) * 60, q(keep), E.time_s, 'linear', 'extrap');
q_python = E.outlet_m3_s;
nse = 1 - sum((q_matlab - q_python).^2) / sum((q_python - mean(q_python)).^2);
peak_error = abs(max(q_matlab) - max(q_python)) / max(q_python);
volume_error = abs(ledger.cumulative_outlet_m3 - E.cumulative_outflow_m3(end)) / E.cumulative_outflow_m3(end);
maps_ok = isfield(Maps.Hydro, 'd8_channel_depth') && all(isfinite(depths.d_t(~idx_nan)));
writetable(table(E.time_s, q_python, q_matlab, 'VariableNames', ...
    {'time_s', 'python_outlet_m3_s', 'matlab_outlet_m3_s'}), ...
    fullfile(summary_dir, "Hydrograph_" + mode + ".csv"));
Row = result_rows("VAL-D8LI-002-" + upper(mode), ...
    "Full HydroPol2D run, D8 " + mode + ", vs Python run_case", ...
    ["ledger_residual_rel", "outlet_nse", "outlet_peak_rel_error", "outflow_volume_rel_error", ...
     "maps_written", "diagonal_transfer_m3"], ...
    [residual, nse, peak_error, volume_error, double(maps_ok), ...
     D8_States.cumulative_diagonal_transfer_m3], ...
    [1e-9, 0.99, 0.05, 0.02, 1, 0], ["le", "ge", "le", "le", "eq", "gt"]);
end

function d8_write_static_rasters(static_dir)
if ~exist(static_dir, 'dir'); mkdir(static_dir); end
[z, ~] = vtilted_dem(0.04);
[ny, nx] = size(z);
R = maprefcells([0, nx * 20], [0, ny * 20], [ny, nx], 'ColumnsStartFrom', 'north');
geotiffwrite(fullfile(static_dir, 'DEM.tif'), z, R, 'CoordRefSysCode', 3857);
geotiffwrite(fullfile(static_dir, 'LULC.tif'), ones(ny, nx), R, 'CoordRefSysCode', 3857);
geotiffwrite(fullfile(static_dir, 'SOIL.tif'), ones(ny, nx), R, 'CoordRefSysCode', 3857);
end

function InputPaths = d8_input_paths(static_dir)
InputPaths = struct();
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
names = {'DTB_path', 'GW_table_path', 'Initial_Soil_Moisture_path', 'LAI_path', 'Albedo_path', ...
    'Initial_SWE_path', 'Initial_Snow_Depth_path', 'ETP_input_spreadsheet', 'Subgrid_DEM_path', ...
    'RiverWidths_path', 'RiverDepths_path', 'Warmup_Depth_path', 'Initial_Buildup_path', ...
    'B1_path', 'B2_path', 'W1_path', 'W2_path', 'Rainfall_Rasters_Folder', ...
    'Transpiration_Rasters_Folder', 'Evaporation_Rasters_Folder', 'Inflow_Hydrograph_CSV', ...
    'Stage_Hydrograph_CSV', 'Observed_Gauges_CSV', 'Rainfall_Timeseries_File', 'Outlet_Cells_CSV'};
for k = 1:numel(names)
    InputPaths.(names{k}) = '';
end
end

function Paths = d8_output_paths(root_dir)
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

function [z, n_land] = vtilted_dem(manning_n)
[col, row] = meshgrid(0:10, 0:8);
z = 0.004 * abs(col - 5) * 20 + 0.004 * (8 - row) * 20;
n_land = manning_n * ones(size(z));
end

function Rows = result_rows(case_id, case_name, metrics, values, thresholds, rules)
passed = true(numel(metrics), 1);
for k = 1:numel(metrics)
    switch rules(k)
        case "eq"
            passed(k) = values(k) == thresholds(k);
        case "le"
            passed(k) = values(k) <= thresholds(k);
        case "gt"
            passed(k) = values(k) > thresholds(k);
        case "ge"
            passed(k) = values(k) >= thresholds(k);
        otherwise
            passed(k) = isfinite(values(k));
    end
end
status = "PASS";
if ~all(passed); status = "FAIL"; end
Rows = table(repmat(case_id, numel(metrics), 1), repmat(case_name, numel(metrics), 1), ...
    repmat(status, numel(metrics), 1), true(numel(metrics), 1), metrics(:), values(:), ...
    thresholds(:), rules(:), passed, ...
    'VariableNames', {'case_id', 'case_name', 'status', 'report_ready', 'metric', ...
    'value', 'threshold', 'rule', 'metric_passed'});
end
