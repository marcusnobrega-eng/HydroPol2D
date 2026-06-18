clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
model_root = fullfile(repo_root, 'HydroPol2D_Model');
functions_dir = fullfile(model_root, 'HydroPol2D_Functions');
topo_path = fullfile(model_root, 'topotoolbox-master');
base_static_dir = fullfile(model_root, 'Validation', 'Phase1_VTilted_Catchment', 'Static');
config_dir = fullfile(case_dir, 'Config');

addpath(genpath(functions_dir));
addpath(genpath(topo_path));
addpath(config_dir);

registry_path = fullfile(config_dir, 'Infiltration_Case_Registry.csv');
Cases = readtable(registry_path, 'TextType', 'string');
case_filter = string(getenv('HYDROPOL2D_VALIDATION_CASE'));
if strlength(case_filter) > 0
    Cases = Cases(Cases.case_id == case_filter | Cases.case_name == case_filter, :);
    if isempty(Cases)
        error('No infiltration validation case matched HYDROPOL2D_VALIDATION_CASE=%s.', case_filter);
    end
end

full_root = fullfile(case_dir, 'FullModelRuns');
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end

RunSummary = table();

for icase = 1:height(Cases)
    ValidationCase = Cases(icase,:);
    case_id = char(ValidationCase.case_id);
    fprintf('\n============================================================\n');
    fprintf('Running full HydroPol2D infiltration case %s\n', case_id);
    fprintf('============================================================\n');

    case_root = fullfile(full_root, case_id);
    static_dir = fullfile(case_root, 'Static');
    output_root = fullfile(case_root, 'Outputs');
    prepare_case_static_rasters(base_static_dir, static_dir, ValidationCase);

    try
        Paths = make_validation_paths(output_root, true);
        InputPaths = make_input_paths(static_dir, topo_path, functions_dir, case_root);
        input_data_bypass_script_path = fullfile(config_dir, 'input_data_bypass_script.m');
        use_inputpaths_bypass = 1;
        use_inputdata_bypass = 1;
        clean_output_folder = true;
        run_postprocessing = false;
        enable_logging = false;
        model_folder = '';
        GD = [];
        resultsDir = Paths.Results;
        export_root_dir = Paths.Root;

        HydroPol2D_preprocessing;
        HydroPol2D_Main_While;

        SummaryRow = summarize_full_model_case(ValidationCase, Paths, ...
            C_a, idx_nan, depths, Soil_Properties, BC_States, Hydro_States, ...
            cumulative_infiltration, cumulative_recharge, outlet_runoff_volume, errors);
        SummaryRow.status = "completed";
        SummaryRow.error_message = "";
        if SummaryRow.max_final_depth_mm > 5000 || ~isfinite(SummaryRow.max_final_depth_mm)
            SummaryRow.status = "unstable";
            SummaryRow.error_message = "Hydraulic sanity check failed: max final depth exceeded 5 m.";
        end
    catch ME
        warning('Full model case %s failed: %s', case_id, ME.message);
        SummaryRow = failed_summary_row(ValidationCase, output_root, ME);
    end

    RunSummary = [RunSummary; SummaryRow]; %#ok<AGROW>
end

RunSummary = evaluate_full_model_summary(RunSummary);
writetable(RunSummary, fullfile(summary_dir, 'VTilted_Infiltration_FullModel_Summary.csv'));
disp(RunSummary);

function InputPaths = make_input_paths(static_dir, topo_path, functions_dir, case_root)
InputPaths = struct();
InputPaths.case_root = case_root;
InputPaths.topo_path = topo_path;
InputPaths.hydropol2d_tools = functions_dir;
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_dir, 'DTB.tif');
InputPaths.GW_table_path = fullfile(static_dir, 'GW_table.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_dir, 'Initial_SM.tif');
InputPaths.LAI_path = fullfile(static_dir, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_dir, 'Albedo.tif');
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
InputPaths.ETP_input_spreadsheet = '';
InputPaths.Rainfall_Timeseries_File = '';
InputPaths.Outlet_Cells_CSV = '';
end

function Paths = make_validation_paths(root_dir, clean_output)
if exist(root_dir, 'dir') && clean_output
    delete_dir_contents(root_dir);
end
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

function prepare_case_static_rasters(base_static_dir, static_dir, Cfg)
if ~exist(static_dir, 'dir'); mkdir(static_dir); end

[DEM, R] = readgeoraster(fullfile(base_static_dir, 'DEM.tif'));
DEM = double(DEM);
[LULC, ~] = readgeoraster(fullfile(base_static_dir, 'LULC.tif'));
[SOIL_base, ~] = readgeoraster(fullfile(base_static_dir, 'SOIL.tif'));
[LAI, ~] = readgeoraster(fullfile(base_static_dir, 'LAI.tif'));
[Albedo, ~] = readgeoraster(fullfile(base_static_dir, 'Albedo.tif'));

valid = isfinite(DEM);

if Cfg.regime == "spatial_soil_contrast"
    SOIL = double(SOIL_base);
else
    SOIL = ones(size(DEM));
    SOIL(~valid) = nan;
end

DTB = Cfg.soil_depth_m * ones(size(DEM));
DTB(~valid) = nan;

GW = DEM - Cfg.groundwater_depth_m;
GW(~valid) = nan;

I0 = zeros(size(DEM));
I0(~valid) = nan;

write_tif(fullfile(static_dir, 'DEM.tif'), DEM, R);
write_tif(fullfile(static_dir, 'LULC.tif'), double(LULC), R);
write_tif(fullfile(static_dir, 'SOIL.tif'), SOIL, R);
write_tif(fullfile(static_dir, 'DTB.tif'), DTB, R);
write_tif(fullfile(static_dir, 'GW_table.tif'), GW, R);
write_tif(fullfile(static_dir, 'Initial_SM.tif'), I0, R);
write_tif(fullfile(static_dir, 'LAI.tif'), double(LAI), R);
write_tif(fullfile(static_dir, 'Albedo.tif'), double(Albedo), R);
end

function write_tif(path, A, R)
A = single(A);
try
    geotiffwrite(path, A, R, 'CoordRefSysCode', 3857);
catch
    geotiffwrite(path, A, R);
end
end

function SummaryRow = summarize_full_model_case(Cfg, Paths, C_a, idx_nan, ...
    depths, Soil_Properties, BC_States, Hydro_States, cumulative_infiltration, ...
    cumulative_recharge, outlet_runoff_volume, errors)

area_m2 = nansum(C_a(~idx_nan), 'all');
rain_cum_mm = Cfg.rainfall_mm_h * Cfg.duration_h;
rain_volume_m3 = rain_cum_mm / 1000 * area_m2;

if ~Cfg.infiltration_enabled
    infiltration_volume_m3 = 0;
elseif exist('cumulative_infiltration','var') && ~isempty(cumulative_infiltration)
    infiltration_volume_m3 = nansum(C_a .* cumulative_infiltration, 'all');
else
    infiltration_volume_m3 = 0;
end

if exist('cumulative_recharge','var') && ~isempty(cumulative_recharge)
    recharge_volume_m3 = nansum(C_a .* cumulative_recharge / 1000, 'all');
else
    recharge_volume_m3 = 0;
end

surface_storage_m3 = nansum(C_a .* depths.d_t / 1000, 'all');
vadose_storage_m3 = nansum(C_a .* Soil_Properties.I_t / 1000, 'all');
outlet_volume_m3 = outlet_runoff_volume;

if isfield(Soil_Properties, 'Layers')
    mean_near_mm = mean(Soil_Properties.Layers.near_surface_storage_mm(:), 'omitnan');
    mean_root_mm = mean(Soil_Properties.Layers.root_zone_storage_mm(:), 'omitnan');
    mean_trans_mm = mean(Soil_Properties.Layers.transmission_storage_mm(:), 'omitnan');
else
    mean_near_mm = nan;
    mean_root_mm = nan;
    mean_trans_mm = nan;
end

SummaryRow = table( ...
    Cfg.case_id, Cfg.case_name, Cfg.regime, string(Paths.Root), ...
    rain_cum_mm, rain_volume_m3, infiltration_volume_m3, recharge_volume_m3, ...
    outlet_volume_m3, surface_storage_m3, vadose_storage_m3, ...
    mean(Hydro_States.f(:), 'omitnan'), max(Hydro_States.f(:), [], 'omitnan'), ...
    mean(depths.d_t(:), 'omitnan'), max(depths.d_t(:), [], 'omitnan'), ...
    mean(BC_States.h_t(:), 'omitnan'), mean(Soil_Properties.I_t(:), 'omitnan'), ...
    mean_near_mm, mean_root_mm, mean_trans_mm, ...
    sum(errors, 'omitnan'), ...
    'VariableNames', {'case_id','case_name','regime','output_root', ...
    'rain_cum_mm','rain_volume_m3','infiltration_volume_m3', ...
    'recharge_volume_m3','outlet_volume_m3','surface_storage_m3', ...
    'vadose_storage_m3','mean_final_infiltration_mm_h', ...
    'max_final_infiltration_mm_h','mean_final_depth_mm','max_final_depth_mm', ...
    'mean_groundwater_head_m','mean_final_soil_storage_mm', ...
    'mean_near_storage_mm','mean_root_storage_mm','mean_trans_storage_mm', ...
    'last_step_module_error_m3'});
end

function SummaryRow = failed_summary_row(Cfg, output_root, ME)
SummaryRow = table( ...
    Cfg.case_id, Cfg.case_name, Cfg.regime, string(output_root), ...
    nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, ...
    nan, nan, nan, nan, ...
    'VariableNames', {'case_id','case_name','regime','output_root', ...
    'rain_cum_mm','rain_volume_m3','infiltration_volume_m3', ...
    'recharge_volume_m3','outlet_volume_m3','surface_storage_m3', ...
    'vadose_storage_m3','mean_final_infiltration_mm_h', ...
    'max_final_infiltration_mm_h','mean_final_depth_mm','max_final_depth_mm', ...
    'mean_groundwater_head_m','mean_final_soil_storage_mm', ...
    'mean_near_storage_mm','mean_root_storage_mm','mean_trans_storage_mm', ...
    'last_step_module_error_m3'});
SummaryRow.status = "failed";
SummaryRow.error_message = string(ME.message);
end

function T = evaluate_full_model_summary(T)
if ~ismember('status', T.Properties.VariableNames)
    T.status = repmat("completed", height(T), 1);
end
if ~ismember('error_message', T.Properties.VariableNames)
    T.error_message = repmat("", height(T), 1);
end

T.infiltration_ratio = T.infiltration_volume_m3 ./ max(T.rain_volume_m3, eps);
T.recharge_ratio = T.recharge_volume_m3 ./ max(T.rain_volume_m3, eps);
T.validation_scope = repmat("coupled_full_model_diagnostic", height(T), 1);
T.report_ready = T.status == "completed" & isfinite(T.infiltration_ratio) & ...
    isfinite(T.max_final_depth_mm) & T.max_final_depth_mm <= 5000;

for i = 1:height(T)
    if T.status(i) ~= "completed"
        T.report_ready(i) = false;
        continue
    end
    switch T.regime(i)
        case "no_infiltration"
            T.report_ready(i) = T.infiltration_volume_m3(i) < 1e-6;
        case "supply_limited"
            T.report_ready(i) = T.infiltration_ratio(i) > 0.90;
        case "capacity_limited"
            T.report_ready(i) = T.infiltration_ratio(i) < 0.80 & T.mean_final_depth_mm(i) > 0;
        case "storage_limited"
            T.report_ready(i) = T.mean_final_soil_storage_mm(i) > 0 & T.mean_final_depth_mm(i) > 0;
        case "layered_percolation"
            T.report_ready(i) = T.recharge_ratio(i) > 0.01;
        case "shallow_groundwater"
            T.report_ready(i) = T.recharge_ratio(i) > 0.01 | T.mean_final_depth_mm(i) > 0;
        case "spatial_soil_contrast"
            T.report_ready(i) = T.mean_final_depth_mm(i) > 0;
    end
end
T.report_ready(:) = false;
end

function delete_dir_contents(d)
if ~exist(d, 'dir'); return; end
L = dir(d);
for i = 1:numel(L)
    if strcmp(L(i).name,'.') || strcmp(L(i).name,'..')
        continue
    end
    target = fullfile(d, L(i).name);
    if L(i).isdir
        rmdir(target, 's');
    else
        delete(target);
    end
end
end
