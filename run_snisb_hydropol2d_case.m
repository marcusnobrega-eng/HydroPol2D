function run_snisb_hydropol2d_case(dam_dir, varargin)
%RUN_SNISB_HYDROPOL2D_CASE Run HydroPol2D for one prepared SNISB dam folder.
%
% Required dam folder structure:
%   dam_dir/
%       rasters/domain_clipped/dem_fabdem_30m_domain.tif
%       rasters/domain_clipped/lulc_mapbiomas_30m_domain.tif
%       rasters/domain_clipped/soil_texture_usda_30m_domain.tif
%       hydrograph/inlet_cells_hydrograph.csv
%       outlet/outlet_cells.csv                 (optional but preferred)
%
% Outputs are written to:
%   dam_dir/hydropol2d/

p = inputParser;
addRequired(p, 'dam_dir', @(x) ischar(x) || isstring(x));
addParameter(p, 'CleanOutputFolder', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'RunPostprocessing', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'EnableLogging', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'DryRun', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SimulationMinutes', 2 * 24 * 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RecordTimeMapsMinutes', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RoutingModel', 'local_inertial', @(x) ischar(x) || isstring(x));
addParameter(p, 'WarmupDepthPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'EnableWarmup', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'EnableInflow', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'EnableRainfall', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'RainfallTimeseriesPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'EnableInfiltration', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'EnablePerimeterOutlet', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Manning', 0.035, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, dam_dir, varargin{:});

dam_dir = char(p.Results.dam_dir);
clean_output_folder = logical(p.Results.CleanOutputFolder);
run_postprocessing = logical(p.Results.RunPostprocessing);
enable_logging = logical(p.Results.EnableLogging);
dry_run = logical(p.Results.DryRun);
simulation_minutes = double(p.Results.SimulationMinutes);
record_time_maps_minutes = double(p.Results.RecordTimeMapsMinutes);
routing_model = lower(strtrim(char(p.Results.RoutingModel)));
warmup_depth_path = char(p.Results.WarmupDepthPath);
enable_warmup = logical(p.Results.EnableWarmup);
enable_inflow = logical(p.Results.EnableInflow);
enable_rainfall = logical(p.Results.EnableRainfall);
rainfall_timeseries_path = char(p.Results.RainfallTimeseriesPath);
enable_infiltration = logical(p.Results.EnableInfiltration);
enable_perimeter_outlet = logical(p.Results.EnablePerimeterOutlet);
manning_value = double(p.Results.Manning);
valid_routing_models = {'local_inertial','full_momentum'};
if ~any(strcmp(routing_model, valid_routing_models))
    error('RoutingModel must be one of: local_inertial, full_momentum.');
end

model_root = fileparts(mfilename('fullpath'));
config_dir = fullfile(model_root, 'Config');
hydropol2d_tools_user = fullfile(model_root, 'HydroPol2D_Functions');
topo_path_user = fullfile(model_root, 'topotoolbox-master');

input_paths_function = fullfile(config_dir, 'input_paths_bypass.m');
input_data_bypass_script_path = fullfile(config_dir, 'input_data_bypass_snisb.m');

domain_rasters = fullfile(dam_dir, 'rasters', 'domain_clipped');
dem_path = fullfile(domain_rasters, 'dem_fabdem_30m_domain.tif');
lulc_path = fullfile(domain_rasters, 'lulc_mapbiomas_30m_domain.tif');
soil_path = fullfile(domain_rasters, 'soil_texture_usda_30m_domain.tif');
inflow_csv = fullfile(dam_dir, 'hydrograph', 'inlet_cells_hydrograph.csv');
if isempty(rainfall_timeseries_path)
    rainfall_timeseries_path = fullfile(dam_dir, 'rainfall', 'rainfall_timeseries.csv');
end
if ~isempty(rainfall_timeseries_path) && ~is_absolute_path(rainfall_timeseries_path)
    rainfall_timeseries_path = fullfile(pwd, rainfall_timeseries_path);
end
outlet_csv = fullfile(dam_dir, 'outlet', 'outlet_cells.csv');
export_root_dir = fullfile(dam_dir, 'hydropol2d');

assert(exist(dem_path, 'file') == 2, 'DEM not found: %s', dem_path);
assert(exist(lulc_path, 'file') == 2, 'LULC not found: %s', lulc_path);
assert(exist(soil_path, 'file') == 2, 'SOIL not found: %s', soil_path);
if enable_inflow
    assert(exist(inflow_csv, 'file') == 2, 'Inflow hydrograph not found: %s', inflow_csv);
end
if enable_rainfall
    assert(exist(rainfall_timeseries_path, 'file') == 2, 'Rainfall hyetograph not found: %s', rainfall_timeseries_path);
end
if enable_warmup
    assert(~isempty(warmup_depth_path), 'WarmupDepthPath is required when EnableWarmup is true.');
    assert(exist(warmup_depth_path, 'file') == 2, 'Warmup depth raster not found: %s', warmup_depth_path);
end
assert(exist(input_paths_function, 'file') == 2, 'input_paths_bypass not found.');
assert(exist(input_data_bypass_script_path, 'file') == 2, 'SNISB input-data bypass not found.');
assert(exist(hydropol2d_tools_user, 'dir') == 7, 'HydroPol2D tools folder not found.');
assert(exist(topo_path_user, 'dir') == 7, 'TopoToolbox folder not found.');

Overrides = struct();
Overrides.DEM_path = dem_path;
Overrides.LULC_path = lulc_path;
Overrides.SOIL_path = soil_path;
Overrides.Inflow_Hydrograph_CSV = inflow_csv;
Overrides.Rainfall_Timeseries_File = rainfall_timeseries_path;
Overrides.Warmup_Depth_path = warmup_depth_path;
Overrides.EnableWarmup = enable_warmup;
Overrides.EnableInflow = enable_inflow;
Overrides.EnableRainfall = enable_rainfall;
Overrides.EnableInfiltration = enable_infiltration;
Overrides.EnablePerimeterOutlet = enable_perimeter_outlet;
Overrides.Manning = manning_value;
if exist(outlet_csv, 'file') == 2
    Overrides.Outlet_Cells_CSV = outlet_csv;
else
    warning('Outlet CSV not found; HydroPol2D will use its fallback outlet placement: %s', outlet_csv);
end
Overrides.export_root_dir = export_root_dir;
Overrides.SimulationMinutes = simulation_minutes;
Overrides.RecordTimeMapsMinutes = record_time_maps_minutes;
Overrides.RoutingModel = routing_model;

run_mode = 'bypass';
addpath(config_dir);
addpath(genpath(hydropol2d_tools_user));
addpath(genpath(topo_path_user));

[~, func_name, ~] = fileparts(input_paths_function);
InputPaths = feval(func_name, topo_path_user, hydropol2d_tools_user, Overrides);

fprintf('\n============================================================\n');
fprintf('SNISB HydroPol2D case\n');
fprintf('Dam folder       : %s\n', dam_dir);
fprintf('Export root      : %s\n', export_root_dir);
fprintf('Simulation min   : %.3f\n', simulation_minutes);
fprintf('Map record min   : %.3f\n', record_time_maps_minutes);
fprintf('Routing model    : %s\n', routing_model);
fprintf('Warmup enabled   : %d\n', enable_warmup);
fprintf('Inflow enabled   : %d\n', enable_inflow);
fprintf('Rainfall enabled : %d\n', enable_rainfall);
fprintf('Infiltration     : %d\n', enable_infiltration);
fprintf('Perimeter outlet : %d\n', enable_perimeter_outlet);
fprintf('Manning n        : %.4f\n', manning_value);
fprintf('Dry run          : %d\n', dry_run);
fprintf('============================================================\n\n');

if dry_run
    write_case_manifest(export_root_dir, dam_dir, dem_path, lulc_path, soil_path, inflow_csv, rainfall_timeseries_path, outlet_csv, warmup_depth_path, enable_warmup, enable_inflow, enable_rainfall, enable_infiltration, enable_perimeter_outlet, dry_run, simulation_minutes, record_time_maps_minutes, routing_model, manning_value);
    fprintf('Dry run requested; paths and InputPaths were validated only.\n');
    return
end

if ~exist(export_root_dir, 'dir')
    mkdir(export_root_dir);
end

% Variables expected by HydroPol2D scripts.
use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
model_folder = '';
GD = [];

Paths = init_results_tree(export_root_dir, clean_output_folder);
ExportRootDir = Paths.Root;
resultsDir = Paths.Results;
write_case_manifest(export_root_dir, dam_dir, dem_path, lulc_path, soil_path, inflow_csv, rainfall_timeseries_path, outlet_csv, warmup_depth_path, enable_warmup, enable_inflow, enable_rainfall, enable_infiltration, enable_perimeter_outlet, dry_run, simulation_minutes, record_time_maps_minutes, routing_model, manning_value);

assignin('base', 'Paths', Paths);
assignin('base', 'ExportRootDir', ExportRootDir);
assignin('base', 'resultsDir', resultsDir);
assignin('base', 'model_folder', model_folder);
assignin('base', 'GD', GD);
assignin('base', 'use_inputpaths_bypass', use_inputpaths_bypass);
assignin('base', 'use_inputdata_bypass', use_inputdata_bypass);
assignin('base', 'InputPaths', InputPaths);
assignin('base', 'Overrides', Overrides);
assignin('base', 'input_data_bypass_script_path', input_data_bypass_script_path);
assignin('base', 'input_paths_function', input_paths_function);

if enable_logging
    if ~exist(Paths.Logs, 'dir')
        mkdir(Paths.Logs);
    end
    diary('off');
    diary(fullfile(Paths.Logs, 'run_log.txt'));
    diary('on');
end

try
    fprintf('STEP 1/3 | HydroPol2D_preprocessing\n');
    HydroPol2D_preprocessing;

    fprintf('STEP 2/3 | HydroPol2D_Main_While\n');
    HydroPol2D_Main_While;

    if run_postprocessing
        fprintf('STEP 3/3 | post_processing\n');
        close all;
        post_processing;
    end

    if enable_logging
        diary('off');
    end
catch ME
    if enable_logging
        diary('off');
    end
    rethrow(ME);
end

fprintf('\nHydroPol2D SNISB case completed: %s\n', export_root_dir);
end

function write_case_manifest(export_root_dir, dam_dir, dem_path, lulc_path, soil_path, inflow_csv, rainfall_timeseries_path, outlet_csv, warmup_depth_path, enable_warmup, enable_inflow, enable_rainfall, enable_infiltration, enable_perimeter_outlet, dry_run, simulation_minutes, record_time_maps_minutes, routing_model, manning_value)
    if ~exist(export_root_dir, 'dir')
        mkdir(export_root_dir);
    end
    manifest_path = fullfile(export_root_dir, 'snisb_hydropol2d_case_manifest.txt');
    fid = fopen(manifest_path, 'w');
    if fid < 0
        warning('Could not write case manifest: %s', manifest_path);
        return
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'SNISB HydroPol2D case manifest\n');
    fprintf(fid, 'Created: %s\n', char(datetime('now')));
    fprintf(fid, 'Dry run: %d\n', dry_run);
    fprintf(fid, 'Dam folder: %s\n', dam_dir);
    fprintf(fid, 'DEM: %s\n', dem_path);
    fprintf(fid, 'LULC: %s\n', lulc_path);
    fprintf(fid, 'SOIL: %s\n', soil_path);
    fprintf(fid, 'Inflow hydrograph: %s\n', inflow_csv);
    fprintf(fid, 'Rainfall hyetograph: %s\n', rainfall_timeseries_path);
    fprintf(fid, 'Outlet cells: %s\n', outlet_csv);
    fprintf(fid, 'Warmup depth raster: %s\n', warmup_depth_path);
    fprintf(fid, 'Warmup enabled: %d\n', enable_warmup);
    fprintf(fid, 'Inflow enabled: %d\n', enable_inflow);
    fprintf(fid, 'Rainfall enabled: %d\n', enable_rainfall);
    fprintf(fid, 'Infiltration enabled: %d\n', enable_infiltration);
    fprintf(fid, 'Perimeter outlet enabled: %d\n', enable_perimeter_outlet);
    fprintf(fid, 'Simulation minutes: %.6g\n', simulation_minutes);
    fprintf(fid, 'Map output interval minutes: %.6g\n', record_time_maps_minutes);
    fprintf(fid, 'Routing model: %s\n', routing_model);
    fprintf(fid, 'Manning n: %.6g\n', manning_value);
    fprintf(fid, 'Configuration: rainfall=%d, infiltration=%d, inflow=%d, warmup=%d, perimeter_outlet=%d, %s.\n', enable_rainfall, enable_infiltration, enable_inflow, enable_warmup, enable_perimeter_outlet, routing_model);
end

function Paths = init_results_tree(exportRootDir, cleanOutputFolder)
    if ~exist(exportRootDir,'dir')
        mkdir(exportRootDir);
    elseif cleanOutputFolder && ~is_dir_empty(exportRootDir)
        delete_dir_contents(exportRootDir);
    end

    Paths = struct();
    Paths.Root = exportRootDir;
    Paths.Results = fullfile(exportRootDir, 'Modeling_Results');
    Paths.Temp = fullfile(exportRootDir, 'Temporary_Files');
    Paths.Logs = fullfile(exportRootDir, 'Logs');
    Paths.FigPDF = fullfile(Paths.Results, 'Figures_PDF');
    Paths.FigFIG = fullfile(Paths.Results, 'Figures_FIG');
    Paths.Tables = fullfile(Paths.Results, 'Tables_CSV');
    Paths.RastersWD = fullfile(Paths.Results, 'Rasters_Water_Depths');
    Paths.RastersWSE = fullfile(Paths.Results, 'Rasters_WSE');
    Paths.RastersStatic = fullfile(Paths.Results, 'Rasters_Static');
    Paths.RastersVelocity     = fullfile(Paths.Results, 'Rasters_Velocity');
    Paths.RastersHazard       = fullfile(Paths.Results, 'Rasters_Hazard');
    Paths.RastersInfiltration = fullfile(Paths.Results, 'Rasters_Infiltration');
    Paths.WQMaps = fullfile(Paths.Results, 'Rasters_WQ');
    Paths.HRMaps = fullfile(Paths.Results, 'Rasters_Human_Risk');
    Paths.Anim = fullfile(Paths.Results, 'GIFs_MP4');
    Paths.Shapes = fullfile(Paths.Results, 'Shapefiles');

    fields = fieldnames(Paths);
    for i = 1:numel(fields)
        d = Paths.(fields{i});
        if ischar(d) && ~exist(d, 'dir')
            mkdir(d);
        end
    end
end

function tf = is_dir_empty(d)
    L = dir(d);
    names = {L.name};
    tf = all(ismember(names, {'.','..'}));
end

function delete_dir_contents(d)
    L = dir(d);
    for i = 1:numel(L)
        name = L(i).name;
        if strcmp(name,'.') || strcmp(name,'..')
            continue;
        end
        target = fullfile(d, name);
        if L(i).isdir
            rmdir(target, 's');
        else
            delete(target);
        end
    end
end

function tf = is_absolute_path(path_value)
    path_value = char(path_value);
    tf = startsWith(path_value, filesep) || ...
         (~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')));
end
