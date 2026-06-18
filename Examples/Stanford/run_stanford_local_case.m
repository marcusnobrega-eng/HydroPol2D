function run_stanford_local_case(varargin)
%RUN_STANFORD_LOCAL_CASE Run the local Stanford example with current code.
%
% This runner is intentionally separate from HydroPol2D_V115.m, which still
% documents the original cluster-oriented workflow. It uses the local
% HydroPol2D_Model tool folders and case inputs under Examples/Stanford.

p = inputParser;
addParameter(p, 'Resolution', 30, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'SimulationMinutes', 180, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RecordTimeMapsMinutes', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RecordTimeHydrographsMinutes', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'StaticFolder', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'RainfallTimeseriesPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputTag', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'RunPostprocessing', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'CleanOutputFolder', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'EnableLogging', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'UseGPU', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'UseSingle', false, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

case_root = fileparts(mfilename('fullpath'));
model_root = fileparts(fileparts(case_root));
config_dir = fullfile(case_root, 'Config');
hydropol2d_tools_user = fullfile(model_root, 'HydroPol2D_Functions');
topo_path_user = fullfile(model_root, 'topotoolbox-master');
input_paths_function = fullfile(config_dir, 'input_paths_bypass.m');
input_data_bypass_script_path = fullfile(config_dir, 'input_data_bypass_script.m');

resolution = double(p.Results.Resolution);
simulation_minutes = double(p.Results.SimulationMinutes);
record_maps_min = double(p.Results.RecordTimeMapsMinutes);
record_hydro_min = double(p.Results.RecordTimeHydrographsMinutes);
run_postprocessing = logical(p.Results.RunPostprocessing);
clean_output_folder = logical(p.Results.CleanOutputFolder);
enable_logging = logical(p.Results.EnableLogging);

static_folder = char(p.Results.StaticFolder);
if isempty(static_folder)
    static_folder = fullfile(case_root, 'Static');
end
rainfall_path = char(p.Results.RainfallTimeseriesPath);
if isempty(rainfall_path)
    rainfall_path = fullfile(case_root, 'Forcing', 'Rainfall', 'Rainfall_Intensity_Event_Clean.csv');
end
output_tag = char(p.Results.OutputTag);
if isempty(output_tag)
    output_tag = sprintf('Stanford_%dm_Event', round(resolution));
end
export_root_dir = fullfile(case_root, 'Outputs', 'Reruns_CurrentModel', output_tag);

assert(exist(hydropol2d_tools_user, 'dir') == 7, 'HydroPol2D tools not found: %s', hydropol2d_tools_user);
assert(exist(topo_path_user, 'dir') == 7, 'TopoToolbox not found: %s', topo_path_user);
assert(exist(input_paths_function, 'file') == 2, 'input_paths_bypass not found: %s', input_paths_function);
assert(exist(input_data_bypass_script_path, 'file') == 2, 'input_data_bypass_script not found: %s', input_data_bypass_script_path);
assert(exist(rainfall_path, 'file') == 2, 'Rainfall timeseries not found: %s', rainfall_path);

dem_path = fullfile(static_folder, 'DEM.tif');
lulc_path = fullfile(static_folder, 'LULC.tif');
soil_path = fullfile(static_folder, 'SOIL.tif');
dtb_path = fullfile(static_folder, 'DTB.tif');
lai_path = fullfile(static_folder, 'LAI.tif');
albedo_path = fullfile(static_folder, 'Albedo.tif');
initial_sm_path = fullfile(static_folder, 'Initial_Soil_Moisture.tif');

assert(exist(dem_path, 'file') == 2, 'DEM not found: %s', dem_path);
assert(exist(lulc_path, 'file') == 2, 'LULC not found: %s', lulc_path);
assert(exist(soil_path, 'file') == 2, 'SOIL not found: %s', soil_path);

addpath(config_dir);
addpath(genpath(hydropol2d_tools_user));
addpath(genpath(topo_path_user));
set(0, 'DefaultFigureVisible', 'off');

Overrides = struct();
Overrides.DEM_path = dem_path;
Overrides.LULC_path = lulc_path;
Overrides.SOIL_path = soil_path;
Overrides.DTB_path = dtb_path;
Overrides.LAI_path = lai_path;
Overrides.Albedo_path = albedo_path;
Overrides.Initial_Soil_Moisture_path = initial_sm_path;
Overrides.Rainfall_Timeseries_File = rainfall_path;
Overrides.export_root_dir = export_root_dir;
Overrides.resolution_resample = resolution;
Overrides.RecordTimeMapsMinutes = record_maps_min;
Overrides.RecordTimeHydrographsMinutes = record_hydro_min;
Overrides.DateBegin = datetime(2025, 5, 1, 0, 0, 0);
Overrides.DateEnd = Overrides.DateBegin + minutes(simulation_minutes);
Overrides.UseGPU = logical(p.Results.UseGPU);
Overrides.UseSingle = logical(p.Results.UseSingle);

[~, func_name, ~] = fileparts(input_paths_function);
InputPaths = feval(func_name, topo_path_user, hydropol2d_tools_user, Overrides);

run_mode = 'bypass'; %#ok<NASGU>
use_inputpaths_bypass = 1;
use_inputdata_bypass = 1;
model_folder = '';
GD = [];

if ~exist(export_root_dir, 'dir')
    mkdir(export_root_dir);
end
Paths = init_results_tree(export_root_dir, clean_output_folder);
ExportRootDir = Paths.Root; %#ok<NASGU>
resultsDir = Paths.Results; %#ok<NASGU>

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

fprintf('\n============================================================\n');
fprintf('Stanford current-model rerun\n');
fprintf('Case root       : %s\n', case_root);
fprintf('Static folder   : %s\n', static_folder);
fprintf('Resolution      : %.3f m\n', resolution);
fprintf('Simulation      : %.3f min\n', simulation_minutes);
fprintf('Rainfall file   : %s\n', rainfall_path);
fprintf('Output root     : %s\n', export_root_dir);
fprintf('GPU/single      : %d / %d\n', Overrides.UseGPU, Overrides.UseSingle);
fprintf('============================================================\n\n');

if enable_logging
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

fprintf('Stanford current-model rerun complete: %s\n', export_root_dir);
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
    Paths.RastersVelocity = fullfile(Paths.Results, 'Rasters_Velocity');
    Paths.RastersHazard = fullfile(Paths.Results, 'Rasters_Hazard');
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
