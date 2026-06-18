function run_snisb_hydropol2d_from_list(list_path, array_index, varargin)
%RUN_SNISB_HYDROPOL2D_FROM_LIST Run one HydroPol2D dam selected from a CSV list.
%
% The CSV is created by scripts/make_hydropol2d_dam_list.py and is designed
% for SLURM arrays. array_index is 1-based and normally equals
% SLURM_ARRAY_TASK_ID.

p = inputParser;
addRequired(p, 'list_path', @(x) ischar(x) || isstring(x));
addRequired(p, 'array_index', @(x) isnumeric(x) || ischar(x) || isstring(x));
addParameter(p, 'SimulationMinutes', 2 * 24 * 60, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RecordTimeMapsMinutes', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'RoutingModel', 'full_momentum', @(x) ischar(x) || isstring(x));
addParameter(p, 'SkipCompleted', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'CleanOutputFolder', true, @(x) islogical(x) || isnumeric(x));
parse(p, list_path, array_index, varargin{:});

list_path = char(p.Results.list_path);
if ischar(array_index) || isstring(array_index)
    array_index = str2double(char(array_index));
end
array_index = double(array_index);

if exist(list_path, 'file') ~= 2
    error('Dam list not found: %s', list_path);
end

T = readtable(list_path, 'Delimiter', ',', 'TextType', 'char');
if array_index < 1 || array_index > height(T)
    error('Array index %.0f is outside dam list range 1..%d.', array_index, height(T));
end

% dam_dir is always a string column
if iscell(T.dam_dir)
    dam_dir = char(T.dam_dir{array_index});
else
    dam_dir = char(T.dam_dir(array_index));
end

% snisb_code may be read as numeric (leading zeros stripped by MATLAB)
if isnumeric(T.snisb_code)
    snisb_code = sprintf('%06d', T.snisb_code(array_index));
elseif iscell(T.snisb_code)
    snisb_code = char(T.snisb_code{array_index});
else
    snisb_code = char(T.snisb_code(array_index));
end
fprintf('Sherlock HydroPol2D task %.0f/%d | SNISB %s\n', array_index, height(T), snisb_code);
fprintf('Dam folder: %s\n', dam_dir);

if logical(p.Results.SkipCompleted) && case_completed(dam_dir, p.Results.SimulationMinutes)
    fprintf('Skipping SNISB %s because a completed %.0f-minute run already exists.\n', ...
        snisb_code, p.Results.SimulationMinutes);
    return
end

run_snisb_hydropol2d_case( ...
    dam_dir, ...
    'CleanOutputFolder', logical(p.Results.CleanOutputFolder), ...
    'RunPostprocessing', true, ...
    'EnableLogging', true, ...
    'DryRun', false, ...
    'SimulationMinutes', double(p.Results.SimulationMinutes), ...
    'RecordTimeMapsMinutes', double(p.Results.RecordTimeMapsMinutes), ...
    'RoutingModel', char(p.Results.RoutingModel));
end

function tf = case_completed(dam_dir, min_minutes)
    summary_file = fullfile(dam_dir, 'hydropol2d', 'Modeling_Results', 'Tables_CSV', 'Summary_Table.csv');
    max_depth_file = fullfile(dam_dir, 'hydropol2d', 'Modeling_Results', 'Rasters_Static', 'Maximum_Depths.tif');
    manifest_file = fullfile(dam_dir, 'hydropol2d', 'snisb_hydropol2d_case_manifest.txt');
    tf = false;
    if exist(summary_file, 'file') ~= 2 || exist(max_depth_file, 'file') ~= 2 || exist(manifest_file, 'file') ~= 2
        return
    end
    text = fileread(manifest_file);
    token = regexp(text, 'Simulation minutes:\s*([0-9.]+)', 'tokens', 'once');
    if isempty(token)
        return
    end
    tf = str2double(token{1}) >= min_minutes;
end
