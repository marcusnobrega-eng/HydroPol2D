function run_code(mesh_mode)
%RUN_CODE Run the V-Tilted example without reading model inputs from Excel.

if nargin < 1 || strlength(string(mesh_mode)) == 0
    mesh_mode = 'regular';
end
mesh_mode = lower(char(string(mesh_mode)));
if ~ismember(mesh_mode, {'regular', 'voronoi'})
    error('HydroPol2D:VTiltedMeshMode', ...
        'mesh_mode must be ''regular'' or ''voronoi''.');
end

example_root = fileparts(mfilename('fullpath'));
model_root = fileparts(fileparts(example_root));
config_root = fullfile(example_root, 'Config');

environment_names = { ...
    'HYDROPOL_RUN_MODE', 'HYDROPOL_VTILTED_MESH', ...
    'HYDROPOL_INPUT_PATHS_FUNCTION', 'HYDROPOL_INPUT_DATA_BYPASS_SCRIPT', ...
    'HYDROPOL_EXPORT_ROOT_DIR', 'HYDROPOL_SKIP_POSTPROCESS'};
previous_environment = cellfun(@getenv, environment_names, 'UniformOutput', false);
restore_environment = onCleanup(@() restore_values( ...
    environment_names, previous_environment));

setenv('HYDROPOL_RUN_MODE', 'bypass');
setenv('HYDROPOL_VTILTED_MESH', mesh_mode);
setenv('HYDROPOL_INPUT_PATHS_FUNCTION', ...
    fullfile(config_root, 'vtilted_input_paths.m'));
setenv('HYDROPOL_INPUT_DATA_BYPASS_SCRIPT', ...
    fullfile(config_root, 'vtilted_input_data.m'));
setenv('HYDROPOL_EXPORT_ROOT_DIR', ...
    fullfile(example_root, 'Outputs', ['code_' mesh_mode]));
setenv('HYDROPOL_SKIP_POSTPROCESS', '');

launcher = strrep(fullfile(model_root, 'HydroPol2D_V115.m'), '''', '''''');
evalin('base', sprintf('run(''%s'')', launcher));
end

function restore_values(names, values)
for index = 1:numel(names)
    setenv(names{index}, values{index});
end
end
