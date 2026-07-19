function runtime = hydropol2d_add_runtime_paths(model_root)
%HYDROPOL2D_ADD_RUNTIME_PATHS Register the bundled HydroPol2D runtime.
%   HYDROPOL2D_ADD_RUNTIME_PATHS() derives the repository root from this
%   file. HYDROPOL2D_ADD_RUNTIME_PATHS(MODEL_ROOT) registers the model
%   functions and the curated TopoToolbox v2.4 runtime shipped with
%   HydroPol2D. No external TopoToolbox installation is used.

if nargin < 1 || isempty(model_root)
    model_root = fileparts(fileparts(mfilename('fullpath')));
end
model_root = char(model_root);

functions_root = fullfile(model_root, 'HydroPol2D_Functions');
topotoolbox_root = fullfile(model_root, 'third_party', 'topotoolbox_lite');

if ~isfolder(functions_root)
    error('HydroPol2D:Runtime:MissingFunctions', ...
        'HydroPol2D functions folder not found: %s', functions_root);
end
if ~isfolder(topotoolbox_root)
    error('HydroPol2D:Runtime:MissingTopoToolboxLite', ...
        'Bundled TopoToolbox Lite runtime not found: %s', topotoolbox_root);
end

% Class folders below topotoolbox_root are resolved by MATLAB from this one
% parent directory; private helpers resolve from their owning class folders.
addpath(functions_root, '-begin');
addpath(topotoolbox_root, '-begin');

required_symbols = {'GRIDobj', 'FLOWobj', 'STREAMobj'};
for i = 1:numel(required_symbols)
    resolved_file = which(required_symbols{i});
    if isempty(resolved_file) || ~startsWith(resolved_file, topotoolbox_root)
        error('HydroPol2D:Runtime:BundledTopoToolboxNotResolved', ...
            'Expected %s to resolve inside %s, but MATLAB resolved: %s', ...
            required_symbols{i}, topotoolbox_root, resolved_file);
    end
end

runtime = struct();
runtime.model_root = model_root;
runtime.functions_root = functions_root;
runtime.topotoolbox_lite_root = topotoolbox_root;
end
