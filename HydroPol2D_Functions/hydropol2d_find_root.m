function model_root = hydropol2d_find_root(start_path)
%HYDROPOL2D_FIND_ROOT Locate a HydroPol2D repository from a descendant path.

if nargin < 1 || isempty(start_path)
    start_path = pwd;
end

model_root = char(start_path);
if isfile(model_root)
    model_root = fileparts(model_root);
end
model_root = canonical_path(model_root);

while ~isfolder(fullfile(model_root, 'HydroPol2D_Functions'))
    parent = fileparts(model_root);
    if strcmp(parent, model_root)
        error('HydroPol2D:Runtime:ModelRootNotFound', ...
            'Could not locate a HydroPol2D repository above: %s', start_path);
    end
    model_root = parent;
end

function path_value = canonical_path(path_value)
% Resolve relative components so downstream path-prefix checks are stable.
path_value = char(java.io.File(path_value).getCanonicalPath());
end
end
