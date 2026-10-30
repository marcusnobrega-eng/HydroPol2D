function hp2d_refuse_mesh_modules(module_flags, case_config, source_label)
%HP2D_REFUSE_MESH_MODULES Stop a mesh run that asks for a module it would skip.
%   module_flags is a struct of 0/1 switches named as in the grid inputs
%   (flag_snow_modeling, flag_reservoir, ...). A mesh run reads its modules
%   from the prepared case (case_config.<field>), not from the grid inputs, so
%   a switch that is on while the case carries no data for that module would
%   be silently dropped. This stops the run instead, as HydroPol2D-Python does
%   (MeshFeatureNotSupportedError). The deep aquifer never runs on meshes in
%   either code.

if nargin < 2 || isempty(case_config), case_config = struct(); end
if nargin < 3, source_label = 'the run inputs'; end

% flag name, label, prepared-case field ('' = not available on meshes)
modules = { ...
    'flag_snow_modeling',        'snow',                   ''; ...
    'flag_reservoir',            'reservoir controls',     'reservoirs'; ...
    'flag_waterquality',         'water quality',          'water_quality'; ...
    'flag_human_instability',    'human risk',             'human_risk'; ...
    'flag_deep_aquifer',         'the deep aquifer',       ''; ...
    'flag_dam_break',            'dam break',              ''; ...
    'flag_detention_ponds',      'detention ponds',        ''; ...
    'flag_urban_drainage',       'urban drainage',         ''; ...
    'flag_green_infrastructure', 'green infrastructure',   ''};

refused = {};
for k = 1:size(modules,1)
    if ~isfield(module_flags, modules{k,1}) || ~is_on(module_flags.(modules{k,1}))
        continue
    end
    field = modules{k,3};
    if ~isempty(field) && isfield(case_config, field) && ~isempty(case_config.(field))
        continue
    end
    if isempty(field)
        refused{end+1} = sprintf('%s (%s = 1): not available on meshes', modules{k,2}, modules{k,1}); %#ok<AGROW>
    else
        refused{end+1} = sprintf('%s (%s = 1): the prepared case has no config.%s', ...
            modules{k,2}, modules{k,1}, field); %#ok<AGROW>
    end
end

if ~isempty(refused)
    error('HydroPol2D:MeshModuleNotSupported', ...
        ['%s switch on modules that this mesh run would skip:\n  - %s\n' ...
         'Switch them off, add them to the prepared case, or run on the regular (D4) grid.'], ...
        source_label, strjoin(refused, sprintf('\n  - ')));
end
end

function tf = is_on(value)
if ischar(value) || isstring(value)
    value = str2double(value);
end
tf = ~isempty(value) && isscalar(value) && isfinite(double(value)) && double(value) ~= 0;
end
