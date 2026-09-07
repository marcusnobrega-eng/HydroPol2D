function options = HydroPol2D_Voronoi_Options(flag_voronoi, input_options)
%HYDROPOL2D_VORONOI_OPTIONS Normalize the compact Voronoi user interface.

if nargin < 1 || isempty(flag_voronoi)
    flag_voronoi = 0;
end
if nargin < 2 || isempty(input_options)
    input_options = struct();
end
if ~isscalar(flag_voronoi) || ~isnumeric(flag_voronoi) || ~ismember(double(flag_voronoi), [0 1])
    error('HydroPol2D:InvalidVoronoiFlag', 'flag_voronoi must be 0 or 1.');
end
if ~isstruct(input_options)
    error('HydroPol2D:InvalidVoronoiOptions', 'Voronoi options must be a structure.');
end

defaults = struct( ...
    'background_target_width_m', 2000, ...
    'minimum_cell_width_m', 100, ...
    'maximum_adjacent_size_ratio', 2, ...
    'urban_target_width_m', 200, ...
    'urban_transition_buffer_m', 1000, ...
    'unresolved_river_policy', 'neal_subgrid', ...
    'voronoi_subgrid_enabled', [], ...
    'subgrid_table_path', '');

% ``subgrid_table_path`` points at the HEC-RAS style property tables that
% HydroBathyDEM writes next to a mesh (cell elevation-volume, face conveyance).
% Leave it empty for the flat-prism closure.  Which cells use the tables is
% decided at table-build time and carried in ``cell_is_subgrid``; a face uses
% them only where BOTH its cells do.  This is a different method from
% Subgrid_Properties_Lookup.m, which builds SFINCS-style tables for structured
% grids and cannot be applied to an unstructured mesh.

% ``river_preferred_cells_across`` is retired.  It was a mesh-generator hint and
% cannot influence a run that loads a finished mesh; the mesh handoff dropped it
% because it forced narrow rivers into artificially wide corridors.  Accept and
% ignore it so existing case files keep loading.
if isfield(input_options, 'river_preferred_cells_across')
    warning('HydroPol2D:RetiredVoronoiOption', ...
        'river_preferred_cells_across is retired and ignored; mesh geometry comes from the mesh file.');
    input_options = rmfield(input_options, 'river_preferred_cells_across');
end

options = defaults;
names = fieldnames(defaults);
for k = 1:numel(names)
    if isfield(input_options, names{k}) && ~isempty(input_options.(names{k}))
        options.(names{k}) = input_options.(names{k});
    end
end

subgrid_choice_supplied = isfield(input_options,'voronoi_subgrid_enabled') && ...
    ~isempty(input_options.voronoi_subgrid_enabled);
if isempty(options.voronoi_subgrid_enabled)
    options.voronoi_subgrid_enabled = strlength(strtrim(string(options.subgrid_table_path))) > 0;
else
    value = options.voronoi_subgrid_enabled;
    if islogical(value), value=double(value); end
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ~ismember(double(value),[0 1])
        error('HydroPol2D:InvalidVoronoiOptions', ...
            'voronoi_subgrid_enabled must be 0 or 1.');
    end
    options.voronoi_subgrid_enabled = logical(value);
end
if subgrid_choice_supplied && ~options.voronoi_subgrid_enabled && ...
        strlength(strtrim(string(options.subgrid_table_path))) > 0
    error('HydroPol2D:VoronoiSubgridConfigurationConflict', ...
        'A subgrid table path was provided while Voronoi subgrid is disabled.');
end

if double(flag_voronoi) == 0
    return
end

positive = {'background_target_width_m','minimum_cell_width_m', ...
    'maximum_adjacent_size_ratio','urban_target_width_m'};
for k = 1:numel(positive)
    value = options.(positive{k});
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || value <= 0
        error('HydroPol2D:InvalidVoronoiOptions', '%s must be a positive finite scalar.', positive{k});
    end
end
if ~isscalar(options.urban_transition_buffer_m) || ~isnumeric(options.urban_transition_buffer_m) || ...
        ~isfinite(options.urban_transition_buffer_m) || options.urban_transition_buffer_m < 0
    error('HydroPol2D:InvalidVoronoiOptions', ...
        'urban_transition_buffer_m must be a nonnegative finite scalar.');
end
if options.background_target_width_m < options.minimum_cell_width_m || ...
        options.urban_target_width_m < options.minimum_cell_width_m
    error('HydroPol2D:InvalidVoronoiOptions', ...
        'Background and urban target widths cannot be below minimum_cell_width_m.');
end
if options.maximum_adjacent_size_ratio < 1
    error('HydroPol2D:InvalidVoronoiOptions', ...
        'maximum_adjacent_size_ratio must be at least 1.');
end
options.unresolved_river_policy = char(lower(strtrim(string(options.unresolved_river_policy))));
if ~any(strcmp(options.unresolved_river_policy, {'neal_subgrid','none'}))
    error('HydroPol2D:InvalidVoronoiOptions', ...
        'unresolved_river_policy must be ''neal_subgrid'' or ''none''.');
end
end
