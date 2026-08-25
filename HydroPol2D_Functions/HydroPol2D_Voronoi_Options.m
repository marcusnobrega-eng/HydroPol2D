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
    'river_preferred_cells_across', 3, ...
    'unresolved_river_policy', 'neal_subgrid');

options = defaults;
names = fieldnames(defaults);
for k = 1:numel(names)
    if isfield(input_options, names{k}) && ~isempty(input_options.(names{k}))
        options.(names{k}) = input_options.(names{k});
    end
end

if double(flag_voronoi) == 0
    return
end

positive = {'background_target_width_m','minimum_cell_width_m', ...
    'maximum_adjacent_size_ratio','urban_target_width_m','river_preferred_cells_across'};
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
if options.river_preferred_cells_across ~= round(options.river_preferred_cells_across)
    error('HydroPol2D:InvalidVoronoiOptions', ...
        'river_preferred_cells_across must be a positive integer.');
end
options.unresolved_river_policy = char(lower(strtrim(string(options.unresolved_river_policy))));
if ~any(strcmp(options.unresolved_river_policy, {'neal_subgrid','none'}))
    error('HydroPol2D:InvalidVoronoiOptions', ...
        'unresolved_river_policy must be ''neal_subgrid'' or ''none''.');
end
end
