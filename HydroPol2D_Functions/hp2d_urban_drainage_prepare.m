function N = hp2d_urban_drainage_prepare(Network, terrain_m, valid, cell_area_m2, hydraulics, inlet_rc)
%HP2D_URBAN_DRAINAGE_PREPARE Check a drainage network and build its arrays.
%   N = hp2d_urban_drainage_prepare(Network, terrain_m, valid, cell_area_m2,
%   hydraulics, inlet_rc) is prepare_urban_drainage of HydroPol2D-Python
%   (modules/urban_drainage.py): the same checks, the same arrays.
%
%   Network:   tables nodes, pipes, inlets, outfalls (hp2d_urban_drainage_read).
%   terrain_m: surface elevation of the grid [m]; valid: active cells.
%   inlet_rc:  [row col] of each inlet in this grid (1-based). Empty uses the
%              inlet layer's row and column (0-based, as Python writes them) + 1.
%   hydraulics: struct, missing fields take these defaults:
%       courant 0.5, min_dt_s 0.01, max_dt_s 5, maximum_substeps 10000,
%       maximum_velocity_m_s 10 (a numerical guard, not a design limit),
%       inlet_defaults.opening_area_m2 0.10, .discharge_coefficient 0.67,
%       .reverse_flow_allowed true, free_outfall_depth "critical_normal".
%   free_outfall_depth: "critical_normal" (default) puts a free outfall at
%   min(critical, normal) depth of its pipe's flow; "invert" keeps it at the
%   node invert as HydroPol2D-Python does, which throttles the pipe to half
%   its upstream depth (Python code review 2026-10-26, F2). Everything else
%   is HydroPol2D-Python's rule.
%
%   Boundary types: -1 dynamic node, 0 free, 1 fixed stage, 2 normal depth.

if nargin < 5 || isempty(hydraulics)
    hydraulics = struct();
end
if nargin < 6
    inlet_rc = [];
end
H = hp2d_urban_drainage_options(hydraulics);
terrain_m = double(terrain_m);
valid = logical(valid);
if ~isequal(size(terrain_m), size(valid))
    fail('Urban drainage requires grid-aligned terrain and mask arrays.');
end
if ~(isfinite(cell_area_m2) && cell_area_m2 > 0)
    fail('Urban drainage requires a positive surface-cell area.');
end

nodes = Network.nodes;
pipes = Network.pipes;
inlets = Network.inlets;
outfalls = Network.outfalls;
if height(nodes) == 0 || height(pipes) == 0 || height(outfalls) == 0
    fail('Urban drainage requires nodes, pipes, and at least one outfall.');
end

node_ids = text_column(nodes, 'node_id');
pipe_ids = text_column(pipes, 'pipe_id');
if height(inlets) > 0
    inlet_ids = text_column(inlets, 'identifier');
else
    inlet_ids = strings(0, 1);
end
require_unique(node_ids, 'node');
require_unique(pipe_ids, 'pipe');
require_unique(inlet_ids, 'inlet');
n_nodes = numel(node_ids);

invert = number_column(nodes, 'invert_elevation_m');
rim = number_column(nodes, 'rim_elevation_m');
shaft = number_column(nodes, 'shaft_area_m2');
if any(~isfinite(invert)) || any(~isfinite(rim))
    fail('Drainage node elevations must be finite.');
end
if any(rim < invert)
    fail('Drainage node rim elevations cannot be below their inverts.');
end
if any(~isfinite(shaft)) || any(shaft <= 0)
    fail('Drainage node shaft areas must be finite and positive.');
end

boundary_type = -ones(n_nodes, 1);
boundary_head = invert;
outfall_nodes = text_column(outfalls, 'node_id');
require_unique(outfall_nodes, 'outfall node');
kinds = lower(text_column(outfalls, 'boundary_type'));
stages = NaN(height(outfalls), 1);
if ismember('fixed_stage_m', outfalls.Properties.VariableNames)
    stages = number_column(outfalls, 'fixed_stage_m');
end
outfall_index = zeros(height(outfalls), 1);
for i = 1:height(outfalls)
    index = find(node_ids == outfall_nodes(i), 1);
    if isempty(index)
        fail(sprintf('Outfall references missing node ''%s''.', outfall_nodes(i)));
    end
    outfall_index(i) = index;
    switch kinds(i)
        case "free"
            boundary_type(index) = 0;
        case "fixed"
            boundary_type(index) = 1;
            if ~isfinite(stages(i))
                fail(sprintf('Fixed outfall ''%s'' requires a finite stage.', outfall_nodes(i)));
            end
            boundary_head(index) = stages(i);
        case "normal"
            boundary_type(index) = 2;
        otherwise
            fail(sprintf('Unsupported outfall boundary type ''%s''.', kinds(i)));
    end
end

start_ids = text_column(pipes, 'upstream_node');
end_ids = text_column(pipes, 'downstream_node');
[found_start, pipe_start] = ismember(start_ids, node_ids);
[found_end, pipe_end] = ismember(end_ids, node_ids);
if any(~found_start) || any(~found_end)
    missing = unique([start_ids(~found_start); end_ids(~found_end)]);
    fail(sprintf('Pipes reference missing drainage nodes: %s.', strjoin(missing(1:min(5, end)), ', ')));
end
if any(pipe_start == pipe_end)
    fail('Drainage pipes cannot begin and end at the same node.');
end
L = number_column(pipes, 'length_m');
D = number_column(pipes, 'diameter_m');
n_manning = number_column(pipes, 'manning_n');
start_invert = number_column(pipes, 'upstream_invert_m');
end_invert = number_column(pipes, 'downstream_invert_m');
check_positive(L, 'length');
check_positive(D, 'diameter');
check_positive(n_manning, 'Manning roughness');
if any(~isfinite(start_invert)) || any(~isfinite(end_invert))
    fail('Drainage pipe endpoint inverts must be finite.');
end
if any(start_invert < invert(pipe_start) - 1e-6) || any(end_invert < invert(pipe_end) - 1e-6)
    fail('A pipe endpoint invert cannot be below its manhole invert.');
end
if any(boundary_type(pipe_start) ~= -1 & boundary_type(pipe_end) ~= -1)
    fail('A drainage pipe cannot connect two outfall boundary nodes.');
end

% Every connected component needs an outfall.
G = graph(pipe_start, pipe_end, [], n_nodes);
component = conncomp(G);
has_outfall = false(max(component), 1);
has_outfall(component(outfall_index)) = true;
if any(~has_outfall)
    orphan = find(~has_outfall, 1);
    members = sort(node_ids(component == orphan));
    fail(sprintf('Drainage component containing ''%s'' has no outfall.', members(1)));
end

normal_slope = zeros(numel(L), 1);
at_start = boundary_type(pipe_start) == 2;
at_end = boundary_type(pipe_end) == 2;
normal_slope(at_start) = (end_invert(at_start) - start_invert(at_start)) ./ L(at_start);
normal_slope(at_end & ~at_start) = (start_invert(at_end & ~at_start) - end_invert(at_end & ~at_start)) ./ L(at_end & ~at_start);
normal_links = at_start | at_end;
if any(normal_links & normal_slope <= 0)
    bad = find(normal_links & normal_slope <= 0, 1);
    fail(sprintf('Normal-depth outfall pipe ''%s'' requires a positive slope toward the outfall.', pipe_ids(bad)));
end
normal_index = find(normal_links);

% Inlets.
n_inlets = numel(inlet_ids);
if isempty(inlet_rc)
    if n_inlets > 0
        inlet_rc = [integer_column(inlets, 'row'), integer_column(inlets, 'column')] + 1;
    else
        inlet_rc = zeros(0, 2);
    end
end
[ny, nx] = size(terrain_m);
if n_inlets > 0
    if any(inlet_rc(:, 1) < 1 | inlet_rc(:, 1) > ny | inlet_rc(:, 2) < 1 | inlet_rc(:, 2) > nx)
        fail('Drainage inlet row/column indices must lie inside the grid.');
    end
    inlet_cell = sub2ind([ny nx], inlet_rc(:, 1), inlet_rc(:, 2));
    if any(~valid(inlet_cell))
        fail('Every drainage inlet must map to an active surface cell.');
    end
    [found, inlet_node] = ismember(text_column(inlets, 'node_id'), node_ids);
    if any(~found)
        fail('Inlets reference missing drainage nodes.');
    end
    if ismember('crest_elevation_m', inlets.Properties.VariableNames)
        crest = number_column(inlets, 'crest_elevation_m');
        if ismember('surface_elevation_m', inlets.Properties.VariableNames)
            fallback = number_column(inlets, 'surface_elevation_m');
            crest(~isfinite(crest)) = fallback(~isfinite(crest));
        end
    else
        crest = number_column(inlets, 'surface_elevation_m');
    end
    opening = with_fallback(inlets, 'opening_area_m2', H.inlet_defaults.opening_area_m2, n_inlets);
    coefficient = with_fallback(inlets, 'discharge_coefficient', H.inlet_defaults.discharge_coefficient, n_inlets);
    reverse = with_fallback(inlets, 'reverse_flow_allowed', double(H.inlet_defaults.reverse_flow_allowed), n_inlets) ~= 0;
    if any(~isfinite(crest))
        fail('Drainage inlet crest elevations must be finite.');
    end
    if any(~isfinite(opening)) || any(opening <= 0)
        fail('Drainage inlet opening areas must be finite and positive.');
    end
    if any(~isfinite(coefficient)) || any(coefficient <= 0 | coefficient > 1)
        fail('Drainage inlet discharge coefficients must be in the interval (0, 1].');
    end
else
    inlet_cell = zeros(0, 1);
    inlet_node = zeros(0, 1);
    crest = zeros(0, 1);
    opening = zeros(0, 1);
    coefficient = zeros(0, 1);
    reverse = false(0, 1);
end

if ismember('initial_head_m', nodes.Properties.VariableNames)
    initial_head = number_column(nodes, 'initial_head_m');
else
    initial_head = NaN(n_nodes, 1);
end
initial_head(~isfinite(initial_head)) = invert(~isfinite(initial_head));
if any(initial_head < invert - 1e-9)
    fail('Initial drainage heads cannot be below node inverts.');
end
initial_head(boundary_type ~= -1) = boundary_head(boundary_type ~= -1);

crown = max(rim, invert);
crown = max(crown, accumarray(pipe_start, start_invert + D, [n_nodes 1], @max, -Inf));
crown = max(crown, accumarray(pipe_end, end_invert + D, [n_nodes 1], @max, -Inf));
[unique_cell, ~, cell_group] = unique(inlet_cell);

N = struct();
N.grid_size = [ny nx];
N.cell_area_m2 = double(cell_area_m2);
N.node_ids = node_ids;
N.pipe_ids = pipe_ids;
N.inlet_ids = inlet_ids;
N.node_invert_m = invert;
N.node_rim_m = rim;
N.node_shaft_area_m2 = shaft;
N.node_initial_head_m = initial_head;
N.node_boundary_type = boundary_type;
N.node_boundary_head_m = boundary_head;
N.node_crown_head_m = crown;
N.dynamic = boundary_type == -1;
N.pipe_start = pipe_start;
N.pipe_end = pipe_end;
N.pipe_length_m = L;
N.pipe_half_length_m = 0.5 * L;
N.pipe_diameter_m = D;
N.pipe_manning_n = n_manning;
N.pipe_friction_factor = 9.80665 * n_manning.^2;
N.pipe_start_invert_m = start_invert;
N.pipe_end_invert_m = end_invert;
N.pipe_normal_slope = normal_slope;
N.pipe_normal_factor = sqrt(max(normal_slope, 0)) ./ n_manning;
N.normal_pipe_index = normal_index;
N.normal_pipe_at_start = boundary_type(pipe_start(normal_index)) == 2;
% Pipes that end at a free outfall, and their slope toward it.
N.free_at_start = boundary_type(pipe_start) == 0;
N.free_at_end = boundary_type(pipe_end) == 0;
free_slope = zeros(numel(L), 1);
free_slope(N.free_at_start) = (end_invert(N.free_at_start) - start_invert(N.free_at_start)) ./ L(N.free_at_start);
free_slope(N.free_at_end) = (start_invert(N.free_at_end) - end_invert(N.free_at_end)) ./ L(N.free_at_end);
N.pipe_free_slope = free_slope;
N.inlet_node = inlet_node(:);
N.inlet_cell = inlet_cell(:);
N.inlet_unique_cell = unique_cell(:);
N.inlet_cell_group = cell_group(:);
N.inlet_terrain_m = terrain_m(inlet_cell(:));
N.inlet_crest_m = crest(:);
N.inlet_opening_area_m2 = opening(:);
N.inlet_discharge_coefficient = coefficient(:);
N.inlet_reverse_allowed = reverse(:);
N.hydraulics = H;
end

function H = hp2d_urban_drainage_options(H)
defaults = struct('courant', 0.5, 'min_dt_s', 0.01, 'max_dt_s', 5.0, ...
    'maximum_substeps', 10000, 'maximum_velocity_m_s', 10.0, ...
    'free_outfall_depth', "critical_normal");
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(H, names{i}) || isempty(H.(names{i}))
        H.(names{i}) = defaults.(names{i});
    end
end
if isfield(H, 'maximum_substeps_per_surface_step') && ~isempty(H.maximum_substeps_per_surface_step)
    H.maximum_substeps = H.maximum_substeps_per_surface_step;
end
inlet = struct('opening_area_m2', 0.10, 'discharge_coefficient', 0.67, 'reverse_flow_allowed', true);
if ~isfield(H, 'inlet_defaults') || ~isstruct(H.inlet_defaults)
    H.inlet_defaults = struct();
end
names = fieldnames(inlet);
for i = 1:numel(names)
    if ~isfield(H.inlet_defaults, names{i}) || isempty(H.inlet_defaults.(names{i}))
        H.inlet_defaults.(names{i}) = inlet.(names{i});
    end
end
H.free_outfall_depth = lower(string(H.free_outfall_depth));
if ~any(H.free_outfall_depth == ["critical_normal", "invert"])
    fail('urban_drainage hydraulics free_outfall_depth must be "critical_normal" or "invert".');
end
if ~(H.courant > 0 && H.courant <= 1)
    fail('urban_drainage hydraulics courant must be in (0, 1].');
end
if ~(H.min_dt_s > 0 && H.max_dt_s > 0 && H.min_dt_s <= H.max_dt_s)
    fail('urban_drainage hydraulics need 0 < min_dt_s <= max_dt_s.');
end
if ~(H.maximum_substeps >= 1) || H.maximum_substeps ~= floor(H.maximum_substeps)
    fail('urban_drainage hydraulics maximum_substeps must be a positive integer.');
end
if ~(H.maximum_velocity_m_s > 0)
    fail('Urban drainage requires a positive maximum pipe velocity.');
end
if ~(H.inlet_defaults.opening_area_m2 > 0)
    fail('The default inlet opening area must be positive.');
end
if ~(H.inlet_defaults.discharge_coefficient > 0 && H.inlet_defaults.discharge_coefficient <= 1)
    fail('The default inlet discharge coefficient must be in (0, 1].');
end
end

function v = text_column(T, name)
if ~ismember(name, T.Properties.VariableNames)
    fail(sprintf('Urban-drainage layer is missing required field ''%s''.', name));
end
v = string(T.(name));
v = v(:);
end

function v = number_column(T, name)
if ~ismember(name, T.Properties.VariableNames)
    fail(sprintf('Urban-drainage layer is missing required field ''%s''.', name));
end
v = T.(name);
if isstring(v) || iscellstr(v) || ischar(v)
    v = str2double(string(v));
end
v = double(v(:));
end

function v = integer_column(T, name)
v = number_column(T, name);
if any(~isfinite(v)) || any(v ~= floor(v))
    fail(sprintf('Urban-drainage field ''%s'' must contain integers.', name));
end
end

function v = with_fallback(T, name, default, n)
if ismember(name, T.Properties.VariableNames)
    v = number_column(T, name);
    v(~isfinite(v)) = default;
else
    v = default * ones(n, 1);
end
end

function check_positive(v, label)
if any(~isfinite(v)) || any(v <= 0)
    fail(sprintf('Drainage pipe %s values must be finite and positive.', label));
end
end

function require_unique(v, label)
if numel(unique(v)) ~= numel(v)
    fail(sprintf('Urban-drainage %s identifiers must be unique.', label));
end
end

function fail(message)
error('HydroPol2D:UrbanDrainage', '%s', message);
end
