function channel_boundary = hp2d_voronoi_auto_channel_boundary(mesh, surface_boundary)
%HP2D_VORONOI_AUTO_CHANNEL_BOUNDARY Give channel nodes their cell's surface outlet.
%   HydroPol2D-Python applies the surface stage, normal-flow and
%   critical-flow edges to the Neal channel node of each boundary cell
%   (voronoi.py, _channel_boundary_step), so a river reaching the domain edge
%   drains through it. MATLAB needed an explicit forcing.channel_boundary.
%   This builds one from surface_boundary (edge_id, type, value; a stage
%   given as a depth above the invert carries stage_depth_m and is measured
%   from the channel bed). As in Python, the node takes the width of its
%   widest link (at least sqrt(plan area)), the mean roughness of its links
%   (0.035 without links) and a length of max(sqrt(plan area), 1 m). A later
%   edge of the same cell overrides an earlier one.

channel_boundary = struct('node_id', zeros(0, 1), 'type', strings(0, 1), 'value', zeros(0, 1));
if mesh.channel.n_nodes == 0 || ~isstruct(surface_boundary) || ...
        ~isfield(surface_boundary, 'edge_id') || isempty(surface_boundary.edge_id)
    return
end
edge_id = double(surface_boundary.edge_id(:));
types = string(surface_boundary.type(:));
if isscalar(types), types = repmat(types, numel(edge_id), 1); end
values = zeros(numel(edge_id), 1);
if isfield(surface_boundary, 'value')
    values = double(surface_boundary.value(:));
    if isscalar(values), values = repmat(values, numel(edge_id), 1); end
end
stage_depth = nan(numel(edge_id), 1);
if isfield(surface_boundary, 'stage_depth_m'), stage_depth = double(surface_boundary.stage_depth_m(:)); end

node_of_cell = zeros(mesh.n_cells, 1);
node_of_cell(mesh.channel.host_cell) = 1:mesh.channel.n_nodes;
control_type = strings(mesh.channel.n_nodes, 1);
control_value = zeros(mesh.channel.n_nodes, 1);
for k = 1:numel(edge_id)
    if ~any(types(k) == ["stage", "normal_flow", "critical_flow"]), continue; end
    node = node_of_cell(mesh.edge_owner(edge_id(k)));
    if node == 0, continue; end
    value = values(k);
    if types(k) == "stage" && isfinite(stage_depth(k))
        value = mesh.channel.bed(node) + stage_depth(k);
    end
    control_type(node) = types(k);
    control_value(node) = value;
end
nodes = find(control_type ~= "");
if isempty(nodes), return; end
plan = mesh.channel.plan_area(nodes);
width = sqrt(plan); roughness = 0.035 * ones(numel(nodes), 1);
for k = 1:numel(nodes)
    if mesh.channel.n_links > 0
        touching = mesh.channel.link_up == nodes(k) | mesh.channel.link_down == nodes(k);
        if any(touching)
            width(k) = max(width(k), max(mesh.channel.link_width(touching)));
            roughness(k) = mean(mesh.channel.link_roughness(touching));
        end
    end
end
channel_boundary = struct('node_id', nodes, 'type', control_type(nodes), ...
    'value', control_value(nodes), 'width_m', width, 'roughness', roughness, ...
    'length_m', max(sqrt(plan), 1));
end
