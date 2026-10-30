function R = hp2d_urban_drainage_step(N, node_volume_m3, pipe_q_m3_s, surface_m3, node_head_m, dt_s)
%HP2D_URBAN_DRAINAGE_STEP One conservative inlet-and-network update
%   (urban_drainage_step in HydroPol2D-Python, same order and limits).
%
%   surface_m3: water volume of each inlet cell, N.inlet_unique_cell order.
%   node_head_m: node heads at the start (empty: solved from the volumes).
%
%   1. Inlets: two-way orifices Q = Cd A sqrt(2 g dh). Capture when the
%      surface head exceeds max(sewer head, crest); surcharge back to the
%      surface when the sewer head exceeds max(surface head, crest) and the
%      inlet allows reverse flow. Capture is capped at the cell's water,
%      surcharge at the node's water.
%   2. Pipes: local inertia with semi-implicit Manning friction,
%      Q = (Q0 - g A dt dH/L) / (1 + g n^2 dt |Q0| / (A R^(4/3))),
%      A and R at the mean end depth. A pipe at a normal-depth outfall
%      carries Manning normal flow at its other end's depth. A free outfall
%      takes no inflow back. A dry donor end (<= 1e-6 m) stops the pipe, and
%      |Q| <= A vmax. Outflow from a node is capped at its water.
%   3. Water entering or leaving boundary nodes is outfall flow.
%
%   R: surface_m3, node_volume_m3, pipe_q_m3_s, inlet_q_m3_s (positive into
%   the network), node_head_m (after the update), outfall_outflow_m3,
%   outfall_inflow_m3, and node_outflow_m3 / node_inflow_m3 (the same per
%   boundary node).

g = 9.80665;
tiny = 1e-12;
dt = dt_s;
dynamic = N.dynamic;
n_nodes = numel(N.node_ids);
node_volume = max(node_volume_m3(:), 0);
if isempty(node_head_m)
    heads = hp2d_urban_drainage_head(N, node_volume);
else
    heads = node_head_m(:);
end

% 1. Inlet exchange.
node = N.inlet_node;
group = N.inlet_cell_group;
n_cells = numel(N.inlet_unique_cell);
surface_at_inlets = max(surface_m3(:), 0);
inlet_q = zeros(numel(node), 1);
if ~isempty(node)
    surface_head = N.inlet_terrain_m + surface_at_inlets(group) / N.cell_area_m2;
    sewer_head = heads(node);
    crest = N.inlet_crest_m;
    capacity = N.inlet_discharge_coefficient .* N.inlet_opening_area_m2;
    capture_head = max(surface_head - max(sewer_head, crest), 0);
    surcharge_head = max(sewer_head - max(surface_head, crest), 0);
    capture = capacity .* sqrt(2 * g * capture_head);
    surcharge = capacity .* sqrt(2 * g * surcharge_head);
    surcharge(~N.inlet_reverse_allowed) = 0;

    cell_demand = accumarray(group, capture * dt, [n_cells 1]);
    cell_scale = ones(n_cells, 1);
    asked = cell_demand > 0;
    cell_scale(asked) = min(1, surface_at_inlets(asked) ./ max(cell_demand(asked), tiny));
    capture = capture .* cell_scale(group);
    node_demand = accumarray(node, surcharge * dt, [n_nodes 1]);
    node_scale = ones(n_nodes, 1);
    asked = node_demand > 0;
    node_scale(asked) = min(1, node_volume(asked) ./ max(node_demand(asked), tiny));
    scale = ones(numel(node), 1);
    scale(dynamic(node)) = node_scale(node(dynamic(node)));
    surcharge = surcharge .* scale;
    inlet_q = capture - surcharge;

    cell_exchange = accumarray(group, -inlet_q * dt, [n_cells 1]);
    surface_at_inlets = max(surface_at_inlets + cell_exchange, 0);
    node_gain = inlet_q * dt;
    node_gain(~dynamic(node)) = 0;
    node_volume = node_volume + accumarray(node, node_gain, [n_nodes 1]);
    node_volume = max(node_volume, 0);
end
heads = hp2d_urban_drainage_head(N, node_volume, heads);

% 2. Pipe momentum.
start = N.pipe_start;
finish = N.pipe_end;
D = N.pipe_diameter_m;
previous_q = pipe_q_m3_s(:);
[head_start, head_end, start_depth, end_depth] = hp2d_urban_drainage_pipe_ends(N, heads, previous_q);
[area, radius] = hp2d_urban_drainage_section(0.5 * (start_depth + end_depth), D);
numerator = previous_q - g * area * dt .* (head_end - head_start) ./ N.pipe_length_m;
friction = N.pipe_friction_factor * dt .* abs(previous_q) ...
    ./ max(area .* radius.^(4/3), tiny);
pipe_q = zeros(numel(D), 1);
wet = area > tiny;
pipe_q(wet) = numerator(wet) ./ (1 + friction(wet));

normal = N.normal_pipe_index;
flow_area = area;
if ~isempty(normal)
    at_start = N.normal_pipe_at_start;
    normal_depth = start_depth(normal);
    normal_depth(at_start) = end_depth(normal(at_start));
    [normal_area, normal_radius] = hp2d_urban_drainage_section(normal_depth, D(normal));
    magnitude = normal_area .* normal_radius.^(2/3) .* N.pipe_normal_factor(normal);
    magnitude(at_start) = -magnitude(at_start);
    pipe_q(normal) = magnitude;
    flow_area(normal) = normal_area;
end
pipe_q(N.free_at_start) = min(pipe_q(N.free_at_start), 0);
pipe_q(N.free_at_end) = max(pipe_q(N.free_at_end), 0);
donor_depth = end_depth;
forward = pipe_q >= 0;
donor_depth(forward) = start_depth(forward);
pipe_q(~(donor_depth > 1e-6)) = 0;
capacity = flow_area * N.hydraulics.maximum_velocity_m_s;
pipe_q = min(max(pipe_q, -capacity), capacity);

donor = finish;
forward = pipe_q >= 0;
donor(forward) = start(forward);
donor_dynamic = dynamic(donor);
demand = abs(pipe_q) * dt;
demand(~donor_dynamic) = 0;
pipe_demand = accumarray(donor, demand, [n_nodes 1]);
pipe_scale = ones(n_nodes, 1);
asked = pipe_demand > 0;
pipe_scale(asked) = min(1, node_volume(asked) ./ max(pipe_demand(asked), tiny));
scale = ones(numel(pipe_q), 1);
scale(donor_dynamic) = pipe_scale(donor(donor_dynamic));
pipe_q = pipe_q .* scale;
transfer = pipe_q * dt;
delta = accumarray(finish, transfer, [n_nodes 1]) - accumarray(start, transfer, [n_nodes 1]);
delta(~dynamic) = 0;
node_volume = max(node_volume + delta, 0);

% 3. Outfall exchange.
boundary_pipe_rate = zeros(numel(pipe_q), 1);
boundary_pipe_rate(~dynamic(start)) = -pipe_q(~dynamic(start));
boundary_pipe_rate(~dynamic(finish)) = boundary_pipe_rate(~dynamic(finish)) + pipe_q(~dynamic(finish));
boundary_inlet_rate = zeros(numel(inlet_q), 1);
if ~isempty(node)
    boundary_inlet_rate(~dynamic(node)) = inlet_q(~dynamic(node));
end
R = struct();
R.surface_m3 = surface_at_inlets;
R.node_volume_m3 = node_volume;
R.pipe_q_m3_s = pipe_q;
R.inlet_q_m3_s = inlet_q;
R.node_head_m = hp2d_urban_drainage_head(N, node_volume, heads);
R.outfall_outflow_m3 = (sum(max(boundary_pipe_rate, 0)) + sum(max(boundary_inlet_rate, 0))) * dt;
R.outfall_inflow_m3 = (sum(max(-boundary_pipe_rate, 0)) + sum(max(-boundary_inlet_rate, 0))) * dt;
% The same volumes per outfall node, to book them where they leave.
boundary_node = start;
boundary_node(~dynamic(finish)) = finish(~dynamic(finish));
on_boundary = ~dynamic(start) | ~dynamic(finish);
R.node_outflow_m3 = accumarray(boundary_node(on_boundary), max(boundary_pipe_rate(on_boundary), 0) * dt, [n_nodes 1]);
R.node_inflow_m3 = accumarray(boundary_node(on_boundary), max(-boundary_pipe_rate(on_boundary), 0) * dt, [n_nodes 1]);
if ~isempty(node)
    R.node_outflow_m3 = R.node_outflow_m3 + accumarray(node, max(boundary_inlet_rate, 0) * dt, [n_nodes 1]);
    R.node_inflow_m3 = R.node_inflow_m3 + accumarray(node, max(-boundary_inlet_rate, 0) * dt, [n_nodes 1]);
end
end
