function [surface_volume, channel_volume, transition_q, diagnostics] = Voronoi_Channel_Transition_Step(mesh, surface_volume, channel_volume, transition_q, dt, options)
%VORONOI_CHANNEL_TRANSITION_STEP Conservative subgrid/resolved connection.

arguments
    mesh struct
    surface_volume (:,1) double
    channel_volume (:,1) double
    transition_q (:,1) double
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.critical_flow (1,1) logical = false
end

channel = mesh.channel;
if channel.n_transitions == 0
    diagnostics = struct('max_velocity_m_s', 0, 'mass_change_m3', 0); return
end
node = channel.transition_node(:); cell_id = channel.transition_cell(:);
channel_depth = max(channel_volume(:) ./ channel.plan_area(:), 0);
surface_depth = max(surface_volume(:) ./ mesh.surface_area(:), 0);
node_wse = channel.bed(node) + channel_depth(node);
cell_wse = mesh.surface_bed(cell_id) + surface_depth(cell_id);
hflow = max(max(node_wse, cell_wse) - max(channel.bed(node), mesh.surface_bed(cell_id)), 0);
area = channel.transition_width(:) .* hflow;
radius = area ./ max(channel.transition_width(:) + 2 .* hflow, eps);
positive_node_to_cell = logical(channel.transition_positive_node_to_cell(:));
slope = (cell_wse - node_wse) ./ channel.transition_length(:);
slope(~positive_node_to_cell) = -slope(~positive_node_to_cell);
wet = hflow > options.dry_tolerance_m;
qnew = zeros(size(transition_q));
denominator = 1 + options.gravity .* dt .* channel.transition_roughness(wet).^2 .* abs(transition_q(wet)) ./ max(area(wet) .* radius(wet).^(4/3), eps);
qnew(wet) = (transition_q(wet) - options.gravity .* area(wet) .* dt .* slope(wet)) ./ denominator;
if options.critical_flow
    qcrit = area .* sqrt(options.gravity .* hflow);
    qnew = sign(qnew) .* min(abs(qnew), qcrit);
end

% Convert the link-oriented discharge to a physical node-to-cell sign for
% conservative donor limiting and storage updates.
node_to_cell_q = qnew;
node_to_cell_q(~positive_node_to_cell) = -node_to_cell_q(~positive_node_to_cell);
node_requested = zeros(channel.n_nodes, 1);
cell_requested = zeros(mesh.n_cells, 1);
if any(node_to_cell_q > 0)
    node_requested = accumarray(node(node_to_cell_q > 0), node_to_cell_q(node_to_cell_q > 0) .* dt, [channel.n_nodes 1], @sum, 0);
end
if any(node_to_cell_q < 0)
    cell_requested = accumarray(cell_id(node_to_cell_q < 0), -node_to_cell_q(node_to_cell_q < 0) .* dt, [mesh.n_cells 1], @sum, 0);
end
node_scale = min(1, channel_volume ./ max(node_requested, eps));
cell_scale = min(1, surface_volume ./ max(cell_requested, eps));
from_node = node_to_cell_q >= 0;
qnew(from_node) = qnew(from_node) .* node_scale(node(from_node));
qnew(~from_node) = qnew(~from_node) .* cell_scale(cell_id(~from_node));
node_to_cell_q = qnew;
node_to_cell_q(~positive_node_to_cell) = -node_to_cell_q(~positive_node_to_cell);

volume = node_to_cell_q .* dt;
channel_volume = channel_volume + accumarray(node, -volume, [channel.n_nodes 1], @sum, 0);
surface_volume = surface_volume + accumarray(cell_id, volume, [mesh.n_cells 1], @sum, 0);
if any(channel_volume < -1e-10) || any(surface_volume < -1e-10)
    error('HydroPol2D:TransitionDrainingLimiter', 'Transition flux over-drew storage.');
end
channel_volume = max(channel_volume, 0); surface_volume = max(surface_volume, 0); transition_q = qnew;
diagnostics.max_velocity_m_s = max(abs(qnew) ./ max(area, eps), [], 'omitnan');
diagnostics.mass_change_m3 = sum(volume) - sum(volume);
end
