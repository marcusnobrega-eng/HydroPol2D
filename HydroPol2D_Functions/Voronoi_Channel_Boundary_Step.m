function [channel_volume, boundary_q, diagnostics] = Voronoi_Channel_Boundary_Step(channel, channel_volume, boundary_q, boundary, dt, options)
%VORONOI_CHANNEL_BOUNDARY_STEP Endpoint inflow/stage/outlet conditions.

arguments
    channel struct
    channel_volume (:,1) double
    boundary_q (:,1) double
    boundary struct
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.evaluation_volume (:,1) double = channel_volume
end
if isempty(boundary) || ~isfield(boundary, 'node_id') || isempty(boundary.node_id)
    diagnostics = struct('net_inflow_volume_m3', 0, 'max_discharge_m3_s', 0); return
end
node = double(boundary.node_id(:)); types = string(boundary.type(:));
if numel(boundary_q) ~= numel(node), boundary_q = zeros(numel(node),1); end
if isscalar(types), types = repmat(types, numel(node), 1); end
values = zeros(numel(node),1); if isfield(boundary,'value'), values = expand(boundary.value,numel(node)); end
lengths = ones(numel(node),1); if isfield(boundary,'length_m'), lengths = expand(boundary.length_m,numel(node)); end
widths = sqrt(channel.plan_area(node)); if isfield(boundary,'width_m'), widths = expand(boundary.width_m,numel(node)); end
roughness = 0.035 * ones(numel(node),1); if isfield(boundary,'roughness'), roughness = expand(boundary.roughness,numel(node)); end
assert(all(node >= 1 & node <= channel.n_nodes));
assert(numel(options.evaluation_volume) == channel.n_nodes);
depth = max(options.evaluation_volume(node) ./ channel.plan_area(node), 0);
wse = channel.bed(node) + depth;
qout = zeros(numel(node),1);
for k = 1:numel(node)
    switch types(k)
        case "inflow"
            qout(k) = -values(k);
        case "stage"
            hflow = max(max(wse(k), values(k)) - channel.bed(node(k)), 0);
            area = widths(k) * hflow; radius = area / max(widths(k) + 2*hflow, eps);
            slope = (values(k) - wse(k)) / lengths(k);
            qout(k) = (boundary_q(k) - options.gravity * area * dt * slope) / ...
                (1 + options.gravity * dt * roughness(k)^2 * abs(boundary_q(k)) / max(area * radius^(4/3), eps));
        case "normal_flow"
            area = widths(k) * depth(k); radius = area / max(widths(k) + 2*depth(k), eps);
            qout(k) = area / roughness(k) * radius^(2/3) * sqrt(max(values(k),0));
        case "critical_flow"
            qout(k) = widths(k) * depth(k) * sqrt(options.gravity * depth(k));
        case "wall"
            qout(k) = 0;
        otherwise
            error('HydroPol2D:UnknownChannelBoundaryType', 'Unknown channel boundary type %s.', types(k));
    end
end
outgoing = zeros(channel.n_nodes, 1);
if any(qout > 0)
    outgoing = accumarray(node(qout > 0), qout(qout > 0) .* dt, [channel.n_nodes 1], @sum, 0);
end
scale = min(1, channel_volume ./ max(outgoing, eps)); qout(qout > 0) = qout(qout > 0) .* scale(node(qout > 0));
channel_volume = max(channel_volume + accumarray(node, -qout .* dt, [channel.n_nodes 1], @sum, 0), 0);
boundary_q = qout;
diagnostics.net_inflow_volume_m3 = -sum(qout) * dt;
diagnostics.max_discharge_m3_s = max(abs(qout), [], 'omitnan');
end

function value = expand(value, count)
if isscalar(value), value = repmat(double(value), count, 1); else, value = double(value(:)); end
assert(numel(value) == count && all(isfinite(value)));
end
