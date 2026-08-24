function [channel_volume, link_q, diagnostics] = Voronoi_Neal_Channel_Step(channel, channel_volume, link_q, dt, options)
%VORONOI_NEAL_CHANNEL_STEP Local-inertial fluxes on an arbitrary river graph.

arguments
    channel struct
    channel_volume (:,1) double
    link_q (:,1) double
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.critical_flow (1,1) logical = false
end

if channel.n_nodes == 0 || channel.n_links == 0
    diagnostics = struct('max_depth_m', 0, 'max_velocity_m_s', 0, 'mass_change_m3', 0);
    return
end
depth = max(channel_volume(:) ./ channel.plan_area(:), 0);
wse = channel.bed(:) + depth;
u = channel.link_up(:); d = channel.link_down(:);
hflow = max(max(wse(u), wse(d)) - max(channel.bed(u), channel.bed(d)), 0);
width = channel.link_width(:);
area = width .* hflow;
radius = area ./ max(width + 2 .* hflow, eps);
slope = (wse(d) - wse(u)) ./ channel.link_length(:);
qold = link_q(:);
wet = hflow > options.dry_tolerance_m;
qnew = zeros(size(qold));
denominator = 1 + options.gravity .* dt .* channel.link_roughness(wet).^2 .* abs(qold(wet)) ./ max(area(wet) .* radius(wet).^(4/3), eps);
qnew(wet) = (qold(wet) - options.gravity .* area(wet) .* dt .* slope(wet)) ./ denominator;
if options.critical_flow
    qcrit = area .* sqrt(options.gravity .* hflow);
    qnew = sign(qnew) .* min(abs(qnew), qcrit);
end
out_node = u; out_node(qnew < 0) = d(qnew < 0);
requested = accumarray(out_node, abs(qnew) .* dt, [channel.n_nodes 1], @sum, 0);
scale = min(1, channel_volume(:) ./ max(requested, eps));
qnew = qnew .* scale(out_node);
delta = accumarray(u, -qnew .* dt, [channel.n_nodes 1], @sum, 0) + accumarray(d, qnew .* dt, [channel.n_nodes 1], @sum, 0);
channel_volume = channel_volume(:) + delta;
if any(channel_volume < -1e-10)
    error('HydroPol2D:NegativeChannelVolume', 'Channel draining limiter failed.');
end
channel_volume = max(channel_volume, 0); link_q = qnew;
diagnostics.max_depth_m = max(depth);
diagnostics.max_velocity_m_s = max(abs(qnew) ./ max(area, eps), [], 'omitnan');
diagnostics.mass_change_m3 = sum(delta);
end
