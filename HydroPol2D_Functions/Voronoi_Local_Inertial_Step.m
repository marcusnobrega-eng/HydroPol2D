function [surface_volume, edge_q, diagnostics] = Voronoi_Local_Inertial_Step(mesh, surface_volume, edge_q, roughness, dt, options)
%VORONOI_LOCAL_INERTIAL_STEP Conservative edge-based local-inertial update.

arguments
    mesh struct
    surface_volume (:,1) double
    edge_q (:,1) double
    roughness double
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.critical_flow (1,1) logical = false
end

n = mesh.n_cells;
if isscalar(roughness), roughness = repmat(roughness, n, 1); end
assert(numel(surface_volume) == n && numel(roughness) == n);
owner = mesh.edge_owner(:); neighbor = mesh.edge_neighbor(:);
internal = neighbor > 0;
o = owner(internal); d = neighbor(internal);
depth = max(surface_volume(:) ./ mesh.surface_area(:), 0);
wse = mesh.surface_bed(:) + depth;
hflow = max(max(wse(o), wse(d)) - max(mesh.surface_bed(o), mesh.surface_bed(d)), 0);
slope = (wse(d) - wse(o)) ./ mesh.edge_distance(internal);
nface = 0.5 .* (roughness(o) + roughness(d));
qold = edge_q(internal);
wet = hflow > options.dry_tolerance_m;
qnew = zeros(size(qold));
denominator = 1 + options.gravity .* dt .* nface(wet).^2 .* abs(qold(wet)) ./ max(hflow(wet).^(7/3), eps);
qnew(wet) = (qold(wet) - options.gravity .* hflow(wet) .* dt .* slope(wet)) ./ denominator;
if options.critical_flow
    qcrit = hflow .* sqrt(options.gravity .* hflow);
    qnew = sign(qnew) .* min(abs(qnew), qcrit);
end

Q = qnew .* mesh.edge_length(internal);
out_cell = o; out_cell(Q < 0) = d(Q < 0);
requested = accumarray(out_cell, abs(Q) .* dt, [n 1], @sum, 0);
scale = min(1, surface_volume(:) ./ max(requested, eps));
Q = Q .* scale(out_cell);
qnew = Q ./ mesh.edge_length(internal);

delta = accumarray(o, -Q .* dt, [n 1], @sum, 0) + accumarray(d, Q .* dt, [n 1], @sum, 0);
surface_volume = surface_volume(:) + delta;
surface_volume(abs(surface_volume) < 10 * eps(max(1, max(surface_volume)))) = 0;
if any(surface_volume < -1e-10)
    error('HydroPol2D:NegativeSurfaceVolume', 'Surface draining limiter failed.');
end
surface_volume = max(surface_volume, 0);
edge_q(internal) = qnew;
edge_q(~internal) = 0; % Boundary fluxes are supplied separately by the runner.
% The face flow depth PAIRED with the discharge stored in edge_q. Diagnostics must
% divide by this rather than recomputing a depth from the post-step state; see
% HydroPol2D_Voronoi_Cell_Velocity.
diagnostics.face_flow_depth_m = zeros(mesh.n_edges,1);
diagnostics.face_flow_depth_m(internal) = hflow;
diagnostics.max_depth_m = max(depth);
diagnostics.max_velocity_m_s = max(abs(qnew) ./ max(hflow, options.dry_tolerance_m), [], 'omitnan');
diagnostics.internal_flux_volume_m3 = sum(abs(Q)) * dt;
diagnostics.mass_change_m3 = sum(delta);
end
