function [surface_volume, edge_q, diagnostics] = Voronoi_Kinematic_Step(mesh, surface_volume, roughness, dt, options)
%VORONOI_KINEMATIC_STEP Conservative explicit kinematic-wave update.

arguments
    mesh struct
    surface_volume (:,1) double
    roughness double
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.critical_flow (1,1) logical = false
end

[edge_q, face_depth] = Voronoi_Kinematic_Face_Flux(mesh, surface_volume, roughness, ...
    gravity=options.gravity, dry_tolerance_m=options.dry_tolerance_m, critical_flow=options.critical_flow);
internal = mesh.edge_neighbor > 0;
owner = mesh.edge_owner(internal); neighbor = mesh.edge_neighbor(internal);
Q = edge_q(internal) .* mesh.edge_length(internal);
out_cell = owner; out_cell(Q < 0) = neighbor(Q < 0);
requested = accumarray(out_cell, abs(Q) .* dt, [mesh.n_cells 1], @sum, 0);
scale = min(1, surface_volume(:) ./ max(requested, eps));
Q = Q .* scale(out_cell);
edge_q(internal) = Q ./ mesh.edge_length(internal);
delta = accumarray(owner, -Q .* dt, [mesh.n_cells 1], @sum, 0) + ...
    accumarray(neighbor, Q .* dt, [mesh.n_cells 1], @sum, 0);
surface_volume = surface_volume(:) + delta;
if any(surface_volume < -1e-10)
    error('HydroPol2D:NegativeSurfaceVolume', 'Kinematic surface draining limiter failed.');
end
surface_volume = max(surface_volume, 0);
velocity = abs(edge_q(internal)) ./ max(face_depth(internal), options.dry_tolerance_m);
diagnostics.max_depth_m = max(surface_volume ./ mesh.surface_area(:));
diagnostics.max_velocity_m_s = max(velocity, [], 'omitnan');
diagnostics.internal_flux_volume_m3 = sum(abs(Q)) * dt;
diagnostics.mass_change_m3 = sum(delta);
end
