function [surface_volume, edge_q, diagnostics] = Voronoi_Diffusive_Step(mesh, surface_volume, roughness, dt, options)
%VORONOI_DIFFUSIVE_STEP Conservative explicit diffusive-wave update.

arguments
    mesh struct
    surface_volume (:,1) double
    roughness double
    dt (1,1) double {mustBePositive}
    options.dry_tolerance_m (1,1) double = 1e-6
    options.slope_regularization (1,1) double {mustBePositive} = 1e-4
end

[edge_q, face_depth] = Voronoi_Diffusive_Face_Flux(mesh, surface_volume, roughness, ...
    dry_tolerance_m=options.dry_tolerance_m, slope_regularization=options.slope_regularization);
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
    error('HydroPol2D:NegativeSurfaceVolume', 'Diffusive surface draining limiter failed.');
end
surface_volume = max(surface_volume, 0);
velocity = abs(edge_q(internal)) ./ max(face_depth(internal), options.dry_tolerance_m);
diagnostics.max_depth_m = max(surface_volume ./ mesh.surface_area(:));
diagnostics.max_velocity_m_s = max(velocity, [], 'omitnan');
diagnostics.internal_flux_volume_m3 = sum(abs(Q)) * dt;
diagnostics.mass_change_m3 = sum(delta);
end
