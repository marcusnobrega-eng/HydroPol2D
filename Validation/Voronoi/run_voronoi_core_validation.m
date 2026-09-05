function summary = run_voronoi_core_validation(mesh_file)
%RUN_VORONOI_CORE_VALIDATION Conservation smoke test for the new runner.

arguments
    mesh_file (1,:) char
end
config = struct('duration_s', 600, 'initial_surface_depth_m', 0.05, ...
    'initial_channel_depth_m', 0.5, 'min_dt_s', 1e-4, 'max_dt_s', 2, ...
    'output_interval_s', 60, 'surface_roughness', 0.04);
mesh = HydroPol2D_Read_UGRID(mesh_file, allow_legacy=true);
config.allow_legacy_mesh = true;
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct());
initial = sum(mesh.surface_area) * 0.05;
if mesh.channel.n_nodes > 0
    initial = initial + sum(mesh.channel.plan_area) * 0.5;
end
final = sum(results.final_surface_volume_m3) + sum(results.final_channel_volume_m3);
relative_error = abs(final - initial) / max(initial, eps);
assert(relative_error <= 1e-10, 'Voronoi mass residual %.3g exceeds tolerance.', relative_error);
assert(all(results.final_surface_volume_m3 >= 0));
assert(all(results.final_channel_volume_m3 >= 0));
summary = struct('relative_mass_error', relative_error, 'steps', numel(results.diagnostics), ...
    'max_depth_m', max(results.surface_depth_m, [], 'all'));
end
