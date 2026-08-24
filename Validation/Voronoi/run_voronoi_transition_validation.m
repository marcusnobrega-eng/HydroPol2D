function summary = run_voronoi_transition_validation(mesh_file)
%RUN_VORONOI_TRANSITION_VALIDATION Subgrid-to-resolved mass and flow test.

mesh = HydroPol2D_Read_UGRID(mesh_file);
assert(mesh.channel.n_transitions >= 1, 'Fixture must contain at least one transition face.');
config = struct('duration_s', 120, 'initial_surface_depth_m', 0, ...
    'initial_channel_depth_m', 1.5, 'min_dt_s', 1e-4, 'max_dt_s', 0.5, ...
    'output_interval_s', 10, 'surface_roughness', 0.04, 'critical_flow', true);
initial = sum(mesh.channel.plan_area) * 1.5;
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct());
final = sum(results.final_surface_volume_m3) + sum(results.final_channel_volume_m3);
relative_error = abs(final - initial) / max(initial, eps);
assert(relative_error <= 1e-10, 'Transition mass error %.3g.', relative_error);
assert(any(abs(results.channel_transition_discharge_m3_s) > 0), 'Transition did not transmit water.');
summary = struct('relative_mass_error', relative_error, ...
    'transition_discharge_m3_s', results.channel_transition_discharge_m3_s, ...
    'resolved_surface_volume_m3', sum(results.final_surface_volume_m3));
end
