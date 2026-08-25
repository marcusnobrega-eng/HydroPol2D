function summary = run_voronoi_channel_uniform_flow_validation(mesh_file)
%RUN_VORONOI_CHANNEL_UNIFORM_FLOW_VALIDATION Manning steady-state channel.

mesh = HydroPol2D_Read_UGRID(mesh_file);
channel = mesh.channel;
assert(channel.n_nodes > 2 && channel.n_links > 1 && channel.n_transitions == 0);
inlet = setdiff((1:channel.n_nodes)', channel.link_down);
outlet = setdiff((1:channel.n_nodes)', channel.link_up);
assert(numel(inlet) == 1 && numel(outlet) == 1);
depth = 1.5; width = channel.link_width(1); roughness = channel.link_roughness(1);
slope = mean((channel.bed(channel.link_up) - channel.bed(channel.link_down)) ./ channel.link_length);
area = width * depth; radius = area / (width + 2 * depth);
Q = area / roughness * radius^(2/3) * sqrt(slope);
config = struct('duration_s', 300, 'initial_channel_depth_m', depth, ...
    'initial_channel_discharge_m3_s', Q, 'min_dt_s', 1e-4, 'max_dt_s', 1, ...
    'output_interval_s', 30, 'surface_roughness', 0.04);
forcing.channel_boundary = struct('node_id', [inlet; outlet], ...
    'type', ["inflow"; "normal_flow"], 'value', [Q; slope], ...
    'width_m', [width; width], 'roughness', [roughness; roughness]);
results = HydroPol2D_Voronoi_Run(mesh_file, config, forcing);
q_error = max(abs(results.channel_discharge_m3_s - Q)) / Q;
depth_error = max(abs(results.final_channel_volume_m3 ./ channel.plan_area - depth));
initial = sum(channel.plan_area) * depth;
final = sum(results.final_channel_volume_m3) + sum(results.final_surface_volume_m3);
mass_error = abs(final - initial) / initial;
assert(q_error <= 1e-10, 'Uniform-flow discharge error %.3g.', q_error);
assert(depth_error <= 1e-10, 'Uniform-flow depth error %.3g.', depth_error);
assert(mass_error <= 1e-10, 'Uniform-flow mass error %.3g.', mass_error);
summary = struct('analytical_discharge_m3_s', Q, 'relative_discharge_error', q_error, ...
    'maximum_depth_error_m', depth_error, 'relative_mass_error', mass_error);
end
