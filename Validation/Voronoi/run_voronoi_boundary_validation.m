function summary = run_voronoi_boundary_validation(surface_mesh_file, channel_mesh_file)
%RUN_VORONOI_BOUNDARY_VALIDATION Prescribed inflow volume checks.

surface_mesh = HydroPol2D_Read_UGRID(surface_mesh_file);
edge_id = find(surface_mesh.edge_neighbor == 0, 1, 'first');
duration = 100; inflow = 2;
config = struct('duration_s', duration, 'min_dt_s', 1e-4, 'max_dt_s', 2, 'output_interval_s', 20);
forcing.surface_boundary = struct('edge_id', edge_id, 'type', "inflow", 'value', inflow);
surface_results = HydroPol2D_Voronoi_Run(surface_mesh_file, config, forcing);
surface_error = abs(sum(surface_results.final_surface_volume_m3) - inflow * duration) / (inflow * duration);
assert(surface_error <= 1e-10, 'Surface inflow volume error %.3g.', surface_error);

channel_mesh = HydroPol2D_Read_UGRID(channel_mesh_file);
assert(channel_mesh.channel.n_nodes > 0);
channel_inflow = 1;
forcing = struct('channel_boundary', struct('node_id', 1, 'type', "inflow", 'value', channel_inflow));
channel_results = HydroPol2D_Voronoi_Run(channel_mesh_file, config, forcing);
channel_total = sum(channel_results.final_surface_volume_m3) + sum(channel_results.final_channel_volume_m3);
channel_error = abs(channel_total - channel_inflow * duration) / (channel_inflow * duration);
assert(channel_error <= 1e-10, 'Channel inflow volume error %.3g.', channel_error);
summary = struct('surface_relative_mass_error', surface_error, ...
    'channel_relative_mass_error', channel_error);
end
