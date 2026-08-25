function summary = run_voronoi_channel_cfl_validation(output_directory)
%RUN_VORONOI_CHANNEL_CFL_VALIDATION Exercise hybrid channel graph CFL limits.
%
% A three-node Neal reach connects to one resolved polygon through a
% transition.  The shared middle/terminal-node signal makes the graph CFL
% limit stricter than an individual link wave-travel limit.

arguments
    output_directory (1,:) char = fullfile(fileparts(mfilename('fullpath')), 'outputs_channel_cfl')
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'HydroPol2D_Functions'));
if exist(output_directory, 'dir') ~= 7, mkdir(output_directory); end
mesh_file = fullfile(output_directory, 'channel-cfl-fixture.nc');
write_fixture(mesh_file);

config = struct('routing_solver', 'local_inertial', 'duration_s', 60, ...
    'initial_surface_depth_m', 0, 'initial_channel_depth_m', 1.5, ...
    'initial_channel_discharge_m3_s', 0, 'surface_roughness', 0.04, ...
    'min_dt_s', 1e-3, 'max_dt_s', 100, 'output_interval_s', 10, ...
    'courant', 0.6, 'critical_flow', true, 'minimum_cell_width_m', 100);
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct());
dt = [results.diagnostics.dt_s]';
initial_volume = 3 * 2000 * 1.5;
final_volume = sum(results.final_surface_volume_m3) + sum(results.final_channel_volume_m3);
mass_error = abs(final_volume - initial_volume) / initial_volume;

% At t=0, each 100 m link alone permits 15.64 s, whereas the two-link
% node plus transition signal permits 7.82 s.  This distinguishes the
% graph/storage CFL from the former shortest-link-only criterion.
single_link_limit_s = config.courant * 100 / sqrt(9.81 * 1.5);
assert(dt(1) < 0.75 * single_link_limit_s, ...
    'Channel graph CFL was not more restrictive than a single-link limit.');
assert(all(isfinite(dt) & dt >= config.min_dt_s));
assert(all(isfinite(results.final_surface_volume_m3) & results.final_surface_volume_m3 >= 0));
assert(all(isfinite(results.final_channel_volume_m3) & results.final_channel_volume_m3 >= 0));
assert(any(abs(results.channel_discharge_m3_s) > 0), 'Channel links did not transmit flow.');
assert(any(abs(results.channel_transition_discharge_m3_s) > 0), 'Transition did not transmit flow.');
assert(mass_error <= 1e-10, 'Hybrid graph mass error %.3g.', mass_error);

summary = struct('mesh_file', mesh_file, 'first_timestep_s', dt(1), ...
    'minimum_timestep_s', min(dt), 'single_link_limit_s', single_link_limit_s, ...
    'relative_mass_error', mass_error, 'final_link_discharge_m3_s', ...
    results.channel_discharge_m3_s, 'final_transition_discharge_m3_s', ...
    results.channel_transition_discharge_m3_s);
end

function write_fixture(mesh_file)
if exist(mesh_file, 'file') == 2, delete(mesh_file); end
write_vector(mesh_file, 'cell_area_m2', 'cell', 10000 * ones(4,1));
write_vector(mesh_file, 'cell_bed_elevation_m', 'cell', [13; 12; 11; 7]);
write_vector(mesh_file, 'mesh2d_face_x', 'cell', [50; 150; 250; 350]);
write_vector(mesh_file, 'mesh2d_face_y', 'cell', 50 * ones(4,1));
write_vector(mesh_file, 'cell_target_width_m', 'cell', 100 * ones(4,1));
write_vector(mesh_file, 'cell_refinement_source', 'cell', zeros(4,1));
write_vector(mesh_file, 'cell_hydraulic_roughness', 'cell', [nan; nan; nan; 0.04]);
write_vector(mesh_file, 'cell_cfl_width_m', 'cell', 100 * ones(4,1));

owner = [0; 1; 2; 0; 1; 2; 3];
neighbor = [1; 2; 3; -1; -1; -1; -1];
write_vector(mesh_file, 'edge_owner', 'edge', owner);
write_vector(mesh_file, 'edge_neighbor', 'edge', neighbor);
write_vector(mesh_file, 'edge_length_m', 'edge', 100 * ones(7,1));
write_vector(mesh_file, 'edge_center_distance_m', 'edge', 100 * ones(7,1));
write_vector(mesh_file, 'edge_midpoint_x', 'edge', [100; 200; 300; 50; 150; 250; 350]);
write_vector(mesh_file, 'edge_midpoint_y', 'edge', 50 * ones(7,1));
write_vector(mesh_file, 'edge_normal_x', 'edge', [1; 1; 1; 0; 0; 0; 0]);
write_vector(mesh_file, 'edge_normal_y', 'edge', [0; 0; 0; 1; 1; 1; 1]);
write_vector(mesh_file, 'edge_boundary_type', 'edge', zeros(7,1));

write_vector(mesh_file, 'channel_node_x', 'channel_node', [50; 150; 250]);
write_vector(mesh_file, 'channel_node_y', 'channel_node', 50 * ones(3,1));
write_vector(mesh_file, 'channel_node_host_cell', 'channel_node', [0; 1; 2]);
write_vector(mesh_file, 'channel_node_bed_elevation_m', 'channel_node', [10; 9; 8]);
write_vector(mesh_file, 'channel_node_bank_elevation_m', 'channel_node', [13; 12; 11]);
write_vector(mesh_file, 'channel_node_plan_area_m2', 'channel_node', 2000 * ones(3,1));

write_vector(mesh_file, 'channel_link_upstream_node', 'channel_link', [0; 1]);
write_vector(mesh_file, 'channel_link_downstream_node', 'channel_link', [1; 2]);
write_vector(mesh_file, 'channel_link_length_m', 'channel_link', 100 * ones(2,1));
write_vector(mesh_file, 'channel_link_width_m', 'channel_link', 20 * ones(2,1));
write_vector(mesh_file, 'channel_link_bankfull_depth_m', 'channel_link', 3 * ones(2,1));
write_vector(mesh_file, 'channel_link_roughness', 'channel_link', 0.03 * ones(2,1));
write_vector(mesh_file, 'channel_link_reach_id', 'channel_link', ones(2,1));

write_vector(mesh_file, 'channel_transition_node', 'channel_transition', 2);
write_vector(mesh_file, 'channel_transition_resolved_cell', 'channel_transition', 3);
write_vector(mesh_file, 'channel_transition_positive_node_to_cell', 'channel_transition', 1);
write_vector(mesh_file, 'channel_transition_length_m', 'channel_transition', 100);
write_vector(mesh_file, 'channel_transition_width_m', 'channel_transition', 20);
write_vector(mesh_file, 'channel_transition_bankfull_depth_m', 'channel_transition', 3);
write_vector(mesh_file, 'channel_transition_roughness', 'channel_transition', 0.03);
end

function write_vector(file, name, dimension, value)
nccreate(file, name, 'Dimensions', {dimension, numel(value)}, 'Datatype', 'double');
ncwrite(file, name, double(value(:)));
end
