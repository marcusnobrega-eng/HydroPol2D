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
overlap_file = fullfile(output_directory, 'channel-cfl-overlap.nc');
write_fixture(mesh_file);
write_overlap_fixture(overlap_file);

config = struct('routing_solver', 'local_inertial', 'duration_s', 60, ...
    'initial_surface_depth_m', 0, 'initial_channel_depth_m', 1.5, ...
    'initial_channel_discharge_m3_s', 0, 'surface_roughness', 0.04, ...
    'min_dt_s', 1e-3, 'max_dt_s', 100, 'output_interval_s', 10, ...
    'courant', 0.2, 'critical_flow', true, 'minimum_cell_width_m', 100, ...
    'allow_legacy_mesh', true);
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
    'overlap_file', overlap_file, ...
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
write_contract_metadata(mesh_file);
end

function write_contract_metadata(mesh_file)
crs = 'LOCAL_CS["Channel CFL fixture",UNIT["metre",1],AXIS["Easting",EAST],AXIS["Northing",NORTH]]';
ncwriteatt(mesh_file,'/','Conventions','CF-1.13, UGRID-1.0');
ncwriteatt(mesh_file,'/','title','HydroPol2D strict channel-CFL validation fixture');
ncwriteatt(mesh_file,'/','mesh_contract_version','1.0');
ncwriteatt(mesh_file,'/','schema_version','hydrobathydem-mesh-1.0');
ncwriteatt(mesh_file,'/','product_type','hydraulic_mesh');
ncwriteatt(mesh_file,'/','file_index_base',int32(0));
ncwriteatt(mesh_file,'/','boundary_neighbor_sentinel',int32(-1));
ncwriteatt(mesh_file,'/','coordinate_units','m');
ncwriteatt(mesh_file,'/','edge_normal_convention','unit normal points outward from edge_owner');
ncwriteatt(mesh_file,'/','crs_wkt',crs);

nccreate(mesh_file,'mesh2d','Datatype','int32');
ncwrite(mesh_file,'mesh2d',int32(0));
ncwriteatt(mesh_file,'mesh2d','cf_role','mesh_topology');
ncwriteatt(mesh_file,'mesh2d','topology_dimension',int32(2));
ncwriteatt(mesh_file,'mesh2d','node_coordinates','mesh2d_node_x mesh2d_node_y');
ncwriteatt(mesh_file,'mesh2d','face_node_connectivity','mesh2d_face_nodes');
ncwriteatt(mesh_file,'mesh2d','face_coordinates','mesh2d_face_x mesh2d_face_y');

node_x=repelem((0:100:400)',2); node_y=repmat([0;100],5,1);
write_vector(mesh_file,'mesh2d_node_x','node',node_x);
write_vector(mesh_file,'mesh2d_node_y','node',node_y);
nccreate(mesh_file,'mesh2d_face_nodes','Dimensions',{'cell',4,'max_face_nodes',4}, ...
    'Datatype','int32','FillValue',int32(-1));
connectivity=int32([0 2 4 6; 2 4 6 8; 3 5 7 9; 1 3 5 7]);
ncwrite(mesh_file,'mesh2d_face_nodes',connectivity);
ncwriteatt(mesh_file,'mesh2d_face_nodes','start_index',int32(0));

for item = {'cell_area_m2','m2'; 'cell_bed_elevation_m','m'; ...
        'cell_cfl_width_m','m'; 'mesh2d_node_x','m'; 'mesh2d_node_y','m'; ...
        'mesh2d_face_x','m'; 'mesh2d_face_y','m'; 'edge_length_m','m'; ...
        'edge_center_distance_m','m'; 'edge_midpoint_x','m'; 'edge_midpoint_y','m'}'
    ncwriteatt(mesh_file,item{1},'units',item{2});
end
end

function write_overlap_fixture(overlap_file)
if exist(overlap_file,'file')==2, delete(overlap_file); end
write_vector(overlap_file,'overlap_mesh_index','overlap',int32((0:3)'));
write_vector(overlap_file,'overlap_raster_index','overlap',int32((0:3)'));
write_vector(overlap_file,'overlap_area_m2','overlap',10000*ones(4,1));
write_vector(overlap_file,'mesh_area_m2','mesh_cell',10000*ones(4,1));
write_vector(overlap_file,'raster_area_m2','raster_cell',10000*ones(4,1));
write_vector(overlap_file,'x_edges','x_edge',(0:100:400)');
write_vector(overlap_file,'y_edges','y_edge',[0;100]);
crs = 'LOCAL_CS["Channel CFL fixture",UNIT["metre",1],AXIS["Easting",EAST],AXIS["Northing",NORTH]]';
ncwriteatt(overlap_file,'/','Conventions','CF-1.13');
ncwriteatt(overlap_file,'/','title','HydroPol2D strict channel-CFL overlap fixture');
ncwriteatt(overlap_file,'/','mesh_contract_version','1.0');
ncwriteatt(overlap_file,'/','schema_version','hydrobathydem-mesh-1.0');
ncwriteatt(overlap_file,'/','product_type','conservative_raster_overlap');
ncwriteatt(overlap_file,'/','file_index_base',int32(0));
ncwriteatt(overlap_file,'/','raster_index_order','south_up_row_major');
ncwriteatt(overlap_file,'/','coordinate_units','m');
ncwriteatt(overlap_file,'/','crs_wkt',crs);
ncwriteatt(overlap_file,'/','raster_rows',int32(1));
ncwriteatt(overlap_file,'/','raster_cols',int32(4));
for item = {'overlap_area_m2','m2'; 'mesh_area_m2','m2'; ...
        'raster_area_m2','m2'; 'x_edges','m'; 'y_edges','m'}'
    ncwriteatt(overlap_file,item{1},'units',item{2});
end
end

function write_vector(file, name, dimension, value)
nccreate(file, name, 'Dimensions', {dimension, numel(value)}, 'Datatype', 'double');
ncwrite(file, name, double(value(:)));
end
