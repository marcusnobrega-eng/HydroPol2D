function summary = run_voronoi_dry_face_cfl_validation(output_directory)
%RUN_VORONOI_DRY_FACE_CFL_VALIDATION Reject stale wet-face CFL velocities.

arguments
    output_directory (1,:) char = fullfile(fileparts(mfilename('fullpath')), 'outputs_dry_face_cfl')
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'HydroPol2D_Functions'));
if exist(output_directory, 'dir') ~= 7, mkdir(output_directory); end
mesh_file = fullfile(output_directory, 'dry-face-cfl-fixture.nc');
write_fixture(mesh_file);

config = struct('routing_solver','local_inertial','duration_s',1, ...
    'initial_surface_depth_m',[1e-4;0], ...
    'initial_surface_discharge_per_width_m2_s',[1;zeros(6,1)], ...
    'surface_roughness',0.04,'min_dt_s',0.1,'max_dt_s',1, ...
    'output_interval_s',1,'courant',0.2,'critical_flow',true, ...
    'minimum_cell_width_m',100,'allow_legacy_mesh',true);
results = HydroPol2D_Voronoi_Run(mesh_file,config,struct());
dt = [results.diagnostics.dt_s]';
initial_volume = 1e-4 * 10000;
final_volume = sum(results.final_surface_volume_m3);
relative_error = abs(final_volume-initial_volume)/initial_volume;

assert(min(dt) >= config.min_dt_s, 'A stale wet-face flux collapsed the timestep.');
assert(relative_error <= 1e-10, 'Dry-face CFL test mass error %.3g.', relative_error);
assert(all(results.final_surface_volume_m3 >= 0));
summary = struct('minimum_timestep_s',min(dt), ...
    'relative_mass_error',relative_error,'passed',true);
end

function write_fixture(path)
if exist(path,'file') == 2, delete(path); end
write_vector(path,'cell_area_m2','cell',10000*ones(2,1));
write_vector(path,'cell_bed_elevation_m','cell',zeros(2,1));
write_vector(path,'mesh2d_face_x','cell',[50;150]);
write_vector(path,'mesh2d_face_y','cell',[50;50]);
write_vector(path,'cell_target_width_m','cell',100*ones(2,1));
write_vector(path,'cell_refinement_source','cell',zeros(2,1));
write_vector(path,'cell_hydraulic_roughness','cell',nan(2,1));
write_vector(path,'cell_cfl_width_m','cell',100*ones(2,1));
write_vector(path,'edge_owner','edge',[0;0;0;0;1;1;1]);
write_vector(path,'edge_neighbor','edge',[1;-1;-1;-1;-1;-1;-1]);
write_vector(path,'edge_length_m','edge',100*ones(7,1));
write_vector(path,'edge_center_distance_m','edge',100*ones(7,1));
write_vector(path,'edge_midpoint_x','edge',[100;0;50;50;200;150;150]);
write_vector(path,'edge_midpoint_y','edge',[50;50;0;100;50;0;100]);
write_vector(path,'edge_normal_x','edge',[1;-1;0;0;1;0;0]);
write_vector(path,'edge_normal_y','edge',[0;0;-1;1;0;-1;1]);
write_vector(path,'edge_boundary_type','edge',zeros(7,1));
end

function write_vector(path,name,dimension,value)
nccreate(path,name,'Dimensions',{dimension,numel(value)},'Datatype','double');
ncwrite(path,name,double(value(:)));
end
