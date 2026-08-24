function summary = run_voronoi_output_validation(mesh_file, output_file)
%RUN_VORONOI_OUTPUT_VALIDATION Native NetCDF state/checkpoint contract.

config = struct('duration_s',10,'initial_surface_depth_m',0.02, ...
    'initial_channel_depth_m',0.1,'max_dt_s',1,'min_dt_s',1e-4, ...
    'output_interval_s',2,'output_netcdf',output_file,'overwrite_output',true);
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct());
saved_depth = double(ncread(output_file,'surface_depth_m'));
saved_volume = double(ncread(output_file,'final_surface_volume_m3'));
assert(isequal(size(saved_depth),size(results.surface_depth_m)));
assert(max(abs(saved_depth-results.surface_depth_m),[],'all') <= 1e-14);
assert(max(abs(saved_volume-results.final_surface_volume_m3)) <= 1e-14);
summary = struct('time_count',numel(results.time_s),'variable_count',numel(ncinfo(output_file).Variables));
end
