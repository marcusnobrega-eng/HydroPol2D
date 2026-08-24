function summary = run_voronoi_lake_at_rest_validation(mesh_file)
%RUN_VORONOI_LAKE_AT_REST_VALIDATION Constant WSE must remain stationary.

mesh = HydroPol2D_Read_UGRID(mesh_file);
initial_depth = 2 - mesh.cell_bed;
config = struct('duration_s', 300, 'initial_surface_depth_m', initial_depth, ...
    'min_dt_s', 1e-4, 'max_dt_s', 5, 'output_interval_s', 30, ...
    'surface_roughness', 0.04);
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct());
final_depth = results.final_surface_volume_m3 ./ mesh.cell_area;
maximum_depth_error = max(abs(final_depth - initial_depth));
assert(maximum_depth_error <= 1e-12, 'Lake-at-rest depth error %.3g.', maximum_depth_error);
assert(max(abs(results.edge_discharge_per_width_m2_s)) <= 1e-12);
summary = struct('maximum_depth_error_m', maximum_depth_error, ...
    'maximum_flux_per_width_m2_s', max(abs(results.edge_discharge_per_width_m2_s)));
end
