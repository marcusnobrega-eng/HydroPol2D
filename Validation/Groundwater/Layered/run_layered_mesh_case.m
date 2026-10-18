function Summary = run_layered_mesh_case()
%RUN_LAYERED_MESH_CASE Layered subsurface on the 20 m V-tilted Voronoi mesh (2 HRUs per cell).
%   Same storm, HRUs and aquifer as Validation/Voronoi/run_voronoi_coupled_hydrology_validation
%   (30 mm/h for 3 h, 6 h run, lateral groundwater on). The fixture mesh is copied to
%   tempdir with zero riverbed area and conductivity (the native-mesh preflight requires
%   them). Runs legacy, layered, and layered + Feddes; every step must close mass within
%   1e-8 of the rain volume.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'));
source = fullfile(root, 'Validation', 'Voronoi', 'fixtures_vtilted_20m', 'vtilted-uniform20-mesh.nc');
work = fullfile(tempdir, 'hp2d_layered_mesh');
if isfolder(work), rmdir(work, 's'); end
mkdir(work);
mesh_file = fullfile(work, 'vtilted-uniform20-mesh.nc');
copyfile(source, mesh_file);
info = ncinfo(mesh_file, 'cell_area_m2');
face_dim = info.Dimensions(1).Name;
n = info.Dimensions(1).Length;
nccreate(mesh_file, 'cell_riverbed_area_m2', 'Dimensions', {face_dim, n});
ncwrite(mesh_file, 'cell_riverbed_area_m2', zeros(n, 1));
nccreate(mesh_file, 'cell_riverbed_hydraulic_conductivity_m_s', 'Dimensions', {face_dim, n});
ncwrite(mesh_file, 'cell_riverbed_hydraulic_conductivity_m_s', 1e-6*ones(n, 1));
mesh = HydroPol2D_Read_UGRID(mesh_file);
outlet = find(mesh.edge_neighbor == 0 & abs(mesh.edge_midpoint_y) < 1e-8 & ...
    mesh.edge_midpoint_x >= 380 & mesh.edge_midpoint_x <= 440);
hydrology = struct('hru_fraction', repmat([0.65 0.35], mesh.n_cells, 1), ...
    'pervious_fraction', [1 0.25], 'lai', [3 0], 'crop_coefficient', [1 0.4], ...
    'soil_depth_m', [2 1.2], 'root_depth_m', [1 0.3], ...
    'theta_r', [0.05 0.07], 'theta_sat', [0.45 0.42], ...
    'alpha_vg_1_m', [3.6 1.9], 'n_vg', [1.56 1.31], ...
    'ksat_m_s', [5e-6 5e-7], 'initial_soil_saturation', [0.35 0.45], ...
    'wilting_saturation', [0.10 0.15], 'capillary_rise_enabled', true, ...
    'soil_base_depth_m', 1.5);
base = struct('duration_s', 6*3600, 'min_dt_s', 0.01, 'max_dt_s', 30, ...
    'output_interval_s', 600, 'surface_roughness', 0.04, 'hydrology_enabled', true, ...
    'hydrology', hydrology, 'groundwater_enabled', true, ...
    'initial_groundwater_head_m', mesh.surface_bed - 2, ...
    'aquifer_bottom_m', mesh.surface_bed - 10, 'hydraulic_conductivity_m_s', 1e-5, ...
    'specific_yield', 0.2, 'groundwater_relative_tolerance', 1e-3, ...
    'groundwater_head_tolerance_m', 1e-3, ...
    'river_groundwater_exchange', false, 'riverbed_thickness_m', 0.5, 'overwrite_output', true);
forcing = struct('surface_source_m_s', @(time_s, state, current_mesh) double(time_s < 3*3600)*30e-3/3600, ...
    'meteorology', struct('temperature_c', 25, 'maximum_temperature_c', 30, ...
    'minimum_temperature_c', 20, 'day_of_year', 200, 'latitude_deg', 20, ...
    'wind_speed_m_s', 2, 'relative_humidity_pct', 60, 'krs', 0.16, 'albedo', 0.23), ...
    'surface_boundary', struct('edge_id', outlet, 'type', "normal_flow", 'value', 0.004));
rain_volume = 3*3600*30e-3/3600*sum(mesh.cell_area);
variants = {'legacy', 'layered', 'layered_feddes'};
residual = zeros(3, 1); recharge = residual; head_min = residual; head_max = residual;
outlet_volume = residual; actual_et = residual;
setenv('HYDROPOL_SUBSURFACE_SCHEME', '');
setenv('HYDROPOL_ENABLE_FEDDES', '');
for i = 1:numel(variants)
    config = base;
    config.output_netcdf = fullfile(work, [variants{i} '.nc']);
    flags = struct('flag_layered_subsurface', double(i > 1), 'flag_feddes', double(i == 3), ...
        'groundwater_target_dt_s', 600, 'subsurface_max_depth_m', 10);
    config.layered_subsurface = hp2d_layered_subsurface_config(flags);
    run = HydroPol2D_Voronoi_Run(mesh_file, config, forcing);
    residual(i) = max(abs([run.diagnostics.step_mass_residual_m3]))/rain_volume;
    recharge(i) = sum([run.diagnostics.recharge_volume_m3]);
    relative_head = double(run.groundwater_head_m) - mesh.surface_bed(:);
    head_min(i) = min(relative_head, [], 'all');
    head_max(i) = max(relative_head, [], 'all');
    outlet_volume(i) = -sum([run.diagnostics.boundary_net_inflow_volume_m3]);
    actual_et(i) = sum([run.diagnostics.actual_et_volume_m3]);
end
Summary = table(variants(:), residual, recharge, actual_et, outlet_volume, head_min, head_max, ...
    'VariableNames', {'variant', 'max_step_residual_fraction', 'recharge_m3', 'actual_et_m3', ...
    'outlet_volume_m3', 'min_head_minus_ground_m', 'max_head_minus_ground_m'});
format long g
disp(Summary);
assert(all(residual < 1e-8), 'HydroPol2D:LayeredMesh', 'A mesh step did not close mass.');
end
