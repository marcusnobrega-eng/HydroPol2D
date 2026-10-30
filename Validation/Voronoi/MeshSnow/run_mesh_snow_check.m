% VAL-MESH-SNOW-001: snow on mesh HRUs.
%
% HydroPol2D-Python runs snow_step on every HRU piece of a mesh: the
% precipitation is split into rain and snow before the canopy, only the rain
% is intercepted, and the melt leaves the pack on the ground and goes
% straight to the surface. MATLAB now does the same with Snow_Model_Function
% on every mesh HRU (hp2d_voronoi_snow_step).
%
% A 2 x 2 mesh of 40 m cells, two HRUs per cell (60 % and 40 %) with
% different degree-day factors and LAI, flat ground, 12 hours: 2 mm/h of
% precipitation for 4 h at -5 C, then +6 C. The check:
%   1. The mesh water ledger closes every step, with the snowpack as a
%      store and sublimation as a loss.
%   2. The recorded SWE equals Snow_Model_Function replayed offline on each
%      HRU with the run's own steps and forcing.
%   3. While it snows the canopy stays empty (snowfall is not intercepted),
%      and all the snow melts to the surface once it is warm.
%   4. The same run on the layered subsurface closes its ledger and gives
%      the same snowpack (the pack does not depend on the soil scheme).

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

nm = 2; dm = 40; n = nm * nm;
mesh_file = fullfile(out_dir, 'square-mesh-2x2.nc');
write_square_mesh(mesh_file, nm, nm, dm);

snow = struct('alpha', 0.8, 'epsilon', 0.97, 'C_e', 0.05, 'DDF', [3 5], ...
    'T_snow_all', 0, 'T_rain_all', 2, 'rho_snow_init', 100, 'rho_max', 450, ...
    'k_t', 1, 'k_swe', 0.05, 'k_D', 0.01);
hydrology = struct('hru_fraction', repmat([0.6 0.4], n, 1), 'pervious_fraction', [1 0.3], ...
    'lai', [3 1], 'soil_depth_m', 2, 'root_depth_m', 0.5, 'theta_r', 0.05, ...
    'theta_sat', 0.45, 'alpha_vg_1_m', 3.6, 'n_vg', 1.56, 'ksat_m_s', 5e-6, ...
    'initial_soil_saturation', 0.3, 'capillary_rise_enabled', false);
config = struct('duration_s', 12 * 3600, 'min_dt_s', 0.01, 'max_dt_s', 60, ...
    'output_interval_s', 1800, 'surface_roughness', 0.03, 'courant', 0.5, ...
    'critical_flow', true, 'hydrology_enabled', true, 'hydrology', hydrology, ...
    'groundwater_enabled', false, 'minimum_cell_width_m', dm, 'allow_legacy_mesh', true, ...
    'snow', snow);
rain_rate = @(t) 2e-3 / 3600 * double(t < 4 * 3600);
temperature = @(t) -5 + 11 * double(t >= 4 * 3600);
meteorology = @(t, state, mesh) struct('temperature_c', temperature(t), ...
    'minimum_temperature_c', temperature(t) - 3, 'wind_speed_m_s', 2, ...
    'latitude_deg', 45, 'day_of_year', 60);
forcing = struct('surface_source_m_s', @(t, state, mesh) rain_rate(t), ...
    'potential_et_m_s', 0, 'meteorology', meteorology);
results = HydroPol2D_Voronoi_Run(mesh_file, config, forcing);
d = results.diagnostics;
M = results.modules;

% 1. water ledger
rain_volume = sum([d.surface_source_volume_m3]);
max_step_residual = max(abs([d.step_mass_residual_m3])) / rain_volume;

% 2. offline replay of every HRU
modules = hp2d_voronoi_modules_initialize(HydroPol2D_Read_UGRID(mesh_file, allow_legacy=true), ...
    mesh_file, config, '', false, false, repmat([0.6 0.4], n, 1));
pack = modules.snow;
replay_swe = zeros(n, numel(results.time_s));
output = 2;
for k = 1:numel(d)
    t0 = d(k).time_s - d(k).dt_s;
    met = meteorology(t0, [], []);
    met.temperature_c = repmat(met.temperature_c, n, 1);
    met.minimum_temperature_c = repmat(met.minimum_temperature_c, n, 1);
    met.wind_speed_m_s = repmat(met.wind_speed_m_s, n, 1);
    met.latitude_deg = repmat(met.latitude_deg, n, 1);
    pack = hp2d_voronoi_snow_step(pack, repmat(rain_rate(t0), n, 1), met, d(k).dt_s);
    if output <= numel(results.time_s) && abs(d(k).time_s - results.time_s(output)) < 1e-6
        replay_swe(:, output) = sum(pack.swe_mm .* pack.fraction, 2);
        output = output + 1;
    end
end
replay_error = max(abs(replay_swe(:) - M.snow_swe_mm(:)));

% 3. no interception of snowfall; full melt when warm
snowing = results.time_s > 0 & results.time_s <= 4 * 3600;
canopy_while_snowing = max(results.hydrology.canopy_storage_m(:, snowing), [], 'all');
peak_swe = max(M.snow_swe_mm, [], 'all');
final_swe = max(M.snow_swe_mm(:, end));
snow_sublimation = sum([d.snow_sublimation_volume_m3]);

% 4. the layered subsurface
layered_config = config;
layered_config.layered_subsurface = hp2d_layered_subsurface_config(struct('flag_layered_subsurface', 1));
layered = HydroPol2D_Voronoi_Run(mesh_file, layered_config, forcing);
layered_residual = max(abs([layered.diagnostics.step_mass_residual_m3])) / rain_volume;
layered_swe_error = max(abs(layered.modules.snow_swe_mm(:, end) - M.snow_swe_mm(:, end)));
layered_peak_error = abs(max(layered.modules.snow_swe_mm, [], 'all') - peak_swe);

passed = max_step_residual < 1e-9 && replay_error < 1e-9 && canopy_while_snowing == 0 && ...
    peak_swe > 7 && final_swe < 1e-9 && layered_residual < 1e-9 && ...
    layered_swe_error < 1e-9 && layered_peak_error < 1e-6;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-SNOW-001", "Snow on mesh HRUs", max_step_residual, replay_error, ...
    canopy_while_snowing, peak_swe, final_swe, snow_sublimation, layered_residual, ...
    layered_peak_error, status, true, ...
    'VariableNames', {'case_id','case_name','max_step_residual_fraction','replay_swe_error_mm', ...
    'canopy_storage_while_snowing_m','peak_swe_mm','final_swe_mm','sublimation_m3', ...
    'layered_step_residual_fraction','layered_peak_swe_error_mm','status','report_ready'});
writetable(table(results.time_s / 3600, M.snow_swe_mm(1, :)', replay_swe(1, :)', ...
    'VariableNames', {'time_h','swe_cell1_mm','replay_swe_cell1_mm'}), ...
    fullfile(out_dir, 'Mesh_Snow_SWE.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-SNOW-001 failed.');
end
fprintf('VAL-MESH-SNOW-001 passed.\n');

function write_square_mesh(path, ny, nx, dx)
if exist(path, 'file') == 2, delete(path); end
n = ny * nx;
id = @(r, c) (r - 1) * nx + c;
[cols, rows] = meshgrid(1:nx, 1:ny);
cell_x = (reshape(cols', [], 1) - 0.5) * dx;
cell_y = (reshape(rows', [], 1) - 0.5) * dx;
owner = []; neighbor = []; mx = []; my = []; normal_x = []; normal_y = [];
for r = 1:ny
    for c = 1:nx
        if c < nx
            owner(end+1) = id(r, c); neighbor(end+1) = id(r, c+1); %#ok<AGROW>
            mx(end+1) = c * dx; my(end+1) = (r - 0.5) * dx; normal_x(end+1) = 1; normal_y(end+1) = 0; %#ok<AGROW>
        end
        if r < ny
            owner(end+1) = id(r, c); neighbor(end+1) = id(r+1, c); %#ok<AGROW>
            mx(end+1) = (c - 0.5) * dx; my(end+1) = r * dx; normal_x(end+1) = 0; normal_y(end+1) = 1; %#ok<AGROW>
        end
    end
end
for r = 1:ny
    owner(end+1) = id(r, 1); neighbor(end+1) = 0; mx(end+1) = 0; my(end+1) = (r - 0.5) * dx; normal_x(end+1) = -1; normal_y(end+1) = 0; %#ok<AGROW>
    owner(end+1) = id(r, nx); neighbor(end+1) = 0; mx(end+1) = nx * dx; my(end+1) = (r - 0.5) * dx; normal_x(end+1) = 1; normal_y(end+1) = 0; %#ok<AGROW>
end
for c = 1:nx
    owner(end+1) = id(1, c); neighbor(end+1) = 0; mx(end+1) = (c - 0.5) * dx; my(end+1) = 0; normal_x(end+1) = 0; normal_y(end+1) = -1; %#ok<AGROW>
    owner(end+1) = id(ny, c); neighbor(end+1) = 0; mx(end+1) = (c - 0.5) * dx; my(end+1) = ny * dx; normal_x(end+1) = 0; normal_y(end+1) = 1; %#ok<AGROW>
end
n_edges = numel(owner);
write_vector(path, 'cell_area_m2', 'cell', dx^2 * ones(n, 1));
write_vector(path, 'cell_bed_elevation_m', 'cell', zeros(n, 1));
write_vector(path, 'mesh2d_face_x', 'cell', cell_x);
write_vector(path, 'mesh2d_face_y', 'cell', cell_y);
write_vector(path, 'cell_target_width_m', 'cell', dx * ones(n, 1));
write_vector(path, 'cell_refinement_source', 'cell', zeros(n, 1));
write_vector(path, 'cell_hydraulic_roughness', 'cell', nan(n, 1));
write_vector(path, 'cell_riverbed_area_m2', 'cell', zeros(n, 1));
write_vector(path, 'cell_riverbed_hydraulic_conductivity_m_s', 'cell', zeros(n, 1));
write_vector(path, 'cell_cfl_width_m', 'cell', dx * ones(n, 1));
write_vector(path, 'edge_owner', 'edge', owner - 1);
write_vector(path, 'edge_neighbor', 'edge', neighbor - 1);
write_vector(path, 'edge_length_m', 'edge', dx * ones(n_edges, 1));
write_vector(path, 'edge_center_distance_m', 'edge', dx * ones(n_edges, 1));
write_vector(path, 'edge_midpoint_x', 'edge', mx);
write_vector(path, 'edge_midpoint_y', 'edge', my);
write_vector(path, 'edge_normal_x', 'edge', normal_x);
write_vector(path, 'edge_normal_y', 'edge', normal_y);
write_vector(path, 'edge_boundary_type', 'edge', zeros(n_edges, 1));
end

function write_vector(path, name, dimension, value)
nccreate(path, name, 'Dimensions', {dimension, numel(value)}, 'Datatype', 'double');
ncwrite(path, name, double(value(:)));
end
