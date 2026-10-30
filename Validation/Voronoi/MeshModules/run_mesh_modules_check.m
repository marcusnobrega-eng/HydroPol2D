% VAL-MESH-MODULES-001: reservoirs, water quality and human risk on meshes.
%
% HydroPol2D-Python runs the grid's optional modules on mesh cells
% (mesh_modules/). MATLAB now does the same for reservoir controls, the
% mass-based wash-off and the human risk index (hp2d_voronoi_modules_*).
% The check has two parts on a 3 x 4 mesh of 10 m squares:
%   1. One module step on fixed volumes and edge discharges, against the
%      Python kernels on the same numbers (reservoir_outlet_transfers,
%      mass_based_washoff_edge_step, cell_discharge_operators and
%      force_balance_risk). The reference is python_reference_*.csv, made by
%      make_python_reference.py from the kernel_*.csv inputs this check
%      writes.
%   2. A one-hour mesh run with rain, a critical-flow outlet and all three
%      modules: the water ledger closes with the reservoir outflow out of
%      the model, the pollutant mass closes with its export, and the risk
%      index stays in [0, 1].
% MATLAB's wash-off coefficient C_3 is per hour (W in kg/h); Python reads
% it per second (audit section 4), so the reference gets C_3/3600.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

ny = 3; nx = 4; dx = 10;
mesh_file = fullfile(out_dir, 'square-mesh-3x4.nc');
write_square_mesh(mesh_file, ny, nx, dx);
mesh = HydroPol2D_Read_UGRID(mesh_file, allow_legacy=true);
n = mesh.n_cells;

% ---- 1. one module step against the Python kernels ----------------------
cell_id = (1:n)';
volume = 10 + 7 * mod(cell_id * 37, 11);                   % m3 (0.1 to 0.8 m)
volume([6 11]) = [65; 48];                                  % reservoirs above every sill
edge_id = (1:mesh.n_edges)';
edge_q = 0.02 * sin(1.3 * edge_id) .* (1 + 0.5 * mod(edge_id, 3));   % m2/s
reservoirs = struct('source_cell', [6; 11], 'ds1_cell', [7; 0], 'ds2_cell', [0; NaN], ...
    'k1', [1.5; 2.0], 'h1', [0.2; 0.1], 'k2', [1.5; 1.0], ...
    'k3', [0.8; 0], 'h2', [0.4; 0], 'k4', [2.0; 1]);
quality = struct('C_3', 0.5 + 0.1 * mod(cell_id, 4), 'C_4', 1.2, 'min_Bt', 0.5, ...
    'initial_mass_kg', 2 + mod(cell_id * 13, 5));
quality.initial_mass_kg(3) = 0.01;                          % below the 0.5 g/m2 threshold
risk = struct('mu', 0.6, 'Cd', 1.0, 'weight_person', 75, 'width1_person', 0.35, ...
    'width2_person', 0.20, 'ro_water', 1000, 'gravity', 9.81);
config = struct('reservoirs', reservoirs, 'water_quality', quality, 'human_risk', risk);
modules = hp2d_voronoi_modules_initialize(mesh, mesh_file, config, '', false, false);
dt = 30;
[volume_after, modules_after, step_diag] = hp2d_voronoi_modules_step(mesh, modules, volume, edge_q, dt);

writetable(table(cell_id, mesh.cell_area, mesh.cell_x, mesh.cell_y, volume, ...
    modules.water_quality.mass_kg, modules.water_quality.coefficient_kg_s, ...
    modules.water_quality.exponent, modules.water_quality.threshold_kg, ...
    'VariableNames', {'cell','area_m2','x_m','y_m','volume_m3','mass_kg', ...
    'coefficient_kg_s','exponent','threshold_kg'}), fullfile(case_dir, 'kernel_cells.csv'));
writetable(table(edge_id, mesh.edge_owner, mesh.edge_neighbor, mesh.edge_length, ...
    mesh.edge_midpoint_x, mesh.edge_midpoint_y, edge_q, ...
    'VariableNames', {'edge','owner','neighbor','length_m','midpoint_x_m', ...
    'midpoint_y_m','q_m2_s'}), fullfile(case_dir, 'kernel_edges.csv'));
outlet_source = [6; 6; 11]; outlet_destination = [7; 0; 0];
outlet_external = [false; true; true]; outlet_control = [1; 1; 2];
writetable(table(outlet_source, outlet_destination, outlet_external, ...
    [1.5; 0.8; 2.0], [0.2; 0.4; 0.1], [1.5; 2.0; 1.0], outlet_control, ...
    'VariableNames', {'source','destination','external','coefficient', ...
    'sill_depth_m','exponent','control'}), fullfile(case_dir, 'kernel_outlets.csv'));

reference_cells = readtable(fullfile(case_dir, 'python_reference_cells.csv'));
reference_scalars = readtable(fullfile(case_dir, 'python_reference_scalars.csv'));
volume_error = max(abs(volume_after - reference_cells.volume_m3));
mass_error = max(abs(modules_after.water_quality.mass_kg - reference_cells.mass_kg));
risk_error = max(abs(modules_after.human_risk.risk - reference_cells.risk_index));
concentration_error = max(abs(modules_after.water_quality.cell_concentration_mg_l - ...
    reference_cells.cell_concentration_mg_l));
external_error = abs(step_diag.reservoir_outflow_m3 - reference_scalars.external_m3);
release_error = max(abs(modules_after.reservoirs.discharge_m3_s * dt - ...
    [reference_scalars.control_1_release_m3; reference_scalars.control_2_release_m3]));
export_error = abs(step_diag.pollutant_export_kg - reference_scalars.exported_kg);
outlet_concentration_error = abs(modules_after.water_quality.outlet_concentration_mg_l - ...
    reference_scalars.outlet_concentration_mg_l);
step_water_residual = abs(sum(volume_after) + step_diag.reservoir_outflow_m3 - sum(volume));
step_pollutant_residual = abs(sum(modules_after.water_quality.mass_kg) + ...
    step_diag.pollutant_export_kg - sum(modules.water_quality.mass_kg));
kernel_errors = [volume_error, mass_error, risk_error, concentration_error, ...
    external_error, release_error, export_error, outlet_concentration_error];
kernel_passed = all(kernel_errors < 1e-9) && step_water_residual < 1e-9 && ...
    step_pollutant_residual < 1e-12;

% ---- 2. a one-hour mesh run with the three modules ----------------------
west = find(mesh.edge_neighbor == 0 & abs(mesh.edge_midpoint_x) < 1e-9);
run_reservoirs = struct('source_cell', 6, 'ds1_cell', 0, 'ds2_cell', NaN, ...
    'k1', 0.05, 'h1', 0, 'k2', 1.5, 'k3', 0, 'h2', 0, 'k4', 1);
run_config = struct('duration_s', 3600, 'min_dt_s', 0.01, 'max_dt_s', 10, ...
    'output_interval_s', 600, 'surface_roughness', 0.03, 'courant', 0.5, ...
    'critical_flow', true, 'hydrology_enabled', false, 'groundwater_enabled', false, ...
    'minimum_cell_width_m', dx, 'allow_legacy_mesh', true, ...
    'reservoirs', run_reservoirs, 'water_quality', quality, 'human_risk', risk);
rain = 50e-3 / 3600;
forcing = struct('surface_source_m_s', @(time_s, state, current_mesh) rain * double(time_s < 1800), ...
    'surface_boundary', struct('edge_id', west, 'type', "critical_flow", 'value', 0));
results = HydroPol2D_Voronoi_Run(mesh_file, run_config, forcing);
d = results.diagnostics;
rain_volume = sum([d.surface_source_volume_m3]);
reservoir_outflow = sum([d.reservoir_outflow_volume_m3]);
boundary_outflow = -sum([d.boundary_net_inflow_volume_m3]);
final_water = sum(results.final_surface_volume_m3);
run_water_residual = abs(rain_volume - reservoir_outflow - boundary_outflow - final_water) / rain_volume;
max_step_residual = max(abs([d.step_mass_residual_m3])) / rain_volume;
M = results.modules;
run_pollutant_residual = abs(sum(quality.initial_mass_kg) - M.total_pollutant_mass_kg(end) - ...
    M.cumulative_pollutant_export_kg(end)) / sum(quality.initial_mass_kg);
risk_in_range = all(M.human_risk_index(:) >= 0 & M.human_risk_index(:) <= 1) && ...
    all(M.maximum_human_risk_index >= max(M.human_risk_index, [], 2) - 1e-15);
run_passed = run_water_residual < 1e-9 && max_step_residual < 1e-9 && ...
    run_pollutant_residual < 1e-9 && reservoir_outflow > 0 && ...
    M.cumulative_pollutant_export_kg(end) > 0 && risk_in_range && ...
    abs(M.cumulative_reservoir_outflow_m3(end) - reservoir_outflow) < 1e-9;

passed = kernel_passed && run_passed;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-MODULES-001", "Reservoirs, water quality and human risk on meshes", ...
    max(kernel_errors), step_water_residual, step_pollutant_residual, run_water_residual, ...
    max_step_residual, run_pollutant_residual, reservoir_outflow, ...
    M.cumulative_pollutant_export_kg(end), max(M.maximum_human_risk_index), status, true, ...
    'VariableNames', {'case_id','case_name','max_kernel_error_vs_python', ...
    'step_water_residual_m3','step_pollutant_residual_kg','run_water_residual_fraction', ...
    'max_step_residual_fraction','run_pollutant_residual_fraction','reservoir_outflow_m3', ...
    'pollutant_export_kg','maximum_risk_index','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(table(kernel_errors', 'RowNames', {'volume','mass','risk','cell_concentration', ...
    'external','release','export','outlet_concentration'}, 'VariableNames', {'max_abs_error'}));
disp(PassFail);
if ~passed
    error('VAL-MESH-MODULES-001 failed.');
end
fprintf('VAL-MESH-MODULES-001 passed.\n');

function write_square_mesh(path, ny, nx, dx)
% Cells row-major from the south-west corner; edges owner to neighbour.
if exist(path, 'file') == 2, delete(path); end
n = ny * nx;
id = @(r, c) (r - 1) * nx + c;
[cols, rows] = meshgrid(1:nx, 1:ny);
cell_x = (reshape(cols', [], 1) - 0.5) * dx;
cell_y = (reshape(rows', [], 1) - 0.5) * dx;
owner = []; neighbor = []; mx = []; my = []; normal_x = []; normal_y = [];
for r = 1:ny
    for c = 1:nx
        if c < nx   % east face
            owner(end+1) = id(r, c); neighbor(end+1) = id(r, c+1); %#ok<AGROW>
            mx(end+1) = c * dx; my(end+1) = (r - 0.5) * dx; normal_x(end+1) = 1; normal_y(end+1) = 0; %#ok<AGROW>
        end
        if r < ny   % north face
            owner(end+1) = id(r, c); neighbor(end+1) = id(r+1, c); %#ok<AGROW>
            mx(end+1) = (c - 0.5) * dx; my(end+1) = r * dx; normal_x(end+1) = 0; normal_y(end+1) = 1; %#ok<AGROW>
        end
    end
end
for r = 1:ny   % west and east boundaries
    owner(end+1) = id(r, 1); neighbor(end+1) = 0; mx(end+1) = 0; my(end+1) = (r - 0.5) * dx; normal_x(end+1) = -1; normal_y(end+1) = 0; %#ok<AGROW>
    owner(end+1) = id(r, nx); neighbor(end+1) = 0; mx(end+1) = nx * dx; my(end+1) = (r - 0.5) * dx; normal_x(end+1) = 1; normal_y(end+1) = 0; %#ok<AGROW>
end
for c = 1:nx   % south and north boundaries
    owner(end+1) = id(1, c); neighbor(end+1) = 0; mx(end+1) = (c - 0.5) * dx; my(end+1) = 0; normal_x(end+1) = 0; normal_y(end+1) = -1; %#ok<AGROW>
    owner(end+1) = id(ny, c); neighbor(end+1) = 0; mx(end+1) = (c - 0.5) * dx; my(end+1) = ny * dx; normal_x(end+1) = 0; normal_y(end+1) = 1; %#ok<AGROW>
end
n_edges = numel(owner);
write_vector(path, 'cell_area_m2', 'cell', dx^2 * ones(n, 1));
write_vector(path, 'cell_bed_elevation_m', 'cell', 0.01 * cell_x);   % falls to the west
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
