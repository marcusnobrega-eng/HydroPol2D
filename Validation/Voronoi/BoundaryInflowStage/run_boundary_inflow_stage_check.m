% VAL-MESH-BC-002: mesh inflow as cell volume; stage edges under every scheme.
%
% HydroPol2D-Python (voronoi.py):
%   - adds an inflow hydrograph to its host cells as volume before the
%     vertical step, so it can infiltrate in the step it arrives (audit L22;
%     MATLAB applied it as an edge flux after routing);
%   - applies the local-inertial face formula at stage edges under every
%     routing scheme (L23; MATLAB's diffusive solver used Manning there), and
%     can give a stage as a depth above the edge cell's bed (12f494b).
% On a 3 x 4 mesh of 10 m cells:
%   A. A 0.05 m3/s inflow onto a dry, infiltrating west cell for 10 min with
%      walls elsewhere: the cell infiltrates in the first step, and the
%      ledger closes with the inflow volume.
%   B. One 1 s diffusive step with a 0.30 m stage above a dry west edge: the
%      edge discharge equals the local-inertial value -g h dt slope.
%   C. The same stage given as a 0.30 m depth above the local invert gives
%      the same discharge.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

mesh_file = fullfile(out_dir, 'square-mesh-3x4.nc');
write_square_mesh(mesh_file, 3, 4, 10);
mesh = HydroPol2D_Read_UGRID(mesh_file, allow_legacy=true);
west = find(mesh.edge_neighbor == 0 & abs(mesh.edge_midpoint_x) < 1e-9);
inflow_edge = west(2);
inflow_cell = mesh.edge_owner(inflow_edge);

% A. inflow before the vertical step
hydrology = struct('pervious_fraction', 1, 'soil_depth_m', 2, 'root_depth_m', 0.5, ...
    'theta_r', 0.05, 'theta_sat', 0.45, 'alpha_vg_1_m', 3.6, 'n_vg', 1.56, 'ksat_m_s', 2e-5, ...
    'initial_soil_saturation', 0.1, 'capillary_rise_enabled', false);
config = struct('duration_s', 600, 'min_dt_s', 0.01, 'max_dt_s', 5, 'output_interval_s', 600, ...
    'surface_roughness', 0.03, 'courant', 0.5, 'critical_flow', true, 'hydrology_enabled', true, ...
    'hydrology', hydrology, 'groundwater_enabled', false, 'minimum_cell_width_m', 10, ...
    'allow_legacy_mesh', true);
boundary = struct('edge_id', inflow_edge, 'type', "inflow", 'value', 0.05);
forcing = struct('surface_source_m_s', 0, 'potential_et_m_s', 0, 'surface_boundary', boundary);
results = HydroPol2D_Voronoi_Run(mesh_file, config, forcing);
d = results.diagnostics;
first_step_infiltration = d(1).infiltration_volume_m3;
inflow_volume = 0.05 * 600;
recorded_inflow = sum([d.boundary_net_inflow_volume_m3]);
max_residual = max(abs([d.step_mass_residual_m3])) / inflow_volume;
inflow_ok = first_step_infiltration > 0 && abs(recorded_inflow - inflow_volume) < 1e-9 && max_residual < 1e-9;

% B and C. a diffusive stage edge
stage_edge = west(1); stage_cell = mesh.edge_owner(stage_edge);
bed = mesh.surface_bed(stage_cell);
one_step = struct('duration_s', 1, 'min_dt_s', 1, 'max_dt_s', 1, 'output_interval_s', 1, ...
    'surface_roughness', 0.03, 'courant', 0.5, 'critical_flow', true, 'hydrology_enabled', false, ...
    'groundwater_enabled', false, 'minimum_cell_width_m', 10, 'allow_legacy_mesh', true, ...
    'routing_solver', 'diffusive', 'timestep_rule', 'unstructured');
absolute = struct('edge_id', stage_edge, 'type', "stage", 'value', bed + 0.30);
relative = struct('edge_id', stage_edge, 'type', "stage", 'value', 0.30, ...
    'stage_reference', "depth_above_local_invert");
run_absolute = HydroPol2D_Voronoi_Run(mesh_file, one_step, struct('surface_boundary', absolute));
run_relative = HydroPol2D_Voronoi_Run(mesh_file, one_step, struct('surface_boundary', relative));
q_absolute = run_absolute.edge_discharge_per_width_m2_s(stage_edge);
q_relative = run_relative.edge_discharge_per_width_m2_s(stage_edge);
h = 0.30; slope = (bed + 0.30 - bed) / mesh.edge_distance(stage_edge);
q_expected = -(-9.81 * h * 1 * slope);                % into the cell: outward discharge negative
stage_error = abs(abs(q_absolute) - q_expected);
reference_error = abs(q_absolute - q_relative);

passed = inflow_ok && stage_error < 1e-12 && reference_error < 1e-15;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-BC-002", "Mesh inflow as cell volume; stage under every scheme", ...
    first_step_infiltration, recorded_inflow, max_residual, abs(q_absolute), q_expected, ...
    stage_error, reference_error, status, true, 'VariableNames', {'case_id','case_name', ...
    'first_step_infiltration_m3','recorded_inflow_m3','max_step_residual_fraction', ...
    'stage_q_m2_s','expected_stage_q_m2_s','stage_error','stage_reference_error', ...
    'status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-BC-002 failed.');
end
fprintf('VAL-MESH-BC-002 passed.\n');

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
