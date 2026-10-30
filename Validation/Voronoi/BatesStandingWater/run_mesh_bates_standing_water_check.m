% VAL-MESH-BATES-001: the mesh time step counts standing water.
%
% On axis-aligned square meshes MATLAB uses the structured Bates rule cell by
% cell (bates_cartesian_timestep). It kept only cells moving at 0.01 m/s or
% more, so a deep lake at rest set no wave limit and the step jumped to
% max_dt_s. HydroPol2D-Python's rule (voronoi.py) takes
% alfa_min * width / (|q/h| + sqrt(g h)) over every wet cell, as the grid
% does since 20f2ded (audit M22).
% A 4 x 4 mesh of 10 m cells holds a 1 m lake at rest for 60 s: every step
% must be alfa_min * 10 / sqrt(9.81) floored to the 1 ms increment, and the
% lake must stay at rest.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

n_side = 4; dx = 10; depth = 1; courant = 0.5;
mesh_file = fullfile(out_dir, 'square-mesh-4x4.nc');
write_square_mesh(mesh_file, n_side, n_side, dx);
config = struct('duration_s', 60, 'min_dt_s', 0.01, 'max_dt_s', 30, ...
    'output_interval_s', 60, 'surface_roughness', 0.03, 'courant', courant, ...
    'critical_flow', true, 'hydrology_enabled', false, 'groundwater_enabled', false, ...
    'minimum_cell_width_m', dx, 'allow_legacy_mesh', true, 'initial_surface_depth_m', depth);
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct('surface_source_m_s', 0));
dt = [results.diagnostics.dt_s];
expected = 0.001 * floor(courant * dx / sqrt(9.81 * depth) / 0.001);
dt_error = max(abs(dt(1:end-1) - expected));        % the last step ends the run
depth_error = max(abs(results.surface_depth_m(:, end) - depth));
passed = dt_error < 1e-12 && depth_error < 1e-12;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-BATES-001", "Mesh time step counts standing water", expected, ...
    max(dt), dt_error, depth_error, status, true, 'VariableNames', {'case_id','case_name', ...
    'expected_dt_s','largest_dt_s','dt_error_s','depth_error_m','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-BATES-001 failed.');
end
fprintf('VAL-MESH-BATES-001 passed.\n');

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
