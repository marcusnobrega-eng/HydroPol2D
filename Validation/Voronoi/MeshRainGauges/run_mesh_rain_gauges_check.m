% VAL-MESH-GAUGE-001: rain gauges on meshes.
%
% HydroPol2D-Python spreads rain-gauge series over a mesh with normalised
% inverse-distance weights (power 2): over the gauges that report, a point
% on a reporting gauge takes its value, and the rain stops one interval
% after the last row (forcing/gauge_rainfall.py, voronoi.py). MATLAB meshes
% could not take gauges at all (audit M17). hp2d_gauge_rainfall now does
% this, and HydroPol2D_Voronoi_Run reads forcing.rain_gauges.
%   A. Weights: three gauges, four points (one on gauge 1), two times (gauge
%      3 not reporting at the second), and the rate after the last row.
%   B. A 2 x 2 mesh of 40 m cells with walls, two gauges and three 30 min
%      rows over 2 h: the stored water must equal the gauge rain exactly.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

% A. weights
gauges = struct('x', [0 100 0], 'y', [0 0 100], 'time_s', [0; 600], ...
    'intensity_mm_h', [10 20 40; 5 15 NaN]);
px = [0; 50; 30; 100]; py = [0; 50; 80; 100];
W = hp2d_gauge_rainfall('weights', gauges, px, py);
idw = @(v, use) arrayfun(@(k) sum(v(use) ./ hypot(px(k) - gauges.x(use), py(k) - gauges.y(use)).^2) / ...
    sum(1 ./ hypot(px(k) - gauges.x(use), py(k) - gauges.y(use)).^2), (1:4)');
expected1 = idw([10 20 40], [true true true]); expected1(1) = 10;
expected2 = idw([5 15 NaN], [true true false]); expected2(1) = 5;
rate1 = hp2d_gauge_rainfall('rate', W, 300);
rate2 = hp2d_gauge_rainfall('rate', W, 900);
rate_after = hp2d_gauge_rainfall('rate', W, 1200);
rate_before = hp2d_gauge_rainfall('rate', W, -1);
weights_error = max(abs([rate1 - expected1; rate2 - expected2]));
stops = all(rate_after == 0) && all(rate_before == 0);

% B. a mesh run
mesh_file = fullfile(out_dir, 'square-mesh-2x2.nc');
write_square_mesh(mesh_file, 2, 2, 40);
mesh = HydroPol2D_Read_UGRID(mesh_file, allow_legacy=true);
run_gauges = struct('x', [0 80], 'y', [0 80], 'time_s', [0; 1800; 3600], ...
    'intensity_mm_h', [20 40; 10 0; 30 30]);
config = struct('duration_s', 7200, 'min_dt_s', 0.01, 'max_dt_s', 30, 'output_interval_s', 1800, ...
    'surface_roughness', 0.03, 'courant', 0.5, 'critical_flow', true, 'hydrology_enabled', false, ...
    'groundwater_enabled', false, 'minimum_cell_width_m', 40, 'allow_legacy_mesh', true, ...
    'forcing_interval_s', 1800);
results = HydroPol2D_Voronoi_Run(mesh_file, config, struct('rain_gauges', run_gauges));
RW = hp2d_gauge_rainfall('weights', run_gauges, mesh.cell_x(:), mesh.cell_y(:));
expected_volume = 0;
for t = [0 1800 3600 5400]
    expected_volume = expected_volume + sum(hp2d_gauge_rainfall('rate', RW, t) / 1000 / 3600 .* mesh.cell_area(:)) * 1800;
end
stored = sum(results.final_surface_volume_m3);
volume_error = abs(stored - expected_volume) / expected_volume;

passed = weights_error < 1e-12 && stops && volume_error < 1e-12;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-GAUGE-001", "Rain gauges on meshes", weights_error, stops, ...
    expected_volume, stored, volume_error, status, true, 'VariableNames', {'case_id','case_name', ...
    'idw_error_mm_h','rain_stops_outside_rows','expected_rain_m3','stored_m3', ...
    'volume_error_fraction','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-GAUGE-001 failed.');
end
fprintf('VAL-MESH-GAUGE-001 passed.\n');

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
