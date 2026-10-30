% VAL-MESH-CHANBC-001: Neal channel nodes take their boundary cell's outlet.
%
% HydroPol2D-Python applies the surface stage, normal-flow and critical-flow
% edges to the Neal channel node of each boundary cell (voronoi.py,
% _channel_boundary_step). A MATLAB mesh needed an explicit
% forcing.channel_boundary, and without one the node was closed (audit M21).
% hp2d_voronoi_auto_channel_boundary now builds it from the surface boundary
% when the case gives none.
%   A. The helper on a 1 x 3 mesh with a channel node in every cell: a
%      critical edge on cell 3 and a stage edge (0.4 m above the invert) on
%      cell 1 give node 3 critical flow and node 1 a stage 0.4 m above its
%      channel bed, with the widest link width and the mean link roughness.
%   B. The same mesh with 0.5 m of water in the channel and a critical
%      outlet on cell 3's east edge: the channel drains through node 3 and
%      the ledger closes. Without the outlet nothing leaves.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

dx = 10; channel_bed = [2.2; 2.1; 2.0]; bank = channel_bed + 1;
mesh_file = fullfile(out_dir, 'channel-mesh-1x3.nc');
write_square_mesh(mesh_file, 1, 3, dx, bank);
write_vector(mesh_file, 'channel_node_x', 'channel_node', [5; 15; 25]);
write_vector(mesh_file, 'channel_node_y', 'channel_node', [5; 5; 5]);
write_vector(mesh_file, 'channel_node_host_cell', 'channel_node', [0; 1; 2]);
write_vector(mesh_file, 'channel_node_bed_elevation_m', 'channel_node', channel_bed);
write_vector(mesh_file, 'channel_node_bank_elevation_m', 'channel_node', bank);
write_vector(mesh_file, 'channel_node_plan_area_m2', 'channel_node', [20; 20; 20]);
write_vector(mesh_file, 'channel_link_upstream_node', 'channel_link', [0; 1]);
write_vector(mesh_file, 'channel_link_downstream_node', 'channel_link', [1; 2]);
write_vector(mesh_file, 'channel_link_length_m', 'channel_link', [10; 10]);
write_vector(mesh_file, 'channel_link_width_m', 'channel_link', [3; 5]);
write_vector(mesh_file, 'channel_link_bankfull_depth_m', 'channel_link', [1; 1]);
write_vector(mesh_file, 'channel_link_roughness', 'channel_link', [0.03; 0.05]);
write_vector(mesh_file, 'channel_link_reach_id', 'channel_link', [1; 1]);
mesh = HydroPol2D_Read_UGRID(mesh_file, allow_legacy=true);
east = find(mesh.edge_neighbor == 0 & abs(mesh.edge_midpoint_x - 3 * dx) < 1e-9);
west = find(mesh.edge_neighbor == 0 & abs(mesh.edge_midpoint_x) < 1e-9);

% A. the helper
surface = struct('edge_id', [west; east], 'type', ["stage"; "critical_flow"], 'value', [bank(1) + 0.4; 0], ...
    'stage_depth_m', [0.4; NaN]);
auto = hp2d_voronoi_auto_channel_boundary(mesh, surface);
helper_ok = isequal(auto.node_id(:), [1; 3]) && auto.type(1) == "stage" && ...
    abs(auto.value(1) - (channel_bed(1) + 0.4)) < 1e-12 && auto.type(2) == "critical_flow" && ...
    abs(auto.width_m(2) - 5) < 1e-12 && abs(auto.roughness(1) - 0.03) < 1e-12 && ...
    abs(auto.roughness(2) - 0.05) < 1e-12 && abs(auto.length_m(1) - sqrt(20)) < 1e-12;

% B. a run
config = struct('duration_s', 600, 'min_dt_s', 0.01, 'max_dt_s', 5, 'output_interval_s', 600, ...
    'surface_roughness', 0.03, 'courant', 0.5, 'critical_flow', true, 'hydrology_enabled', false, ...
    'groundwater_enabled', false, 'minimum_cell_width_m', dx, 'allow_legacy_mesh', true, ...
    'initial_channel_depth_m', 0.5);
outlet = struct('edge_id', east, 'type', "critical_flow", 'value', 0);
open_run = HydroPol2D_Voronoi_Run(mesh_file, config, struct('surface_boundary', outlet));
closed_run = HydroPol2D_Voronoi_Run(mesh_file, config, struct());
initial = 0.5 * 60;
open_out = -sum([open_run.diagnostics.boundary_net_inflow_volume_m3]);
open_left = sum(open_run.final_channel_volume_m3) + sum(open_run.final_surface_volume_m3);
closed_left = sum(closed_run.final_channel_volume_m3) + sum(closed_run.final_surface_volume_m3);
open_residual = abs(initial - open_out - open_left) / initial;
closed_residual = abs(initial - closed_left) / initial;
passed = helper_ok && open_out > 0.05 * initial && open_residual < 1e-9 && closed_residual < 1e-9;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-CHANBC-001", "Channel nodes take their boundary cell's outlet", ...
    helper_ok, initial, open_out, open_residual, closed_residual, status, true, ...
    'VariableNames', {'case_id','case_name','helper_ok','initial_m3','channel_outflow_m3', ...
    'open_residual_fraction','closed_residual_fraction','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-CHANBC-001 failed.');
end
fprintf('VAL-MESH-CHANBC-001 passed.\n');

function write_square_mesh(path, ny, nx, dx, bed)
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
write_vector(path, 'cell_bed_elevation_m', 'cell', bed(:));
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
