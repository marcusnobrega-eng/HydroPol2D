% VAL-MESH-FILLBC-001: unassigned mesh boundary edges can get one condition.
%
% HydroPol2D-Python can give every exterior edge that no boundary lists
% (and that is not on an inflow cell) normal flow at a set slope
% (normal_flow_all_unassigned_boundary_edges, voronoi.py). MATLAB mesh
% cases had to list every outlet edge (audit M20). hp2d_voronoi_fill_boundary
% now reads boundary.fill_unassigned_type / fill_unassigned_value.
%   A. On a 3 x 4 mesh: with the fill and nothing listed, all 14 exterior
%      edges become normal flow; with one inflow edge listed, the inflow
%      cell's other exterior edges are left out; without the fill nothing
%      changes.
%   B. Rain on a plane falling to the west: with the fill water leaves
%      through the boundary and the ledger closes; without it the walls
%      keep all of it.

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
exterior = find(mesh.edge_neighbor == 0);
fill = struct('edge_id', zeros(0, 1), 'type', strings(0, 1), 'value', zeros(0, 1), ...
    'fill_unassigned_type', "normal_flow", 'fill_unassigned_value', 0.01);
all_filled = hp2d_voronoi_fill_boundary(mesh, fill);
inflow_edge = exterior(1);
with_inflow = fill; with_inflow.edge_id = inflow_edge; with_inflow.type = "inflow"; with_inflow.value = 0.1;
with_inflow = hp2d_voronoi_fill_boundary(mesh, with_inflow);
inflow_cell = mesh.edge_owner(inflow_edge);
skipped = exterior(mesh.edge_owner(exterior) == inflow_cell & exterior ~= inflow_edge);
plain = struct('edge_id', exterior(2), 'type', "critical_flow", 'value', 0);
unchanged = isequal(hp2d_voronoi_fill_boundary(mesh, plain), plain);
helper_ok = numel(all_filled.edge_id) == numel(exterior) && all(all_filled.type == "normal_flow") && ...
    ~any(ismember(skipped, with_inflow.edge_id)) && ...
    numel(with_inflow.edge_id) == numel(exterior) - numel(skipped) && unchanged;

config = struct('duration_s', 1800, 'min_dt_s', 0.01, 'max_dt_s', 10, 'output_interval_s', 1800, ...
    'surface_roughness', 0.03, 'courant', 0.5, 'critical_flow', true, 'hydrology_enabled', false, ...
    'groundwater_enabled', false, 'minimum_cell_width_m', 10, 'allow_legacy_mesh', true);
rain = 50e-3 / 3600;
open_run = HydroPol2D_Voronoi_Run(mesh_file, config, struct('surface_source_m_s', rain, 'surface_boundary', fill));
walled = HydroPol2D_Voronoi_Run(mesh_file, config, struct('surface_source_m_s', rain));
rain_volume = rain * 1800 * sum(mesh.cell_area);
outflow = -sum([open_run.diagnostics.boundary_net_inflow_volume_m3]);
open_residual = abs(rain_volume - outflow - sum(open_run.final_surface_volume_m3)) / rain_volume;
walled_residual = abs(rain_volume - sum(walled.final_surface_volume_m3)) / rain_volume;
passed = helper_ok && outflow > 0.01 * rain_volume && open_residual < 1e-9 && walled_residual < 1e-9;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-FILLBC-001", "Unassigned mesh boundary edges get one condition", ...
    helper_ok, rain_volume, outflow, open_residual, walled_residual, status, true, ...
    'VariableNames', {'case_id','case_name','helper_ok','rain_m3','outflow_m3', ...
    'open_residual_fraction','walled_residual_fraction','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-FILLBC-001 failed.');
end
fprintf('VAL-MESH-FILLBC-001 passed.\n');

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
