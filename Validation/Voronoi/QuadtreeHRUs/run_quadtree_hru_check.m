% VAL-MESH-HRU-001: quadtree cells keep their land-use and soil mix.
%
% HydroPol2D-Python gives every mesh-cell x raster-cell overlap its own HRU
% (voronoi.py). HydroPol2D_Run_Prepared_Quadtree used to give each mesh cell
% its majority land-use and soil class and a 0/1 pervious fraction, so a
% coarse cell that is 40 % roof and 60 % lawn infiltrated like a lawn.
% hp2d_quadtree_hrus now groups the overlaps of a cell by land-use x soil
% pair: one HRU per pair, its share of the cell, and the area means of the
% continuous rasters.
%
% A 2 x 2 mesh of 40 m cells over an 8 x 8 raster of 10 m cells:
%   A. The HRU shares and means equal the exact counts and means.
%   B. A two-hour storm on flat ground with infiltration and canopy, run
%      three ways: one HRU per overlap (Python's construction), the class
%      HRUs, and the old majority class. The class HRUs reproduce the
%      per-overlap run (the class parameters are uniform within a pair, and
%      LAI and initial moisture here too) to the soil solver's precision;
%      the majority run does not.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

% ---- raster (south-up, row-major) and its overlap with the mesh ---------
nr = 8; nc = 8; dr = 10; nm = 2; dm = 40;
[col, row] = meshgrid(1:nc, 1:nr);                       % row 1 = south
row = reshape(row', [], 1); col = reshape(col', [], 1);
n_raster = nr * nc;
host = (ceil(row / 4) - 1) * nm + ceil(col / 4);         % mesh cell of each raster cell
raster_to_mesh = sparse(host, (1:n_raster)', dr^2 / dm^2, nm * nm, n_raster);
lulc = 30 * ones(n_raster, 1);                            % grass
lulc(mod(row + 2 * col, 5) < 2) = 50;                     % roofs, about 40 %
lulc(host == 4) = 30;                                     % one cell all grass
soil = 2 + (col > 5);
soil(host == 4) = 3;                                      % and one soil
lulc(1) = NaN;                                            % outside the domain
lai = 2.5 * (lulc == 30);
theta = 0.20 + 0.05 * (soil == 3);
hrus = hp2d_quadtree_hrus(raster_to_mesh, lulc, soil, struct('lai', lai, 'theta', theta));

% A. exact shares and means
share_error = 0; mean_error = 0;
for cell = 1:nm * nm
    inside = host == cell & isfinite(lulc);
    pairs = unique([lulc(inside) soil(inside)], 'rows');
    for p = 1:size(pairs, 1)
        members = inside & lulc == pairs(p, 1) & soil == pairs(p, 2);
        slot = find(hrus.lulc(cell, :) == pairs(p, 1) & hrus.soil(cell, :) == pairs(p, 2) & ...
            hrus.fraction(cell, :) > 0);
        assert(numel(slot) == 1, 'Pair missing from the HRUs.');
        share_error = max(share_error, abs(hrus.fraction(cell, slot) - nnz(members) / nnz(inside)));
        mean_error = max(mean_error, abs(hrus.lai(cell, slot) - mean(lai(members))));
    end
    share_error = max(share_error, abs(sum(hrus.fraction(cell, :)) - 1));
end
single_class_cell_hrus = nnz(hrus.fraction(4, :) > 0);

% B. a storm run three ways
mesh_file = fullfile(out_dir, 'square-mesh-2x2.nc');
write_square_mesh(mesh_file, nm, nm, dm);
classes = struct('code', [30 50], 'pervious', [1 0], 'root', [0.6 0.1]);
soils = struct('code', [2 3], 'ksat_m_s', [8e-6 2e-6], 'theta_r', [0.05 0.07], ...
    'theta_sat', [0.45 0.40], 'alpha', [3.6 2.0], 'n_vg', [1.56 1.35]);
active = isfinite(lulc);
piece_fraction = zeros(nm * nm, 16); piece_lulc = 30 * ones(nm * nm, 16);
piece_soil = 2 * ones(nm * nm, 16); piece_lai = zeros(nm * nm, 16); piece_theta = 0.2 * ones(nm * nm, 16);
for cell = 1:nm * nm
    members = find(host == cell & active);
    piece_fraction(cell, 1:numel(members)) = 1 / numel(members);
    piece_lulc(cell, 1:numel(members)) = lulc(members);
    piece_soil(cell, 1:numel(members)) = soil(members);
    piece_lai(cell, 1:numel(members)) = lai(members);
    piece_theta(cell, 1:numel(members)) = theta(members);
end
majority_lulc = zeros(nm * nm, 1); majority_soil = zeros(nm * nm, 1);
majority_lai = zeros(nm * nm, 1); majority_theta = zeros(nm * nm, 1);
for cell = 1:nm * nm
    inside = host == cell & active;
    majority_lulc(cell) = mode(lulc(inside)); majority_soil(cell) = mode(soil(inside));
    majority_lai(cell) = mean(lai(inside)); majority_theta(cell) = mean(theta(inside));
end
per_piece = run_storm(mesh_file, piece_fraction, piece_lulc, piece_soil, piece_lai, piece_theta, classes, soils);
class_hrus = run_storm(mesh_file, hrus.fraction, hrus.lulc, hrus.soil, hrus.lai, hrus.theta, classes, soils);
majority = run_storm(mesh_file, ones(nm * nm, 1), majority_lulc, majority_soil, majority_lai, majority_theta, classes, soils);
class_error = max(abs(class_hrus.infiltration_m3 - per_piece.infiltration_m3)) / max(per_piece.infiltration_m3);
class_depth_error = max(abs(class_hrus.final_depth_m - per_piece.final_depth_m));
majority_error = max(abs(majority.infiltration_m3 - per_piece.infiltration_m3)) / max(per_piece.infiltration_m3);

Cells = table((1:nm * nm)', per_piece.infiltration_m3, class_hrus.infiltration_m3, ...
    majority.infiltration_m3, per_piece.final_depth_m, class_hrus.final_depth_m, majority.final_depth_m, ...
    'VariableNames', {'cell','infiltration_per_overlap_m3','infiltration_class_hru_m3', ...
    'infiltration_majority_m3','depth_per_overlap_m','depth_class_hru_m','depth_majority_m'});
passed = share_error < 1e-12 && mean_error < 1e-12 && single_class_cell_hrus == 1 && ...
    class_error < 1e-6 && class_depth_error < 1e-8 && majority_error > 0.05;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-HRU-001", "Quadtree cells keep their land-use and soil mix", ...
    share_error, mean_error, class_error, class_depth_error, majority_error, status, true, ...
    'VariableNames', {'case_id','case_name','max_share_error','max_mean_error', ...
    'class_hru_infiltration_error','class_hru_depth_error_m', ...
    'majority_infiltration_error','status','report_ready'});
writetable(Cells, fullfile(out_dir, 'Quadtree_HRU_Cells.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(Cells); disp(PassFail);
if ~passed
    error('VAL-MESH-HRU-001 failed.');
end
fprintf('VAL-MESH-HRU-001 passed.\n');

function out = run_storm(mesh_file, fraction, lulc, soil, lai, theta, classes, soils)
lookup = @(table_codes, values, codes) reshape(values(arrayfun(@(c) find(table_codes == c), codes(:))), size(codes));
tr = lookup(soils.code, soils.theta_r, soil); ts = lookup(soils.code, soils.theta_sat, soil);
hydrology = struct('hru_fraction', fraction, ...
    'pervious_fraction', lookup(classes.code, classes.pervious, lulc), ...
    'lai', lai, 'soil_depth_m', 2, 'root_depth_m', lookup(classes.code, classes.root, lulc), ...
    'theta_r', tr, 'theta_sat', ts, 'alpha_vg_1_m', lookup(soils.code, soils.alpha, soil), ...
    'n_vg', lookup(soils.code, soils.n_vg, soil), 'ksat_m_s', lookup(soils.code, soils.ksat_m_s, soil), ...
    'initial_soil_saturation', min(max((theta - tr) ./ (ts - tr), 0), 1), ...
    'capillary_rise_enabled', false);
config = struct('duration_s', 7200, 'min_dt_s', 0.01, 'max_dt_s', 30, ...
    'output_interval_s', 3600, 'surface_roughness', 0.03, 'courant', 0.5, ...
    'critical_flow', true, 'hydrology_enabled', true, 'hydrology', hydrology, ...
    'groundwater_enabled', false, 'minimum_cell_width_m', 40, 'allow_legacy_mesh', true);
forcing = struct('surface_source_m_s', 40e-3 / 3600, 'potential_et_m_s', 0);
results = HydroPol2D_Voronoi_Run(mesh_file, config, forcing);
area = 1600;
out.infiltration_m3 = results.hydrology.cumulative_infiltration_m(:, end) * area;
out.final_depth_m = results.surface_depth_m(:, end);
end

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
