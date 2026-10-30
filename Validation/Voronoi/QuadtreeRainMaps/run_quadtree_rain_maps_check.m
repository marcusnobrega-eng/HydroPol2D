% VAL-MESH-RAINMAP-001: rain maps on a quadtree mesh.
%
% HydroPol2D-Python samples each rain map on every HRU piece (one mesh cell x
% raster cell overlap) and stops the rain one interval after the last map
% (voronoi.py). HydroPol2D_Run_Prepared_Quadtree took the nearest map value
% at each mesh cell's centre and held the last map to the end of the run
% (audit M16). hp2d_quadtree_rainfall_rate now samples the map at every
% land-use raster cell and averages it over each mesh cell with the
% overlap weights; after the last map the rate is 0.
% A 4 x 6 raster (south-up) under two mesh cells (columns 1-3 and 4-6), one
% raster cell outside the domain, and two 30 min rain maps that vary by
% column.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

rows = 4; cols = 6; d = 0.01;                          % degrees
mapping = struct('x_edges', (0:cols)' * d, 'y_edges', (0:rows)' * d, ...
    'raster_rows', rows, 'raster_cols', cols);
[lon, lat] = hp2d_overlap_raster_centres(mapping);
[col_index, row_index] = meshgrid(1:cols, 1:rows);     % row 1 = south
col_index = reshape(col_index.', [], 1); row_index = reshape(row_index.', [], 1);
host = 1 + (col_index > 3);
raster_to_mesh = sparse(host, (1:rows * cols)', 1 / (rows * 3), 2, rows * cols);
mapping.raster_to_mesh = raster_to_mesh;
valid = true(rows * cols, 1); valid(1) = false;        % south-west cell outside
order_ok = all(abs(lat - (row_index - 0.5) * d) < 1e-12) && all(abs(lon - (col_index - 0.5) * d) < 1e-12);

forcing_root = tempname; mkdir(fullfile(forcing_root, 'Rainfall'));
R = georefcells([0 rows * d], [0 cols * d], [rows cols]);  % north-up raster
map1 = repmat(10 * (1:cols), rows, 1);                 % mm/h by column
map2 = 2 * map1;
geotiffwrite(fullfile(forcing_root, 'Rainfall', 'rain_0000.tif'), map1, R);
geotiffwrite(fullfile(forcing_root, 'Rainfall', 'rain_0001.tif'), map2, R);
files = dir(fullfile(forcing_root, 'Rainfall', '*.tif'));
[~, sort_index] = sort({files.name}); files = files(sort_index);
clear hp2d_quadtree_rainfall_rate
rate1 = hp2d_quadtree_rainfall_rate(600, files, forcing_root, 1800, lat, lon, mapping, valid);
rate2 = hp2d_quadtree_rainfall_rate(2400, files, forcing_root, 1800, lat, lon, mapping, valid);
rate3 = hp2d_quadtree_rainfall_rate(3700, files, forcing_root, 1800, lat, lon, mapping, valid);
rmdir(forcing_root, 's');

column_rain = 10 * col_index;
expected1 = [mean(column_rain(host == 1 & valid)); mean(column_rain(host == 2 & valid))] / 3.6e6;
error1 = max(abs(rate1 - expected1)); error2 = max(abs(rate2 - 2 * expected1));
after_last = max(abs(rate3));
passed = order_ok && error1 < 1e-15 && error2 < 1e-15 && after_last == 0;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-RAINMAP-001", "Rain maps on a quadtree mesh", order_ok, ...
    rate1(1) * 3.6e6, expected1(1) * 3.6e6, error1, error2, after_last, status, true, ...
    'VariableNames', {'case_id','case_name','raster_order_ok','cell1_rain_mm_h', ...
    'expected_cell1_mm_h','map1_error_m_s','map2_error_m_s','rate_after_last_map', ...
    'status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-RAINMAP-001 failed.');
end
fprintf('VAL-MESH-RAINMAP-001 passed.\n');
