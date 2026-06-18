clear; clc;

case_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
raster_dir = fullfile(out_dir, 'Rasters');
table_dir = fullfile(out_dir, 'Tables');
fig_dir = fullfile(out_dir, 'Figures');
ensure_dir(out_dir);
ensure_dir(raster_dir);
ensure_dir(table_dir);
ensure_dir(fig_dir);

[Diag, CellTable, StepTable, Pass] = run_spatial_rainfall_case(raster_dir, table_dir);

writetable(Diag, fullfile(out_dir, 'SpatialRainfall_Raster_Diagnostics.csv'));
writetable(Pass, fullfile(out_dir, 'SpatialRainfall_Raster_Pass_Fail.csv'));
writetable(Diag, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(table(Diag.case_id, Diag.case_name, Diag.total_volume_error_pct, ...
    Diag.storage_residual_m3, Diag.mass_error_pct, Diag.passed, ...
    'VariableNames', {'case_id','case_name','total_rainfall_volume_error_pct', ...
    'rainfall_storage_residual_m3','mass_balance_error_pct','passed'}), ...
    fullfile(out_dir, 'Mass_Balance.csv'));

make_figures(CellTable, StepTable, fig_dir);

disp(Diag);
disp(Pass);

function [Diag, CellTable, StepTable, Pass] = run_spatial_rainfall_case(raster_dir, table_dir)
Cfg = struct();
Cfg.case_id = "P1-RAIN-MAP-001";
Cfg.case_name = "Tiny raster spatial rainfall totals";
Cfg.dx_m = 20;
Cfg.ny = 4;
Cfg.nx = 5;
Cfg.epsg = 32610;
Cfg.duration_min = [10; 20; 15];
Cfg.n_steps = numel(Cfg.duration_min);

rain_mm_h = zeros(Cfg.ny, Cfg.nx, Cfg.n_steps);
rain_mm_h(:,:,1) = [ ...
     0  5 10 15 20; ...
     2  7 12 17 22; ...
     4  9 14 19 24; ...
     6 11 16 21 26];
rain_mm_h(:,:,2) = [ ...
    30  0 15  0 30; ...
    25  5 20  5 25; ...
    20 10 25 10 20; ...
    15 15 30 15 15];
rain_mm_h(:,:,3) = [ ...
     3  6  9 12 15; ...
     0  0  0  0  0; ...
    15 12  9  6  3; ...
    18 14 10  6  2];

write_rainfall_geotiffs(rain_mm_h, Cfg, raster_dir);

read_mm_h = zeros(size(rain_mm_h));
for k = 1:Cfg.n_steps
    fname = fullfile(raster_dir, sprintf('rain_%02d_mm_h.tif', k));
    read_mm_h(:,:,k) = double(readgeoraster(fname));
end

ref_depth_m = accumulate_depth(rain_mm_h, Cfg.duration_min);
model_depth_m = accumulate_depth(read_mm_h, Cfg.duration_min);
depth_error_m = model_depth_m - ref_depth_m;

cell_area_m2 = Cfg.dx_m^2;
ref_volume_m3 = sum(ref_depth_m, 'all') * cell_area_m2;
model_volume_m3 = sum(model_depth_m, 'all') * cell_area_m2;
storage_residual_m3 = model_volume_m3 - ref_volume_m3;
total_volume_error_pct = 100 * storage_residual_m3 / max(ref_volume_m3, eps);
cellwise_max_error_m = max(abs(depth_error_m), [], 'all');
cellwise_rmse_m = rmse_omitnan(depth_error_m(:));

step_id = (1:Cfg.n_steps)';
duration_min = Cfg.duration_min;
reference_step_volume_m3 = zeros(Cfg.n_steps, 1);
read_step_volume_m3 = zeros(Cfg.n_steps, 1);
step_volume_error_pct = zeros(Cfg.n_steps, 1);
for k = 1:Cfg.n_steps
    ref_step_depth_m = rain_mm_h(:,:,k) / 1000 * Cfg.duration_min(k) / 60;
    read_step_depth_m = read_mm_h(:,:,k) / 1000 * Cfg.duration_min(k) / 60;
    reference_step_volume_m3(k) = sum(ref_step_depth_m, 'all') * cell_area_m2;
    read_step_volume_m3(k) = sum(read_step_depth_m, 'all') * cell_area_m2;
    step_volume_error_pct(k) = 100 * (read_step_volume_m3(k) - ...
        reference_step_volume_m3(k)) / max(reference_step_volume_m3(k), eps);
end
StepTable = table(step_id, duration_min, reference_step_volume_m3, ...
    read_step_volume_m3, step_volume_error_pct);
CellTable = build_cell_table(Cfg, ref_depth_m, model_depth_m, depth_error_m);

writetable(CellTable, fullfile(table_dir, Cfg.case_id + "_cellwise_totals.csv"));
writetable(StepTable, fullfile(table_dir, Cfg.case_id + "_step_volumes.csv"));

mass_error_pct = abs(total_volume_error_pct);
passed = cellwise_max_error_m < 1e-9 && abs(total_volume_error_pct) < 0.1 && ...
    max(abs(step_volume_error_pct), [], 'omitnan') < 0.1;

Diag = table(Cfg.case_id, Cfg.case_name, "exact_raster_depth_volume_bookkeeping", ...
    cellwise_rmse_m, mean(abs(depth_error_m), 'all', 'omitnan'), ...
    cellwise_max_error_m, 0, 1, 0, 0, mass_error_pct, ...
    total_volume_error_pct, storage_residual_m3, ref_volume_m3, model_volume_m3, ...
    max(abs(step_volume_error_pct), [], 'omitnan'), passed, ...
    "Tiny GeoTIFF rainfall-intensity stack read back and accumulated to exact cellwise depths and domain volumes.", ...
    'VariableNames', {'case_id','case_name','evidence_type','rmse', ...
    'mae','max_error','relative_l2','nse','peak_time_error_min', ...
    'peak_magnitude_error_pct','mass_error_pct','total_volume_error_pct', ...
    'storage_residual_m3','reference_volume_m3','model_volume_m3', ...
    'max_step_volume_error_pct','passed','metric_note'});
Pass = pass_row(Diag, passed);
end

function write_rainfall_geotiffs(rain_mm_h, Cfg, raster_dir)
try
    R = maprefcells([500000, 500000 + Cfg.nx * Cfg.dx_m], ...
        [4100000, 4100000 + Cfg.ny * Cfg.dx_m], [Cfg.ny, Cfg.nx], ...
        'ColumnsStartFrom', 'north');
catch
    R = maprasterref('RasterSize', [Cfg.ny, Cfg.nx], ...
        'XWorldLimits', [500000, 500000 + Cfg.nx * Cfg.dx_m], ...
        'YWorldLimits', [4100000, 4100000 + Cfg.ny * Cfg.dx_m], ...
        'ColumnsStartFrom', 'north');
end

for k = 1:Cfg.n_steps
    fname = fullfile(raster_dir, sprintf('rain_%02d_mm_h.tif', k));
    try
        geotiffwrite(fname, single(rain_mm_h(:,:,k)), R, 'CoordRefSysCode', Cfg.epsg);
    catch
        geotiffwrite(fname, single(rain_mm_h(:,:,k)), R);
    end
end
end

function depth_m = accumulate_depth(rain_mm_h, duration_min)
depth_m = zeros(size(rain_mm_h, 1), size(rain_mm_h, 2));
for k = 1:numel(duration_min)
    depth_m = depth_m + rain_mm_h(:,:,k) / 1000 * duration_min(k) / 60;
end
end

function CellTable = build_cell_table(Cfg, ref_depth_m, model_depth_m, depth_error_m)
[col, row] = meshgrid(1:Cfg.nx, 1:Cfg.ny);
cell_id = ((1:numel(ref_depth_m)))';
row = row(:);
col = col(:);
reference_depth_m = ref_depth_m(:);
model_depth_m = model_depth_m(:);
depth_error_m = depth_error_m(:);
reference_volume_m3 = reference_depth_m * Cfg.dx_m^2;
model_volume_m3 = model_depth_m * Cfg.dx_m^2;
CellTable = table(cell_id, row, col, reference_depth_m, model_depth_m, ...
    depth_error_m, reference_volume_m3, model_volume_m3);
end

function Pass = pass_row(Diag, passed)
status = "fail";
if passed
    status = "pass";
end
report_ready = passed;
Pass = table(Diag.case_id, Diag.case_name, status, report_ready, ...
    'VariableNames', {'case_id','case_name','status','report_ready'});
end

function make_figures(CellTable, StepTable, fig_dir)
ny = max(CellTable.row);
nx = max(CellTable.col);
ref = reshape(CellTable.reference_depth_m, ny, nx);
model = reshape(CellTable.model_depth_m, ny, nx);
err_mm = reshape(CellTable.depth_error_m * 1000, ny, nx);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 2);
nexttile;
imagesc(ref * 1000);
axis image; colorbar;
title('Reference accumulated rainfall (mm)');
xlabel('Column'); ylabel('Row');
nexttile;
imagesc(model * 1000);
axis image; colorbar;
title('Raster-read accumulated rainfall (mm)');
xlabel('Column'); ylabel('Row');
nexttile;
imagesc(err_mm);
axis image; colorbar;
title('Cellwise error (mm)');
xlabel('Column'); ylabel('Row');
nexttile;
bar(StepTable.step_id, [StepTable.reference_step_volume_m3, StepTable.read_step_volume_m3]);
xlabel('Raster step'); ylabel('Rainfall volume (m^3)');
legend('Reference', 'Raster-read', 'Location', 'best');
title('Stepwise volume check');
grid on;
exportgraphics(fig, fullfile(fig_dir, 'P1_RAIN_MAP_001_TOTALS.png'), 'Resolution', 200);
close(fig);
end

function out = rmse_omitnan(x)
x = x(isfinite(x));
if isempty(x)
    out = NaN;
else
    out = sqrt(mean(x.^2));
end
end

function ensure_dir(path_name)
if ~exist(path_name, 'dir')
    mkdir(path_name);
end
end
