clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

domain_dir = fullfile(repo_root, 'HydroPol2D_Model', 'Validation', 'Phase1_VTilted_Catchment');
dem_path = fullfile(domain_dir, 'Static', 'DEM.tif');

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
map_dir = fullfile(out_dir, 'Maps');
ref_dir = fullfile(case_dir, 'Reference');
ensure_dir(out_dir);
ensure_dir(ts_dir);
ensure_dir(fig_dir);
ensure_dir(map_dir);
ensure_dir(ref_dir);

[t_bench, q_bench] = benchmark_hydrograph();
Benchmark = table(t_bench, q_bench, 'VariableNames', {'Time_min','Benchmark_m3s'});
writetable(Benchmark, fullfile(ref_dir, 'VTilted_Benchmark_Hydrograph.csv'));

dem_raw = double(readgeoraster(dem_path));
dem_raw(~isfinite(dem_raw)) = NaN;
dem_full = flipud(dem_raw);
dx_fine = 20;
n_manning = 0.015;
routing_duration_min = 180;
dt_s = 1.0;
record_dt_s = 60;
Rain = readtable(fullfile(ref_dir, 'Rainfall_Intensity_Data.csv'));
rain_time_min = double(Rain.Time_min(:));
rain_intensity_mm_h = double(Rain.Rainfall_Intensity_mm_h(:));
writetable(Rain, fullfile(ref_dir, 'Rainfall_Intensity_Data.csv'));

channel_cols = find(abs(dem_full(end,:) - min(dem_full(end,:))) < 1e-9);
if isempty(channel_cols)
    channel_cols = round(size(dem_full,2) / 2);
end
roughness_full = vtilted_roughness(size(dem_full), channel_cols);

Full20 = simulate_ordinary_li(dem_full, dx_fine, roughness_full, rain_time_min, ...
    rain_intensity_mm_h, routing_duration_min, dt_s, record_dt_s, channel_cols);
Hydro = score_against_benchmark(Full20.series, t_bench, q_bench);
Hydro.case_id = "P1-HYDRO-LI-VTILT-BENCH-001";
Hydro.case_name = "Canonical V-tilted local-inertial benchmark hydrograph";
Hydro.method = "local_inertial";
Hydro.evidence_type = "benchmark_hydrograph";
Hydro.passed = Hydro.mass_error_pct < 0.1 && Hydro.nse > 0.95 && ...
    Hydro.outlet_volume_error_pct < 2 && Hydro.peak_time_error_min <= 10;

HydroComparison = table(t_bench, q_bench, Hydro.q_model_i, ...
    'VariableNames', {'Time_min','Benchmark_m3s','Model20mLI_m3s'});
writetable(HydroComparison, fullfile(ts_dir, 'Hydrograph_Comparison.csv'));
writetable(Full20.series, fullfile(ts_dir, 'P1-HYDRO-LI-VTILT-BENCH-001_20m.csv'));

% The 20 m DEM has 50 rows by 81 columns. The subgrid test uses an exact
% 3-by-3 aggregation, so crop two upstream rows after flipping the DEM while
% preserving the downstream outlet row.
dem_crop = dem_full(end-47:end, :);
channel_cols_crop = find(abs(dem_crop(end,:) - min(dem_crop(end,:))) < 1e-9);
roughness_crop = vtilted_roughness(size(dem_crop), channel_cols_crop);
Fine20 = simulate_ordinary_li(dem_crop, dx_fine, roughness_crop, rain_time_min, ...
    rain_intensity_mm_h, routing_duration_min, dt_s, record_dt_s, channel_cols_crop);

ratio = 3;
dx_coarse = dx_fine * ratio;
[SubgridTables, coarse_ny, coarse_nx, roughness_coarse] = build_subgrid_tables(dem_crop, ...
    roughness_crop, dx_fine, ratio);
coarse_channel_col = unique(ceil(channel_cols_crop ./ ratio));
coarse_channel_col = coarse_channel_col(coarse_channel_col >= 1 & coarse_channel_col <= coarse_nx);
Sub60 = simulate_subgrid_li(SubgridTables, coarse_ny, coarse_nx, dx_coarse, ...
    roughness_coarse, rain_time_min, rain_intensity_mm_h, routing_duration_min, dt_s, ...
    record_dt_s, coarse_channel_col);

SubScores = score_subgrid(Fine20.series, Sub60.series, Fine20.max_depth_m, ...
    Sub60.max_eta_m, dem_crop, SubgridTables, ratio, dx_fine, dx_coarse);
SubScores.case_id = "P1-SUBGRID-VTILT-001";
SubScores.case_name = "V-tilted 60 m lookup-subgrid local-inertial projection";
SubScores.method = "local_inertial_subgrid";
SubScores.evidence_type = "numerical_reference";
SubScores.passed = SubScores.mass_error_pct < 0.1 && SubScores.nse > 0.90 && ...
    SubScores.outlet_volume_error_pct < 5 && SubScores.max_depth_rmse_m < 0.10 && ...
    abs(SubScores.wet_area_error_pct_001m) < 10;

writetable(Fine20.series, fullfile(ts_dir, 'P1-SUBGRID-VTILT-001_20m_reference.csv'));
writetable(Sub60.series, fullfile(ts_dir, 'P1-SUBGRID-VTILT-001_60m_subgrid.csv'));
writetable(SubScores.hydrograph_comparison, fullfile(ts_dir, 'Subgrid_Hydrograph_Comparison.csv'));
writematrix(flipud(Full20.max_depth_m), fullfile(map_dir, 'P1_HYDRO_LI_VTILT_BENCH_001_max_depth_20m.csv'));
writematrix(flipud(Fine20.max_depth_m), fullfile(map_dir, 'P1_SUBGRID_VTILT_001_reference_max_depth_20m.csv'));
writematrix(flipud(SubScores.projected_depth_m), fullfile(map_dir, 'P1_SUBGRID_VTILT_001_projected_max_depth_20m.csv'));
writematrix(flipud(SubScores.projected_depth_m - Fine20.max_depth_m), ...
    fullfile(map_dir, 'P1_SUBGRID_VTILT_001_max_depth_difference_20m.csv'));

MetricSummary = metrics_table(Hydro, SubScores);
MassBalance = mass_table(Hydro, SubScores);
PassFail = passfail_table(Hydro, SubScores);
MapMetrics = map_metrics_table(SubScores);

writetable(MetricSummary, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(MassBalance, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(MapMetrics, fullfile(out_dir, 'MaxDepth_MapMetrics.csv'));

make_hydrograph_figure(HydroComparison, Full20.series, Hydro, fig_dir);
make_cumulative_figure(Full20.series, Hydro, fig_dir);
make_subgrid_hydrograph_figure(SubScores.hydrograph_comparison, Fine20.series, Sub60.series, SubScores, fig_dir);
make_depth_map_figure(Fine20.max_depth_m, SubScores.projected_depth_m, fig_dir);
make_depth_scatter_figure(Fine20.max_depth_m, SubScores.projected_depth_m, fig_dir);

disp(MetricSummary);
disp(MassBalance);
disp(PassFail);
disp(MapMetrics);

function [t_bench, q_bench] = benchmark_hydrograph()
t_bench = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 200 220 240]';
q_bench = [0.00 0.02 0.35 1.80 3.60 4.60 4.85 4.86 4.86 4.75 ...
           4.20 3.20 2.20 1.45 0.95 0.65 0.45 0.32 0.23 0.13 0.07 0.03]';
end

function Run = simulate_ordinary_li(z, dx, roughness, rain_time_min, rain_intensity_mm_h, ...
    routing_duration_min, dt_s, record_dt_s, outlet_cols)
[ny, nx] = size(z);
h = zeros(ny, nx);
state.outflow = zeros(ny, nx, 3);
state.q_out_m3_s = 0;
Outlet = outlet_down(ny, nx, outlet_cols);
area = dx^2;
n_steps = round(routing_duration_min * 60 / dt_s);
record_every = round(record_dt_s / dt_s);
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
rain_volume_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);
max_depth_m = zeros(ny, nx);
cum_out = 0;
cum_rain = 0;
rec = 1;
record();

for it = 1:n_steps
    rain_mm_h = rainfall_at_time((it - 1) * dt_s / 60, rain_time_min, rain_intensity_mm_h);
    if rain_mm_h > 0
        rain_m_s = rain_mm_h / 1000 / 3600;
        h = h + rain_m_s * dt_s;
        cum_rain = cum_rain + rain_m_s * dt_s * ny * nx * area;
    end
    [h, state, q_out] = ordinary_li_step(h, z, roughness, dt_s, dx, state, Outlet);
    cum_out = cum_out + q_out * dt_s;
    max_depth_m = max(max_depth_m, h);
    if mod(it, record_every) == 0
        rec = rec + 1;
        record();
    end
end

Run.series = table(t_min, q_out_m3_s, stored_m3, rain_volume_m3, ...
    outlet_volume_m3, mass_error_pct);
Run.max_depth_m = max_depth_m;

    function record()
        t_min(rec) = (rec - 1) * record_dt_s / 60;
        q_out_m3_s(rec) = state.q_out_m3_s;
        stored_m3(rec) = sum(h, 'all', 'omitnan') * area;
        rain_volume_m3(rec) = cum_rain;
        outlet_volume_m3(rec) = cum_out;
        mass_error_pct(rec) = 100 * (stored_m3(rec) + outlet_volume_m3(rec) - ...
            rain_volume_m3(rec)) / max(rain_volume_m3(rec), eps);
    end
end

function [h, state, q_out_m3_s] = ordinary_li_step(h, z, roughness, dt_s, dx, state, Outlet)
idx_nan = ~isfinite(z);
h(idx_nan) = NaN;
cell_area = dx^2;
time_step_min = dt_s / 60;
d_tot = h * 1000;
d_p = d_tot;
roughness_squared = roughness.^2;
River_Width = zeros(size(h));
River_Depth = zeros(size(h));

[~,~,~,~,outlet_flow,d_t,~,state.outflow,~,~,~,~,~,~] = Local_Inertial_Model_D4( ...
    1, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, z, d_tot, d_p, roughness, roughness_squared, cell_area, ...
    time_step_min, dx, Outlet.mask, 1, 0.02, Outlet.row, Outlet.col, ...
    1e-6, state.outflow, idx_nan, 0, 0, [], [], River_Width, River_Depth, ...
    0, 0, 0, 0, cell_area, [], 0, 0, []);

h = max(double(d_t) / 1000, 0);
h(idx_nan) = NaN;
q_out_m3_s = sum(double(outlet_flow), 'all', 'omitnan') / 1000 / 3600 * cell_area;
state.q_out_m3_s = q_out_m3_s;
end

function [SubgridTables, nyc, nxc, coarse_roughness] = build_subgrid_tables(fine_z, fine_roughness, fine_dx, ratio)
[ny_f, nx_f] = size(fine_z);
if mod(ny_f, ratio) ~= 0 || mod(nx_f, ratio) ~= 0
    error('Fine grid must be divisible by the requested subgrid ratio.');
end
nyc = ny_f / ratio;
nxc = nx_f / ratio;
coarse_z = block_reduce(fine_z, ratio, @min);
coarse_roughness = block_reduce(fine_roughness, ratio, @mean);
DEM_raster = struct('Z', fine_z, 'cellsize', fine_dx);
Reference_raster = struct('Z', coarse_z, 'cellsize', fine_dx * ratio);
Roughness_raster = struct('Z', coarse_roughness, 'cellsize', fine_dx * ratio);
[SubgridTables, ~] = Subgrid_Properties_Lookup(DEM_raster, Roughness_raster, ...
    Reference_raster, fine_dx * ratio, 'dz', 0.02, 'maxDepth', 3.0, 'verbose', false);
end

function Run = simulate_subgrid_li(SubgridTables, ny, nx, dx, roughness, rain_time_min, ...
    rain_intensity_mm_h, routing_duration_min, dt_s, record_dt_s, outlet_cols)
state.d_tot = zeros(ny, nx);
state.d_p = state.d_tot;
state.outflow = zeros(ny, nx, 3);
state.C_a = dx^2;
state.q_out_m3_s = 0;
z = SubgridTables.invert_el;
Outlet = outlet_down(ny, nx, outlet_cols);
n_steps = round(routing_duration_min * 60 / dt_s);
record_every = round(record_dt_s / dt_s);
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
rain_volume_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);
max_eta_m = z;
cum_out = 0;
cum_rain = 0;
rec = 1;
record();

for it = 1:n_steps
    rain_mm_h = rainfall_at_time((it - 1) * dt_s / 60, rain_time_min, rain_intensity_mm_h);
    if rain_mm_h > 0
        rain_step_mm = rain_mm_h / 3600 * dt_s;
        drep_new = hp2d_subgrid_apply_volume_change(state.d_tot ./ 1000, ...
            rain_step_mm, SubgridTables, dx^2);
        state.d_tot = 1000 .* drep_new;
        cum_rain = cum_rain + rain_step_mm / 1000 * ny * nx * dx^2;
    end
    [state, q_out] = subgrid_li_step(state, z, roughness, dt_s, dx, Outlet, SubgridTables);
    cum_out = cum_out + q_out * dt_s;
    max_eta_m = max(max_eta_m, z + state.d_tot ./ 1000);
    if mod(it, record_every) == 0
        rec = rec + 1;
        record();
    end
end

Run.series = table(t_min, q_out_m3_s, stored_m3, rain_volume_m3, ...
    outlet_volume_m3, mass_error_pct);
Run.max_eta_m = max_eta_m;

    function record()
        t_min(rec) = (rec - 1) * record_dt_s / 60;
        q_out_m3_s(rec) = state.q_out_m3_s;
        stored_m3(rec) = sum(hp2d_subgrid_lookup_depth(SubgridTables.volume_cell, ...
            max(state.d_tot ./ 1000, 0), SubgridTables.dz, SubgridTables.maxDepth), ...
            'all', 'omitnan');
        rain_volume_m3(rec) = cum_rain;
        outlet_volume_m3(rec) = cum_out;
        mass_error_pct(rec) = 100 * (stored_m3(rec) + outlet_volume_m3(rec) - ...
            rain_volume_m3(rec)) / max(rain_volume_m3(rec), eps);
    end
end

function [state, q_out_m3_s] = subgrid_li_step(state, z, roughness, dt_s, dx, Outlet, SubgridTables)
idx_nan = false(size(z));
cell_area = dx^2;
time_step_min = dt_s / 60;
River_Width = zeros(size(z));
River_Depth = zeros(size(z));

[~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = Local_Inertial_Model_D4( ...
    2, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, z, state.d_tot, state.d_p, roughness, roughness.^2, cell_area, ...
    time_step_min, dx, Outlet.mask, 1, 0.02, Outlet.row, Outlet.col, 1e-6, ...
    state.outflow, idx_nan, 0, 1, roughness, [], River_Width, River_Depth, ...
    [], [], [], [], state.C_a, [], 0, 0, SubgridTables);

state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.C_a = C_a;
q_out_m3_s = sum(double(outlet_flow), 'all', 'omitnan') / 1000 / 3600 * cell_area;
state.q_out_m3_s = q_out_m3_s;
end

function Outlet = outlet_down(ny, nx, outlet_cols)
Outlet.mask = false(ny, nx);
outlet_cols = outlet_cols(outlet_cols >= 1 & outlet_cols <= nx);
Outlet.mask(ny, outlet_cols) = true;
[Outlet.row, Outlet.col] = find(Outlet.mask);
end

function S = score_against_benchmark(Series, t_bench, q_bench)
q_model_i = interp1(Series.t_min, Series.q_out_m3_s, t_bench, 'linear', NaN);
valid = isfinite(q_model_i) & isfinite(q_bench);
err = q_model_i(valid) - q_bench(valid);
S.rmse = sqrt(mean(err.^2, 'omitnan'));
S.mae = mean(abs(err), 'omitnan');
S.max_error = max(abs(err), [], 'omitnan');
S.relative_l2 = norm(err) / max(norm(q_bench(valid)), eps);
S.nse = 1 - sum(err.^2, 'omitnan') / max(sum((q_bench(valid) - mean(q_bench(valid), 'omitnan')).^2, 'omitnan'), eps);
[tp_m, qp_m] = peak_time_q(t_bench(valid), q_model_i(valid));
[tp_b, qp_b] = peak_time_q(t_bench(valid), q_bench(valid));
S.peak_time_error_min = abs(tp_m - tp_b);
S.peak_magnitude_error_pct = abs(qp_m - qp_b) / max(qp_b, eps) * 100;
S.q_model_i = q_model_i;
S.benchmark_outlet_volume_m3 = trapz(t_bench * 60, q_bench);
S.model_outlet_volume_m3 = trapz(Series.t_min * 60, Series.q_out_m3_s);
S.rain_volume_m3 = Series.rain_volume_m3(end);
S.final_stored_volume_m3 = Series.stored_m3(end);
S.mass_residual_m3 = S.final_stored_volume_m3 + S.model_outlet_volume_m3 - S.rain_volume_m3;
S.mass_error_pct = abs(100 * S.mass_residual_m3 / max(S.rain_volume_m3, eps));
S.outlet_volume_error_pct = abs(S.model_outlet_volume_m3 - S.benchmark_outlet_volume_m3) / ...
    max(S.benchmark_outlet_volume_m3, eps) * 100;
end

function S = score_subgrid(FineSeries, SubSeries, fine_max_depth, sub_max_eta, fine_z, SubgridTables, ratio, dx_fine, dx_coarse)
q_sub_i = interp1(SubSeries.t_min, SubSeries.q_out_m3_s, FineSeries.t_min, 'linear', NaN);
valid = isfinite(q_sub_i) & isfinite(FineSeries.q_out_m3_s);
err = q_sub_i(valid) - FineSeries.q_out_m3_s(valid);
S.rmse = sqrt(mean(err.^2, 'omitnan'));
S.mae = mean(abs(err), 'omitnan');
S.max_error = max(abs(err), [], 'omitnan');
S.relative_l2 = norm(err) / max(norm(FineSeries.q_out_m3_s(valid)), eps);
S.nse = 1 - sum(err.^2, 'omitnan') / max(sum((FineSeries.q_out_m3_s(valid) - mean(FineSeries.q_out_m3_s(valid), 'omitnan')).^2, 'omitnan'), eps);
[tp_s, qp_s] = peak_time_q(FineSeries.t_min(valid), q_sub_i(valid));
[tp_f, qp_f] = peak_time_q(FineSeries.t_min(valid), FineSeries.q_out_m3_s(valid));
S.peak_time_error_min = abs(tp_s - tp_f);
S.peak_magnitude_error_pct = abs(qp_s - qp_f) / max(qp_f, eps) * 100;
S.hydrograph_comparison = table(FineSeries.t_min, FineSeries.q_out_m3_s, q_sub_i, ...
    'VariableNames', {'Time_min','Reference20mLI_m3s','Subgrid60mLI_m3s'});

S.reference_outlet_volume_m3 = trapz(FineSeries.t_min * 60, FineSeries.q_out_m3_s);
S.subgrid_outlet_volume_m3 = trapz(SubSeries.t_min * 60, SubSeries.q_out_m3_s);
S.outlet_volume_error_pct = abs(S.subgrid_outlet_volume_m3 - S.reference_outlet_volume_m3) / ...
    max(S.reference_outlet_volume_m3, eps) * 100;
S.rain_volume_m3 = SubSeries.rain_volume_m3(end);
S.final_stored_volume_m3 = SubSeries.stored_m3(end);
S.mass_residual_m3 = S.final_stored_volume_m3 + SubSeries.outlet_volume_m3(end) - S.rain_volume_m3;
S.mass_error_pct = abs(100 * S.mass_residual_m3 / max(S.rain_volume_m3, eps));

S.projected_depth_m = project_subgrid_depth(sub_max_eta, fine_z, ratio);
dref = fine_max_depth;
dsub = S.projected_depth_m;
map_err = dsub - dref;
S.max_depth_rmse_m = sqrt(mean(map_err(:).^2, 'omitnan'));
S.max_depth_mae_m = mean(abs(map_err(:)), 'omitnan');
S.max_depth_bias_m = mean(map_err(:), 'omitnan');
S.max_depth_max_error_m = max(abs(map_err(:)), [], 'omitnan');
S.reference_depth_volume_m3 = sum(dref, 'all', 'omitnan') * dx_fine^2;
S.projected_depth_volume_m3 = sum(dsub, 'all', 'omitnan') * dx_fine^2;
S.depth_volume_error_pct = 100 * (S.projected_depth_volume_m3 - S.reference_depth_volume_m3) / ...
    max(S.reference_depth_volume_m3, eps);
[S.wet_area_error_pct_001m, S.csi_001m] = wet_metrics(dref, dsub, 0.01, dx_fine);
[S.wet_area_error_pct_010m, S.csi_010m] = wet_metrics(dref, dsub, 0.10, dx_fine);
S.coarse_storage_volume_m3 = sum(hp2d_subgrid_lookup_depth(SubgridTables.volume_cell, ...
    max(sub_max_eta - SubgridTables.invert_el, 0), SubgridTables.dz, SubgridTables.maxDepth), ...
    'all', 'omitnan');
S.dx_coarse = dx_coarse;
end

function d_proj = project_subgrid_depth(max_eta, fine_z, ratio)
[nyc, nxc] = size(max_eta);
d_proj = zeros(size(fine_z));
for i = 1:nyc
    rows = (i - 1) * ratio + (1:ratio);
    for j = 1:nxc
        cols = (j - 1) * ratio + (1:ratio);
        d_proj(rows, cols) = max(max_eta(i,j) - fine_z(rows, cols), 0);
    end
end
end

function [wet_area_error_pct, csi] = wet_metrics(ref, pred, threshold, dx)
ref_wet = ref > threshold;
pred_wet = pred > threshold;
area_ref = nnz(ref_wet) * dx^2;
area_pred = nnz(pred_wet) * dx^2;
wet_area_error_pct = 100 * (area_pred - area_ref) / max(area_ref, eps);
tp = nnz(ref_wet & pred_wet);
fp = nnz(~ref_wet & pred_wet);
fn = nnz(ref_wet & ~pred_wet);
csi = tp / max(tp + fp + fn, eps);
end

function A = block_reduce(Z, ratio, fun)
[ny, nx] = size(Z);
A = zeros(ny / ratio, nx / ratio);
for i = 1:size(A,1)
    rows = (i - 1) * ratio + (1:ratio);
    for j = 1:size(A,2)
        cols = (j - 1) * ratio + (1:ratio);
        block = Z(rows, cols);
        A(i,j) = fun(block(:));
    end
end
end

function roughness = vtilted_roughness(sz, channel_cols)
roughness = 0.015 * ones(sz);
roughness(:, channel_cols) = 0.15;
end

function rain_mm_h = rainfall_at_time(t_min, rain_time_min, rain_intensity_mm_h)
if t_min < min(rain_time_min) || t_min > max(rain_time_min)
    rain_mm_h = 0;
    return;
end
idx = find(rain_time_min <= t_min, 1, 'last');
if isempty(idx)
    rain_mm_h = 0;
else
    rain_mm_h = rain_intensity_mm_h(idx);
end
end

function T = metrics_table(H, S)
case_id = [H.case_id; S.case_id];
case_name = [H.case_name; S.case_name];
method = [H.method; S.method];
evidence_type = [H.evidence_type; S.evidence_type];
rmse = [H.rmse; S.rmse];
mae = [H.mae; S.mae];
max_error = [H.max_error; S.max_error];
relative_l2 = [H.relative_l2; S.relative_l2];
nse = [H.nse; S.nse];
peak_time_error_min = [H.peak_time_error_min; S.peak_time_error_min];
peak_magnitude_error_pct = [H.peak_magnitude_error_pct; S.peak_magnitude_error_pct];
mass_error_pct = [H.mass_error_pct; S.mass_error_pct];
outlet_volume_error_pct = [H.outlet_volume_error_pct; S.outlet_volume_error_pct];
max_depth_rmse_m = [NaN; S.max_depth_rmse_m];
wet_area_error_pct_001m = [NaN; S.wet_area_error_pct_001m];
passed = [H.passed; S.passed];
T = table(case_id, case_name, method, evidence_type, rmse, mae, max_error, ...
    relative_l2, nse, peak_time_error_min, peak_magnitude_error_pct, ...
    mass_error_pct, outlet_volume_error_pct, max_depth_rmse_m, ...
    wet_area_error_pct_001m, passed);
end

function T = mass_table(H, S)
case_id = [H.case_id; S.case_id];
rain_volume_m3 = [H.rain_volume_m3; S.rain_volume_m3];
reference_outlet_volume_m3 = [H.benchmark_outlet_volume_m3; S.reference_outlet_volume_m3];
model_outlet_volume_m3 = [H.model_outlet_volume_m3; S.subgrid_outlet_volume_m3];
final_stored_volume_m3 = [H.final_stored_volume_m3; S.final_stored_volume_m3];
mass_residual_m3 = [H.mass_residual_m3; S.mass_residual_m3];
mass_error_pct = [H.mass_error_pct; S.mass_error_pct];
T = table(case_id, rain_volume_m3, reference_outlet_volume_m3, ...
    model_outlet_volume_m3, final_stored_volume_m3, mass_residual_m3, mass_error_pct);
end

function T = passfail_table(H, S)
case_id = [H.case_id; S.case_id];
case_name = [H.case_name; S.case_name];
method = [H.method; S.method];
status = strings(2,1);
status(1) = ternary(H.passed, "pass", "fail");
status(2) = ternary(S.passed, "pass", "fail");
report_ready = [H.passed; S.passed];
T = table(case_id, case_name, method, status, report_ready);
end

function T = map_metrics_table(S)
case_id = S.case_id;
max_depth_rmse_m = S.max_depth_rmse_m;
max_depth_mae_m = S.max_depth_mae_m;
max_depth_bias_m = S.max_depth_bias_m;
max_depth_max_error_m = S.max_depth_max_error_m;
reference_depth_volume_m3 = S.reference_depth_volume_m3;
projected_depth_volume_m3 = S.projected_depth_volume_m3;
depth_volume_error_pct = S.depth_volume_error_pct;
wet_area_error_pct_001m = S.wet_area_error_pct_001m;
csi_001m = S.csi_001m;
wet_area_error_pct_010m = S.wet_area_error_pct_010m;
csi_010m = S.csi_010m;
coarse_storage_volume_m3 = S.coarse_storage_volume_m3;
T = table(case_id, max_depth_rmse_m, max_depth_mae_m, max_depth_bias_m, ...
    max_depth_max_error_m, reference_depth_volume_m3, projected_depth_volume_m3, ...
    depth_volume_error_pct, wet_area_error_pct_001m, csi_001m, ...
    wet_area_error_pct_010m, csi_010m, coarse_storage_volume_m3);
end

function make_hydrograph_figure(HydroComparison, RawSeries, H, fig_dir)
fig = figure('Visible','off','Color','w');
plot(HydroComparison.Time_min, HydroComparison.Benchmark_m3s, 'ko-', 'LineWidth', 1.6); hold on;
plot(HydroComparison.Time_min, HydroComparison.Model20mLI_m3s, 'rs-', 'LineWidth', 1.4);
plot(RawSeries.t_min, RawSeries.q_out_m3_s, 'r:', 'LineWidth', 1.0);
grid on; box on;
xlabel('Time (min)');
ylabel('Outlet discharge (m^3/s)');
title(sprintf('Canonical V-tilted hydrograph, NSE = %.3f', H.nse));
legend('Benchmark', '20 m LI at benchmark times', '20 m LI raw', 'Location', 'best');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_LI_VTILT_BENCH_001_HYDROGRAPH.png'), 'Resolution', 250);
close(fig);
end

function make_cumulative_figure(Series, H, fig_dir)
fig = figure('Visible','off','Color','w');
plot(Series.t_min, Series.outlet_volume_m3, 'r-', 'LineWidth', 1.5); hold on;
plot(Series.t_min, Series.rain_volume_m3, 'b--', 'LineWidth', 1.5);
plot(Series.t_min, Series.stored_m3 + Series.outlet_volume_m3, 'k:', 'LineWidth', 1.5);
yline(H.benchmark_outlet_volume_m3, 'Color', [0.2 0.2 0.2], 'LineStyle', '-.');
grid on; box on;
xlabel('Time (min)');
ylabel('Cumulative volume (m^3)');
title(sprintf('V-tilted event volume balance, residual = %.3g%%', H.mass_error_pct));
legend('Outlet volume', 'Rain volume', 'Stored + outlet', 'Benchmark outlet volume', 'Location', 'best');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_LI_VTILT_BENCH_001_VOLUME.png'), 'Resolution', 250);
close(fig);
end

function make_subgrid_hydrograph_figure(Comparison, FineSeries, SubSeries, S, fig_dir)
fig = figure('Visible','off','Color','w');
tiledlayout(2,1);
nexttile;
plot(Comparison.Time_min, Comparison.Reference20mLI_m3s, 'k-', 'LineWidth', 1.5); hold on;
plot(Comparison.Time_min, Comparison.Subgrid60mLI_m3s, 'r--', 'LineWidth', 1.5);
grid on; box on;
xlabel('Time (min)');
ylabel('Outlet discharge (m^3/s)');
title(sprintf('20 m LI vs 60 m lookup-subgrid LI, NSE = %.3f', S.nse));
legend('20 m LI reference', '60 m subgrid LI', 'Location', 'best');
nexttile;
plot(FineSeries.t_min, FineSeries.mass_error_pct, 'k-', 'LineWidth', 1.2); hold on;
plot(SubSeries.t_min, SubSeries.mass_error_pct, 'r--', 'LineWidth', 1.2);
grid on; box on;
xlabel('Time (min)');
ylabel('Mass error (%)');
legend('20 m LI reference', '60 m subgrid LI', 'Location', 'best');
exportgraphics(fig, fullfile(fig_dir, 'P1_SUBGRID_VTILT_001_HYDROGRAPH.png'), 'Resolution', 250);
close(fig);
end

function make_depth_map_figure(ref_depth, pred_depth, fig_dir)
diff_depth = pred_depth - ref_depth;
lim = max([ref_depth(:); pred_depth(:)], [], 'omitnan');
fig = figure('Visible','off','Color','w');
tiledlayout(1,3);
nexttile; imagesc(flipud(ref_depth)); axis image off; colorbar; clim([0 lim]); title('20 m LI max depth');
nexttile; imagesc(flipud(pred_depth)); axis image off; colorbar; clim([0 lim]); title('60 m projected max depth');
nexttile; imagesc(flipud(diff_depth)); axis image off; colorbar; title('Projected - reference');
exportgraphics(fig, fullfile(fig_dir, 'P1_SUBGRID_VTILT_001_MAX_DEPTH_MAPS.png'), 'Resolution', 250);
close(fig);
end

function make_depth_scatter_figure(ref_depth, pred_depth, fig_dir)
fig = figure('Visible','off','Color','w');
scatter(ref_depth(:), pred_depth(:), 8, 'filled', 'MarkerFaceAlpha', 0.25);
hold on;
mx = max([ref_depth(:); pred_depth(:)], [], 'omitnan');
plot([0 mx], [0 mx], 'k--', 'LineWidth', 1.2);
grid on; box on; axis equal;
xlabel('20 m LI max depth (m)');
ylabel('Projected 60 m subgrid max depth (m)');
title('V-tilted max-depth projection comparison');
exportgraphics(fig, fullfile(fig_dir, 'P1_SUBGRID_VTILT_001_DEPTH_SCATTER.png'), 'Resolution', 250);
close(fig);
end

function [tp, qp] = peak_time_q(t, q)
[qp, idx] = max(q);
tp = t(idx);
end

function s = ternary(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end

function ensure_dir(p)
if ~exist(p, 'dir')
    mkdir(p);
end
end
