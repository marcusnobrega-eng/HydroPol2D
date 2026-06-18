clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
ensure_dir(out_dir);
ensure_dir(ts_dir);
ensure_dir(fig_dir);

[Diag, Series, Pass] = run_inflow_case(ts_dir);

writetable(Diag, fullfile(out_dir, 'InflowHydrograph_Boundary_Diagnostics.csv'));
writetable(Pass, fullfile(out_dir, 'InflowHydrograph_Boundary_Pass_Fail.csv'));
writetable(Diag, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(table(Diag.case_id, Diag.case_name, ...
    Diag.imposed_volume_error_pct, Diag.mass_error_pct, ...
    Series.boundary_volume_m3(end), Series.outlet_volume_m3(end), ...
    Series.stored_m3(end), Diag.passed, ...
    'VariableNames', {'case_id','case_name','imposed_volume_error_pct', ...
    'mass_balance_error_pct','imposed_volume_m3','outlet_volume_m3', ...
    'final_storage_m3','passed'}), fullfile(out_dir, 'Mass_Balance.csv'));

make_figures(Series, fig_dir);

disp(Diag);
disp(Pass);

function [Diag, Series, Pass] = run_inflow_case(ts_dir)
Cfg = struct();
Cfg.case_id = "P1-BC-INFLOW-001";
Cfg.case_name = "Full momentum inflow hydrograph boundary";
Cfg.dx = 10;
Cfg.nx = 36;
Cfg.ny = 5;
Cfg.slope = 0.003;
Cfg.n_manning = 0.03;
Cfg.dt_s = 1.0;
Cfg.duration_s = 3 * 3600;
Cfg.record_dt_s = 60;
Cfg.start_s = 5 * 60;
Cfg.peak_s = 25 * 60;
Cfg.end_s = 45 * 60;
Cfg.peak_q_m3_s = 8.0;

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = repmat(max(x) * Cfg.slope - x * Cfg.slope, Cfg.ny, 1);
h = zeros(Cfg.ny, Cfg.nx);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
outflow = zeros(Cfg.ny, Cfg.nx, 5);
Outlet = outlet_faces(Cfg.ny, Cfg.nx, 'right', true(Cfg.ny, 1));
area = Cfg.dx^2;

n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
q_in_m3_s = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);
reference_cumulative_inflow_m3 = zeros(n_records, 1);
boundary_volume_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
mass_residual_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);
imposed_volume_error_pct = zeros(n_records, 1);
mean_depth_m = zeros(n_records, 1);
downstream_depth_m = zeros(n_records, 1);

cum_boundary = 0;
cum_out = 0;
rec = 1;
record_state(0, 0);

for it = 1:n_steps
    t_mid_s = (it - 0.5) * Cfg.dt_s;
    q_in = triangular_inflow(t_mid_s, Cfg);
    added_m3 = q_in * Cfg.dt_s;
    h(:, 1) = h(:, 1) + added_m3 / (Cfg.ny * area);
    cum_boundary = cum_boundary + added_m3;

    [h, outflow, q_out] = full_momentum_step(h, z, n, Cfg.dt_s, Cfg.dx, ...
        outflow, Outlet, Cfg.slope);
    cum_out = cum_out + q_out * Cfg.dt_s;

    if mod(it, record_every) == 0
        rec = rec + 1;
        record_state(it * Cfg.dt_s, q_out);
    end
end

Series = table(t_min, q_in_m3_s, q_out_m3_s, ...
    reference_cumulative_inflow_m3, boundary_volume_m3, outlet_volume_m3, ...
    stored_m3, mass_residual_m3, mass_error_pct, imposed_volume_error_pct, ...
    mean_depth_m, downstream_depth_m);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));

imposed_volume_error = abs(boundary_volume_m3(end) - ...
    exact_cumulative_inflow(Cfg.duration_s, Cfg)) / ...
    max(exact_cumulative_inflow(Cfg.duration_s, Cfg), eps) * 100;
mass_abs_pct = max(abs(mass_error_pct), [], 'omitnan');
outlet_response = max(q_out_m3_s, [], 'omitnan');
[t_in_peak, q_in_peak] = peak_time_q(t_min, q_in_m3_s);
[t_out_peak, q_out_peak] = peak_time_q(t_min, q_out_m3_s);
peak_time_error_min = max(t_out_peak - t_in_peak, 0);
peak_q_error_pct = abs(q_out_peak - q_in_peak) / max(q_in_peak, eps) * 100;
final_accounted_pct = abs(outlet_volume_m3(end) + stored_m3(end) - ...
    boundary_volume_m3(end)) / max(boundary_volume_m3(end), eps) * 100;

passed = imposed_volume_error < 0.01 && mass_abs_pct < 0.1 && ...
    final_accounted_pct < 0.1 && outlet_response > 1e-4;

Diag = table(Cfg.case_id, Cfg.case_name, "exact_inflow_volume_bookkeeping", ...
    imposed_volume_error, mean(abs(imposed_volume_error_pct), 'omitnan'), ...
    outlet_response, NaN, NaN, peak_time_error_min, peak_q_error_pct, ...
    mass_abs_pct, final_accounted_pct, passed, ...
    "Known triangular inflow volume imposed into the upstream boundary and routed with full momentum; downstream hydrograph metrics are diagnostic.", ...
    'VariableNames', {'case_id','case_name','evidence_type','rmse', ...
    'mae','max_error','relative_l2','nse','peak_time_error_min', ...
    'peak_magnitude_error_pct','mass_error_pct','imposed_volume_error_pct', ...
    'passed','metric_note'});
Pass = pass_row(Diag, passed);

    function record_state(t_s, q_out)
        q_now = triangular_inflow(t_s, Cfg);
        exact_v = exact_cumulative_inflow(t_s, Cfg);
        t_min(rec) = t_s / 60;
        q_in_m3_s(rec) = q_now;
        q_out_m3_s(rec) = q_out;
        reference_cumulative_inflow_m3(rec) = exact_v;
        boundary_volume_m3(rec) = cum_boundary;
        outlet_volume_m3(rec) = cum_out;
        stored_m3(rec) = sum(h, 'all') * area;
        mass_residual_m3(rec) = outlet_volume_m3(rec) + stored_m3(rec) - boundary_volume_m3(rec);
        mass_error_pct(rec) = 100 * mass_residual_m3(rec) / max(boundary_volume_m3(rec), eps);
        imposed_volume_error_pct(rec) = 100 * (boundary_volume_m3(rec) - exact_v) / max(exact_v, eps);
        mean_depth_m(rec) = mean(h, 'all', 'omitnan');
        downstream_depth_m(rec) = mean(h(:, end), 'omitnan');
    end
end

function q = triangular_inflow(t_s, Cfg)
if t_s <= Cfg.start_s || t_s >= Cfg.end_s
    q = 0;
elseif t_s <= Cfg.peak_s
    q = Cfg.peak_q_m3_s * (t_s - Cfg.start_s) / (Cfg.peak_s - Cfg.start_s);
else
    q = Cfg.peak_q_m3_s * (Cfg.end_s - t_s) / (Cfg.end_s - Cfg.peak_s);
end
end

function v = exact_cumulative_inflow(t_s, Cfg)
t = min(max(t_s, 0), Cfg.end_s);
if t <= Cfg.start_s
    v = 0;
elseif t <= Cfg.peak_s
    v = 0.5 * Cfg.peak_q_m3_s * (t - Cfg.start_s)^2 / ...
        (Cfg.peak_s - Cfg.start_s);
else
    rise_volume = 0.5 * Cfg.peak_q_m3_s * (Cfg.peak_s - Cfg.start_s);
    tau = t - Cfg.peak_s;
    recession_duration = Cfg.end_s - Cfg.peak_s;
    v = rise_volume + Cfg.peak_q_m3_s * tau - ...
        0.5 * Cfg.peak_q_m3_s * tau^2 / recession_duration;
end
end

function [h, outflow, q_out_m3_s] = full_momentum_step(h, z, roughness, dt_s, dx, outflow, Outlet, outlet_slope)
idx_nan = false(size(h));
outlet_index = false(size(h));
roughness_squared = roughness.^2;
cell_area = dx^2;
time_step_min = dt_s / 60;
d_tot = h * 1000;
d_p = d_tot;
d_tolerance_mm = 1e-6;

[~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,~] = Full_Momentum_Model_D4( ...
    1, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, z, d_tot, d_p, roughness, roughness_squared, cell_area, ...
    time_step_min, dx, outlet_index, 1, outlet_slope, [], [], ...
    d_tolerance_mm, outflow, idx_nan, 0, 0, [], [], [], [], ...
    0, 0, 0, 0, 0, [], 0, 0, [], Outlet);

h = max(double(d_t) / 1000, 0);
q_out_m3_s = sum(double(outlet_flow), 'all') / 1000 / 3600 * cell_area;
end

function Outlet = outlet_faces(ny, nx, side, mask_in)
Outlet = struct();
Outlet.face_right = false(ny, nx);
Outlet.face_left = false(ny, nx);
Outlet.face_up = false(ny, nx);
Outlet.face_down = false(ny, nx);

switch lower(side)
    case 'right'
        Outlet.face_right(:, nx) = logical(mask_in(:));
    case 'left'
        Outlet.face_left(:, 1) = logical(mask_in(:));
    case 'down'
        row_mask = false(1, nx);
        row_mask(:) = logical(mask_in(:));
        Outlet.face_down(ny, :) = row_mask;
    case 'up'
        row_mask = false(1, nx);
        row_mask(:) = logical(mask_in(:));
        Outlet.face_up(1, :) = row_mask;
    case 'none'
        return
    otherwise
        error('Unknown outlet side "%s".', side);
end
end

function [t_peak, q_peak] = peak_time_q(t_min, q)
[q_peak, idx] = max(q);
t_peak = t_min(idx);
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

function make_figures(Series, fig_dir)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(3, 1);
nexttile;
plot(Series.t_min, Series.q_in_m3_s, '-', ...
    Series.t_min, Series.q_out_m3_s, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Discharge (m^3/s)');
legend('Prescribed inflow', 'Downstream outlet', 'Location', 'best');
title('Inflow-hydrograph boundary response');
grid on;
nexttile;
plot(Series.t_min, Series.reference_cumulative_inflow_m3, '-', ...
    Series.t_min, Series.boundary_volume_m3, '--', ...
    Series.t_min, Series.outlet_volume_m3, ':', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Volume (m^3)');
legend('Exact inflow integral', 'Imposed boundary volume', 'Outlet volume', ...
    'Location', 'best');
title('Cumulative volumes');
grid on;
nexttile;
yyaxis left
plot(Series.t_min, Series.stored_m3, 'LineWidth', 1.4);
ylabel('Stored volume (m^3)');
yyaxis right
plot(Series.t_min, Series.mass_error_pct, '--', 'LineWidth', 1.4);
ylabel('Mass residual (%)');
xlabel('Time (min)');
grid on;
exportgraphics(fig, fullfile(fig_dir, 'P1_BC_INFLOW_001_HYDROGRAPH.png'), 'Resolution', 200);
close(fig);
end

function ensure_dir(path_name)
if ~exist(path_name, 'dir')
    mkdir(path_name);
end
end
