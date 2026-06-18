clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
profile_dir = fullfile(out_dir, 'Profiles');
fig_dir = fullfile(out_dir, 'Figures');
ensure_dir(out_dir);
ensure_dir(ts_dir);
ensure_dir(profile_dir);
ensure_dir(fig_dir);

[Diag, Series, Profiles, Pass] = run_nonbreaking_case(ts_dir, profile_dir);
Diag.metric_note = "Boundary-driven non-breaking wave compared against analytical traveling profile; mass metric is stored-volume residual against analytical domain storage.";

writetable(Diag, fullfile(out_dir, 'NonBreakingWave_Hydrodynamics_Diagnostics.csv'));
writetable(Pass, fullfile(out_dir, 'NonBreakingWave_Hydrodynamics_Pass_Fail.csv'));
writetable(Diag, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(table(Diag.case_id, Diag.case_name, Diag.mass_error_pct, Diag.passed, ...
    'VariableNames', {'case_id','case_name','max_storage_residual_pct','passed'}), ...
    fullfile(out_dir, 'Mass_Balance.csv'));

make_figures(Series, Profiles, fig_dir);

disp(Diag);
disp(Pass);

function [Diag, Series, Profiles, Pass] = run_nonbreaking_case(ts_dir, profile_dir)
Cfg = struct();
Cfg.case_id = "P1-HYDRO-FM-NBW-001";
Cfg.case_name = "Full momentum non-breaking wave";
Cfg.dx = 5;
Cfg.nx = 850;
Cfg.ny = 5;
Cfg.n_manning = 0.02;
Cfg.wave_speed_m_s = 1.0;
Cfg.duration_s = 60 * 60;
Cfg.dt_s = 0.50;
Cfg.record_dt_s = 60;
Cfg.compare_times_s = [12; 36; 48; 60] * 60;
Cfg.boundary_cols = 3;
Cfg.depth_threshold_m = 1e-4;

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = zeros(Cfg.ny, Cfg.nx);
h = repmat(analytical_nonbreaking_depth(x, 0, Cfg.n_manning, Cfg.wave_speed_m_s), Cfg.ny, 1);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
outflow = zeros(Cfg.ny, Cfg.nx, 5);
Outlet = outlet_faces(Cfg.ny, Cfg.nx, 'none', false(Cfg.ny, Cfg.nx));

% The analytical wave is a right-moving normal-flow profile; initialize
% momentum consistently where the analytical profile is wet.
outflow(:,:,4) = h * Cfg.wave_speed_m_s;
outflow(:,:,5) = 0;

n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;
compare_steps = round(Cfg.compare_times_s / Cfg.dt_s);
area = Cfg.dx^2;
domain_width = Cfg.ny * Cfg.dx;

t_min = zeros(n_records, 1);
model_storage_m3 = zeros(n_records, 1);
analytical_storage_m3 = zeros(n_records, 1);
storage_residual_m3 = zeros(n_records, 1);
storage_residual_pct = zeros(n_records, 1);
front_model_m = zeros(n_records, 1);
front_analytical_m = zeros(n_records, 1);
front_error_m = zeros(n_records, 1);
max_depth_model_m = zeros(n_records, 1);
max_depth_analytical_m = zeros(n_records, 1);
outlet_q_m3_s = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
boundary_volume_m3 = zeros(n_records, 1);

Profiles = table();
rmse_each = zeros(numel(compare_steps), 1);
mae_each = zeros(numel(compare_steps), 1);
max_each = zeros(numel(compare_steps), 1);
rel_l2_each = zeros(numel(compare_steps), 1);
front_error_each = zeros(numel(compare_steps), 1);
volume_error_each = zeros(numel(compare_steps), 1);

cum_out = 0;
cum_boundary = 0;
rec = 1;
cmp = 0;
record_state(0, 0);

for it = 1:n_steps
    t_s = it * Cfg.dt_s;

    [h, outflow, boundary_added] = impose_upstream_analytical_strip(h, outflow, x, t_s, Cfg, area);
    cum_boundary = cum_boundary + boundary_added;

    [h, outflow, q_out] = full_momentum_step(h, z, n, Cfg.dt_s, Cfg.dx, ...
        outflow, Outlet, 1e-4);
    cum_out = cum_out + q_out * Cfg.dt_s;

    if any(it == compare_steps)
        cmp = cmp + 1;
        add_profile(t_s, cmp);
    end

    if mod(it, record_every) == 0
        rec = rec + 1;
        record_state(t_s, q_out);
    end
end

Series = table(t_min, model_storage_m3, analytical_storage_m3, storage_residual_m3, ...
    storage_residual_pct, front_model_m, front_analytical_m, front_error_m, ...
    max_depth_model_m, max_depth_analytical_m, outlet_q_m3_s, outlet_volume_m3, ...
    boundary_volume_m3);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));
writetable(Profiles, fullfile(profile_dir, Cfg.case_id + "_profiles.csv"));

rmse = max(rmse_each, [], 'omitnan');
mae = max(mae_each, [], 'omitnan');
max_error = max(max_each, [], 'omitnan');
rel_l2 = max(rel_l2_each, [], 'omitnan');
nse = nse_metric(Profiles.model_depth_m, Profiles.analytical_depth_m);
peak_time_error_min = 0;
peak_magnitude_error_pct = max(abs(front_error_each), [], 'omitnan');
mass_abs_pct = max(abs(storage_residual_pct), [], 'omitnan');

passed = rmse < 0.05 && rel_l2 < 0.05 && mass_abs_pct < 5 && ...
    max(abs(front_error_each), [], 'omitnan') <= 8 * Cfg.dx;

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, "analytical_nonbreaking_wave", ...
    rmse, mae, max_error, rel_l2, nse, peak_time_error_min, ...
    peak_magnitude_error_pct, mass_abs_pct, passed);
Pass = pass_row(Diag, passed);

    function record_state(t_s, q_out)
        h_ref = analytical_nonbreaking_depth(x, t_s, Cfg.n_manning, Cfg.wave_speed_m_s);
        t_min(rec) = t_s / 60;
        model_storage_m3(rec) = sum(h, 'all') * area;
        analytical_storage_m3(rec) = sum(h_ref, 'all') * domain_width * Cfg.dx;
        storage_residual_m3(rec) = model_storage_m3(rec) - analytical_storage_m3(rec);
        storage_residual_pct(rec) = 100 * storage_residual_m3(rec) / max(analytical_storage_m3(rec), eps);
        front_model_m(rec) = detect_front(mean(h, 1, 'omitnan'), x, Cfg.depth_threshold_m);
        front_analytical_m(rec) = detect_front(h_ref, x, Cfg.depth_threshold_m);
        front_error_m(rec) = front_model_m(rec) - front_analytical_m(rec);
        max_depth_model_m(rec) = max(h, [], 'all', 'omitnan');
        max_depth_analytical_m(rec) = max(h_ref, [], 'all', 'omitnan');
        outlet_q_m3_s(rec) = q_out;
        outlet_volume_m3(rec) = cum_out;
        boundary_volume_m3(rec) = cum_boundary;
    end

    function add_profile(t_s, icmp)
        h_profile = mean(h, 1, 'omitnan');
        h_ref = analytical_nonbreaking_depth(x, t_s, Cfg.n_manning, Cfg.wave_speed_m_s);
        valid = isfinite(h_profile(:)) & isfinite(h_ref(:));
        err = h_profile(:) - h_ref(:);
        rmse_each(icmp) = rmse_omitnan(err(valid));
        mae_each(icmp) = mean(abs(err(valid)), 'omitnan');
        max_each(icmp) = max(abs(err(valid)), [], 'omitnan');
        rel_l2_each(icmp) = norm(err(valid)) / max(norm(h_ref(valid)), eps);
        front_error_each(icmp) = detect_front(h_profile, x, Cfg.depth_threshold_m) - ...
            detect_front(h_ref, x, Cfg.depth_threshold_m);
        volume_error_each(icmp) = 100 * (sum(h, 'all') * area - ...
            sum(h_ref, 'all') * domain_width * Cfg.dx) / ...
            max(sum(h_ref, 'all') * domain_width * Cfg.dx, eps);

        case_id = repmat(Cfg.case_id, Cfg.nx, 1);
        time_s = repmat(t_s, Cfg.nx, 1);
        Profiles = [Profiles; table(case_id, time_s, x(:), h_profile(:), ...
            h_ref(:), err(:), 'VariableNames', {'case_id','time_s','x_m', ...
            'model_depth_m','analytical_depth_m','depth_error_m'})]; %#ok<AGROW>
    end
end

function [h, outflow, boundary_added_m3] = impose_upstream_analytical_strip(h, outflow, x, t_s, Cfg, area)
cols = 1:Cfg.boundary_cols;
h_old = h(:, cols);
h_ref = analytical_nonbreaking_depth(x(cols), t_s, Cfg.n_manning, Cfg.wave_speed_m_s);
h(:, cols) = repmat(h_ref, size(h,1), 1);
outflow(:, cols, 4) = h(:, cols) * Cfg.wave_speed_m_s;
outflow(:, cols, 5) = 0;
boundary_added_m3 = sum(h(:, cols) - h_old, 'all') * area;
end

function h = analytical_nonbreaking_depth(x, t_s, n_manning, wave_speed)
arg = -(7/3) * n_manning^2 * wave_speed^2 .* (x - wave_speed * t_s);
h = zeros(size(x));
wet = arg > 0;
h(wet) = arg(wet).^(3/7);
h(~isfinite(h)) = 0;
h(h < 0) = 0;
end

function front = detect_front(h_profile, x, threshold)
idx = find(h_profile > threshold, 1, 'last');
if isempty(idx)
    front = 0;
else
    front = x(idx);
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

function Diag = diagnostic_row(case_id, case_name, evidence_type, rmse, mae, max_error, ...
    rel_l2, nse, peak_time_error, peak_q_error_pct, mass_error_pct, passed)
Diag = table(case_id, case_name, evidence_type, rmse, mae, max_error, ...
    rel_l2, nse, peak_time_error, peak_q_error_pct, mass_error_pct, passed, ...
    'VariableNames', {'case_id','case_name','evidence_type','rmse', ...
    'mae','max_error','relative_l2','nse','peak_time_error_min', ...
    'peak_magnitude_error_pct','mass_error_pct','passed'});
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

function make_figures(Series, Profiles, fig_dir)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(Series.t_min, Series.model_storage_m3, 'o-', ...
    Series.t_min, Series.analytical_storage_m3, '--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Stored volume (m^3)');
legend('HydroPol2D full momentum', 'Analytical profile', 'Location', 'best');
title('Non-breaking wave stored volume');
nexttile;
yyaxis left
plot(Series.t_min, Series.storage_residual_pct, 'LineWidth', 1.5);
ylabel('Storage residual (%)');
yyaxis right
plot(Series.t_min, Series.front_error_m, '--', 'LineWidth', 1.5);
ylabel('Front error (m)');
xlabel('Time (min)');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_FM_NBW_001_VOLUME.png'), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
times = unique(Profiles.time_s);
tiledlayout(numel(times), 1);
for i = 1:numel(times)
    nexttile;
    idx = Profiles.time_s == times(i);
    plot(Profiles.x_m(idx), Profiles.model_depth_m(idx), 'o-', ...
        Profiles.x_m(idx), Profiles.analytical_depth_m(idx), '--', ...
        'LineWidth', 1.2);
    ylabel('Depth (m)');
    title(sprintf('Non-breaking wave, t = %.0f min', times(i)/60));
    if i == 1
        legend('HydroPol2D full momentum', 'Analytical profile', 'Location', 'best');
    end
end
xlabel('x (m)');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_FM_NBW_001_PROFILES.png'), 'Resolution', 200);
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

function out = nse_metric(model, ref)
valid = isfinite(model) & isfinite(ref);
model = model(valid);
ref = ref(valid);
den = sum((ref - mean(ref)).^2);
if isempty(ref) || den <= 0
    out = NaN;
else
    out = 1 - sum((model - ref).^2) / den;
end
end

function ensure_dir(path_name)
if ~exist(path_name, 'dir')
    mkdir(path_name);
end
end
