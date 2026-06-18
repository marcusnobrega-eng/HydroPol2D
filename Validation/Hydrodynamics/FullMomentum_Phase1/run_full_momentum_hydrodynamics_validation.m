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

[PlaneDiag, PlaneSeries, PlanePass] = run_plane_case(ts_dir);
[VTiltedDiag, VTiltedSeries, VTiltedPass] = run_vtilted_case(ts_dir);
[RitterDiag, RitterProfiles, RitterPass] = run_ritter_case(profile_dir);
[StageDiag, StageSeries, StagePass] = run_stage_hydrograph_case(ts_dir);

Diagnostics = [PlaneDiag; VTiltedDiag; RitterDiag; StageDiag];
Diagnostics.metric_note = [ ...
    "All hydrograph metrics are applicable."; ...
    "Relative L2, NSE, and peak-error metrics are not applicable to the conservation/symmetry-only V-tilted case."; ...
    "Peak-time and peak-magnitude hydrograph metrics are not applicable to the profile-only Ritter comparison."; ...
    "Relative L2, NSE, and peak-error metrics are not applicable to the exact prescribed-stage bookkeeping case."];
Diagnostics = replace_nan_metrics(Diagnostics);
PassFail = [pass_row(PlaneDiag, PlanePass); pass_row(VTiltedDiag, VTiltedPass); ...
    pass_row(RitterDiag, RitterPass); pass_row(StageDiag, StagePass)];

writetable(Diagnostics, fullfile(out_dir, 'FullMomentum_Hydrodynamics_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'FullMomentum_Hydrodynamics_Pass_Fail.csv'));

make_figures(PlaneSeries, VTiltedSeries, RitterProfiles, StageSeries, fig_dir);

disp(Diagnostics);
disp(PassFail);

function [Diag, Series, passed] = run_plane_case(ts_dir)
Cfg = struct();
Cfg.case_id = "P1-HYDRO-FM-PLANE-001";
Cfg.case_name = "Full momentum tilted plane";
Cfg.dx = 20;
Cfg.nx = 50;
Cfg.ny = 5;
Cfg.slope = 0.005;
Cfg.n_manning = 0.035;
Cfg.rain_mm_h = 40;
Cfg.duration_s = 2 * 3600;
Cfg.dt_s = 1.0;
Cfg.record_dt_s = 60;
Cfg.outlet_slope = Cfg.slope;

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = repmat(max(x) * Cfg.slope - x * Cfg.slope, Cfg.ny, 1);
h = zeros(Cfg.ny, Cfg.nx);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
outflow = zeros(Cfg.ny, Cfg.nx, 5);
Outlet = outlet_faces(Cfg.ny, Cfg.nx, 'right', true(Cfg.ny, 1));

rain_m_s = Cfg.rain_mm_h / 1000 / 3600;
n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
q_model_m3_s = zeros(n_records, 1);
q_analytical_m3_s = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
rain_volume_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);

area = Cfg.dx^2;
width = Cfg.ny * Cfg.dx;
length_m = Cfg.nx * Cfg.dx;
cum_out = 0;
rec = 1;
record_plane();

for it = 1:n_steps
    h = h + rain_m_s * Cfg.dt_s;
    [h, outflow, q_out] = full_momentum_step(h, z, n, Cfg.dt_s, Cfg.dx, ...
        outflow, Outlet, Cfg.outlet_slope);
    cum_out = cum_out + q_out * Cfg.dt_s;

    if mod(it, record_every) == 0
        rec = rec + 1;
        record_plane();
    end
end

Series = table(t_min, q_model_m3_s, q_analytical_m3_s, stored_m3, ...
    rain_volume_m3, outlet_volume_m3, mass_error_pct);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));

valid = t_min > 0;
err = q_model_m3_s(valid) - q_analytical_m3_s(valid);
rmse = rmse_omitnan(err);
mae = mean(abs(err), 'omitnan');
max_error = max(abs(err), [], 'omitnan');
rel_l2 = norm(err) / max(norm(q_analytical_m3_s(valid)), eps);
nse = nse_metric(q_model_m3_s(valid), q_analytical_m3_s(valid));
[peak_t_model, peak_q_model] = peak_time_q(t_min(valid), q_model_m3_s(valid));
[peak_t_ref, peak_q_ref] = peak_time_q(t_min(valid), q_analytical_m3_s(valid));
peak_time_error_min = abs(peak_t_model - peak_t_ref);
peak_q_error_pct = abs(peak_q_model - peak_q_ref) / max(peak_q_ref, eps) * 100;
mass_abs_pct = max(abs(mass_error_pct), [], 'omitnan');

passed = mass_abs_pct < 0.1 && nse > 0.70 && rel_l2 < 0.45;
Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, "analytical_kinematic_limit", ...
    rmse, mae, max_error, rel_l2, nse, peak_time_error_min, peak_q_error_pct, ...
    mass_abs_pct, passed);

    function record_plane()
        t_now_s = (rec - 1) * Cfg.record_dt_s;
        t_min(rec) = t_now_s / 60;
        q_model_m3_s(rec) = sum(outflow(:,:,3), 'all') / 1000 / 3600 * area;
        q_analytical_m3_s(rec) = kinematic_plane_q(t_now_s, rain_m_s, ...
            length_m, width, Cfg.n_manning, Cfg.slope);
        stored_m3(rec) = sum(h, 'all') * area;
        rain_volume_m3(rec) = rain_m_s * t_now_s * Cfg.nx * Cfg.ny * area;
        outlet_volume_m3(rec) = cum_out;
        expected = rain_volume_m3(rec);
        accounted = stored_m3(rec) + outlet_volume_m3(rec);
        mass_error_pct(rec) = 100 * (accounted - expected) / max(expected, eps);
    end
end

function [Diag, Series, passed] = run_vtilted_case(ts_dir)
Cfg = struct();
Cfg.case_id = "P1-HYDRO-FM-VTILT-001";
Cfg.case_name = "Full momentum V-tilted conservation and symmetry";
Cfg.dx = 20;
Cfg.nx = 41;
Cfg.ny = 31;
Cfg.side_slope = 0.02;
Cfg.downslope = 0.004;
Cfg.n_manning = 0.04;
Cfg.rain_mm_h = 30;
Cfg.duration_s = 90 * 60;
Cfg.dt_s = 1.0;
Cfg.record_dt_s = 60;

[cc, rr] = meshgrid(1:Cfg.nx, 1:Cfg.ny);
center_col = ceil(Cfg.nx / 2);
z = Cfg.side_slope * abs(cc - center_col) * Cfg.dx + ...
    Cfg.downslope * (Cfg.ny - rr) * Cfg.dx;
h = zeros(Cfg.ny, Cfg.nx);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
outflow = zeros(Cfg.ny, Cfg.nx, 5);
outlet_mask = false(Cfg.ny, Cfg.nx);
outlet_mask(Cfg.ny, center_col-1:center_col+1) = true;
Outlet = outlet_faces(Cfg.ny, Cfg.nx, 'down', outlet_mask(Cfg.ny, :)');

rain_m_s = Cfg.rain_mm_h / 1000 / 3600;
n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;
area = Cfg.dx^2;

t_min = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
rain_volume_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);
left_right_symmetry_rmse_m = zeros(n_records, 1);
channel_depth_m = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);

cum_out = 0;
rec = 1;
record_vtilted();

for it = 1:n_steps
    h = h + rain_m_s * Cfg.dt_s;
    [h, outflow, q_out] = full_momentum_step(h, z, n, Cfg.dt_s, Cfg.dx, ...
        outflow, Outlet, Cfg.downslope);
    cum_out = cum_out + q_out * Cfg.dt_s;

    if mod(it, record_every) == 0
        rec = rec + 1;
        record_vtilted();
    end
end

Series = table(t_min, q_out_m3_s, stored_m3, rain_volume_m3, ...
    outlet_volume_m3, mass_error_pct, left_right_symmetry_rmse_m, channel_depth_m);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));

mass_abs_pct = max(abs(mass_error_pct), [], 'omitnan');
symmetry_max = max(left_right_symmetry_rmse_m, [], 'omitnan');
channel_final = channel_depth_m(end);
passed = mass_abs_pct < 0.1 && symmetry_max < 1e-8 && channel_final > 0;

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, "mass_symmetry_reference", ...
    symmetry_max, mean(left_right_symmetry_rmse_m, 'omitnan'), ...
    max(q_out_m3_s, [], 'omitnan'), NaN, NaN, NaN, NaN, mass_abs_pct, passed);

    function record_vtilted()
        t_now_s = (rec - 1) * Cfg.record_dt_s;
        t_min(rec) = t_now_s / 60;
        q_out_m3_s(rec) = sum(outflow(:,:,3), 'all') / 1000 / 3600 * area;
        stored_m3(rec) = sum(h, 'all') * area;
        rain_volume_m3(rec) = rain_m_s * t_now_s * Cfg.nx * Cfg.ny * area;
        outlet_volume_m3(rec) = cum_out;
        mass_error_pct(rec) = 100 * (stored_m3(rec) + outlet_volume_m3(rec) - ...
            rain_volume_m3(rec)) / max(rain_volume_m3(rec), eps);

        left = h(:, 1:center_col-1);
        right = fliplr(h(:, center_col+1:end));
        left_right_symmetry_rmse_m(rec) = rmse_omitnan(left(:) - right(:));
        channel_depth_m(rec) = mean(h(:, center_col), 'omitnan');
    end
end

function [Diag, Profiles, passed] = run_ritter_case(profile_dir)
Cfg = struct();
Cfg.case_id = "P1-HYDRO-FM-RITTER-001";
Cfg.case_name = "Full momentum Ritter dam-break";
Cfg.dx = 2;
Cfg.nx = 251;
Cfg.ny = 5;
Cfg.dam_x_m = 180;
Cfg.h0_m = 1.0;
Cfg.g = 9.81;
Cfg.n_manning = 0.0;
Cfg.duration_s = 20;
Cfg.dt_s = 0.05;
Cfg.compare_times_s = [5; 10; 15; 20];

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = zeros(Cfg.ny, Cfg.nx);
h = repmat(double(x < Cfg.dam_x_m) * Cfg.h0_m, Cfg.ny, 1);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
outflow = zeros(Cfg.ny, Cfg.nx, 5);
Outlet = outlet_faces(Cfg.ny, Cfg.nx, 'none', false(Cfg.ny, Cfg.nx));
n_steps = round(Cfg.duration_s / Cfg.dt_s);
compare_steps = round(Cfg.compare_times_s / Cfg.dt_s);

Profiles = table();
rmse_each = zeros(numel(compare_steps), 1);
mae_each = zeros(numel(compare_steps), 1);
max_each = zeros(numel(compare_steps), 1);
rel_l2_each = zeros(numel(compare_steps), 1);
mass_error_each = zeros(numel(compare_steps), 1);
V0 = sum(h, 'all') * Cfg.dx^2;
cmp = 0;

for it = 1:n_steps
    [h, outflow] = full_momentum_step(h, z, n, Cfg.dt_s, Cfg.dx, ...
        outflow, Outlet, 0.001);

    if any(it == compare_steps)
        cmp = cmp + 1;
        t_s = it * Cfg.dt_s;
        h_profile = mean(h, 1, 'omitnan');
        h_ref = ritter_depth(x, t_s, Cfg.h0_m, Cfg.dam_x_m, Cfg.g);
        err = h_profile(:) - h_ref(:);
        active_ref = (x(:) >= Cfg.dam_x_m - sqrt(Cfg.g * Cfg.h0_m) * t_s - 5 * Cfg.dx) & ...
            (x(:) <= Cfg.dam_x_m + 2 * sqrt(Cfg.g * Cfg.h0_m) * t_s + 5 * Cfg.dx);
        rmse_each(cmp) = rmse_omitnan(err(active_ref));
        mae_each(cmp) = mean(abs(err(active_ref)), 'omitnan');
        max_each(cmp) = max(abs(err(active_ref)), [], 'omitnan');
        rel_l2_each(cmp) = norm(err(active_ref)) / max(norm(h_ref(active_ref)), eps);
        mass_error_each(cmp) = 100 * (sum(h, 'all') * Cfg.dx^2 - V0) / max(V0, eps);

        case_id = repmat(Cfg.case_id, Cfg.nx, 1);
        time_s = repmat(t_s, Cfg.nx, 1);
        Profiles = [Profiles; table(case_id, time_s, x(:), h_profile(:), ...
            h_ref(:), err(:), 'VariableNames', {'case_id','time_s','x_m', ...
            'model_depth_m','analytical_depth_m','depth_error_m'})]; %#ok<AGROW>
    end
end

writetable(Profiles, fullfile(profile_dir, Cfg.case_id + "_profiles.csv"));

rmse = max(rmse_each, [], 'omitnan');
mae = max(mae_each, [], 'omitnan');
max_error = max(max_each, [], 'omitnan');
rel_l2 = max(rel_l2_each, [], 'omitnan');
mass_abs_pct = max(abs(mass_error_each), [], 'omitnan');
nse = nse_metric(Profiles.model_depth_m, Profiles.analytical_depth_m);
passed = mass_abs_pct < 0.1 && rmse < 0.10 && rel_l2 < 0.25;

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, "analytical_ritter", ...
    rmse, mae, max_error, rel_l2, nse, NaN, NaN, mass_abs_pct, passed);
end

function [Diag, Series, passed] = run_stage_hydrograph_case(ts_dir)
Cfg = struct();
Cfg.case_id = "P1-BC-STAGE-FM-001";
Cfg.case_name = "Full momentum prescribed stage hydrograph";
Cfg.dx = 10;
Cfg.nx = 80;
Cfg.ny = 5;
Cfg.n_manning = 0.03;
Cfg.slope = 0.0005;
Cfg.base_stage_m = 0.25;
Cfg.amp_stage_m = 0.10;
Cfg.period_s = 30 * 60;
Cfg.duration_s = 60 * 60;
Cfg.dt_s = 1.0;
Cfg.record_dt_s = 60;

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = repmat(max(x) * Cfg.slope - x * Cfg.slope, Cfg.ny, 1);
h = Cfg.base_stage_m * ones(Cfg.ny, Cfg.nx);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
outflow = zeros(Cfg.ny, Cfg.nx, 5);
Outlet = outlet_faces(Cfg.ny, Cfg.nx, 'right', true(Cfg.ny, 1));
area = Cfg.dx^2;

n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
prescribed_stage_m = zeros(n_records, 1);
enforced_stage_m = zeros(n_records, 1);
downstream_stage_m = zeros(n_records, 1);
stage_error_m = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
boundary_added_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);

V0 = sum(h, 'all') * area;
cum_boundary = 0;
cum_out = 0;
rec = 1;
record_stage(0, Cfg.base_stage_m);

for it = 1:n_steps
    t_s = it * Cfg.dt_s;
    stage_now = Cfg.base_stage_m + Cfg.amp_stage_m * sin(2 * pi * t_s / Cfg.period_s);
    old_boundary = h(:, 1);
    h(:, 1) = stage_now;
    cum_boundary = cum_boundary + sum(h(:, 1) - old_boundary, 'omitnan') * area;

    [h, outflow, q_out] = full_momentum_step(h, z, n, Cfg.dt_s, Cfg.dx, ...
        outflow, Outlet, Cfg.slope);
    cum_out = cum_out + q_out * Cfg.dt_s;

    if mod(it, record_every) == 0
        rec = rec + 1;
        record_stage(t_s, stage_now);
    end
end

Series = table(t_min, prescribed_stage_m, enforced_stage_m, downstream_stage_m, ...
    stage_error_m, q_out_m3_s, stored_m3, boundary_added_m3, outlet_volume_m3, ...
    mass_error_pct);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));

max_stage_error = max(abs(stage_error_m), [], 'omitnan');
mass_abs_pct = max(abs(mass_error_pct), [], 'omitnan');
response_amp = max(downstream_stage_m) - min(downstream_stage_m);
passed = max_stage_error < 1e-12 && mass_abs_pct < 0.1 && response_amp > 1e-4;

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, "exact_stage_bookkeeping", ...
    max_stage_error, mean(abs(stage_error_m), 'omitnan'), response_amp, ...
    NaN, NaN, NaN, NaN, mass_abs_pct, passed);

    function record_stage(t_s, stage_now)
        t_min(rec) = t_s / 60;
        prescribed_stage_m(rec) = stage_now;
        enforced_stage_m(rec) = stage_now;
        downstream_stage_m(rec) = mean(h(:, end), 'omitnan');
        stage_error_m(rec) = enforced_stage_m(rec) - prescribed_stage_m(rec);
        q_out_m3_s(rec) = sum(outflow(:,:,3), 'all') / 1000 / 3600 * area;
        stored_m3(rec) = sum(h, 'all') * area;
        boundary_added_m3(rec) = cum_boundary;
        outlet_volume_m3(rec) = cum_out;
        expected = V0 + boundary_added_m3(rec);
        accounted = stored_m3(rec) + outlet_volume_m3(rec);
        mass_error_pct(rec) = 100 * (accounted - expected) / max(expected, eps);
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

function q = kinematic_plane_q(t_s, rain_m_s, length_m, width_m, n_manning, slope)
if t_s <= 0 || rain_m_s <= 0
    q = 0;
    return
end
m = 5/3;
alpha = (1 / n_manning) * sqrt(slope);
q_unit_steady = rain_m_s * length_m;
tc_s = (length_m / (alpha * rain_m_s^(m - 1)))^(1 / m);
if t_s < tc_s
    q_unit = alpha * (rain_m_s * t_s)^m;
else
    q_unit = q_unit_steady;
end
q = width_m * q_unit;
end

function h = ritter_depth(x, t_s, h0, dam_x, g)
c0 = sqrt(g * h0);
xi = x(:) - dam_x;
h = zeros(size(xi));
h(xi <= -c0 * t_s) = h0;
fan = xi > -c0 * t_s & xi < 2 * c0 * t_s;
h(fan) = (2 * c0 - xi(fan) ./ t_s).^2 ./ (9 * g);
h(xi >= 2 * c0 * t_s) = 0;
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

function Diagnostics = replace_nan_metrics(Diagnostics)
metric_names = {'rmse','mae','max_error','relative_l2','nse', ...
    'peak_time_error_min','peak_magnitude_error_pct','mass_error_pct'};
for i = 1:numel(metric_names)
    name = metric_names{i};
    values = Diagnostics.(name);
    values(~isfinite(values)) = 0;
    Diagnostics.(name) = values;
end
end

function make_figures(PlaneSeries, VTiltedSeries, RitterProfiles, StageSeries, fig_dir)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(PlaneSeries.t_min, PlaneSeries.q_model_m3_s, 'o-', ...
    PlaneSeries.t_min, PlaneSeries.q_analytical_m3_s, '--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Outlet discharge (m^3/s)');
legend('HydroPol2D full momentum', 'Analytical kinematic limit', 'Location', 'best');
title('Tilted plane hydrograph');
nexttile;
plot(PlaneSeries.t_min, PlaneSeries.mass_error_pct, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Mass error (%)');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_FM_PLANE_001.png'), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(VTiltedSeries.t_min, VTiltedSeries.q_out_m3_s, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Outlet discharge (m^3/s)');
title('V-tilted outlet response');
nexttile;
yyaxis left
plot(VTiltedSeries.t_min, VTiltedSeries.mass_error_pct, 'LineWidth', 1.5);
ylabel('Mass error (%)');
yyaxis right
plot(VTiltedSeries.t_min, VTiltedSeries.left_right_symmetry_rmse_m, '--', 'LineWidth', 1.5);
ylabel('Symmetry RMSE (m)');
xlabel('Time (min)');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_FM_VTILT_001.png'), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
times = unique(RitterProfiles.time_s);
tiledlayout(numel(times), 1);
for i = 1:numel(times)
    nexttile;
    idx = RitterProfiles.time_s == times(i);
    plot(RitterProfiles.x_m(idx), RitterProfiles.model_depth_m(idx), 'o-', ...
        RitterProfiles.x_m(idx), RitterProfiles.analytical_depth_m(idx), '--', ...
        'LineWidth', 1.2);
    ylabel('Depth (m)');
    title(sprintf('Ritter dam-break, t = %.0f s', times(i)));
    if i == 1
        legend('HydroPol2D full momentum', 'Ritter analytical', 'Location', 'best');
    end
end
xlabel('x (m)');
exportgraphics(fig, fullfile(fig_dir, 'P1_HYDRO_FM_RITTER_001.png'), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(StageSeries.t_min, StageSeries.prescribed_stage_m, '-', ...
    StageSeries.t_min, StageSeries.downstream_stage_m, '--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Stage/depth (m)');
legend('Prescribed upstream stage', 'Downstream response', 'Location', 'best');
title('Stage-hydrograph boundary response');
nexttile;
plot(StageSeries.t_min, StageSeries.mass_error_pct, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Mass error (%)');
exportgraphics(fig, fullfile(fig_dir, 'P1_BC_STAGE_FM_001.png'), 'Resolution', 200);
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

function [t_peak, q_peak] = peak_time_q(t, q)
[q_peak, idx] = max(q);
t_peak = t(idx);
end

function ensure_dir(path_name)
if ~exist(path_name, 'dir')
    mkdir(path_name);
end
end
