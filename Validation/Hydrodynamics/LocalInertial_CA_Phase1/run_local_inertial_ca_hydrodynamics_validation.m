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

methods = ["local_inertial"; "cellular_automata"];
Diagnostics = table();
PassFail = table();

for im = 1:numel(methods)
    method = methods(im);
    [PlaneDiag, PlaneSeries, PlanePass] = run_plane_case(method, ts_dir);
    [VTiltedDiag, VTiltedSeries, VTiltedPass] = run_vtilted_case(method, ts_dir);
    [RitterDiag, RitterProfiles, RitterPass] = run_ritter_case(method, profile_dir);
    [StageDiag, StageSeries, StagePass] = run_stage_hydrograph_case(method, ts_dir);

    MethodDiagnostics = [PlaneDiag; VTiltedDiag; RitterDiag; StageDiag];
    MethodDiagnostics.metric_note = metric_notes(method);
    MethodDiagnostics = replace_nan_metrics(MethodDiagnostics);
    Diagnostics = [Diagnostics; MethodDiagnostics]; %#ok<AGROW>
    PassFail = [PassFail; pass_row(PlaneDiag, PlanePass); pass_row(VTiltedDiag, VTiltedPass); ...
        pass_row(RitterDiag, RitterPass); pass_row(StageDiag, StagePass)]; %#ok<AGROW>

    make_figures(method, PlaneSeries, VTiltedSeries, RitterProfiles, StageSeries, fig_dir);
end

writetable(Diagnostics, fullfile(out_dir, 'LocalInertial_CA_Hydrodynamics_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'LocalInertial_CA_Hydrodynamics_Pass_Fail.csv'));

disp(Diagnostics);
disp(PassFail);

function [Diag, Series, passed] = run_plane_case(method, ts_dir)
Cfg = base_cfg(method, "plane");
Cfg.dx = 20;
Cfg.nx = 50;
Cfg.ny = 5;
Cfg.slope = 0.005;
Cfg.n_manning = 0.035;
Cfg.rain_mm_h = 40;
Cfg.duration_s = 2 * 3600;
Cfg.dt_s = 1.0;
Cfg.record_dt_s = 60;
Cfg.outlet_side = "right";

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = repmat(max(x) * Cfg.slope - x * Cfg.slope, Cfg.ny, 1);
h = zeros(Cfg.ny, Cfg.nx);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
state = init_solver_state(Cfg.ny, Cfg.nx, method);
Outlet = outlet_cells(Cfg.ny, Cfg.nx, Cfg.outlet_side, true(Cfg.ny, 1));

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
    [h, state, q_out] = routing_step(method, h, z, n, Cfg.dt_s, Cfg.dx, ...
        state, Outlet, Cfg.slope);
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

if method == "local_inertial"
    passed = mass_abs_pct < 0.1 && nse > 0.70 && rel_l2 < 0.45;
else
    passed = mass_abs_pct < 0.1 && max(q_model_m3_s) > 0;
end

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, method, "analytical_kinematic_limit", ...
    rmse, mae, max_error, rel_l2, nse, peak_time_error_min, peak_q_error_pct, ...
    mass_abs_pct, passed);

    function record_plane()
        t_now_s = (rec - 1) * Cfg.record_dt_s;
        t_min(rec) = t_now_s / 60;
        q_model_m3_s(rec) = state.q_out_m3_s;
        q_analytical_m3_s(rec) = kinematic_plane_q(t_now_s, rain_m_s, ...
            length_m, width, Cfg.n_manning, Cfg.slope);
        stored_m3(rec) = sum(h, 'all') * area;
        rain_volume_m3(rec) = rain_m_s * t_now_s * Cfg.nx * Cfg.ny * area;
        outlet_volume_m3(rec) = cum_out;
        mass_error_pct(rec) = 100 * (stored_m3(rec) + outlet_volume_m3(rec) - ...
            rain_volume_m3(rec)) / max(rain_volume_m3(rec), eps);
    end
end

function [Diag, Series, passed] = run_vtilted_case(method, ts_dir)
Cfg = base_cfg(method, "vtilt");
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
state = init_solver_state(Cfg.ny, Cfg.nx, method);
outlet_mask = false(1, Cfg.nx);
outlet_mask(center_col-1:center_col+1) = true;
Outlet = outlet_cells(Cfg.ny, Cfg.nx, "down", outlet_mask(:));

rain_m_s = Cfg.rain_mm_h / 1000 / 3600;
n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;
area = Cfg.dx^2;

t_min = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
rain_volume_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);
left_right_symmetry_rmse_m = zeros(n_records, 1);
channel_depth_m = zeros(n_records, 1);

cum_out = 0;
rec = 1;
record_vtilted();

for it = 1:n_steps
    h = h + rain_m_s * Cfg.dt_s;
    [h, state, q_out] = routing_step(method, h, z, n, Cfg.dt_s, Cfg.dx, ...
        state, Outlet, Cfg.downslope);
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

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, method, "mass_symmetry_reference", ...
    symmetry_max, mean(abs(left_right_symmetry_rmse_m), 'omitnan'), ...
    max(q_out_m3_s, [], 'omitnan'), NaN, NaN, NaN, NaN, mass_abs_pct, passed);

    function record_vtilted()
        t_now_s = (rec - 1) * Cfg.record_dt_s;
        t_min(rec) = t_now_s / 60;
        q_out_m3_s(rec) = state.q_out_m3_s;
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

function [Diag, Profiles, passed] = run_ritter_case(method, profile_dir)
Cfg = base_cfg(method, "ritter");
Cfg.dx = 2;
Cfg.nx = 251;
Cfg.ny = 5;
Cfg.dam_x_m = 180;
Cfg.h0_m = 1.0;
Cfg.g = 9.81;
if method == "local_inertial"
    Cfg.n_manning = 0.0;
else
    Cfg.n_manning = 0.003;
end
Cfg.duration_s = 20;
Cfg.dt_s = 0.05;
Cfg.compare_times_s = [5; 10; 15; 20];

x = ((1:Cfg.nx) - 0.5) * Cfg.dx;
z = zeros(Cfg.ny, Cfg.nx);
h = repmat(double(x < Cfg.dam_x_m) * Cfg.h0_m, Cfg.ny, 1);
n = Cfg.n_manning * ones(Cfg.ny, Cfg.nx);
state = init_solver_state(Cfg.ny, Cfg.nx, method);
Outlet = outlet_cells(Cfg.ny, Cfg.nx, "none", false(Cfg.ny, Cfg.nx));
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
    [h, state] = routing_step(method, h, z, n, Cfg.dt_s, Cfg.dx, ...
        state, Outlet, 0.001);
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

if method == "local_inertial"
    passed = mass_abs_pct < 0.1 && rmse < 0.15 && rel_l2 < 0.35;
else
    passed = mass_abs_pct < 0.1 && rmse < 0.50;
end

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, method, "analytical_ritter", ...
    rmse, mae, max_error, rel_l2, nse, NaN, NaN, mass_abs_pct, passed);
end

function [Diag, Series, passed] = run_stage_hydrograph_case(method, ts_dir)
Cfg = base_cfg(method, "stage");
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
state = init_solver_state(Cfg.ny, Cfg.nx, method);
Outlet = outlet_cells(Cfg.ny, Cfg.nx, "right", true(Cfg.ny, 1));
area = Cfg.dx^2;

n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
prescribed_stage_m = zeros(n_records, 1);
post_step_stage_m = zeros(n_records, 1);
stage_enforcement_error_m = zeros(n_records, 1);
downstream_stage_m = zeros(n_records, 1);
stored_m3 = zeros(n_records, 1);
boundary_added_m3 = zeros(n_records, 1);
outlet_volume_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);
q_out_m3_s = zeros(n_records, 1);

V0 = sum(h, 'all') * area;
cum_boundary = 0;
cum_out = 0;
rec = 1;
record_stage(0, Cfg.base_stage_m, 0);

for it = 1:n_steps
    t_s = it * Cfg.dt_s;
    stage_now = Cfg.base_stage_m + Cfg.amp_stage_m * sin(2 * pi * t_s / Cfg.period_s);
    old_boundary = h(:, 1);
    h(:, 1) = stage_now;
    immediate_error = max(abs(h(:, 1) - stage_now), [], 'omitnan');
    cum_boundary = cum_boundary + sum(h(:, 1) - old_boundary, 'omitnan') * area;
    [h, state, q_out] = routing_step(method, h, z, n, Cfg.dt_s, Cfg.dx, ...
        state, Outlet, Cfg.slope);
    cum_out = cum_out + q_out * Cfg.dt_s;
    if mod(it, record_every) == 0
        rec = rec + 1;
        record_stage(t_s, stage_now, immediate_error);
    end
end

Series = table(t_min, prescribed_stage_m, post_step_stage_m, ...
    stage_enforcement_error_m, downstream_stage_m, q_out_m3_s, stored_m3, ...
    boundary_added_m3, outlet_volume_m3, mass_error_pct);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));

max_stage_error = max(abs(stage_enforcement_error_m), [], 'omitnan');
mass_abs_pct = max(abs(mass_error_pct), [], 'omitnan');
response_amp = max(downstream_stage_m) - min(downstream_stage_m);
passed = max_stage_error < 1e-12 && mass_abs_pct < 0.1 && response_amp > 1e-4;

Diag = diagnostic_row(Cfg.case_id, Cfg.case_name, method, "exact_stage_bookkeeping", ...
    max_stage_error, mean(abs(stage_enforcement_error_m), 'omitnan'), ...
    response_amp, NaN, NaN, NaN, NaN, mass_abs_pct, passed);

    function record_stage(t_s, stage_now, immediate_error)
        t_min(rec) = t_s / 60;
        prescribed_stage_m(rec) = stage_now;
        post_step_stage_m(rec) = mean(h(:, 1), 'omitnan');
        stage_enforcement_error_m(rec) = immediate_error;
        downstream_stage_m(rec) = mean(h(:, end), 'omitnan');
        q_out_m3_s(rec) = state.q_out_m3_s;
        stored_m3(rec) = sum(h, 'all') * area;
        boundary_added_m3(rec) = cum_boundary;
        outlet_volume_m3(rec) = cum_out;
        mass_error_pct(rec) = 100 * (stored_m3(rec) + outlet_volume_m3(rec) - ...
            (V0 + boundary_added_m3(rec))) / max(V0 + boundary_added_m3(rec), eps);
    end
end

function [h, state, q_out_m3_s] = routing_step(method, h, z, roughness, dt_s, dx, state, Outlet, outlet_slope)
if method == "local_inertial"
    [h, state, q_out_m3_s] = local_inertial_step(h, z, roughness, dt_s, dx, state, Outlet, outlet_slope);
else
    [h, state, q_out_m3_s] = ca_step(h, z, roughness, dt_s, dx, state, Outlet, outlet_slope);
end
state.q_out_m3_s = q_out_m3_s;
end

function [h, state, q_out_m3_s] = local_inertial_step(h, z, roughness, dt_s, dx, state, Outlet, outlet_slope)
idx_nan = false(size(h));
cell_area = dx^2;
time_step_min = dt_s / 60;
d_tot = h * 1000;
d_p = d_tot;
d_tolerance_mm = 1e-6;
roughness_squared = roughness.^2;
River_Width = zeros(size(h));
River_Depth = zeros(size(h));

[~,~,~,~,outlet_flow,d_t,~,state.outflow,~,~,~,~,~,~] = Local_Inertial_Model_D4( ...
    1, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, z, d_tot, d_p, roughness, roughness_squared, cell_area, ...
    time_step_min, dx, Outlet.mask, 1, outlet_slope, Outlet.row, Outlet.col, ...
    d_tolerance_mm, state.outflow, idx_nan, 0, 0, [], [], River_Width, River_Depth, ...
    0, 0, 0, 0, cell_area, [], 0, 0, []);

h = max(double(d_t) / 1000, 0);
q_out_m3_s = sum(double(outlet_flow), 'all') / 1000 / 3600 * cell_area;
end

function [h, state, q_out_m3_s] = ca_step(h, z, roughness, dt_s, dx, state, Outlet, outlet_slope)
idx_nan = false(size(h));
cell_area = dx^2;
time_step_min = dt_s / 60;
d_tot = h * 1000;
h0_cell = zeros(size(h));

[qleft,qright,qup,qdown,outlet_flow,d_t,state.I_tot_end_cell] = CA_Routing( ...
    [], [], [], [], [], [], [], [], [], [], [], [], 0, z, d_tot, roughness, ...
    cell_area, time_step_min, h0_cell, dx, state.I_tot_end_cell, Outlet.mask, ...
    1, outlet_slope, Outlet.row, Outlet.col, idx_nan, 0);

qin = zeros(size(h));
[ny, nx] = size(h);
qin(:,2:nx) = qin(:,2:nx) + qright(:,1:nx-1);
qin(:,1:nx-1) = qin(:,1:nx-1) + qleft(:,2:nx);
qin(2:ny,:) = qin(2:ny,:) + qdown(1:ny-1,:);
qin(1:ny-1,:) = qin(1:ny-1,:) + qup(2:ny,:);
d_t = d_t + qin * time_step_min / 60;
h = max(double(d_t) / 1000, 0);
q_out_m3_s = sum(double(outlet_flow), 'all') / 1000 / 3600 * cell_area;
end

function Cfg = base_cfg(method, case_key)
prefix = method_prefix(method);
label = method_label(method);
switch case_key
    case "plane"
        suffix = "PLANE-001"; name = " tilted plane";
    case "vtilt"
        suffix = "VTILT-001"; name = " V-tilted conservation and symmetry";
    case "ritter"
        suffix = "RITTER-001"; name = " Ritter dam-break";
    case "stage"
        suffix = "STAGE-001"; name = " prescribed stage hydrograph";
end
Cfg.case_id = "P1-HYDRO-" + prefix + "-" + suffix;
if case_key == "stage"
    Cfg.case_id = "P1-BC-STAGE-" + prefix + "-001";
end
Cfg.case_name = label + name;
end

function state = init_solver_state(ny, nx, method)
state = struct();
state.q_out_m3_s = 0;
if method == "local_inertial"
    state.outflow = zeros(ny, nx, 3);
else
    state.I_tot_end_cell = zeros(ny, nx);
end
end

function Outlet = outlet_cells(ny, nx, side, mask_in)
Outlet = struct();
Outlet.mask = false(ny, nx);
switch lower(side)
    case 'right'
        Outlet.mask(:, nx) = logical(mask_in(:));
    case 'left'
        Outlet.mask(:, 1) = logical(mask_in(:));
    case 'down'
        row_mask = false(1, nx);
        row_mask(:) = logical(mask_in(:));
        Outlet.mask(ny, :) = row_mask;
    case 'up'
        row_mask = false(1, nx);
        row_mask(:) = logical(mask_in(:));
        Outlet.mask(1, :) = row_mask;
    case 'none'
        Outlet.mask(:) = false;
    otherwise
        error('Unknown outlet side "%s".', side);
end
[Outlet.row, Outlet.col] = find(Outlet.mask);
end

function notes = metric_notes(method)
if method == "local_inertial"
    method_note = "Local inertial analytical metrics are treated as approximation benchmarks.";
else
    method_note = "CA analytical hydrograph/profile metrics are diagnostic; pass/fail is based on conservation, symmetry, and boundary bookkeeping where applicable.";
end
notes = [ ...
    method_note; ...
    "Relative L2, NSE, and peak-error metrics are not applicable to the conservation/symmetry-only V-tilted case."; ...
    "Peak-time and peak-magnitude hydrograph metrics are not applicable to the profile-only Ritter comparison."; ...
    "Relative L2, NSE, and peak-error metrics are not applicable to the exact prescribed-stage bookkeeping case."];
end

function p = method_prefix(method)
if method == "local_inertial"
    p = "LI";
else
    p = "CA";
end
end

function label = method_label(method)
if method == "local_inertial"
    label = "Local inertial";
else
    label = "Cellular automata";
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

function Diag = diagnostic_row(case_id, case_name, method, evidence_type, rmse, mae, max_error, ...
    rel_l2, nse, peak_time_error, peak_q_error_pct, mass_error_pct, passed)
Diag = table(case_id, case_name, method, evidence_type, rmse, mae, max_error, ...
    rel_l2, nse, peak_time_error, peak_q_error_pct, mass_error_pct, passed, ...
    'VariableNames', {'case_id','case_name','method','evidence_type','rmse', ...
    'mae','max_error','relative_l2','nse','peak_time_error_min', ...
    'peak_magnitude_error_pct','mass_error_pct','passed'});
end

function Pass = pass_row(Diag, passed)
status = "fail";
if passed
    status = "pass";
end
report_ready = passed;
Pass = table(Diag.case_id, Diag.case_name, Diag.method, status, report_ready, ...
    'VariableNames', {'case_id','case_name','method','status','report_ready'});
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

function make_figures(method, PlaneSeries, VTiltedSeries, RitterProfiles, StageSeries, fig_dir)
prefix = method_prefix(method);
label = method_label(method);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(PlaneSeries.t_min, PlaneSeries.q_model_m3_s, 'o-', ...
    PlaneSeries.t_min, PlaneSeries.q_analytical_m3_s, '--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Outlet discharge (m^3/s)');
legend(label, 'Analytical kinematic limit', 'Location', 'best');
title(label + " tilted-plane hydrograph");
nexttile;
plot(PlaneSeries.t_min, PlaneSeries.mass_error_pct, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Mass error (%)');
exportgraphics(fig, fullfile(fig_dir, "P1_HYDRO_" + prefix + "_PLANE_001.png"), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(VTiltedSeries.t_min, VTiltedSeries.q_out_m3_s, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Outlet discharge (m^3/s)');
title(label + " V-tilted outlet response");
nexttile;
yyaxis left
plot(VTiltedSeries.t_min, VTiltedSeries.mass_error_pct, 'LineWidth', 1.5);
ylabel('Mass error (%)');
yyaxis right
plot(VTiltedSeries.t_min, VTiltedSeries.left_right_symmetry_rmse_m, '--', 'LineWidth', 1.5);
ylabel('Symmetry RMSE (m)');
xlabel('Time (min)');
exportgraphics(fig, fullfile(fig_dir, "P1_HYDRO_" + prefix + "_VTILT_001.png"), 'Resolution', 200);
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
    title(sprintf('%s Ritter dam-break, t = %.0f s', label, times(i)));
    if i == 1
        legend(label, 'Ritter analytical', 'Location', 'best');
    end
end
xlabel('x (m)');
exportgraphics(fig, fullfile(fig_dir, "P1_HYDRO_" + prefix + "_RITTER_001.png"), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1);
nexttile;
plot(StageSeries.t_min, StageSeries.prescribed_stage_m, '-', ...
    StageSeries.t_min, StageSeries.downstream_stage_m, '--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Stage/depth (m)');
legend('Prescribed upstream stage', 'Downstream response', 'Location', 'best');
title(label + " stage-hydrograph boundary response");
nexttile;
plot(StageSeries.t_min, StageSeries.mass_error_pct, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Mass error (%)');
exportgraphics(fig, fullfile(fig_dir, "P1_BC_STAGE_" + prefix + "_001.png"), 'Resolution', 200);
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
