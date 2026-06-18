clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
profile_dir = fullfile(out_dir, 'Profiles');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(profile_dir, 'dir'); mkdir(profile_dir); end
if ~exist(ts_dir, 'dir'); mkdir(ts_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

[SteadyDiag, SteadyProfiles, SteadyPass] = run_steady_dupuit_case(functions_dir);
[TransientDiag, TransientProfiles, TransientSeries, TransientPass] = run_transient_linearized_case(functions_dir);

Diagnostics = [SteadyDiag; TransientDiag];
PassFail = [pass_row(SteadyDiag, SteadyPass); pass_row(TransientDiag, TransientPass)];

writetable(SteadyProfiles, fullfile(profile_dir, 'P1-GW-HILL-STEADY-001_profiles.csv'));
writetable(TransientProfiles, fullfile(profile_dir, 'P1-GW-HILL-TRANSIENT-001_profiles.csv'));
writetable(TransientSeries, fullfile(ts_dir, 'P1-GW-HILL-TRANSIENT-001_hydrograph.csv'));
writetable(Diagnostics, fullfile(out_dir, 'Hillslope_Groundwater_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'Hillslope_Groundwater_Pass_Fail.csv'));

make_figures(SteadyProfiles, TransientProfiles, TransientSeries, fig_dir);

disp(Diagnostics);
disp(PassFail);

function [Diag, Profiles, passed] = run_steady_dupuit_case(functions_dir)
Cfg = steady_cfg();
resolutions = [41; 81; 161];
Profiles = table();
rmse_by_resolution = zeros(numel(resolutions), 1);
max_error_by_resolution = zeros(numel(resolutions), 1);
equilibrium_error_by_resolution = zeros(numel(resolutions), 1);
discharge_error_pct = zeros(numel(resolutions), 1);

for ir = 1:numel(resolutions)
    n_col = resolutions(ir);
    dx = Cfg.length_m / (n_col - 1);
    x = (0:n_col-1) * dx;
    H_continuous = sqrt(Cfg.drain_head_m^2 + (Cfg.recharge_m_s / Cfg.K_m_s) .* ...
        (Cfg.length_m^2 - x.^2));
    H_discrete = steady_discrete_profile(Cfg, n_col, dx);

    h0 = repmat(H_discrete, Cfg.n_rows, 1);
    z_bed = zeros(Cfg.n_rows, n_col);
    Sy = Cfg.Sy * ones(Cfg.n_rows, n_col);
    K = Cfg.K_m_s * ones(Cfg.n_rows, n_col);
    R = Cfg.recharge_m_s * ones(Cfg.n_rows, n_col);
    R(:, end) = 0;
    h_soil = Cfg.soil_depth_m * ones(Cfg.n_rows, n_col);
    catchment_mask = true(Cfg.n_rows, n_col);
    dirichlet_mask = false(Cfg.n_rows, n_col);
    dirichlet_mask(:, end) = true;
    h_dirichlet = zeros(Cfg.n_rows, n_col);
    h_dirichlet(:, end) = Cfg.drain_head_m;
    perimeter = false(Cfg.n_rows, n_col);
    river_mask = false(Cfg.n_rows, n_col);
    h_river = zeros(Cfg.n_rows, n_col);
    z_river = zeros(Cfg.n_rows, n_col);

    [h_model, ~, ~, ~, ~, ~] = Boussinesq_2D_explicit( ...
        Cfg.check_dt_s, dx, Cfg.dy_m, h0, z_bed, Sy, R, K, ...
        river_mask, K, h_river, z_river, Cfg.courant, h_soil, ...
        catchment_mask, dirichlet_mask, h_dirichlet, perimeter);

    H_model = mean(h_model, 1, 'omitnan');
    error_continuous = H_model - H_continuous;
    error_equilibrium = H_model - H_discrete;
    model_toe_q_m2_s = toe_flux_per_width(H_model, Cfg.K_m_s, dx);
    expected_toe_q_m2_s = Cfg.recharge_m_s * (n_col - 1) * dx;

    rmse_by_resolution(ir) = rmse_omitnan(error_continuous);
    max_error_by_resolution(ir) = max_abs_omitnan(error_continuous);
    equilibrium_error_by_resolution(ir) = max_abs_omitnan(error_equilibrium);
    discharge_error_pct(ir) = abs(model_toe_q_m2_s - expected_toe_q_m2_s) / ...
        max(expected_toe_q_m2_s, eps) * 100;

    case_id = repmat("P1-GW-HILL-STEADY-001", n_col, 1);
    resolution_cells = repmat(n_col, n_col, 1);
    Profiles = [Profiles; table(case_id, resolution_cells, x(:), ...
        H_model(:), H_continuous(:), H_discrete(:), error_continuous(:), ...
        error_equilibrium(:), ...
        'VariableNames', {'case_id','n_columns','x_m','model_head_m', ...
        'analytical_head_m','discrete_steady_head_m','analytical_error_m', ...
        'equilibrium_error_m'})]; %#ok<AGROW>
end

fine_idx = numel(resolutions);
convergence_ratio = rmse_by_resolution(1) / max(rmse_by_resolution(end), eps);
passed = rmse_by_resolution(fine_idx) < 1e-2 && ...
    max_error_by_resolution(fine_idx) < 2e-2 && ...
    equilibrium_error_by_resolution(fine_idx) < 1e-5 && ...
    discharge_error_pct(fine_idx) < 0.1 && ...
    convergence_ratio > 2;

Diag = diagnostic_row("P1-GW-HILL-STEADY-001", "Steady Dupuit hillslope", ...
    "steady_dupuit", rmse_by_resolution(fine_idx), ...
    max_error_by_resolution(fine_idx), equilibrium_error_by_resolution(fine_idx), ...
    discharge_error_pct(fine_idx), nan, nan, convergence_ratio, passed);
end

function [Diag, Profiles, Series, passed] = run_transient_linearized_case(functions_dir)
Cfg = transient_cfg();
n_col = Cfg.n_col;
dx = Cfg.length_m / (n_col - 1);
x = (0:n_col-1) * dx;
lambda = pi / (2 * Cfg.length_m);
D = Cfg.K_m_s * Cfg.base_head_m / Cfg.Sy;

z_bed = zeros(Cfg.n_rows, n_col);
Sy = Cfg.Sy * ones(Cfg.n_rows, n_col);
K = Cfg.K_m_s * ones(Cfg.n_rows, n_col);
R = zeros(Cfg.n_rows, n_col);
h_soil = Cfg.soil_depth_m * ones(Cfg.n_rows, n_col);
catchment_mask = true(Cfg.n_rows, n_col);
dirichlet_mask = false(Cfg.n_rows, n_col);
dirichlet_mask(:, end) = true;
h_dirichlet = zeros(Cfg.n_rows, n_col);
h_dirichlet(:, end) = Cfg.base_head_m;
perimeter = false(Cfg.n_rows, n_col);
river_mask = false(Cfg.n_rows, n_col);
h_river = zeros(Cfg.n_rows, n_col);
z_river = zeros(Cfg.n_rows, n_col);

u0 = Cfg.amplitude_m * cos(lambda .* x);
h_model = repmat(Cfg.base_head_m + u0, Cfg.n_rows, 1);
h_model(:, end) = Cfg.base_head_m;

n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = round(Cfg.record_dt_s / Cfg.dt_s);
n_records = floor(n_steps / record_every) + 1;

Profiles = table();
time_days = zeros(n_records, 1);
model_toe_q_m2_s = zeros(n_records, 1);
analytical_toe_q_m2_s = zeros(n_records, 1);
head_rmse_m = zeros(n_records, 1);
max_head_error_m = zeros(n_records, 1);

record_idx = 1;
[Profiles, time_days, model_toe_q_m2_s, analytical_toe_q_m2_s, ...
    head_rmse_m, max_head_error_m] = record_transient_state( ...
    Profiles, time_days, model_toe_q_m2_s, analytical_toe_q_m2_s, ...
    head_rmse_m, max_head_error_m, record_idx, 0, h_model, x, Cfg, D, lambda, dx);

for it = 1:n_steps
    [h_model, ~, ~, ~, ~, ~] = Boussinesq_2D_explicit( ...
        Cfg.dt_s, dx, Cfg.dy_m, h_model, z_bed, Sy, R, K, ...
        river_mask, K, h_river, z_river, Cfg.courant, h_soil, ...
        catchment_mask, dirichlet_mask, h_dirichlet, perimeter);

    if mod(it, record_every) == 0
        record_idx = record_idx + 1;
        t_s = it * Cfg.dt_s;
        [Profiles, time_days, model_toe_q_m2_s, analytical_toe_q_m2_s, ...
            head_rmse_m, max_head_error_m] = record_transient_state( ...
            Profiles, time_days, model_toe_q_m2_s, analytical_toe_q_m2_s, ...
            head_rmse_m, max_head_error_m, record_idx, t_s, h_model, x, Cfg, D, lambda, dx);
    end
end

Series = table(time_days, model_toe_q_m2_s, analytical_toe_q_m2_s, ...
    head_rmse_m, max_head_error_m);

valid = isfinite(analytical_toe_q_m2_s) & isfinite(model_toe_q_m2_s);
nse = 1 - sum((model_toe_q_m2_s(valid) - analytical_toe_q_m2_s(valid)).^2) / ...
    max(sum((analytical_toe_q_m2_s(valid) - mean(analytical_toe_q_m2_s(valid))).^2), eps);
final_rmse = head_rmse_m(end);
max_rmse = max(head_rmse_m, [], 'omitnan');
max_profile_error = max(max_head_error_m, [], 'omitnan');
peak_q_error_pct = abs(max(model_toe_q_m2_s) - max(analytical_toe_q_m2_s)) / ...
    max(max(analytical_toe_q_m2_s), eps) * 100;

passed = max_rmse < 0.02 && max_profile_error < 0.05 && nse > 0.99 && ...
    peak_q_error_pct < 5;

Diag = diagnostic_row("P1-GW-HILL-TRANSIENT-001", ...
    "Transient linearized Boussinesq hillslope", ...
    "transient_linearized_boussinesq", final_rmse, max_profile_error, ...
    nan, nan, nse, peak_q_error_pct, nan, passed);
end

function H = steady_discrete_profile(Cfg, n_col, dx)
H2 = zeros(1, n_col);
H2(end) = Cfg.drain_head_m^2;
for j = n_col-1:-1:1
    H2(j) = H2(j+1) + 2 * Cfg.recharge_m_s * j * dx^2 / Cfg.K_m_s;
end
H = sqrt(H2);
end

function q = toe_flux_per_width(H, K, dx)
H_face = 0.5 * (H(end-1) + H(end));
q = -K * H_face * (H(end) - H(end-1)) / dx;
end

function [Profiles, time_days, model_q, analytical_q, head_rmse, max_head_error] = ...
    record_transient_state(Profiles, time_days, model_q, analytical_q, ...
    head_rmse, max_head_error, record_idx, t_s, h_model, x, Cfg, D, lambda, dx)

H_model = mean(h_model, 1, 'omitnan');
decay = exp(-D * lambda^2 * t_s);
H_analytical = Cfg.base_head_m + Cfg.amplitude_m * cos(lambda .* x) .* decay;
H_analytical(end) = Cfg.base_head_m;
err = H_model - H_analytical;

time_days(record_idx) = t_s / 86400;
model_q(record_idx) = toe_flux_per_width(H_model, Cfg.K_m_s, dx);
analytical_q(record_idx) = Cfg.K_m_s * Cfg.base_head_m * Cfg.amplitude_m * lambda * decay;
head_rmse(record_idx) = rmse_omitnan(err);
max_head_error(record_idx) = max_abs_omitnan(err);

case_id = repmat("P1-GW-HILL-TRANSIENT-001", numel(x), 1);
time_days_profile = repmat(time_days(record_idx), numel(x), 1);
Profiles = [Profiles; table(case_id, time_days_profile, x(:), H_model(:), ...
    H_analytical(:), err(:), ...
    'VariableNames', {'case_id','time_days','x_m','model_head_m', ...
    'analytical_head_m','head_error_m'})]; %#ok<AGROW>
end

function make_figures(SteadyProfiles, TransientProfiles, TransientSeries, fig_dir)
fig = figure('Visible', 'off');
tiledlayout(2, 1);
nexttile;
idx = SteadyProfiles.n_columns == max(SteadyProfiles.n_columns);
plot(SteadyProfiles.x_m(idx), SteadyProfiles.model_head_m(idx), 'o-', ...
    SteadyProfiles.x_m(idx), SteadyProfiles.analytical_head_m(idx), '--');
ylabel('Head above bedrock (m)');
legend('HydroPol2D', 'Continuous Dupuit', 'Location', 'best');
title('Steady Dupuit hillslope');
nexttile;
semilogy(SteadyProfiles.x_m(idx), abs(SteadyProfiles.analytical_error_m(idx)), 'o-');
xlabel('Distance from divide (m)');
ylabel('|Error| (m)');
exportgraphics(fig, fullfile(fig_dir, 'P1_GW_HILL_STEADY_001_profile.png'), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off');
tiledlayout(2, 1);
nexttile;
times_to_plot = unique([0; max(TransientProfiles.time_days) / 2; max(TransientProfiles.time_days)]);
for i = 1:numel(times_to_plot)
    [~, idx_time] = min(abs(TransientProfiles.time_days - times_to_plot(i)));
    t_plot = TransientProfiles.time_days(idx_time);
    idx = abs(TransientProfiles.time_days - t_plot) < 1e-12;
    plot(TransientProfiles.x_m(idx), TransientProfiles.model_head_m(idx), 'o-'); hold on
    plot(TransientProfiles.x_m(idx), TransientProfiles.analytical_head_m(idx), '--');
end
ylabel('Head above bedrock (m)');
title('Transient linearized Boussinesq profiles');
nexttile;
plot(TransientSeries.time_days, TransientSeries.model_toe_q_m2_s, 'o-', ...
    TransientSeries.time_days, TransientSeries.analytical_toe_q_m2_s, '--');
xlabel('Time (days)');
ylabel('Toe flux per width (m^2/s)');
legend('HydroPol2D', 'Analytical', 'Location', 'best');
exportgraphics(fig, fullfile(fig_dir, 'P1_GW_HILL_TRANSIENT_001_profiles_hydrograph.png'), 'Resolution', 200);
close(fig);
end

function Cfg = steady_cfg()
Cfg = struct();
Cfg.length_m = 100;
Cfg.n_rows = 3;
Cfg.dy_m = 20;
Cfg.K_m_s = 1e-4;
Cfg.Sy = 0.25;
Cfg.recharge_m_s = 1e-7;
Cfg.drain_head_m = 0.50;
Cfg.soil_depth_m = 5.0;
Cfg.check_dt_s = 60;
Cfg.courant = 0.25;
end

function Cfg = transient_cfg()
Cfg = struct();
Cfg.length_m = 100;
Cfg.n_col = 81;
Cfg.n_rows = 3;
Cfg.dy_m = 20;
Cfg.K_m_s = 1e-4;
Cfg.Sy = 0.25;
Cfg.base_head_m = 2.0;
Cfg.amplitude_m = 0.01;
Cfg.soil_depth_m = 5.0;
Cfg.dt_s = 600;
Cfg.record_dt_s = 86400;
Cfg.duration_s = 20 * 86400;
Cfg.courant = 0.25;
end

function row = diagnostic_row(case_id, case_name, regime, head_rmse_m, ...
    max_head_error_m, equilibrium_max_error_m, discharge_error_pct, ...
    hydrograph_nse, peak_q_error_pct, convergence_ratio, passed)
row = table(case_id, case_name, regime, head_rmse_m, max_head_error_m, ...
    equilibrium_max_error_m, discharge_error_pct, hydrograph_nse, ...
    peak_q_error_pct, convergence_ratio, passed, ...
    'VariableNames', {'case_id','case_name','regime','head_rmse_m', ...
    'max_head_error_m','equilibrium_max_error_m','discharge_error_pct', ...
    'hydrograph_nse','peak_q_error_pct','convergence_ratio','passed'});
end

function row = pass_row(Diag, passed)
status = "fail";
if passed
    status = "pass";
end
row = table(Diag.case_id, Diag.case_name, Diag.regime, status, passed, ...
    'VariableNames', {'case_id','case_name','regime','status','report_ready'});
end

function value = rmse_omitnan(x)
value = sqrt(mean(x(:).^2, 'omitnan'));
end

function value = max_abs_omitnan(x)
value = max(abs(x(:)), [], 'omitnan');
end
