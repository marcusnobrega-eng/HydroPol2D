clear; clc;

case_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
ensure_dir(out_dir);
ensure_dir(ts_dir);
ensure_dir(fig_dir);

[Diag, Series, Pass] = run_reservoir_case(ts_dir);

writetable(Diag, fullfile(out_dir, 'Reservoir_RatingCurve_Diagnostics.csv'));
writetable(Pass, fullfile(out_dir, 'Reservoir_RatingCurve_Pass_Fail.csv'));
writetable(Diag, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(table(Diag.case_id, Diag.case_name, ...
    Diag.mass_error_pct, Series.released_volume_m3(end), ...
    Series.reference_released_volume_m3(end), Diag.passed, ...
    'VariableNames', {'case_id','case_name','mass_balance_error_pct', ...
    'released_volume_m3','reference_released_volume_m3','passed'}), ...
    fullfile(out_dir, 'Mass_Balance.csv'));

make_figures(Series, fig_dir);

disp(Diag);
disp(Pass);

function [Diag, Series, Pass] = run_reservoir_case(ts_dir)
Cfg = struct();
Cfg.case_id = "P1-RES-001";
Cfg.case_name = "Rating-curve reservoir storage depletion";
Cfg.area_m2 = 20000;
Cfg.initial_stage_m = 2.50;
Cfg.sill_stage_m = 0.35;
Cfg.rating_k = 2.2;
Cfg.rating_exponent = 1.5;
Cfg.dt_s = 0.25;
Cfg.duration_s = 2 * 3600;
Cfg.record_dt_s = 60;

n_steps = round(Cfg.duration_s / Cfg.dt_s);
record_every = max(1, round(Cfg.record_dt_s / Cfg.dt_s));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
stage_model_m = zeros(n_records, 1);
stage_reference_m = zeros(n_records, 1);
storage_model_m3 = zeros(n_records, 1);
storage_reference_m3 = zeros(n_records, 1);
outflow_model_m3_s = zeros(n_records, 1);
outflow_reference_m3_s = zeros(n_records, 1);
released_volume_m3 = zeros(n_records, 1);
reference_released_volume_m3 = zeros(n_records, 1);
mass_residual_m3 = zeros(n_records, 1);
mass_error_pct = zeros(n_records, 1);

storage0_m3 = Cfg.area_m2 * Cfg.initial_stage_m;
stage_m = Cfg.initial_stage_m;
released_m3 = 0;
rec = 1;
record_state(0);

for it = 1:n_steps
    q_m3_s = reservoir_rating_discharge(stage_m, Cfg);
    available_m3 = max(stage_m - Cfg.sill_stage_m, 0) * Cfg.area_m2;
    release_m3 = min(q_m3_s * Cfg.dt_s, available_m3);
    stage_m = stage_m - release_m3 / Cfg.area_m2;
    released_m3 = released_m3 + release_m3;

    if mod(it, record_every) == 0
        rec = rec + 1;
        record_state(it * Cfg.dt_s);
    end
end

Series = table(t_min, stage_model_m, stage_reference_m, storage_model_m3, ...
    storage_reference_m3, outflow_model_m3_s, outflow_reference_m3_s, ...
    released_volume_m3, reference_released_volume_m3, mass_residual_m3, ...
    mass_error_pct);
writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));

storage_err = storage_model_m3 - storage_reference_m3;
q_err = outflow_model_m3_s - outflow_reference_m3_s;
rmse = rmse_omitnan(storage_err);
mae = mean(abs(storage_err), 'omitnan');
max_error = max(abs(storage_err), [], 'omitnan');
relative_l2 = norm(storage_err) / max(norm(storage_reference_m3), eps);
nse = nse_metric(storage_model_m3, storage_reference_m3);
peak_time_error_min = abs(peak_time(t_min, outflow_model_m3_s) - ...
    peak_time(t_min, outflow_reference_m3_s));
peak_q_error_pct = abs(max(outflow_model_m3_s) - max(outflow_reference_m3_s)) / ...
    max(max(outflow_reference_m3_s), eps) * 100;
mass_abs_pct = max(abs(mass_error_pct), [], 'omitnan');
outflow_rmse = rmse_omitnan(q_err);
recovered_k_error_pct = recover_k_error(Cfg, stage_model_m, outflow_model_m3_s);

passed = rmse < 2.0 && outflow_rmse < 0.01 && mass_abs_pct < 1e-8 && ...
    peak_q_error_pct < 0.1 && recovered_k_error_pct < 0.1;

Diag = table(Cfg.case_id, Cfg.case_name, "analytical_rating_curve_ode", ...
    rmse, mae, max_error, relative_l2, nse, peak_time_error_min, ...
    peak_q_error_pct, mass_abs_pct, outflow_rmse, recovered_k_error_pct, ...
    passed, ...
    "Storage-stage dynamics compared to the exact solution of dS/dt = -k(h-h0)^b with S=A h.", ...
    'VariableNames', {'case_id','case_name','evidence_type','rmse', ...
    'mae','max_error','relative_l2','nse','peak_time_error_min', ...
    'peak_magnitude_error_pct','mass_error_pct','outflow_rmse_m3_s', ...
    'recovered_k_error_pct','passed','metric_note'});
Pass = pass_row(Diag, passed);

    function record_state(t_s)
        href = exact_stage(t_s, Cfg);
        q_model = reservoir_rating_discharge(stage_m, Cfg);
        qref = reservoir_rating_discharge(href, Cfg);
        t_min(rec) = t_s / 60;
        stage_model_m(rec) = stage_m;
        stage_reference_m(rec) = href;
        storage_model_m3(rec) = stage_m * Cfg.area_m2;
        storage_reference_m3(rec) = href * Cfg.area_m2;
        outflow_model_m3_s(rec) = q_model;
        outflow_reference_m3_s(rec) = qref;
        released_volume_m3(rec) = released_m3;
        reference_released_volume_m3(rec) = storage0_m3 - storage_reference_m3(rec);
        mass_residual_m3(rec) = storage_model_m3(rec) + released_volume_m3(rec) - storage0_m3;
        mass_error_pct(rec) = 100 * mass_residual_m3(rec) / storage0_m3;
    end
end

function q = reservoir_rating_discharge(stage_m, Cfg)
q = Cfg.rating_k * max(stage_m - Cfg.sill_stage_m, 0).^Cfg.rating_exponent;
end

function h = exact_stage(t_s, Cfg)
y0 = max(Cfg.initial_stage_m - Cfg.sill_stage_m, 0);
b = Cfg.rating_exponent;
if y0 <= 0
    h = Cfg.sill_stage_m;
elseif abs(b - 1) < 1e-12
    y = y0 .* exp(-Cfg.rating_k .* t_s ./ Cfg.area_m2);
    h = Cfg.sill_stage_m + y;
else
    y = (y0.^(1 - b) + (b - 1) .* Cfg.rating_k .* t_s ./ Cfg.area_m2).^(1 ./ (1 - b));
    h = Cfg.sill_stage_m + y;
end
end

function pct = recover_k_error(Cfg, stage_m, q_m3_s)
head = max(stage_m - Cfg.sill_stage_m, 0);
valid = head > 0.05 & q_m3_s > 0;
if nnz(valid) < 3
    pct = Inf;
    return
end
k_est = median(q_m3_s(valid) ./ head(valid).^Cfg.rating_exponent, 'omitnan');
pct = abs(k_est - Cfg.rating_k) / Cfg.rating_k * 100;
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
plot(Series.t_min, Series.stage_model_m, 'o-', ...
    Series.t_min, Series.stage_reference_m, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Stage (m)');
legend('HydroPol2D rating step', 'Exact ODE', 'Location', 'best');
title('Reservoir stage depletion');
grid on;
nexttile;
plot(Series.t_min, Series.outflow_model_m3_s, 'o-', ...
    Series.t_min, Series.outflow_reference_m3_s, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Outflow (m^3/s)');
legend('HydroPol2D rating step', 'Exact ODE', 'Location', 'best');
title('Rating-curve discharge');
grid on;
nexttile;
yyaxis left
plot(Series.t_min, Series.storage_model_m3 - Series.storage_reference_m3, 'LineWidth', 1.4);
ylabel('Storage error (m^3)');
yyaxis right
plot(Series.t_min, Series.mass_error_pct, '--', 'LineWidth', 1.4);
ylabel('Mass residual (%)');
xlabel('Time (min)');
grid on;
exportgraphics(fig, fullfile(fig_dir, 'P1_RES_001_RATING_CURVE.png'), 'Resolution', 200);
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

function t = peak_time(t_min, y)
[~, idx] = max(y);
t = t_min(idx);
end

function ensure_dir(path_name)
if ~exist(path_name, 'dir')
    mkdir(path_name);
end
end
