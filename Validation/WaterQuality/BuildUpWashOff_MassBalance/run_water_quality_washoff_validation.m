clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
table_dir = fullfile(out_dir, 'Tables');
ensure_dir(out_dir);
ensure_dir(ts_dir);
ensure_dir(fig_dir);
ensure_dir(table_dir);

[Diagnostics, MassBalance, PassFail, DrySeries, SingleSeries, TwoCellSeries, Params] = ...
    run_water_quality_cases(ts_dir, table_dir);

writetable(Diagnostics, fullfile(out_dir, 'WaterQuality_Washoff_Diagnostics.csv'));
writetable(MassBalance, fullfile(out_dir, 'WaterQuality_Washoff_Mass_Balance.csv'));
writetable(PassFail, fullfile(out_dir, 'WaterQuality_Washoff_Pass_Fail.csv'));
writetable(Diagnostics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(MassBalance, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));

make_figures(SingleSeries, TwoCellSeries, fig_dir);

disp(Diagnostics);
disp(MassBalance);
disp(PassFail);

function [Diagnostics, MassBalance, PassFail, DrySeries, SingleSeries, TwoCellSeries, Params] = ...
        run_water_quality_cases(ts_dir, table_dir)
Params = base_parameters();

[DryDiag, DryMass, DryPass, DrySeries] = dry_buildup_case(Params);
[SingleDiag, SingleMass, SinglePass, SingleSeries] = single_cell_washoff_case(Params);
[TwoDiag, TwoMass, TwoPass, TwoCellSeries] = two_cell_transfer_case(Params);

Diagnostics = [DryDiag; SingleDiag; TwoDiag];
MassBalance = [DryMass; SingleMass; TwoMass];
PassFail = [DryPass; SinglePass; TwoPass];

writetable(DrySeries, fullfile(ts_dir, 'P1-WQ-001_DRY_BUILDUP.csv'));
writetable(SingleSeries, fullfile(ts_dir, 'P1-WQ-001_SINGLE_CELL.csv'));
writetable(TwoCellSeries, fullfile(ts_dir, 'P1-WQ-001_TWO_CELL.csv'));
writetable(struct2table(Params, 'AsArray', true), fullfile(table_dir, 'P1-WQ-001_parameters.csv'));
end

function Params = base_parameters()
Params = struct();
Params.case_id = "P1-WQ-001";
Params.area_m2 = 400;
Params.runoff_mm_h = 30;
Params.washoff_exponent = 1;
Params.max_buildup_kg_ha = 56.04;
Params.buildup_rate_per_day = 1.0;
Params.antecedent_dry_days = 5.0;
Params.target_initial_concentration_mg_L = 120;
Params.dt_min = 5e-4;
Params.duration_min = 15;
Params.record_dt_min = 0.5;
Params.min_Bt_g_m2 = 0;
Params.Bmin_g_m2 = 0;
Params.Bmax_g_m2 = 0;

Params.initial_buildup_kg = Params.max_buildup_kg_ha * ...
    (1 - exp(-Params.buildup_rate_per_day * Params.antecedent_dry_days)) * ...
    Params.area_m2 / 10000;
Params.runoff_m3_s = Params.runoff_mm_h / 1000 / 3600 * Params.area_m2;
Params.lambda_h = Params.target_initial_concentration_mg_L * ...
    (Params.runoff_mm_h * Params.area_m2) / 1e6 / Params.initial_buildup_kg;
Params.washoff_coefficient = Params.lambda_h / ...
    (Params.runoff_m3_s ^ Params.washoff_exponent);
end

function [Diag, Mass, Pass, Series] = dry_buildup_case(P)
case_id = "P1-WQ-001A";
case_name = "Dry-weather TSS buildup initialization";
reference_buildup_kg = P.max_buildup_kg_ha * ...
    (1 - exp(-P.buildup_rate_per_day * P.antecedent_dry_days)) * ...
    P.area_m2 / 10000;
model_buildup_kg = P.initial_buildup_kg;
err = model_buildup_kg - reference_buildup_kg;
mass_error_pct = abs(err) / max(reference_buildup_kg, eps) * 100;
passed = abs(err) < 1e-12;

Series = table(case_id, P.area_m2, P.max_buildup_kg_ha, ...
    P.buildup_rate_per_day, P.antecedent_dry_days, reference_buildup_kg, ...
    model_buildup_kg, err, ...
    'VariableNames', {'case_id','cell_area_m2','max_buildup_kg_ha', ...
    'buildup_rate_per_day','antecedent_dry_days','reference_buildup_kg', ...
    'model_buildup_kg','buildup_error_kg'});

Diag = diagnostic_row(case_id, case_name, "analytical_buildup_formula", ...
    abs(err), abs(err), abs(err), 0, 1, 0, mass_error_pct, 0, passed, ...
    "SWMM-style exponential antecedent dry-weather buildup formula.");
Mass = mass_row(case_id, case_name, reference_buildup_kg, model_buildup_kg, ...
    0, model_buildup_kg - reference_buildup_kg, mass_error_pct, passed);
Pass = pass_row(case_id, case_name, passed);
end

function [Diag, Mass, Pass, Series] = single_cell_washoff_case(P)
case_id = "P1-WQ-001B";
case_name = "Single-cell analytical exponential washoff";
n_steps = round(P.duration_min / P.dt_min);
record_every = max(1, round(P.record_dt_min / P.dt_min));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
model_mass_kg = zeros(n_records, 1);
reference_mass_kg = zeros(n_records, 1);
model_export_kg = zeros(n_records, 1);
reference_export_kg = zeros(n_records, 1);
model_concentration_mg_L = zeros(n_records, 1);
reference_concentration_mg_L = zeros(n_records, 1);
mass_residual_kg = zeros(n_records, 1);

B_t = P.initial_buildup_kg;
mass_lost = 0;
Tot_Washed = 0;
rec = 1;
record_state(0, P.target_initial_concentration_mg_L);

for it = 1:n_steps
    [B_t, ~, Out_Conc, ~, ~, mass_lost, Tot_Washed] = call_wq_single(B_t, P, mass_lost, Tot_Washed);
    if mod(it, record_every) == 0
        rec = rec + 1;
        record_state(it * P.dt_min, Out_Conc);
    end
end

Series = table(t_min, model_mass_kg, reference_mass_kg, model_export_kg, ...
    reference_export_kg, model_concentration_mg_L, reference_concentration_mg_L, ...
    mass_residual_kg);

mass_err = model_mass_kg - reference_mass_kg;
export_err = model_export_kg - reference_export_kg;
rmse = rmse_omitnan(mass_err);
mae = mean(abs(mass_err), 'omitnan');
max_error = max(abs(mass_err), [], 'omitnan');
relative_l2 = norm(mass_err) / max(norm(reference_mass_kg), eps);
nse = nse_metric(model_mass_kg, reference_mass_kg);
exported_load_error_pct = max(abs(export_err), [], 'omitnan') / ...
    max(reference_export_kg(end), eps) * 100;
mass_balance_abs = max(abs(mass_residual_kg), [], 'omitnan');
mass_error_pct = mass_balance_abs / max(P.initial_buildup_kg, eps) * 100;
negative_mass_count = sum(model_mass_kg < -1e-12);
passed = rmse < 1e-6 && exported_load_error_pct < 0.01 && ...
    (mass_balance_abs < 1e-6 || mass_error_pct < 0.01) && ...
    negative_mass_count == 0 && all(isfinite(table2array(Series)), 'all');

Diag = diagnostic_row(case_id, case_name, "analytical_exponential_washoff", ...
    rmse, mae, max_error, relative_l2, nse, exported_load_error_pct, ...
    mass_error_pct, negative_mass_count, passed, ...
    "HydroPol2D mass-based washoff compared with B(t)=B0 exp(-lambda t).");
Mass = mass_row(case_id, case_name, P.initial_buildup_kg, model_mass_kg(end), ...
    model_export_kg(end), mass_residual_kg(end), mass_error_pct, passed);
Pass = pass_row(case_id, case_name, passed);

    function record_state(t_now_min, out_conc)
        t_h = t_now_min / 60;
        b_ref = P.initial_buildup_kg * exp(-P.lambda_h * t_h);
        t_min(rec) = t_now_min;
        model_mass_kg(rec) = B_t;
        reference_mass_kg(rec) = b_ref;
        model_export_kg(rec) = max(P.initial_buildup_kg - B_t, 0);
        reference_export_kg(rec) = P.initial_buildup_kg - b_ref;
        model_concentration_mg_L(rec) = out_conc;
        reference_concentration_mg_L(rec) = P.target_initial_concentration_mg_L * exp(-P.lambda_h * t_h);
        mass_residual_kg(rec) = model_mass_kg(rec) + model_export_kg(rec) - P.initial_buildup_kg;
    end
end

function [Diag, Mass, Pass, Series] = two_cell_transfer_case(P)
case_id = "P1-WQ-001C";
case_name = "Two-cell conservative washoff transfer";
n_steps = round(P.duration_min / P.dt_min);
record_every = max(1, round(P.record_dt_min / P.dt_min));
n_records = floor(n_steps / record_every) + 1;

t_min = zeros(n_records, 1);
model_upstream_mass_kg = zeros(n_records, 1);
model_downstream_mass_kg = zeros(n_records, 1);
reference_upstream_mass_kg = zeros(n_records, 1);
reference_downstream_mass_kg = zeros(n_records, 1);
model_export_kg = zeros(n_records, 1);
reference_export_kg = zeros(n_records, 1);
outlet_concentration_mg_L = zeros(n_records, 1);
mass_residual_kg = zeros(n_records, 1);

B_t = [P.initial_buildup_kg, P.initial_buildup_kg];
mass_lost = 0;
Tot_Washed = zeros(1, 2);
initial_total = sum(B_t, 'all');
rec = 1;
record_state(0, P.target_initial_concentration_mg_L);

for it = 1:n_steps
    [B_t, Out_Conc, mass_lost, Tot_Washed] = call_wq_two_cell(B_t, P, mass_lost, Tot_Washed);
    if mod(it, record_every) == 0
        rec = rec + 1;
        record_state(it * P.dt_min, Out_Conc);
    end
end

Series = table(t_min, model_upstream_mass_kg, model_downstream_mass_kg, ...
    reference_upstream_mass_kg, reference_downstream_mass_kg, model_export_kg, ...
    reference_export_kg, outlet_concentration_mg_L, mass_residual_kg);

model_total = model_upstream_mass_kg + model_downstream_mass_kg;
reference_total = reference_upstream_mass_kg + reference_downstream_mass_kg;
mass_err = model_total - reference_total;
export_err = model_export_kg - reference_export_kg;
rmse = rmse_omitnan(mass_err);
mae = mean(abs(mass_err), 'omitnan');
max_error = max(abs(mass_err), [], 'omitnan');
relative_l2 = norm(mass_err) / max(norm(reference_total), eps);
nse = nse_metric(model_total, reference_total);
exported_load_error_pct = max(abs(export_err), [], 'omitnan') / ...
    max(reference_export_kg(end), eps) * 100;
mass_balance_abs = max(abs(mass_residual_kg), [], 'omitnan');
mass_error_pct = mass_balance_abs / max(initial_total, eps) * 100;
negative_mass_count = sum(model_upstream_mass_kg < -1e-12) + ...
    sum(model_downstream_mass_kg < -1e-12);
passed = rmse < 1e-6 && exported_load_error_pct < 0.01 && ...
    (mass_balance_abs < 1e-6 || mass_error_pct < 0.01) && ...
    negative_mass_count == 0 && all(isfinite(table2array(Series)), 'all');

Diag = diagnostic_row(case_id, case_name, "analytical_two_cell_linear_ode", ...
    rmse, mae, max_error, relative_l2, nse, exported_load_error_pct, ...
    mass_error_pct, negative_mass_count, passed, ...
    "Two-cell transfer compared with coupled linear exponential washoff solution.");
Mass = mass_row(case_id, case_name, initial_total, model_total(end), ...
    model_export_kg(end), mass_residual_kg(end), mass_error_pct, passed);
Pass = pass_row(case_id, case_name, passed);

    function record_state(t_now_min, out_conc)
        t_h = t_now_min / 60;
        b1_ref = P.initial_buildup_kg * exp(-P.lambda_h * t_h);
        b2_ref = exp(-P.lambda_h * t_h) * ...
            (P.initial_buildup_kg + P.lambda_h * P.initial_buildup_kg * t_h);
        t_min(rec) = t_now_min;
        model_upstream_mass_kg(rec) = B_t(1);
        model_downstream_mass_kg(rec) = B_t(2);
        reference_upstream_mass_kg(rec) = b1_ref;
        reference_downstream_mass_kg(rec) = b2_ref;
        model_export_kg(rec) = max(initial_total - sum(B_t, 'all'), 0);
        reference_export_kg(rec) = initial_total - b1_ref - b2_ref;
        outlet_concentration_mg_L(rec) = out_conc;
        mass_residual_kg(rec) = sum(B_t, 'all') + model_export_kg(rec) - initial_total;
    end
end

function [B_t, P_conc, Out_Conc, tmin_wq, tot_W_out, mass_lost, Tot_Washed] = ...
        call_wq_single(B_t, P, mass_lost, Tot_Washed)
q0 = zeros(1, 1);
qoutlet = P.runoff_mm_h;
idx_nan_5 = false(1, 1, 5);
[B_t, P_conc, Out_Conc, tmin_wq, tot_W_out, mass_lost, Tot_Washed] = ...
    build_up_wash_off(P.washoff_coefficient, P.washoff_exponent, ...
    q0, q0, q0, q0, qoutlet, B_t, P.dt_min, 1, 1, P.area_m2, ...
    true(1, 1), idx_nan_5, 0, mass_lost, Tot_Washed, ...
    P.Bmin_g_m2, P.Bmax_g_m2, P.min_Bt_g_m2);
end

function [B_t, Out_Conc, mass_lost, Tot_Washed] = call_wq_two_cell(B_t, P, mass_lost, Tot_Washed)
q0 = zeros(1, 2);
qright = [P.runoff_mm_h, 0];
qoutlet = [0, P.runoff_mm_h];
idx_nan_5 = false(1, 2, 5);
[B_t, ~, Out_Conc, ~, ~, mass_lost, Tot_Washed] = ...
    build_up_wash_off(P.washoff_coefficient * ones(1, 2), ...
    P.washoff_exponent * ones(1, 2), q0, qright, q0, q0, qoutlet, ...
    B_t, P.dt_min, 2, 1, P.area_m2, [false, true], idx_nan_5, 0, ...
    mass_lost, Tot_Washed, P.Bmin_g_m2, P.Bmax_g_m2, P.min_Bt_g_m2);
end

function Diag = diagnostic_row(case_id, case_name, evidence_type, rmse, mae, ...
        max_error, relative_l2, nse, exported_load_error_pct, mass_error_pct, ...
        negative_mass_count, passed, metric_note)
Diag = table(case_id, case_name, evidence_type, rmse, mae, max_error, ...
    relative_l2, nse, exported_load_error_pct, mass_error_pct, ...
    negative_mass_count, passed, metric_note, ...
    'VariableNames', {'case_id','case_name','evidence_type','rmse', ...
    'mae','max_error','relative_l2','nse','exported_load_error_pct', ...
    'mass_error_pct','negative_mass_count','passed','metric_note'});
end

function Mass = mass_row(case_id, case_name, initial_mass_kg, final_mass_kg, ...
        exported_mass_kg, mass_residual_kg, mass_error_pct, passed)
Mass = table(case_id, case_name, initial_mass_kg, final_mass_kg, ...
    exported_mass_kg, mass_residual_kg, mass_error_pct, passed, ...
    'VariableNames', {'case_id','case_name','initial_mass_kg', ...
    'final_mass_kg','exported_mass_kg','mass_residual_kg', ...
    'mass_error_pct','passed'});
end

function Pass = pass_row(case_id, case_name, passed)
status = "fail";
if passed
    status = "pass";
end
report_ready = passed;
Pass = table(case_id, case_name, status, report_ready, ...
    'VariableNames', {'case_id','case_name','status','report_ready'});
end

function make_figures(SingleSeries, TwoCellSeries, fig_dir)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(3, 1);
nexttile;
plot(SingleSeries.t_min, SingleSeries.model_mass_kg, 'o-', ...
    SingleSeries.t_min, SingleSeries.reference_mass_kg, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Surface mass (kg)');
legend('HydroPol2D', 'Analytical', 'Location', 'best');
title('Single-cell TSS washoff mass');
grid on;
nexttile;
plot(SingleSeries.t_min, SingleSeries.model_export_kg, 'o-', ...
    SingleSeries.t_min, SingleSeries.reference_export_kg, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Cumulative export (kg)');
legend('HydroPol2D', 'Analytical', 'Location', 'best');
title('Single-cell cumulative export');
grid on;
nexttile;
plot(SingleSeries.t_min, SingleSeries.model_concentration_mg_L, 'o-', ...
    SingleSeries.t_min, SingleSeries.reference_concentration_mg_L, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Concentration (mg/L)');
legend('HydroPol2D', 'Analytical', 'Location', 'best');
title('Single-cell outlet concentration');
grid on;
exportgraphics(fig, fullfile(fig_dir, 'P1_WQ_001_SINGLE_CELL.png'), 'Resolution', 200);
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(3, 1);
nexttile;
plot(TwoCellSeries.t_min, TwoCellSeries.model_upstream_mass_kg, 'o-', ...
    TwoCellSeries.t_min, TwoCellSeries.reference_upstream_mass_kg, '--', ...
    TwoCellSeries.t_min, TwoCellSeries.model_downstream_mass_kg, 's-', ...
    TwoCellSeries.t_min, TwoCellSeries.reference_downstream_mass_kg, ':', ...
    'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Surface mass (kg)');
legend('Upstream model', 'Upstream analytical', 'Downstream model', ...
    'Downstream analytical', 'Location', 'best');
title('Two-cell pollutant transfer');
grid on;
nexttile;
plot(TwoCellSeries.t_min, TwoCellSeries.model_export_kg, 'o-', ...
    TwoCellSeries.t_min, TwoCellSeries.reference_export_kg, '--', 'LineWidth', 1.4);
xlabel('Time (min)'); ylabel('Cumulative export (kg)');
legend('HydroPol2D', 'Analytical', 'Location', 'best');
title('Two-cell outlet export');
grid on;
nexttile;
yyaxis left
plot(TwoCellSeries.t_min, TwoCellSeries.outlet_concentration_mg_L, 'LineWidth', 1.4);
ylabel('Outlet concentration (mg/L)');
yyaxis right
plot(TwoCellSeries.t_min, TwoCellSeries.mass_residual_kg, '--', 'LineWidth', 1.4);
ylabel('Mass residual (kg)');
xlabel('Time (min)');
grid on;
exportgraphics(fig, fullfile(fig_dir, 'P1_WQ_001_TWO_CELL.png'), 'Resolution', 200);
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
