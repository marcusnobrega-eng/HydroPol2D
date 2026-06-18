clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

registry_path = fullfile(case_dir, 'Config', 'Infiltration_Case_Registry.csv');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(ts_dir, 'dir'); mkdir(ts_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

Cases = readtable(registry_path, 'TextType', 'string');

Diagnostics = table();
PassFail = table();

for i = 1:height(Cases)
    Cfg = Cases(i,:);

    if Cfg.regime == "spatial_soil_contrast"
        [DiagRow, passed] = run_spatial_contrast_case(Cfg, ts_dir);
    else
        [Result, Series] = run_column_case(Cfg);
        writetable(Series, fullfile(ts_dir, Cfg.case_id + ".csv"));
        [DiagRow, passed] = evaluate_case(Cfg, Result);
    end

    Diagnostics = [Diagnostics; DiagRow]; %#ok<AGROW>

    status = "fail";
    if passed
        status = "pass";
    end

    PassFail = [PassFail; table( ...
        Cfg.case_id, Cfg.case_name, Cfg.regime, status, passed, ...
        'VariableNames', {'case_id','case_name','regime','status','report_ready'})]; %#ok<AGROW>
end

writetable(Diagnostics, fullfile(out_dir, 'VTilted_Infiltration_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'VTilted_Infiltration_Pass_Fail.csv'));
make_infiltration_figures(Cases, ts_dir, fig_dir);

disp(Diagnostics);
disp(PassFail);

function [Result, Series] = run_column_case(Cfg)

idx_nan = false(1,1);
elevation = 100;
dt_h = Cfg.dt_min / 60;
n_steps = round(Cfg.duration_h / dt_h);
rain_depth_step_mm = Cfg.rainfall_mm_h * dt_h;

[Soil_Properties, LULC_Properties, BC_States] = make_soil_profile(Cfg, Cfg.ksat_mm_h, elevation);
Soil_Properties = initialize_layered_soil_storage(Soil_Properties, idx_nan);

surface_runoff_cum_mm = 0;
rain_cum_mm = 0;
infiltration_cum_mm = 0;
recharge_cum_mm = 0;
gw_excess_cum_mm = 0;
capacity_limited_steps = 0;
supply_limited_steps = 0;
storage_limited_steps = 0;

t_h = zeros(n_steps,1);
rain_mm = zeros(n_steps,1);
f_mm_h = zeros(n_steps,1);
capacity_mm_h = zeros(n_steps,1);
runoff_mm = zeros(n_steps,1);
storage_mm = zeros(n_steps,1);
near_mm = zeros(n_steps,1);
root_mm = zeros(n_steps,1);
trans_mm = zeros(n_steps,1);
recharge_mm = zeros(n_steps,1);
gw_depth_m = zeros(n_steps,1);
excess_mm = zeros(n_steps,1);

for it = 1:n_steps
    t_h(it) = it * dt_h;
    rain_mm(it) = rain_depth_step_mm;
    rain_cum_mm = rain_cum_mm + rain_depth_step_mm;

    if Cfg.infiltration_enabled == 0
        accepted_mm = 0;
        cap_mm_h = 0;
        surface_excess_mm = rain_depth_step_mm;
    else
        cap_mm_h = reference_infiltration_capacity_mm_h(Soil_Properties, rain_depth_step_mm);
        remaining_capacity_mm = layered_soil_remaining_capacity_mm(Soil_Properties, idx_nan);
        storage_limit_mm_h = remaining_capacity_mm / max(dt_h, eps);

        requested_f_mm_h = min([Cfg.rainfall_mm_h, cap_mm_h, storage_limit_mm_h]);
        [Soil_Properties, accepted_mm] = add_layered_infiltration( ...
            Soil_Properties, requested_f_mm_h * dt_h, idx_nan);
        surface_excess_mm = max(rain_depth_step_mm - accepted_mm, 0);

        if cap_mm_h < Cfg.rainfall_mm_h && cap_mm_h <= storage_limit_mm_h
            capacity_limited_steps = capacity_limited_steps + 1;
        elseif storage_limit_mm_h < Cfg.rainfall_mm_h && storage_limit_mm_h < cap_mm_h
            storage_limited_steps = storage_limited_steps + 1;
        else
            supply_limited_steps = supply_limited_steps + 1;
        end
    end

    infiltration_cum_mm = infiltration_cum_mm + accepted_mm;
    step_recharge_mm = 0;
    step_gw_excess_mm = 0;
    step_saturation_excess_mm = 0;

    if Cfg.enable_recharge == 1
        [recharge_rate, Soil_Properties, ~] = simulate_layered_groundwater_recharge( ...
            Soil_Properties, Cfg.dt_min * 60, idx_nan, zeros(1,1));

        step_recharge_mm = recharge_rate * Cfg.dt_min * 60 * 1000;
        recharge_cum_mm = recharge_cum_mm + step_recharge_mm;

        Sy_local = max(Cfg.theta_sat - Cfg.theta_r, 1e-6);
        h_candidate = BC_States.h_t + (step_recharge_mm / 1000) / Sy_local;
        step_gw_excess_mm = max(h_candidate - elevation, 0) * Sy_local * 1000;
        BC_States.h_t = min(h_candidate, elevation);
        gw_excess_cum_mm = gw_excess_cum_mm + step_gw_excess_mm;

        [Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
            Soil_Properties, LULC_Properties, BC_States, elevation, idx_nan, ...
            struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005));
        [Soil_Properties, step_saturation_excess_mm] = clamp_layered_soil_storage(Soil_Properties, idx_nan);
    end

    surface_runoff_cum_mm = surface_runoff_cum_mm + surface_excess_mm + ...
        step_gw_excess_mm + step_saturation_excess_mm;

    capacity_mm_h(it) = cap_mm_h;
    f_mm_h(it) = accepted_mm / max(dt_h, eps);
    runoff_mm(it) = surface_excess_mm + step_gw_excess_mm + step_saturation_excess_mm;
    storage_mm(it) = Soil_Properties.I_t;
    near_mm(it) = Soil_Properties.Layers.near_surface_storage_mm;
    root_mm(it) = Soil_Properties.Layers.root_zone_storage_mm;
    trans_mm(it) = Soil_Properties.Layers.transmission_storage_mm;
    recharge_mm(it) = step_recharge_mm;
    gw_depth_m(it) = elevation - BC_States.h_t;
    excess_mm(it) = step_gw_excess_mm + step_saturation_excess_mm;
end

total_capacity_mm = Soil_Properties.Layers.total_vadose_capacity_mm;
storage_closure_error_mm = rain_cum_mm - surface_runoff_cum_mm - ...
    recharge_cum_mm - Soil_Properties.I_t + Soil_Properties.I_0 - gw_excess_cum_mm;

Result = struct();
Result.rain_cum_mm = rain_cum_mm;
Result.infiltration_cum_mm = infiltration_cum_mm;
Result.runoff_cum_mm = surface_runoff_cum_mm;
Result.recharge_cum_mm = recharge_cum_mm;
Result.final_storage_mm = Soil_Properties.I_t;
Result.initial_storage_mm = Soil_Properties.I_0;
Result.total_capacity_mm = total_capacity_mm;
Result.storage_fraction = Soil_Properties.I_t / max(total_capacity_mm, eps);
Result.capacity_limited_fraction = capacity_limited_steps / n_steps;
Result.supply_limited_fraction = supply_limited_steps / n_steps;
Result.storage_limited_fraction = storage_limited_steps / n_steps;
Result.max_excess_mm = max(excess_mm);
Result.storage_closure_error_mm = storage_closure_error_mm;
Result.first_near_h = first_positive_time(t_h, near_mm);
Result.first_root_h = first_positive_time(t_h, root_mm);
Result.first_trans_h = first_positive_time(t_h, trans_mm);
Result.first_recharge_h = first_positive_time(t_h, recharge_mm);

Series = table(t_h, rain_mm, capacity_mm_h, f_mm_h, runoff_mm, storage_mm, ...
    near_mm, root_mm, trans_mm, recharge_mm, gw_depth_m, excess_mm);
end

function [DiagRow, passed] = evaluate_case(Cfg, Result)

rain = max(Result.rain_cum_mm, eps);
infiltration_ratio = Result.infiltration_cum_mm / rain;
runoff_ratio = Result.runoff_cum_mm / rain;
recharge_ratio = Result.recharge_cum_mm / rain;
closure_error_mm = Result.storage_closure_error_mm;

switch Cfg.regime
    case "no_infiltration"
        passed = abs(Result.infiltration_cum_mm) < 1e-9 && ...
            abs(runoff_ratio - 1) < 1e-9;
    case "supply_limited"
        passed = infiltration_ratio > 0.95 && runoff_ratio < 0.05;
    case "capacity_limited"
        passed = runoff_ratio > 0.15 && Result.capacity_limited_fraction > 0.50;
    case "storage_limited"
        passed = Result.storage_fraction > 0.98 && runoff_ratio > 0.30 && ...
            Result.storage_limited_fraction > 0;
    case "layered_percolation"
        passed = Result.first_near_h <= Result.first_root_h && ...
            Result.first_root_h <= Result.first_trans_h && ...
            isfinite(Result.first_trans_h) && isfinite(Result.first_recharge_h) && ...
            Result.recharge_cum_mm > 1.0 && recharge_ratio > 0.01;
    case "shallow_groundwater"
        passed = runoff_ratio > 0.30 && Result.max_excess_mm > 0;
    otherwise
        passed = false;
end

DiagRow = table( ...
    Cfg.case_id, Cfg.case_name, Cfg.regime, ...
    Result.rain_cum_mm, Result.infiltration_cum_mm, Result.runoff_cum_mm, ...
    Result.recharge_cum_mm, Result.final_storage_mm, Result.total_capacity_mm, ...
    infiltration_ratio, runoff_ratio, recharge_ratio, ...
    Result.capacity_limited_fraction, Result.storage_limited_fraction, ...
    Result.first_near_h, Result.first_root_h, Result.first_trans_h, ...
    Result.first_recharge_h, Result.max_excess_mm, closure_error_mm, passed, ...
    'VariableNames', {'case_id','case_name','regime','rain_cum_mm', ...
    'infiltration_cum_mm','runoff_cum_mm','recharge_cum_mm', ...
    'final_storage_mm','total_capacity_mm','infiltration_ratio', ...
    'runoff_ratio','recharge_ratio','capacity_limited_fraction', ...
    'storage_limited_fraction','first_near_h','first_root_h', ...
    'first_trans_h','first_recharge_h','max_excess_mm', ...
    'closure_error_mm','passed'});
end

function [DiagRow, passed] = run_spatial_contrast_case(Cfg, ts_dir)

LowCfg = Cfg;
HighCfg = Cfg;
LowCfg.case_id = Cfg.case_id + "_lowKs";
HighCfg.case_id = Cfg.case_id + "_highKs";
HighCfg.ksat_mm_h = Cfg.ksat_alt_mm_h;

[Low, LowSeries] = run_column_case(LowCfg);
[High, HighSeries] = run_column_case(HighCfg);

writetable(LowSeries, fullfile(ts_dir, Cfg.case_id + "_lowKs.csv"));
writetable(HighSeries, fullfile(ts_dir, Cfg.case_id + "_highKs.csv"));

rain = max(Low.rain_cum_mm, eps);
low_runoff_ratio = Low.runoff_cum_mm / rain;
high_runoff_ratio = High.runoff_cum_mm / rain;
low_infiltration_ratio = Low.infiltration_cum_mm / rain;
high_infiltration_ratio = High.infiltration_cum_mm / rain;

passed = low_runoff_ratio > high_runoff_ratio && ...
    high_infiltration_ratio > low_infiltration_ratio;

DiagRow = table( ...
    Cfg.case_id, Cfg.case_name, Cfg.regime, ...
    Low.rain_cum_mm, 0.5 * (Low.infiltration_cum_mm + High.infiltration_cum_mm), ...
    0.5 * (Low.runoff_cum_mm + High.runoff_cum_mm), ...
    0.5 * (Low.recharge_cum_mm + High.recharge_cum_mm), ...
    0.5 * (Low.final_storage_mm + High.final_storage_mm), ...
    0.5 * (Low.total_capacity_mm + High.total_capacity_mm), ...
    0.5 * (low_infiltration_ratio + high_infiltration_ratio), ...
    0.5 * (low_runoff_ratio + high_runoff_ratio), 0, ...
    Low.capacity_limited_fraction, 0, ...
    Low.first_near_h, Low.first_root_h, Low.first_trans_h, Low.first_recharge_h, ...
    low_runoff_ratio - high_runoff_ratio, ...
    (Low.storage_closure_error_mm + High.storage_closure_error_mm) / 2, passed, ...
    'VariableNames', {'case_id','case_name','regime','rain_cum_mm', ...
    'infiltration_cum_mm','runoff_cum_mm','recharge_cum_mm', ...
    'final_storage_mm','total_capacity_mm','infiltration_ratio', ...
    'runoff_ratio','recharge_ratio','capacity_limited_fraction', ...
    'storage_limited_fraction','first_near_h','first_root_h', ...
    'first_trans_h','first_recharge_h','max_excess_mm', ...
    'closure_error_mm','passed'});
end

function [Soil_Properties, LULC_Properties, BC_States] = make_soil_profile(Cfg, ksat_mm_h, elevation)

Soil_Properties = struct();
Soil_Properties.Soil_Depth = Cfg.soil_depth_m;
Soil_Properties.theta_sat = Cfg.theta_sat;
Soil_Properties.theta_r = Cfg.theta_r;
Soil_Properties.theta_i = Cfg.theta_i;
Soil_Properties.alpha_vg = Cfg.alpha_vg;
Soil_Properties.n_vg = Cfg.n_vg;
Soil_Properties.ksat = ksat_mm_h;
Soil_Properties.Ks_multiplier_near_surface = Cfg.ks_mult_near;
Soil_Properties.Ks_multiplier_root_zone = Cfg.ks_mult_root;
Soil_Properties.Ks_multiplier_transmission = Cfg.ks_mult_trans;
Soil_Properties.I_t = 0;
Soil_Properties.I_p = 0;
Soil_Properties.I_0 = 0;

LULC_Properties = struct();
LULC_Properties.root_depth_m = Cfg.root_depth_m;
LULC_Properties.idx_imp = false(1,1);

BC_States = struct();
BC_States.h_t = elevation - Cfg.groundwater_depth_m;

[Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
    Soil_Properties, LULC_Properties, BC_States, elevation, false(1,1), ...
    struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005));
end

function cap_mm_h = reference_infiltration_capacity_mm_h(Soil_Properties, available_depth_mm)

layers = Soil_Properties.Layers;

if layers.near_surface_capacity_mm <= 0
    cap_mm_h = 0;
    return
end

min_layer = 0.005;
if isfield(layers, 'min_layer_thickness_m') && ~isempty(layers.min_layer_thickness_m)
    min_layer = layers.min_layer_thickness_m;
end
if layers.near_surface_thickness_m < min_layer
    cap_mm_h = 0;
    return
end

theta_r = layers.theta_r_near_surface;
theta_s = layers.theta_sat_near_surface;
n = layers.n_vg_near_surface;
m = 1 - 1 / n;
a = layers.alpha_vg_near_surface;
Ks = layers.ksat_near_surface / 1000 / 3600;

near_thick = max(layers.near_surface_thickness_m, 1e-6);
Ltop = min(0.05, near_thick);
Ltop = max(Ltop, 1e-6);

S_near_m = max(layers.near_surface_storage_mm, 0) / 1000;
theta_bucket = theta_r + S_near_m / near_thick;
theta_bucket = min(max(theta_bucket, theta_r), theta_s);

Se_bucket = (theta_bucket - theta_r) / max(theta_s - theta_r, 1e-12);
Se_bucket = min(max(Se_bucket, 1e-6), 1);
h_bucket = -(1 / a) * ((Se_bucket ^ (-1 / m) - 1) ^ (1 / n));
if Se_bucket >= 0.999999
    h_bucket = 0;
end

w_avail = max(available_depth_mm, 0) / 1000;
w_top_deficit = max(layers.near_surface_capacity_mm - layers.near_surface_storage_mm, 0) / 1000;
w_wet = min(w_avail, w_top_deficit);

theta_top = theta_bucket + w_wet / Ltop;
theta_top = min(max(theta_top, theta_r), theta_s);
Se_top = (theta_top - theta_r) / max(theta_s - theta_r, 1e-12);
Se_top = min(max(Se_top, 1e-6), 1);

term_top = 1 - Se_top ^ (1 / m);
term_top = min(max(term_top, 0), 1);
Kr_top = Se_top ^ 0.5 * (1 - term_top ^ m) ^ 2;
Kr_top = min(max(Kr_top, 0), 1);

Ktop = Ks * Kr_top;
Ksoil = 0.5 * Ks + 0.5 * Ktop;
dh = max((available_depth_mm / 1000) - h_bucket, 0);
dh = min(dh, 1.0);
grad = dh / Ltop + 1;

cap_mm_h = Ksoil * grad * 1000 * 3600;
cap_mm_h = max(cap_mm_h, 0);
end

function t_first = first_positive_time(t_h, y)

idx = find(y > 1e-9, 1, 'first');
if isempty(idx)
    t_first = nan;
else
    t_first = t_h(idx);
end
end

function make_infiltration_figures(Cases, ts_dir, fig_dir)

close all force
for i = 1:height(Cases)
    Cfg = Cases(i,:);
    make_case_figure(Cfg, ts_dir, fig_dir);
end
make_overview_hydrograph_figure(Cases, ts_dir, fig_dir);
make_overview_infiltration_figure(Cases, ts_dir, fig_dir);
end

function make_case_figure(Cfg, ts_dir, fig_dir)

[Series, label_suffix] = load_plot_series(Cfg, ts_dir);

f = figure('Visible','off', 'Color','w', 'Units','pixels', ...
    'Position', [100 100 1200 780]);
tiledlayout(f, 3, 1, 'TileSpacing','compact', 'Padding','compact');

nexttile
plot(Series.t_h, Series.rain_rate_mm_h, 'k-', 'LineWidth', 1.5); hold on
plot(Series.t_h, Series.runoff_rate_mm_h, 'b-', 'LineWidth', 1.5);
plot(Series.t_h, Series.recharge_rate_mm_h, 'Color', [0.1 0.55 0.1], 'LineWidth', 1.5);
if ismember('runoff_rate_high_mm_h', Series.Properties.VariableNames)
    plot(Series.t_h, Series.runoff_rate_high_mm_h, 'b--', 'LineWidth', 1.2);
end
grid on
ylabel('Hydrograph (mm/h)')
title(sprintf('%s: %s%s', Cfg.case_id, Cfg.case_name, label_suffix), ...
    'Interpreter','none')
legend_entries = {'Rainfall','Runoff/excess','Recharge'};
if ismember('runoff_rate_high_mm_h', Series.Properties.VariableNames)
    legend_entries = {'Rainfall','Runoff/excess low Ks','Recharge low Ks', ...
        'Runoff/excess high Ks'};
end
legend(legend_entries, 'Location','best')
set_rate_ylim(gca, [Series.rain_rate_mm_h; Series.runoff_rate_mm_h; ...
    Series.recharge_rate_mm_h])

nexttile
plot(Series.t_h, Series.f_mm_h, 'r-', 'LineWidth', 1.6); hold on
plot(Series.t_h, Series.capacity_mm_h, 'Color', [0.45 0.2 0.7], ...
    'LineWidth', 1.3);
if ismember('f_high_mm_h', Series.Properties.VariableNames)
    plot(Series.t_h, Series.f_high_mm_h, 'r--', 'LineWidth', 1.2);
    plot(Series.t_h, Series.capacity_high_mm_h, '--', ...
        'Color', [0.45 0.2 0.7], 'LineWidth', 1.1);
    legend({'Infiltration low Ks','Capacity low Ks', ...
        'Infiltration high Ks','Capacity high Ks'}, 'Location','best')
else
    legend({'Accepted infiltration','Capacity'}, 'Location','best')
end
grid on
ylabel('Infiltration (mm/h)')
set_rate_ylim(gca, [Series.f_mm_h; Series.capacity_mm_h])

nexttile
plot(Series.t_h, Series.near_mm, 'LineWidth', 1.4); hold on
plot(Series.t_h, Series.root_mm, 'LineWidth', 1.4);
plot(Series.t_h, Series.trans_mm, 'LineWidth', 1.4);
plot(Series.t_h, Series.storage_mm, 'k-', 'LineWidth', 1.2);
grid on
xlabel('Time (h)')
ylabel('Storage (mm)')
legend({'Near surface','Root zone','Transmission','Total'}, 'Location','best')

safe_id = matlab.lang.makeValidName(char(Cfg.case_id));
exportgraphics(f, fullfile(fig_dir, [safe_id '_hydrograph_infiltration.png']), ...
    'Resolution', 220);
exportgraphics(f, fullfile(fig_dir, [safe_id '_hydrograph_infiltration.pdf']), ...
    'ContentType','vector');
close(f)
end

function make_overview_hydrograph_figure(Cases, ts_dir, fig_dir)

f = figure('Visible','off', 'Color','w', 'Units','pixels', ...
    'Position', [100 100 1500 1050]);
tiledlayout(f, 4, 2, 'TileSpacing','compact', 'Padding','compact');

for i = 1:height(Cases)
    Cfg = Cases(i,:);
    Series = load_plot_series(Cfg, ts_dir);
    nexttile
    plot(Series.t_h, Series.rain_rate_mm_h, 'k-', 'LineWidth', 1.0); hold on
    plot(Series.t_h, Series.runoff_rate_mm_h, 'b-', 'LineWidth', 1.2);
    plot(Series.t_h, Series.recharge_rate_mm_h, 'Color', [0.1 0.55 0.1], ...
        'LineWidth', 1.1);
    if ismember('runoff_rate_high_mm_h', Series.Properties.VariableNames)
        plot(Series.t_h, Series.runoff_rate_high_mm_h, 'b--', 'LineWidth', 1.0);
    end
    grid on
    title(sprintf('%s %s', Cfg.case_id, Cfg.regime), 'Interpreter','none')
    xlabel('Time (h)')
    ylabel('Rate (mm/h)')
    set_rate_ylim(gca, [Series.rain_rate_mm_h; Series.runoff_rate_mm_h; ...
        Series.recharge_rate_mm_h])
end

nexttile
axis off
text(0, 0.7, sprintf(['Hydrograph overview: black = rainfall, blue = ' ...
    'runoff/excess, green = recharge. Dashed blue is high-Ks runoff ' ...
    'for the spatial contrast case.']), 'FontSize', 11)

exportgraphics(f, fullfile(fig_dir, 'VTilted_Infiltration_Hydrograph_Overview.png'), ...
    'Resolution', 220);
exportgraphics(f, fullfile(fig_dir, 'VTilted_Infiltration_Hydrograph_Overview.pdf'), ...
    'ContentType','vector');
close(f)
end

function make_overview_infiltration_figure(Cases, ts_dir, fig_dir)

f = figure('Visible','off', 'Color','w', 'Units','pixels', ...
    'Position', [100 100 1500 1050]);
tiledlayout(f, 4, 2, 'TileSpacing','compact', 'Padding','compact');

for i = 1:height(Cases)
    Cfg = Cases(i,:);
    Series = load_plot_series(Cfg, ts_dir);
    nexttile
plot(Series.t_h, Series.rain_rate_mm_h, 'k-', 'LineWidth', 1.0); hold on
plot(Series.t_h, Series.f_mm_h, 'r-', 'LineWidth', 1.2);
if ismember('f_high_mm_h', Series.Properties.VariableNames)
    plot(Series.t_h, Series.f_high_mm_h, 'r--', 'LineWidth', 1.0);
end
    grid on
    title(sprintf('%s %s', Cfg.case_id, Cfg.regime), 'Interpreter','none')
    xlabel('Time (h)')
    ylabel('Rate (mm/h)')
    if ismember('f_high_mm_h', Series.Properties.VariableNames)
        y_values = [Series.rain_rate_mm_h; Series.f_mm_h; Series.f_high_mm_h];
    else
        y_values = [Series.rain_rate_mm_h; Series.f_mm_h];
    end
    set_rate_ylim(gca, y_values)
end

nexttile
axis off
text(0, 0.75, {'Infiltration-rate overview:', ...
    'black = rainfall', ...
    'red = accepted infiltration', ...
    'dashed red = high-Ks accepted infiltration in spatial contrast case'}, ...
    'FontSize', 11)

exportgraphics(f, fullfile(fig_dir, 'VTilted_Infiltration_Rate_Overview.png'), ...
    'Resolution', 220);
exportgraphics(f, fullfile(fig_dir, 'VTilted_Infiltration_Rate_Overview.pdf'), ...
    'ContentType','vector');
close(f)
end

function [Series, label_suffix] = load_plot_series(Cfg, ts_dir)

label_suffix = "";
if Cfg.regime == "spatial_soil_contrast"
    Low = readtable(fullfile(ts_dir, Cfg.case_id + "_lowKs.csv"));
    High = readtable(fullfile(ts_dir, Cfg.case_id + "_highKs.csv"));
    Series = convert_depths_to_rates(Low);
    High = convert_depths_to_rates(High);
    Series.f_high_mm_h = High.f_mm_h;
    Series.capacity_high_mm_h = High.capacity_mm_h;
    Series.runoff_rate_high_mm_h = High.runoff_rate_mm_h;
    label_suffix = " (low Ks solid, high Ks dashed)";
else
    Series = readtable(fullfile(ts_dir, Cfg.case_id + ".csv"));
    Series = convert_depths_to_rates(Series);
end
end

function Series = convert_depths_to_rates(Series)

dt_h = [Series.t_h(1); diff(Series.t_h)];
dt_h = max(dt_h, eps);
Series.rain_rate_mm_h = Series.rain_mm ./ dt_h;
Series.runoff_rate_mm_h = Series.runoff_mm ./ dt_h;
Series.recharge_rate_mm_h = Series.recharge_mm ./ dt_h;
end

function set_rate_ylim(ax, values)

values = values(isfinite(values));
if isempty(values)
    ylim(ax, [0 1])
    return
end
upper = max(values);
if upper <= 0
    ylim(ax, [0 1])
else
    ylim(ax, [0, upper * 1.10])
end
end
