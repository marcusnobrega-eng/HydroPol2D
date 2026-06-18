clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
ts_dir = fullfile(out_dir, 'TimeSeries');
fig_dir = fullfile(out_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(ts_dir, 'dir'); mkdir(ts_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

Diagnostics = table();
PassFail = table();

[Diag, Series, passed] = run_constant_recharge_case(functions_dir, true);
writetable(Series, fullfile(ts_dir, 'P1-GW-ASYNC-001.csv'));
Diagnostics = [Diagnostics; Diag]; %#ok<AGROW>
PassFail = [PassFail; pass_row(Diag, passed)]; %#ok<AGROW>

[Diag, Series, passed] = run_scheduler_equivalence_case(functions_dir);
writetable(Series, fullfile(ts_dir, 'P1-GW-ASYNC-002.csv'));
Diagnostics = [Diagnostics; Diag]; %#ok<AGROW>
PassFail = [PassFail; pass_row(Diag, passed)]; %#ok<AGROW>

[Diag, Series, passed] = run_capillary_case(functions_dir);
writetable(Series, fullfile(ts_dir, 'P1-GW-CAP-001.csv'));
Diagnostics = [Diagnostics; Diag]; %#ok<AGROW>
PassFail = [PassFail; pass_row(Diag, passed)]; %#ok<AGROW>

[Diag, Series, passed] = run_boussinesq_mound_case(functions_dir);
writetable(Series, fullfile(ts_dir, 'P1-GW-MOUND-001.csv'));
Diagnostics = [Diagnostics; Diag]; %#ok<AGROW>
PassFail = [PassFail; pass_row(Diag, passed)]; %#ok<AGROW>

writetable(Diagnostics, fullfile(out_dir, 'VTilted_Groundwater_Diagnostics.csv'));
writetable(PassFail, fullfile(out_dir, 'VTilted_Groundwater_Pass_Fail.csv'));
make_figures(ts_dir, fig_dir);

disp(Diagnostics);
disp(PassFail);

function [Diag, Series, passed] = run_constant_recharge_case(functions_dir, async_enabled)
Cfg = base_cfg();
Cfg.case_id = "P1-GW-ASYNC-001";
Cfg.case_name = "Constant recharge water-table rise";
Cfg.regime = "local_recharge";
Cfg.flag_baseflow = 0;
Cfg.flag_groundwater_async = double(async_enabled);
Cfg.flag_capillary_rise = 0;
Cfg.target_dt_min = 60;
Cfg.dt_min = 10;
Cfg.duration_min = 360;
Cfg.recharge_mm_h = 1.0;
Cfg.ksat_mm_h = 1e10;

[State, Env] = make_model_state(Cfg);
[Series, State] = run_forcing_series(State, Env, Cfg, functions_dir);

total_recharge_m = Cfg.recharge_mm_h * Cfg.duration_min / 60 / 1000;
expected_head_m = Env.initial_head_mean_m + total_recharge_m / Cfg.Sy;
model_head_m = mean_omitnan(State.BC_States.h_t);
head_error_m = model_head_m - expected_head_m;

initial_gw = Env.initial_gw_storage_m3;
final_gw = groundwater_storage_m3(State, Env);
input_volume_m3 = total_recharge_m * Env.active_area_m2;
mass_residual_m3 = input_volume_m3 - (final_gw - initial_gw);

passed = abs(head_error_m) < 1e-10 && abs(mass_residual_m3) < 1e-6 && ...
    State.GW_States.n_updates == Cfg.duration_min / Cfg.target_dt_min;

Diag = diagnostic_row(Cfg, expected_head_m, model_head_m, head_error_m, ...
    0, 0, 0, input_volume_m3, State.GW_States.total_recharge_m3, ...
    State.GW_States.total_capillary_m3, final_gw - initial_gw, ...
    mass_residual_m3, State.GW_States.n_updates, passed);
end

function [Diag, Series, passed] = run_scheduler_equivalence_case(functions_dir)
CfgAsync = base_cfg();
CfgAsync.case_id = "P1-GW-ASYNC-002";
CfgAsync.case_name = "Recharge accumulation equivalence";
CfgAsync.regime = "scheduler_equivalence";
CfgAsync.flag_baseflow = 0;
CfgAsync.flag_groundwater_async = 1;
CfgAsync.flag_capillary_rise = 0;
CfgAsync.target_dt_min = 60;
CfgAsync.dt_min = 10;
CfgAsync.duration_min = 360;
CfgAsync.recharge_mm_h = 1.0;
CfgAsync.ksat_mm_h = 1e10;

CfgLegacy = CfgAsync;
CfgLegacy.flag_groundwater_async = 0;

[AsyncState, AsyncEnv] = make_model_state(CfgAsync);
[AsyncSeries, AsyncState] = run_forcing_series(AsyncState, AsyncEnv, CfgAsync, functions_dir);

[LegacyState, LegacyEnv] = make_model_state(CfgLegacy);
[LegacySeries, LegacyState] = run_forcing_series(LegacyState, LegacyEnv, CfgLegacy, functions_dir);

head_error_grid = AsyncState.BC_States.h_t - LegacyState.BC_States.h_t;
head_rmse_m = rmse_omitnan(head_error_grid);
max_head_error_m = max_abs_omitnan(head_error_grid);

async_storage = groundwater_storage_m3(AsyncState, AsyncEnv);
legacy_storage = groundwater_storage_m3(LegacyState, LegacyEnv);
storage_error_m3 = async_storage - legacy_storage;
relative_storage_error_pct = abs(storage_error_m3) / max(abs(legacy_storage), eps) * 100;

Series = table(AsyncSeries.time_h, AsyncSeries.mean_head_m, LegacySeries.mean_head_m, ...
    AsyncSeries.pending_net_mm, LegacySeries.pending_net_mm, ...
    'VariableNames', {'time_h','async_head_m','legacy_head_m', ...
    'async_pending_net_mm','legacy_pending_net_mm'});

passed = head_rmse_m < 1e-10 && max_head_error_m < 1e-10 && relative_storage_error_pct < 0.1;

Diag = diagnostic_row(CfgAsync, mean_omitnan(LegacyState.BC_States.h_t), ...
    mean_omitnan(AsyncState.BC_States.h_t), mean_omitnan(head_error_grid), ...
    head_rmse_m, max_head_error_m, relative_storage_error_pct, ...
    AsyncState.GW_States.total_recharge_m3, AsyncState.GW_States.total_recharge_m3, ...
    0, async_storage - AsyncEnv.initial_gw_storage_m3, storage_error_m3, ...
    AsyncState.GW_States.n_updates, passed);
end

function [Diag, Series, passed] = run_capillary_case(functions_dir)
Cfg = base_cfg();
Cfg.case_id = "P1-GW-CAP-001";
Cfg.case_name = "Capillary rise drawdown";
Cfg.regime = "capillary_rise";
Cfg.flag_baseflow = 0;
Cfg.flag_groundwater_async = 1;
Cfg.flag_capillary_rise = 1;
Cfg.target_dt_min = 60;
Cfg.dt_min = 60;
Cfg.duration_min = 60;
Cfg.recharge_mm_h = 0;
Cfg.ksat_mm_h = 12;
Cfg.initial_depth_to_gw_m = 0.50;
Cfg.root_depth_m = 1.50;

[State, Env] = make_model_state(Cfg);
[Series, State] = run_forcing_series(State, Env, Cfg, functions_dir);

expected_capillary_mm = Cfg.ksat_mm_h * (1 - Cfg.initial_depth_to_gw_m / 2.0);
expected_capillary_mm = max(expected_capillary_mm, 0);
expected_head_m = Env.initial_head_mean_m - expected_capillary_mm / 1000 / Cfg.Sy;
model_head_m = mean_omitnan(State.BC_States.h_t);
head_error_m = model_head_m - expected_head_m;

soil_storage_gain_m3 = vadose_storage_m3(State, Env);
gw_storage_change_m3 = groundwater_storage_m3(State, Env) - Env.initial_gw_storage_m3;
capillary_volume_m3 = expected_capillary_mm / 1000 * Env.active_area_m2;
mass_residual_m3 = soil_storage_gain_m3 + gw_storage_change_m3;

passed = abs(head_error_m) < 1e-10 && ...
    abs(State.GW_States.total_capillary_m3 - capillary_volume_m3) < 1e-6 && ...
    abs(mass_residual_m3) < 1e-6;

Diag = diagnostic_row(Cfg, expected_head_m, model_head_m, head_error_m, ...
    0, 0, 0, 0, State.GW_States.total_recharge_m3, ...
    State.GW_States.total_capillary_m3, gw_storage_change_m3, ...
    mass_residual_m3, State.GW_States.n_updates, passed);
end

function [Diag, Series, passed] = run_boussinesq_mound_case(functions_dir)
CfgFine = base_cfg();
CfgFine.case_id = "P1-GW-MOUND-001";
CfgFine.case_name = "Boussinesq recharge mound";
CfgFine.regime = "boussinesq_mound";
CfgFine.flag_baseflow = 1;
CfgFine.flag_groundwater_async = 0;
CfgFine.flag_capillary_rise = 0;
CfgFine.target_dt_min = 30;
CfgFine.dt_min = 30;
CfgFine.duration_min = 360;
CfgFine.recharge_mm_h = 2.0;
CfgFine.recharge_mask = "channel";
CfgFine.ksat_mm_h = 1e10;
CfgFine.ksat_gw_mm_h = 5;
CfgFine.grid_n = 7;

CfgAsync = CfgFine;
CfgAsync.flag_groundwater_async = 1;
CfgAsync.target_dt_min = 60;

[FineState, FineEnv] = make_model_state(CfgFine);
[FineSeries, FineState] = run_forcing_series(FineState, FineEnv, CfgFine, functions_dir);

[AsyncState, AsyncEnv] = make_model_state(CfgAsync);
[AsyncSeries, AsyncState] = run_forcing_series(AsyncState, AsyncEnv, CfgAsync, functions_dir);

head_error_grid = AsyncState.BC_States.h_t - FineState.BC_States.h_t;
head_rmse_m = rmse_omitnan(head_error_grid);
max_head_error_m = max_abs_omitnan(head_error_grid);
fine_storage = groundwater_storage_m3(FineState, FineEnv);
async_storage = groundwater_storage_m3(AsyncState, AsyncEnv);
storage_error_m3 = async_storage - fine_storage;
relative_storage_error_pct = abs(storage_error_m3) / max(abs(fine_storage), eps) * 100;

Series = table(AsyncSeries.time_h, AsyncSeries.mean_head_m, FineSeries.mean_head_m, ...
    AsyncSeries.mean_recharge_mm_h, FineSeries.mean_recharge_mm_h, ...
    'VariableNames', {'time_h','async_head_m','reference_head_m', ...
    'async_recharge_mm_h','reference_recharge_mm_h'});

passed = head_rmse_m < 5e-3 && max_head_error_m < 2e-2 && relative_storage_error_pct < 0.1;

Diag = diagnostic_row(CfgFine, mean_omitnan(FineState.BC_States.h_t), ...
    mean_omitnan(AsyncState.BC_States.h_t), mean_omitnan(head_error_grid), ...
    head_rmse_m, max_head_error_m, relative_storage_error_pct, ...
    AsyncState.GW_States.total_recharge_m3, AsyncState.GW_States.total_recharge_m3, ...
    0, async_storage - AsyncEnv.initial_gw_storage_m3, storage_error_m3, ...
    AsyncState.GW_States.n_updates, passed);
end

function [Series, State] = run_forcing_series(State, Env, Cfg, functions_dir)
n_steps = round(Cfg.duration_min / Cfg.dt_min);
time_h = zeros(n_steps, 1);
mean_head_m = zeros(n_steps, 1);
pending_net_mm = zeros(n_steps, 1);
mean_recharge_mm_h = zeros(n_steps, 1);
mean_capillary_mm_h = zeros(n_steps, 1);
gw_updates = zeros(n_steps, 1);

for it = 1:n_steps
    State.k = it;
    State.t = it * Cfg.dt_min;
    State.time_step = Cfg.dt_min;

    if Cfg.recharge_mm_h > 0
        step_recharge_mm = Cfg.recharge_mm_h * Cfg.dt_min / 60;
        State.Soil_Properties = configure_recharge_source_layer( ...
            State.Soil_Properties, State.idx_nan, step_recharge_mm);
        State.Soil_Properties.Layers.transmission_storage_mm(Env.recharge_mask) = ...
            State.Soil_Properties.Layers.transmission_storage_mm(Env.recharge_mask) + step_recharge_mm;
        State.Soil_Properties = sync_layered_soil_storage(State.Soil_Properties, State.idx_nan);
    end

    State = run_groundwater_module(State, functions_dir);

    time_h(it) = State.t / 60;
    mean_head_m(it) = mean_omitnan(State.BC_States.h_t);
    pending_net_mm(it) = mean_omitnan(State.GW_States.pending_net_exchange_m) * 1000;
    mean_recharge_mm_h(it) = mean_omitnan(State.GW_States.last_surface_recharge_rate_m_s) * 1000 * 3600;
    mean_capillary_mm_h(it) = mean_omitnan(State.GW_States.last_capillary_rate_m_s) * 1000 * 3600;
    gw_updates(it) = State.GW_States.n_updates;
end

Series = table(time_h, mean_head_m, pending_net_mm, mean_recharge_mm_h, ...
    mean_capillary_mm_h, gw_updates);
end

function State = run_groundwater_module(State, functions_dir)
fields = fieldnames(State);
for i = 1:numel(fields)
    eval([fields{i} ' = State.(fields{i});']); %#ok<EVLDOT>
end

run(fullfile(functions_dir, 'Groundwater_Module.m'));

for i = 1:numel(fields)
    if exist(fields{i}, 'var')
        State.(fields{i}) = eval(fields{i}); %#ok<EVLDOT>
    end
end
extra_fields = {'GW_States','cumulative_recharge','max_GW_depth','recharge_rate','q_exf','q_river'};
for i = 1:numel(extra_fields)
    if exist(extra_fields{i}, 'var')
        State.(extra_fields{i}) = eval(extra_fields{i}); %#ok<EVLDOT>
    end
end
end

function [State, Env] = make_model_state(Cfg)
[elevation, zone_id, idx_nan] = make_vtilted_domain(Cfg.grid_n);
cell_area = Cfg.resolution_m^2;

State = struct();
State.k = 1;
State.t = Cfg.dt_min;
State.time_step = Cfg.dt_min;
State.idx_nan = idx_nan;
State.idx_rivers = false(size(elevation));
State.elevation = elevation;
State.DEM_raster = struct('Z', elevation);
State.Elevation_Properties = struct('elevation_cell', elevation);
State.Wshed_Properties = struct();
State.Wshed_Properties.Resolution = Cfg.resolution_m;
State.Wshed_Properties.cell_area = cell_area;
State.Wshed_Properties.domain = ~idx_nan;
State.Wshed_Properties.perimeter = false(size(elevation));
State.Wshed_Properties.perimeter([1 end], :) = true;
State.Wshed_Properties.perimeter(:, [1 end]) = true;

State.running_control = struct('routing_time', Cfg.duration_min);

State.flags = struct();
State.flags.flag_groundwater_modeling = 1;
State.flags.flag_baseflow = Cfg.flag_baseflow;
State.flags.flag_groundwater_async = Cfg.flag_groundwater_async;
State.flags.flag_capillary_rise = Cfg.flag_capillary_rise;
State.flags.groundwater_target_dt_min = Cfg.target_dt_min;
State.flags.groundwater_min_dt_min = 1;
State.flags.groundwater_max_head_change_m = 0.25;
State.flags.groundwater_courant = 0.25;

State.LULC_Properties = struct();
State.LULC_Properties.idx_imp = false(size(elevation));
State.LULC_Properties.River_K_coeff = 1;
State.LULC_Properties.root_depth_m = Cfg.root_depth_m * ones(size(elevation));

State.Soil_Properties = make_soil(Cfg, elevation, idx_nan);
z_bed = elevation - State.Soil_Properties.Soil_Depth;
State.BC_States = struct();
State.BC_States.h_t = elevation - Cfg.initial_depth_to_gw_m;
State.BC_States.h_t = max(State.BC_States.h_t, z_bed);
State.BC_States.h_0 = State.BC_States.h_t;

State.depths = struct('d_t', zeros(size(elevation)));
State.Hydro_States = struct('f', zeros(size(elevation)));
State.errors = zeros(5, 1);

layer_options = struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005);
[State.Soil_Properties, State.LULC_Properties] = derive_layered_soil_profile( ...
    State.Soil_Properties, State.LULC_Properties, State.BC_States, ...
    elevation, idx_nan, layer_options);
State.Soil_Properties = initialize_layered_soil_storage(State.Soil_Properties, idx_nan);
State.Soil_Properties.Layers.near_surface_storage_mm(:) = 0;
State.Soil_Properties.Layers.root_zone_storage_mm(:) = 0;
State.Soil_Properties.Layers.transmission_storage_mm(:) = 0;
State.Soil_Properties = sync_layered_soil_storage(State.Soil_Properties, idx_nan);
State.Soil_Properties.I_0 = State.Soil_Properties.I_t;

if Cfg.recharge_mask == "channel"
    recharge_mask = zone_id == 2;
else
    recharge_mask = ~idx_nan;
end

Env = struct();
Env.zone_id = zone_id;
Env.recharge_mask = recharge_mask;
Env.active_area_m2 = sum(~idx_nan, 'all') * cell_area;
Env.initial_head_mean_m = mean_omitnan(State.BC_States.h_t);
Env.initial_gw_storage_m3 = groundwater_storage_m3(State, struct('active_area_m2', Env.active_area_m2));
end

function Soil = make_soil(Cfg, elevation, idx_nan)
Soil = struct();
template = ones(size(elevation));
Soil.theta_sat = Cfg.theta_sat * template;
Soil.theta_r = Cfg.theta_r * template;
Soil.theta_i = Cfg.theta_r * template;
Soil.Sy = Cfg.Sy * template;
Soil.Soil_Depth = Cfg.soil_depth_m * template;
Soil.ksat = Cfg.ksat_mm_h * template;
Soil.ksat_gw = Cfg.ksat_gw_mm_h * template;
Soil.alpha_vg = 3.5 * template;
Soil.n_vg = 1.45 * template;
Soil.Ks_multiplier_near_surface = template;
Soil.Ks_multiplier_root_zone = template;
Soil.Ks_multiplier_transmission = template;

fields = fieldnames(Soil);
for i = 1:numel(fields)
    Soil.(fields{i})(idx_nan) = nan;
end
Soil.I_t = zeros(size(elevation));
Soil.I_0 = Soil.I_t;
Soil.I_p = Soil.I_t;
end

function Soil_Properties = configure_recharge_source_layer(Soil_Properties, idx_nan, source_depth_mm)
layers = Soil_Properties.Layers;
capacity_mm = max(source_depth_mm, 1e-6);
theta_mobile = max(Soil_Properties.theta_sat - Soil_Properties.theta_r, 1e-6);
layers.transmission_capacity_mm(:) = capacity_mm;
layers.transmission_thickness_m = capacity_mm ./ theta_mobile ./ 1000;
layers.transmission_capacity_mm(idx_nan) = nan;
layers.transmission_thickness_m(idx_nan) = nan;
layers.total_vadose_capacity_mm = layers.near_surface_capacity_mm + ...
    layers.root_zone_capacity_mm + layers.transmission_capacity_mm;
Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end

function Cfg = base_cfg()
Cfg = struct();
Cfg.case_id = "";
Cfg.case_name = "";
Cfg.regime = "";
Cfg.grid_n = 3;
Cfg.resolution_m = 20;
Cfg.dt_min = 10;
Cfg.duration_min = 360;
Cfg.target_dt_min = 60;
Cfg.recharge_mm_h = 0;
Cfg.recharge_mask = "all";
Cfg.flag_baseflow = 0;
Cfg.flag_groundwater_async = 1;
Cfg.flag_capillary_rise = 0;
Cfg.initial_depth_to_gw_m = 2.0;
Cfg.soil_depth_m = 3.0;
Cfg.root_depth_m = 1.2;
Cfg.theta_sat = 0.43;
Cfg.theta_r = 0.08;
Cfg.Sy = 0.30;
Cfg.ksat_mm_h = 50;
Cfg.ksat_gw_mm_h = 5;
end

function [DEM, zone_id, idx_nan] = make_vtilted_domain(n)
if n <= 3
    DEM = [112 100 112; 111 100 111; 110 100 110];
    zone_id = [1 2 3; 1 2 3; 1 2 3];
else
    [x, y] = meshgrid(1:n, 1:n);
    center = ceil(n / 2);
    DEM = 100 + abs(x - center) * 2 + (n - y) * 0.2;
    zone_id = ones(n, n);
    zone_id(:, center) = 2;
    zone_id(:, center+1:end) = 3;
end
idx_nan = false(size(DEM));
end

function storage = groundwater_storage_m3(State, Env)
area = State.Wshed_Properties.cell_area;
z_bed = State.Elevation_Properties.elevation_cell - State.Soil_Properties.Soil_Depth;
storage = nansum(nansum(area .* State.Soil_Properties.Sy .* max(State.BC_States.h_t - z_bed, 0)));
if isfield(State, 'GW_States') && isfield(State.GW_States, 'pending_net_exchange_m')
    storage = storage + nansum(nansum(area .* State.GW_States.pending_net_exchange_m));
end
if nargin > 1 && isfield(Env, 'active_area_m2')
    storage = double(storage);
end
end

function storage = vadose_storage_m3(State, ~)
area = State.Wshed_Properties.cell_area;
storage = nansum(nansum(area .* State.Soil_Properties.I_t ./ 1000));
end

function row = diagnostic_row(Cfg, reference_head_m, model_head_m, mean_head_error_m, ...
    head_rmse_m, max_head_error_m, relative_storage_error_pct, input_volume_m3, ...
    model_recharge_volume_m3, capillary_volume_m3, groundwater_storage_change_m3, ...
    mass_residual_m3, groundwater_update_count, passed)
row = table(Cfg.case_id, Cfg.case_name, Cfg.regime, reference_head_m, ...
    model_head_m, mean_head_error_m, head_rmse_m, max_head_error_m, ...
    relative_storage_error_pct, input_volume_m3, model_recharge_volume_m3, ...
    capillary_volume_m3, groundwater_storage_change_m3, mass_residual_m3, ...
    groundwater_update_count, passed, ...
    'VariableNames', {'case_id','case_name','regime','reference_head_m', ...
    'model_head_m','mean_head_error_m','head_rmse_m','max_head_error_m', ...
    'relative_storage_error_pct','input_volume_m3','model_recharge_volume_m3', ...
    'capillary_volume_m3','groundwater_storage_change_m3','mass_residual_m3', ...
    'groundwater_update_count','passed'});
end

function row = pass_row(Diag, passed)
status = "fail";
if passed
    status = "pass";
end
row = table(Diag.case_id, Diag.case_name, Diag.regime, status, passed, ...
    'VariableNames', {'case_id','case_name','regime','status','report_ready'});
end

function make_figures(ts_dir, fig_dir)
plot_case(fullfile(ts_dir, 'P1-GW-ASYNC-001.csv'), fig_dir, ...
    'P1_GW_ASYNC_001_head_recharge.png', ...
    'Constant recharge water-table rise');
plot_case(fullfile(ts_dir, 'P1-GW-CAP-001.csv'), fig_dir, ...
    'P1_GW_CAP_001_capillary.png', ...
    'Capillary rise drawdown');

T = readtable(fullfile(ts_dir, 'P1-GW-ASYNC-002.csv'));
fig = figure('Visible', 'off');
tiledlayout(2, 1);
nexttile;
plot(T.time_h, T.async_head_m, '-o', T.time_h, T.legacy_head_m, '--');
ylabel('Head (m)'); legend('Async', 'Legacy', 'Location', 'best');
nexttile;
plot(T.time_h, T.async_pending_net_mm, '-o');
xlabel('Time (h)'); ylabel('Pending recharge (mm)');
title('Recharge accumulation equivalence');
exportgraphics(fig, fullfile(fig_dir, 'P1_GW_ASYNC_002_equivalence.png'), 'Resolution', 200);
close(fig);

T = readtable(fullfile(ts_dir, 'P1-GW-MOUND-001.csv'));
fig = figure('Visible', 'off');
tiledlayout(2, 1);
nexttile;
plot(T.time_h, T.async_head_m, '-o', T.time_h, T.reference_head_m, '--');
ylabel('Mean head (m)'); legend('Async', 'Fine reference', 'Location', 'best');
nexttile;
plot(T.time_h, T.async_recharge_mm_h, '-o', T.time_h, T.reference_recharge_mm_h, '--');
xlabel('Time (h)'); ylabel('Recharge (mm/h)');
title('Boussinesq recharge mound');
exportgraphics(fig, fullfile(fig_dir, 'P1_GW_MOUND_001_reference.png'), 'Resolution', 200);
close(fig);
end

function plot_case(path, fig_dir, filename, fig_title)
T = readtable(path);
fig = figure('Visible', 'off');
tiledlayout(2, 1);
nexttile;
plot(T.time_h, T.mean_head_m, '-o');
ylabel('Mean head (m)');
nexttile;
plot(T.time_h, T.mean_recharge_mm_h, '-o', T.time_h, T.mean_capillary_mm_h, '--');
xlabel('Time (h)'); ylabel('Flux (mm/h)');
legend('Recharge', 'Capillary rise', 'Location', 'best');
title(fig_title);
exportgraphics(fig, fullfile(fig_dir, filename), 'Resolution', 200);
close(fig);
end

function value = mean_omitnan(x)
value = mean(x(:), 'omitnan');
end

function value = rmse_omitnan(x)
value = sqrt(mean(x(:).^2, 'omitnan'));
end

function value = max_abs_omitnan(x)
value = max(abs(x(:)), [], 'omitnan');
end
