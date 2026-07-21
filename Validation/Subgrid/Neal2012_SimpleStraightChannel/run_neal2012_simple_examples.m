%RUN_NEAL2012_SIMPLE_EXAMPLES
% Simple straight-channel benchmarks for the Neal et al. (2012) channel-
% subgrid method. These examples isolate the subgrid effect using a one-
% cell-wide coarse domain, a resolved fine-channel reference, and a coarse
% Neal-mode run.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(fileparts(fileparts(case_dir)));
func_dir = fullfile(model_dir, 'HydroPol2D_Functions');
addpath(func_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
fig_dir = fullfile(case_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

geom = struct();
geom.channel_width_m = 10.0;
geom.bank_height_m = 1.0;
geom.coarse_dx_m = 30.0;
geom.fine_dx_m = 10.0;
geom.slope_mpm = 0.0010;
geom.n_channel = 0.035;
geom.n_floodplain = 0.035;

dt_min = 0.05;       % 3 s
record_dt_min = 1.0; % 1 min

cases = define_cases();

metric_rows = table();
mass_rows = table();
pass_rows = table();
time_rows = table();
profile_rows = table();

for icase = 1:numel(cases)
    cfg = cases(icase);
    fprintf('Running %s (%s)\n', cfg.id, cfg.name);

    Grid = build_case_grids(geom, cfg.n_cells_coarse);
    Result = run_case(cfg, Grid, dt_min, record_dt_min);

    [metric_case, mass_case, pass_case, profile_case] = summarize_case(Result, cfg);

    metric_rows = [metric_rows; metric_case]; %#ok<AGROW>
    mass_rows = [mass_rows; mass_case]; %#ok<AGROW>
    pass_rows = [pass_rows; pass_case]; %#ok<AGROW>
    time_rows = [time_rows; Result.TimeSeries]; %#ok<AGROW>
    profile_rows = [profile_rows; profile_case]; %#ok<AGROW>

    plot_hydrograph_and_stage(Result, cfg, fig_dir);
    plot_profiles(profile_case, cfg, fig_dir);
end

writetable(metric_rows, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(mass_rows, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(pass_rows, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(time_rows, fullfile(out_dir, 'Hydrograph_Stage_TimeSeries.csv'));
writetable(profile_rows, fullfile(out_dir, 'Profile_Comparisons.csv'));

fprintf('Simple Neal examples complete. %d of %d cases passed.\n', ...
    nnz(pass_rows.Pass), height(pass_rows));

%% ------------------------------------------------------------------------
function cases = define_cases()
cases(1).id = "P1-SUBGRID-NEAL-SIMPLE-001";
cases(1).name = "FiveCell_WithinBank";
cases(1).n_cells_coarse = 5;
cases(1).warmup_q_m3s = 0.0;
cases(1).warmup_duration_min = 0.0;
cases(1).duration_min = 120.0;
cases(1).qfun = @(t) 5.0 + 0.*t;
cases(1).comparison_state = "final";

cases(2).id = "P1-SUBGRID-NEAL-SIMPLE-002";
cases(2).name = "ThirtyCell_WithinBank";
cases(2).n_cells_coarse = 30;
cases(2).warmup_q_m3s = 0.0;
cases(2).warmup_duration_min = 0.0;
cases(2).duration_min = 420.0;
cases(2).qfun = @(t) 5.0 + 0.*t;
cases(2).comparison_state = "final";

cases(3).id = "P1-SUBGRID-NEAL-SIMPLE-003";
cases(3).name = "TenCell_WetStart_BankfullPulse";
cases(3).n_cells_coarse = 10;
cases(3).warmup_q_m3s = 5.0;
cases(3).warmup_duration_min = 120.0;
cases(3).duration_min = 220.0;
cases(3).qfun = @bankfull_pulse;
cases(3).comparison_state = "peak";
end

function q = bankfull_pulse(t_min)
t = t_min;
q = 5.0 .* ones(size(t));

idx = t <= 40;
q(idx) = 5.0 + (12.0 - 5.0) .* (t(idx) ./ 40.0);

idx = t > 40 & t <= 100;
q(idx) = 12.0;

idx = t > 100 & t <= 160;
q(idx) = 12.0 - (12.0 - 5.0) .* ((t(idx) - 100.0) ./ 60.0);
end

function Grid = build_case_grids(geom, n_cells_coarse)
ratio = round(geom.coarse_dx_m / geom.fine_dx_m);
n_cells_fine = n_cells_coarse * ratio;

Grid = struct();
Grid.geom = geom;
Grid.n_cells_coarse = n_cells_coarse;
Grid.n_cells_fine = n_cells_fine;
Grid.gauge_x_m = 0.5 * n_cells_coarse * geom.coarse_dx_m;

x_fine = ((1:n_cells_fine) - 0.5) .* geom.fine_dx_m;
x_coarse = ((1:n_cells_coarse) - 0.5) .* geom.coarse_dx_m;

z_fine_channel = -geom.bank_height_m - geom.slope_mpm .* x_fine;
z_coarse_plain = -geom.slope_mpm .* x_coarse;

Grid.fine = make_grid("fine_channel", geom.fine_dx_m, x_fine, z_fine_channel, ...
    geom.n_channel, zeros(1, n_cells_fine), zeros(1, n_cells_fine), ...
    geom.n_channel .* ones(1, n_cells_fine), geom.n_floodplain .* ones(1, n_cells_fine));

Grid.coarse = make_grid("ordinary", geom.coarse_dx_m, x_coarse, z_coarse_plain, ...
    geom.n_channel, zeros(1, n_cells_coarse), zeros(1, n_cells_coarse), ...
    geom.n_channel .* ones(1, n_cells_coarse), geom.n_floodplain .* ones(1, n_cells_coarse));

Grid.neal = make_grid("neal", geom.coarse_dx_m, x_coarse, z_coarse_plain, ...
    geom.n_channel, geom.channel_width_m .* ones(1, n_cells_coarse), ...
    geom.bank_height_m .* ones(1, n_cells_coarse), ...
    geom.n_channel .* ones(1, n_cells_coarse), geom.n_floodplain .* ones(1, n_cells_coarse));
end

function grid = make_grid(mode, dx, x, z_row, roughness_value, River_Width, River_Depth, n_channel, n_flood)
grid = struct();
grid.mode = mode;
grid.dx = dx;
grid.area = dx^2;
grid.x = x(:);
grid.nx = numel(x);
grid.ny = 1;
grid.z = reshape(z_row, 1, []);
grid.z_floodplain = grid.z;
if mode == "fine_channel"
    grid.z_floodplain = grid.z + abs(min(River_Depth(:))) + 1.0;
end
grid.roughness = roughness_value .* ones(1, grid.nx);
grid.River_Width = reshape(River_Width, 1, []);
grid.River_Depth = reshape(River_Depth, 1, []);
grid.n_channel = reshape(n_channel, 1, []);
grid.n_flood = reshape(n_flood, 1, []);
grid.row_outlet = 1;
grid.col_outlet = grid.nx;
grid.outlet_index = false(1, grid.nx);
grid.outlet_type = 1;
grid.slope_outlet = abs((grid.z(end-1) - grid.z(end)) / dx);
grid.inlet_row = 1;
grid.inlet_col = 1;
end

function Result = run_case(cfg, Grid, dt_min, record_dt_min)
dt_s = dt_min * 60;
nt = round(cfg.duration_min / dt_min);
record_every = max(1, round(record_dt_min / dt_min));
record_count = floor(nt / record_every) + 1;

fine_state = initialize_state(Grid.fine);
coarse_state = initialize_state(Grid.coarse);
neal_state = initialize_state(Grid.neal);

if cfg.warmup_duration_min > 0 && cfg.warmup_q_m3s > 0
    nw = round(cfg.warmup_duration_min / dt_min);
    for iw = 1:nw
        fine_state = apply_inflow(fine_state, Grid.fine, cfg.warmup_q_m3s, dt_s);
        coarse_state = apply_inflow(coarse_state, Grid.coarse, cfg.warmup_q_m3s, dt_s);
        neal_state = apply_inflow(neal_state, Grid.neal, cfg.warmup_q_m3s, dt_s);

        fine_state = step_model(fine_state, Grid.fine, dt_min);
        coarse_state = step_model(coarse_state, Grid.coarse, dt_min);
        neal_state = step_model(neal_state, Grid.neal, dt_min);
    end
end

records = initialize_records(record_count, cfg.id);
cum = struct('inflow_m3', 0.0, 'fine_outlet_m3', 0.0, ...
    'coarse_outlet_m3', 0.0, 'neal_outlet_m3', 0.0);
peak = struct('q_m3s', -inf, 'time_min', 0.0, 'fine_state', fine_state, ...
    'coarse_state', coarse_state, 'neal_state', neal_state);

rec = 1;
records = store_record(records, rec, 0.0, cfg.qfun(0.0), fine_state, coarse_state, ...
    neal_state, Grid);

for it = 1:nt
    t0 = (it - 1) * dt_min;
    qin = cfg.qfun(t0);

    fine_state = apply_inflow(fine_state, Grid.fine, qin, dt_s);
    coarse_state = apply_inflow(coarse_state, Grid.coarse, qin, dt_s);
    neal_state = apply_inflow(neal_state, Grid.neal, qin, dt_s);

    fine_state = step_model(fine_state, Grid.fine, dt_min);
    coarse_state = step_model(coarse_state, Grid.coarse, dt_min);
    neal_state = step_model(neal_state, Grid.neal, dt_min);

    qf = outlet_discharge(fine_state, Grid.fine);
    qc = outlet_discharge(coarse_state, Grid.coarse);
    qn = outlet_discharge(neal_state, Grid.neal);

    cum.inflow_m3 = cum.inflow_m3 + qin * dt_s;
    cum.fine_outlet_m3 = cum.fine_outlet_m3 + qf * dt_s;
    cum.coarse_outlet_m3 = cum.coarse_outlet_m3 + qc * dt_s;
    cum.neal_outlet_m3 = cum.neal_outlet_m3 + qn * dt_s;

    if qf > peak.q_m3s
        peak.q_m3s = qf;
        peak.time_min = it * dt_min;
        peak.fine_state = fine_state;
        peak.coarse_state = coarse_state;
        peak.neal_state = neal_state;
    end

    if mod(it, record_every) == 0 || it == nt
        rec = rec + 1;
        records = store_record(records, rec, it * dt_min, cfg.qfun(it * dt_min), ...
            fine_state, coarse_state, neal_state, Grid);
    end
end

records = trim_records(records, rec);

Result = struct();
Result.case_id = cfg.id;
Result.case_name = cfg.name;
Result.TimeSeries = make_time_table(records);
Result.cumulative = cum;
Result.Grid = Grid;
Result.final = struct('fine', fine_state, 'coarse', coarse_state, 'neal', neal_state);
Result.peak = peak;
end

function state = initialize_state(grid)
state = struct();
state.d_tot = zeros(1, grid.nx);
state.d_p = zeros(1, grid.nx);
state.outflow = zeros(1, grid.nx, 5);
state.outlet_flow = zeros(1, grid.nx);
state.Qc = zeros(1, grid.nx, 2);
state.Qf = zeros(1, grid.nx, 2);
state.Qci = zeros(1, grid.nx, 2);
state.Qfi = zeros(1, grid.nx, 2);
if grid.mode == "neal"
    state.C_a = hp2d_neal_cell_area(zeros(1, grid.nx), ...
        grid.River_Width, grid.River_Depth, grid.dx);
else
    state.C_a = grid.area .* ones(1, grid.nx);
end
end

function state = apply_inflow(state, grid, qin_m3s, dt_s)
if qin_m3s <= 0
    return;
end

if isscalar(state.C_a)
    active_area = max(state.C_a, eps);
else
    active_area = max(state.C_a(grid.inlet_row, grid.inlet_col), eps);
end
vol = qin_m3s * dt_s;
if grid.mode == "neal"
    row = grid.inlet_row;
    col = grid.inlet_col;
    h = hp2d_neal_apply_volume_change( ...
        max(state.d_tot(row, col) / 1000, 0), vol, ...
        grid.River_Width(row, col), grid.River_Depth(row, col), grid.dx);
    state.d_tot(row, col) = 1000 * h;
else
    state.d_tot(grid.inlet_row, grid.inlet_col) = state.d_tot(grid.inlet_row, grid.inlet_col) + ...
        1000 * vol / active_area;
end
end

function state = step_model(state, grid, dt_min)
idx_nan = false(1, grid.nx);

if grid.mode == "neal"
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,Qc,Qf,Qci,Qfi,C_a] = Local_Inertial_Model_D4( ...
        2, [], [], [], [], [], [], [], [], [], [], [], [], ...
        0, grid.z, state.d_tot, state.d_p, grid.roughness, grid.roughness.^2, grid.area, ...
        dt_min, grid.dx, grid.outlet_index, grid.outlet_type, grid.slope_outlet, ...
        grid.row_outlet, grid.col_outlet, 1e-6, state.outflow, idx_nan, 0, ...
        1, grid.n_channel, grid.n_flood, grid.River_Width, grid.River_Depth, ...
        state.Qc, state.Qf, state.Qci, state.Qfi, state.C_a, [], 1, 0, []);
    state.Qc = Qc;
    state.Qf = Qf;
    state.Qci = Qci;
    state.Qfi = Qfi;
else
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = Local_Inertial_Model_D4( ...
        2, [], [], [], [], [], [], [], [], [], [], [], [], ...
        0, grid.z, state.d_tot, state.d_p, grid.roughness, grid.roughness.^2, grid.area, ...
        dt_min, grid.dx, grid.outlet_index, grid.outlet_type, grid.slope_outlet, ...
        grid.row_outlet, grid.col_outlet, 1e-6, state.outflow, idx_nan, 0, ...
        0, grid.n_channel, grid.n_flood, grid.River_Width, grid.River_Depth, ...
        [], [], [], [], state.C_a, [], 0, 0, []);
end

state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function records = initialize_records(nrec, case_id)
records.case_id = repmat(string(case_id), nrec, 1);
records.time_min = nan(nrec, 1);
records.inflow_m3s = nan(nrec, 1);
records.outlet_fine_m3s = nan(nrec, 1);
records.outlet_coarse_m3s = nan(nrec, 1);
records.outlet_neal_m3s = nan(nrec, 1);
records.stage_fine_m = nan(nrec, 1);
records.stage_coarse_m = nan(nrec, 1);
records.stage_neal_m = nan(nrec, 1);
records.storage_fine_m3 = nan(nrec, 1);
records.storage_coarse_m3 = nan(nrec, 1);
records.storage_neal_m3 = nan(nrec, 1);
end

function records = store_record(records, rec, t_min, qin, fine_state, coarse_state, neal_state, Grid)
records.time_min(rec) = t_min;
records.inflow_m3s(rec) = qin;
records.outlet_fine_m3s(rec) = outlet_discharge(fine_state, Grid.fine);
records.outlet_coarse_m3s(rec) = outlet_discharge(coarse_state, Grid.coarse);
records.outlet_neal_m3s(rec) = outlet_discharge(neal_state, Grid.neal);
records.stage_fine_m(rec) = gauge_stage(fine_state, Grid.fine, Grid.gauge_x_m);
records.stage_coarse_m(rec) = gauge_stage(coarse_state, Grid.coarse, Grid.gauge_x_m);
records.stage_neal_m(rec) = gauge_stage(neal_state, Grid.neal, Grid.gauge_x_m);
records.storage_fine_m3(rec) = compute_storage(fine_state, Grid.fine);
records.storage_coarse_m3(rec) = compute_storage(coarse_state, Grid.coarse);
records.storage_neal_m3(rec) = compute_storage(neal_state, Grid.neal);
end

function records = trim_records(records, nrec)
fn = fieldnames(records);
for k = 1:numel(fn)
    v = records.(fn{k});
    if size(v, 1) == numel(records.time_min)
        records.(fn{k}) = v(1:nrec, :);
    end
end
end

function T = make_time_table(records)
T = table(records.case_id, records.time_min, records.inflow_m3s, ...
    records.outlet_fine_m3s, records.outlet_coarse_m3s, records.outlet_neal_m3s, ...
    records.stage_fine_m, records.stage_coarse_m, records.stage_neal_m, ...
    records.storage_fine_m3, records.storage_coarse_m3, records.storage_neal_m3, ...
    'VariableNames', {'CaseID','Time_min','Inflow_m3s','Outlet_Fine_m3s', ...
    'Outlet_Coarse_m3s','Outlet_Neal_m3s','Stage_Fine_m','Stage_Coarse_m', ...
    'Stage_Neal_m','Storage_Fine_m3','Storage_Coarse_m3','Storage_Neal_m3'});
end

function q = outlet_discharge(state, grid)
q = sum(state.outlet_flow(:), 'omitnan') * grid.area / 1000 / 3600;
end

function storage_m3 = compute_storage(state, grid)
d = max(state.d_tot ./ 1000, 0);
if grid.mode == "neal"
    V = hp2d_neal_cell_volume(d, grid.River_Width, grid.River_Depth, grid.dx);
else
    V = d .* grid.area;
end
storage_m3 = sum(V(:), 'omitnan');
end

function wse = model_wse(state, grid)
d = max(state.d_tot ./ 1000, 0);
if grid.mode == "neal"
    wse = grid.z - grid.River_Depth + d;
else
    wse = grid.z + d;
end
end

function stage = gauge_stage(state, grid, x_gauge_m)
[~, col] = min(abs(grid.x - x_gauge_m));
wse = model_wse(state, grid);
stage = wse(1, col);
end

function profile = x_profile(state, grid)
wse = model_wse(state, grid);
profile = table(grid.x, wse(:), max(state.d_tot(:) ./ 1000, 0), ...
    'VariableNames', {'X_m','WaterSurface_m','Depth_m'});
end

function [metric_case, mass_case, pass_case, profile_case] = summarize_case(Result, cfg)
T = Result.TimeSeries;
fine = compare_to_reference(T.Time_min, T.Outlet_Fine_m3s, T.Outlet_Fine_m3s, ...
    T.Stage_Fine_m, T.Stage_Fine_m);
coarse = compare_to_reference(T.Time_min, T.Outlet_Fine_m3s, T.Outlet_Coarse_m3s, ...
    T.Stage_Fine_m, T.Stage_Coarse_m);
neal = compare_to_reference(T.Time_min, T.Outlet_Fine_m3s, T.Outlet_Neal_m3s, ...
    T.Stage_Fine_m, T.Stage_Neal_m);

if cfg.comparison_state == "peak"
    fine_state = Result.peak.fine_state;
    coarse_state = Result.peak.coarse_state;
    neal_state = Result.peak.neal_state;
    profile_time = Result.peak.time_min;
else
    fine_state = Result.final.fine;
    coarse_state = Result.final.coarse;
    neal_state = Result.final.neal;
    profile_time = T.Time_min(end);
end

fine_profile = x_profile(fine_state, Result.Grid.fine);
coarse_profile = x_profile(coarse_state, Result.Grid.coarse);
neal_profile = x_profile(neal_state, Result.Grid.neal);

fine_cols = arrayfun(@(x) nearest_idx(Result.Grid.fine.x, x), Result.Grid.coarse.x);
coarse.ProfileRMSE_m = rmse(coarse_profile.WaterSurface_m, fine_profile.WaterSurface_m(fine_cols));
neal.ProfileRMSE_m = rmse(neal_profile.WaterSurface_m, fine_profile.WaterSurface_m(fine_cols));
fine.ProfileRMSE_m = 0.0;

fine.MassResidual_pct = mass_residual_pct(Result.cumulative.inflow_m3, ...
    Result.cumulative.fine_outlet_m3, T.Storage_Fine_m3);
coarse.MassResidual_pct = mass_residual_pct(Result.cumulative.inflow_m3, ...
    Result.cumulative.coarse_outlet_m3, T.Storage_Coarse_m3);
neal.MassResidual_pct = mass_residual_pct(Result.cumulative.inflow_m3, ...
    Result.cumulative.neal_outlet_m3, T.Storage_Neal_m3);

metric_case = table( ...
    repmat(string(cfg.id), 3, 1), ...
    ["FineChannel"; "CoarseOrdinary"; "CoarseNeal"], ...
    [fine.HydrographRMSE_m3s; coarse.HydrographRMSE_m3s; neal.HydrographRMSE_m3s], ...
    [fine.HydrographNSE; coarse.HydrographNSE; neal.HydrographNSE], ...
    [fine.PeakTimingError_min; coarse.PeakTimingError_min; neal.PeakTimingError_min], ...
    [fine.PeakError_pct; coarse.PeakError_pct; neal.PeakError_pct], ...
    [fine.StageRMSE_m; coarse.StageRMSE_m; neal.StageRMSE_m], ...
    [fine.ProfileRMSE_m; coarse.ProfileRMSE_m; neal.ProfileRMSE_m], ...
    [fine.MassResidual_pct; coarse.MassResidual_pct; neal.MassResidual_pct], ...
    'VariableNames', {'CaseID','Model','HydrographRMSE_m3s','HydrographNSE', ...
    'PeakTimingError_min','PeakError_pct','StageRMSE_m','ProfileRMSE_m','MassResidual_pct'});

initial_storage = [T.Storage_Fine_m3(1); T.Storage_Coarse_m3(1); T.Storage_Neal_m3(1)];
final_storage = [T.Storage_Fine_m3(end); T.Storage_Coarse_m3(end); T.Storage_Neal_m3(end)];
outlet_volume = [Result.cumulative.fine_outlet_m3; Result.cumulative.coarse_outlet_m3; ...
    Result.cumulative.neal_outlet_m3];
storage_change = final_storage - initial_storage;
residual_m3 = Result.cumulative.inflow_m3 - outlet_volume - storage_change;
available_volume = initial_storage + Result.cumulative.inflow_m3;
accounted_volume = outlet_volume + final_storage;

mass_case = table( ...
    repmat(string(cfg.id), 3, 1), ...
    ["FineChannel"; "CoarseOrdinary"; "CoarseNeal"], ...
    [Result.cumulative.inflow_m3; Result.cumulative.inflow_m3; Result.cumulative.inflow_m3], ...
    initial_storage, available_volume, outlet_volume, final_storage, accounted_volume, ...
    storage_change, residual_m3, ...
    [fine.MassResidual_pct; coarse.MassResidual_pct; neal.MassResidual_pct], ...
    'VariableNames', {'CaseID','Model','EventInputVolume_m3','InitialStorage_m3', ...
    'AvailableVolume_m3','OutletVolume_m3','FinalStorage_m3','AccountedVolume_m3', ...
    'StorageChange_m3','Residual_m3','MassResidual_pct'});

pass_flag = abs(neal.MassResidual_pct) < 0.1 && ...
    neal.HydrographRMSE_m3s <= coarse.HydrographRMSE_m3s && ...
    neal.HydrographNSE >= coarse.HydrographNSE && ...
    neal.StageRMSE_m <= coarse.StageRMSE_m;

if cfg.id ~= "P1-SUBGRID-NEAL-SIMPLE-003"
    pass_flag = pass_flag && neal.HydrographNSE > 0.95 && neal.StageRMSE_m < 0.10;
end

pass_case = table(string(cfg.id), logical(pass_flag), ...
    string(sprintf(['Neal vs coarse: RMSE %.3f vs %.3f m3/s; NSE %.3f vs %.3f; ', ...
    'stage RMSE %.3f vs %.3f m; profile RMSE %.3f vs %.3f m'], ...
    neal.HydrographRMSE_m3s, coarse.HydrographRMSE_m3s, neal.HydrographNSE, ...
    coarse.HydrographNSE, neal.StageRMSE_m, coarse.StageRMSE_m, ...
    neal.ProfileRMSE_m, coarse.ProfileRMSE_m)), ...
    'VariableNames', {'CaseID','Pass','Notes'});

profile_case = table( ...
    repmat(string(cfg.id), numel(Result.Grid.coarse.x) * 3, 1), ...
    repmat(profile_time, numel(Result.Grid.coarse.x) * 3, 1), ...
    [repmat("FineChannel", numel(Result.Grid.coarse.x), 1); ...
     repmat("CoarseOrdinary", numel(Result.Grid.coarse.x), 1); ...
     repmat("CoarseNeal", numel(Result.Grid.coarse.x), 1)], ...
    repmat(Result.Grid.coarse.x, 3, 1), ...
    [fine_profile.WaterSurface_m(fine_cols); coarse_profile.WaterSurface_m; neal_profile.WaterSurface_m], ...
    'VariableNames', {'CaseID','Time_min','Model','X_m','WaterSurface_m'});
end

function M = compare_to_reference(time_min, q_ref, q_sim, stage_ref, stage_sim)
M = struct();
M.HydrographRMSE_m3s = rmse(q_ref, q_sim);
M.HydrographNSE = nse(q_ref, q_sim);
M.PeakTimingError_min = peak_timing_error(time_min, q_ref, q_sim);
M.PeakError_pct = rel_pct(max(q_sim) - max(q_ref), max(q_ref));
M.StageRMSE_m = rmse(stage_ref, stage_sim);
end

function residual_pct = mass_residual_pct(input_vol, out_vol, storage_m3)
delta_s = storage_m3(end) - storage_m3(1);
residual_pct = rel_pct(input_vol - out_vol - delta_s, max(input_vol, 1.0));
end

function err = peak_timing_error(time_min, ref, sim)
[~, ir] = max(ref);
[~, is] = max(sim);
err = time_min(is) - time_min(ir);
end

function idx = nearest_idx(xvec, x)
[~, idx] = min(abs(xvec - x));
end

function v = rmse(a, b)
a = double(a(:));
b = double(b(:));
mask = isfinite(a) & isfinite(b);
if ~any(mask)
    v = NaN;
    return;
end
v = sqrt(mean((a(mask) - b(mask)).^2));
end

function v = nse(obs, sim)
obs = double(obs(:));
sim = double(sim(:));
mask = isfinite(obs) & isfinite(sim);
if nnz(mask) < 2
    v = NaN;
    return;
end
den = sum((obs(mask) - mean(obs(mask))).^2);
if den <= eps
    v = NaN;
    return;
end
v = 1 - sum((sim(mask) - obs(mask)).^2) / den;
end

function pct = rel_pct(delta, ref)
pct = 100 * delta / max(abs(ref), eps);
end

function plot_hydrograph_and_stage(Result, cfg, fig_dir)
T = Result.TimeSeries;
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 750]);
tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(T.Time_min, T.Inflow_m3s, 'k-', 'LineWidth', 1.5); hold on;
plot(T.Time_min, T.Outlet_Fine_m3s, 'b-', 'LineWidth', 1.8);
plot(T.Time_min, T.Outlet_Coarse_m3s, 'Color', [0.82 0.37 0.10], 'LineWidth', 1.5);
plot(T.Time_min, T.Outlet_Neal_m3s, 'Color', [0.06 0.47 0.18], 'LineWidth', 1.8);
grid on;
xlabel('Time (min)');
ylabel('Discharge (m^3 s^{-1})');
title(sprintf('%s outlet hydrograph', cfg.name), 'Interpreter', 'none');
legend({'Inflow','Fine channel','Coarse ordinary','Coarse Neal'}, 'Location', 'best');

nexttile;
plot(T.Time_min, T.Stage_Fine_m, 'b-', 'LineWidth', 1.8); hold on;
plot(T.Time_min, T.Stage_Coarse_m, 'Color', [0.82 0.37 0.10], 'LineWidth', 1.5);
plot(T.Time_min, T.Stage_Neal_m, 'Color', [0.06 0.47 0.18], 'LineWidth', 1.8);
grid on;
xlabel('Time (min)');
ylabel('Water-surface elevation (m)');
title(sprintf('Gauge stage at x = %.0f m', Result.Grid.gauge_x_m));
legend({'Fine channel','Coarse ordinary','Coarse Neal'}, 'Location', 'best');

exportgraphics(f, fullfile(fig_dir, sprintf('%s_HydrographStage.png', cfg.id)), 'Resolution', 200);
close(f);
end

function plot_profiles(profile_case, cfg, fig_dir)
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 950 500]);
fine = profile_case(strcmp(profile_case.Model, "FineChannel"), :);
coarse = profile_case(strcmp(profile_case.Model, "CoarseOrdinary"), :);
neal = profile_case(strcmp(profile_case.Model, "CoarseNeal"), :);

plot(fine.X_m, fine.WaterSurface_m, 'b-', 'LineWidth', 1.8); hold on;
plot(coarse.X_m, coarse.WaterSurface_m, 'Color', [0.82 0.37 0.10], 'LineWidth', 1.5);
plot(neal.X_m, neal.WaterSurface_m, 'Color', [0.06 0.47 0.18], 'LineWidth', 1.8);
grid on;
xlabel('Distance downstream (m)');
ylabel('Water-surface elevation (m)');
title(sprintf('%s profile comparison', cfg.name), 'Interpreter', 'none');
legend({'Fine channel','Coarse ordinary','Coarse Neal'}, 'Location', 'best');

exportgraphics(f, fullfile(fig_dir, sprintf('%s_Profile.png', cfg.id)), 'Resolution', 200);
close(f);
end
