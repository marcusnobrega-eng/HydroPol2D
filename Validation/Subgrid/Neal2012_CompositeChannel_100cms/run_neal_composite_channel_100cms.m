function run_neal_composite_channel_100cms(forcing_type)
%RUN_NEAL_COMPOSITE_CHANNEL_100CMS Compare aligned fine, coarse, and Neal grids.

if nargin < 1
    forcing_type = "constant";
end
forcing_type = lower(string(forcing_type));

case_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(fileparts(fileparts(case_dir)));
addpath(fullfile(model_dir, 'HydroPol2D_Functions'));

Cfg = struct();
Cfg.nominal_length_m = 1000;
Cfg.length_m = 990;             % exact 10 m to 30 m alignment
Cfg.channel_width_m = 10;
Cfg.channel_depth_m = 1;
Cfg.overbank_width_each_m = 90;
Cfg.outer_berm_width_each_m = 10;
Cfg.outer_berm_height_m = 2;
Cfg.slope_mpm = 0.005;
Cfg.manning_n = 0.035;
Cfg.dt_min = 0.01;
Cfg.record_dt_min = 1;
Cfg.fine_dx_m = 10;
Cfg.coarse_dx_m = 30;
Cfg.gauge_x_m = 750;
Cfg.midpoint_x_m = Cfg.length_m / 2;

switch forcing_type
    case "constant"
        Cfg.case_id = "P1-SUBGRID-NEAL-COMPOSITE-100CMS-001";
        Cfg.forcing_name = "Constant 100 m3/s";
        Cfg.peak_inflow_m3s = 100;
        Cfg.duration_min = 120;
        Cfg.qfun = @(t) Cfg.peak_inflow_m3s + 0 .* t;
        out_dir = fullfile(case_dir, 'Outputs', 'Validation');
        fig_dir = fullfile(case_dir, 'Figures');
    case "nash"
        Cfg.case_id = "P1-SUBGRID-NEAL-COMPOSITE-NASH100-001";
        Cfg.forcing_name = "Nash hydrograph, peak 100 m3/s at 30 min";
        Cfg.peak_inflow_m3s = 100;
        Cfg.nash_shape = 4;
        Cfg.nash_time_to_peak_min = 30;
        Cfg.duration_min = 180;
        Cfg.qfun = @(t) nash_hydrograph(t, Cfg.peak_inflow_m3s, ...
            Cfg.nash_time_to_peak_min, Cfg.nash_shape);
        out_dir = fullfile(case_dir, 'Outputs', 'Validation', 'NashTransient');
        fig_dir = fullfile(case_dir, 'Figures', 'NashTransient');
    otherwise
        error('Unknown forcing type: %s. Use constant or nash.', forcing_type);
end
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

G = build_grids(Cfg);
[T, Mass, Final] = run_models(Cfg, G);
Metrics = score_models(Cfg, T, Mass);
Pass = make_pass_table(Cfg, Metrics);

writetable(T, fullfile(out_dir, 'Hydrograph_Comparison.csv'));
writetable(Mass, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(Metrics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(Pass, fullfile(out_dir, 'Pass_Fail.csv'));

plot_hydrographs(Cfg, T, fig_dir);
plot_storage_stage(Cfg, T, fig_dir);
plot_midpoint_depth(Cfg, T, fig_dir);
if forcing_type == "constant"
    plot_cross_section(Cfg, G, Final, fig_dir);
end

fprintf('Composite-channel %s test complete. Neal pass = %d.\n', forcing_type, Pass.Pass);
end

function G = build_grids(Cfg)
ratio = Cfg.coarse_dx_m / Cfg.fine_dx_m;
if ratio ~= round(ratio)
    error('Fine-to-coarse resolution ratio must be an integer.');
end

dx = Cfg.fine_dx_m;
nx = round(Cfg.length_m / dx);
total_width = Cfg.channel_width_m + 2 * ...
    (Cfg.overbank_width_each_m + Cfg.outer_berm_width_each_m);
ny = round(total_width / dx);
x = ((1:nx) - 0.5) .* dx;
y = ((1:ny) - (ny + 1) / 2) .* dx;
[X,Y] = meshgrid(x, y);

base = -Cfg.slope_mpm .* X;
z_reference = base + Cfg.channel_depth_m;
outer_berm = abs(Y) > ...
    (Cfg.channel_width_m / 2 + Cfg.overbank_width_each_m);
z_reference(outer_berm) = base(outer_berm) + Cfg.outer_berm_height_m;

z_fine = z_reference;
channel = abs(Y) < Cfg.channel_width_m / 2 + 10 * eps;
z_fine(channel) = base(channel);

roughness = Cfg.manning_n .* ones(size(z_fine));
G.fine = make_grid("fine", dx, x, y, z_fine, z_reference, roughness, Cfg);
G.fine.channel_rows = find(any(channel, 2));
G.fine.inlet_rows = G.fine.channel_rows;

z_coarse = block_reduce(z_fine, ratio, @mean);
z_reference_coarse = block_reduce(z_reference, ratio, @mean);
roughness_coarse = block_reduce(roughness, ratio, @mean);
dx = Cfg.coarse_dx_m;
nx = size(z_coarse, 2);
ny = size(z_coarse, 1);
x = ((1:nx) - 0.5) .* dx;
y = ((1:ny) - (ny + 1) / 2) .* dx;

G.coarse = make_grid("coarse", dx, x, y, z_coarse, ...
    z_reference_coarse, roughness_coarse, Cfg);
[~, G.coarse.channel_row] = min(z_coarse(:, 1));
G.coarse.channel_rows = G.coarse.channel_row;
G.coarse.inlet_rows = G.coarse.channel_row;

z_neal = z_reference_coarse;
roughness_neal = Cfg.manning_n .* ones(size(z_neal));
G.neal = make_grid("neal", dx, x, y, z_neal, ...
    z_reference_coarse, roughness_neal, Cfg);
[~, G.neal.channel_row] = min(abs(y));
G.neal.channel_rows = G.neal.channel_row;
G.neal.inlet_rows = G.neal.channel_row;
G.neal.River_Width(G.neal.channel_row, :) = Cfg.channel_width_m;
G.neal.River_Depth(G.neal.channel_row, :) = Cfg.channel_depth_m;
end

function grid = make_grid(mode, dx, x, y, z, z_reference, roughness, Cfg)
grid = struct();
grid.mode = mode;
grid.dx = dx;
grid.area = dx^2;
grid.nx = size(z, 2);
grid.ny = size(z, 1);
grid.x = x(:);
grid.y = y(:);
grid.z = z;
grid.z_reference = z_reference;
grid.physical_channel_depth_m = Cfg.channel_depth_m;
grid.roughness = roughness;
grid.n_channel = Cfg.manning_n .* ones(size(z));
grid.n_flood = Cfg.manning_n .* ones(size(z));
grid.River_Width = zeros(size(z));
grid.River_Depth = zeros(size(z));
grid.outlet_index = false(size(z));
grid.outlet_type = 1;
grid.slope_outlet = Cfg.slope_mpm;
grid.row_outlet = (1:grid.ny)';
grid.col_outlet = repmat(grid.nx, grid.ny, 1);
grid.channel_row = ceil(grid.ny / 2);
grid.channel_rows = grid.channel_row;
grid.inlet_rows = grid.channel_row;
end

function [T, Mass, Final] = run_models(Cfg, G)
dt_s = Cfg.dt_min * 60;
nsteps = round(Cfg.duration_min / Cfg.dt_min);
record_every = round(Cfg.record_dt_min / Cfg.dt_min);
nrecords = floor(nsteps / record_every) + 1;

S.fine = initialize_state(G.fine);
S.coarse = initialize_state(G.coarse);
S.neal = initialize_state(G.neal);

initial_storage = [compute_storage(S.fine, G.fine); ...
    compute_storage(S.coarse, G.coarse); compute_storage(S.neal, G.neal)];
cum_out = zeros(3, 1);
cum_input = 0;

R.Time_min = zeros(nrecords, 1);
R.Inflow_m3s = zeros(nrecords, 1);
R.Outlet_Fine_m3s = zeros(nrecords, 1);
R.Outlet_Coarse_m3s = zeros(nrecords, 1);
R.Outlet_Neal_m3s = zeros(nrecords, 1);
R.Gauge_Fine_m3s = zeros(nrecords, 1);
R.Gauge_Coarse_m3s = zeros(nrecords, 1);
R.Gauge_Neal_m3s = zeros(nrecords, 1);
R.Storage_Fine_m3 = zeros(nrecords, 1);
R.Storage_Coarse_m3 = zeros(nrecords, 1);
R.Storage_Neal_m3 = zeros(nrecords, 1);
R.Depth_Fine_m = zeros(nrecords, 1);
R.Depth_Coarse_m = zeros(nrecords, 1);
R.Depth_Neal_m = zeros(nrecords, 1);
R.MidDepth_Fine_m = zeros(nrecords, 1);
R.MidDepth_Coarse_m = zeros(nrecords, 1);
R.MidDepth_Neal_m = zeros(nrecords, 1);

R = store_record(R, 1, 0, S, G, Cfg);
rec = 1;

for it = 1:nsteps
    qin = Cfg.qfun((it - 1) * Cfg.dt_min);
    S.fine = apply_inflow(S.fine, G.fine, qin, dt_s);
    S.coarse = apply_inflow(S.coarse, G.coarse, qin, dt_s);
    S.neal = apply_inflow(S.neal, G.neal, qin, dt_s);

    S.fine = step_model(S.fine, G.fine, Cfg.dt_min);
    S.coarse = step_model(S.coarse, G.coarse, Cfg.dt_min);
    S.neal = step_model(S.neal, G.neal, Cfg.dt_min);

    qout = [outlet_discharge(S.fine, G.fine); ...
        outlet_discharge(S.coarse, G.coarse); ...
        outlet_discharge(S.neal, G.neal)];
    cum_out = cum_out + qout .* dt_s;
    cum_input = cum_input + qin .* dt_s;

    if mod(it, record_every) == 0 || it == nsteps
        rec = rec + 1;
        R = store_record(R, rec, it * Cfg.dt_min, S, G, Cfg);
    end
end

T = struct2table(R);
final_storage = [T.Storage_Fine_m3(end); T.Storage_Coarse_m3(end); ...
    T.Storage_Neal_m3(end)];
input_volume = cum_input;
available_volume = initial_storage + input_volume;
accounted_volume = cum_out + final_storage;
residual = available_volume - accounted_volume;

Mass = table(repmat(Cfg.case_id, 3, 1), ...
    ["FineExplicit10m"; "CoarseOrdinary30m"; "NealSubgrid30m"], ...
    repmat(input_volume, 3, 1), initial_storage, available_volume, ...
    cum_out, final_storage, accounted_volume, residual, ...
    100 .* residual ./ max(input_volume, 1), ...
    'VariableNames', {'CaseID','Model','InputVolume_m3','InitialStorage_m3', ...
    'AvailableVolume_m3','OutletVolume_m3','FinalStorage_m3', ...
    'AccountedVolume_m3','Residual_m3','MassResidual_pct'});

Final = S;
end

function state = initialize_state(grid)
state = struct();
state.d_tot = zeros(grid.ny, grid.nx);
state.d_p = zeros(grid.ny, grid.nx);
state.outflow = zeros(grid.ny, grid.nx, 5);
state.outlet_flow = zeros(grid.ny, grid.nx);
state.Qc = zeros(grid.ny, grid.nx, 2);
state.Qf = zeros(grid.ny, grid.nx, 2);
state.Qci = zeros(grid.ny, grid.nx, 2);
state.Qfi = zeros(grid.ny, grid.nx, 2);
if grid.mode == "neal"
    state.C_a = hp2d_neal_cell_area(zeros(grid.ny, grid.nx), ...
        grid.River_Width, grid.River_Depth, grid.dx);
else
    state.C_a = grid.area .* ones(grid.ny, grid.nx);
end
end

function state = apply_inflow(state, grid, q_m3s, dt_s)
if q_m3s <= 0
    return;
end
rows = grid.inlet_rows(:);
vol_per_cell = q_m3s * dt_s / numel(rows);
if grid.mode == "neal"
    h = max(state.d_tot(rows, 1) ./ 1000, 0);
    h = hp2d_neal_apply_volume_change(h, vol_per_cell .* ones(size(h)), ...
        grid.River_Width(rows, 1), grid.River_Depth(rows, 1), grid.dx);
    state.d_tot(rows, 1) = 1000 .* h;
else
    state.d_tot(rows, 1) = state.d_tot(rows, 1) + ...
        1000 .* vol_per_cell ./ grid.area;
end
end

function state = step_model(state, grid, dt_min)
idx_nan = false(grid.ny, grid.nx);
if grid.mode == "neal"
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,Qc,Qf,Qci,Qfi,C_a] = ...
        Local_Inertial_Model_D4(2, [], [], [], [], [], [], [], [], [], [], [], [], ...
        0, grid.z, state.d_tot, state.d_p, grid.roughness, grid.roughness.^2, ...
        grid.area, dt_min, grid.dx, grid.outlet_index, grid.outlet_type, ...
        grid.slope_outlet, grid.row_outlet, grid.col_outlet, 1e-6, ...
        state.outflow, idx_nan, 0, 1, grid.n_channel, grid.n_flood, ...
        grid.River_Width, grid.River_Depth, state.Qc, state.Qf, state.Qci, ...
        state.Qfi, state.C_a, [], 1, 0, []);
    state.Qc = Qc;
    state.Qf = Qf;
    state.Qci = Qci;
    state.Qfi = Qfi;
else
    [~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = ...
        Local_Inertial_Model_D4(2, [], [], [], [], [], [], [], [], [], [], [], [], ...
        0, grid.z, state.d_tot, state.d_p, grid.roughness, grid.roughness.^2, ...
        grid.area, dt_min, grid.dx, grid.outlet_index, grid.outlet_type, ...
        grid.slope_outlet, grid.row_outlet, grid.col_outlet, 1e-6, ...
        state.outflow, idx_nan, 0, 0, grid.n_channel, grid.n_flood, ...
        grid.River_Width, grid.River_Depth, [], [], [], [], state.C_a, [], 0, 0, []);
end
state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function R = store_record(R, i, time_min, S, G, Cfg)
R.Time_min(i) = time_min;
R.Inflow_m3s(i) = Cfg.qfun(time_min);
R.Outlet_Fine_m3s(i) = outlet_discharge(S.fine, G.fine);
R.Outlet_Coarse_m3s(i) = outlet_discharge(S.coarse, G.coarse);
R.Outlet_Neal_m3s(i) = outlet_discharge(S.neal, G.neal);
R.Gauge_Fine_m3s(i) = section_discharge(S.fine, G.fine, Cfg.gauge_x_m);
R.Gauge_Coarse_m3s(i) = section_discharge(S.coarse, G.coarse, Cfg.gauge_x_m);
R.Gauge_Neal_m3s(i) = section_discharge(S.neal, G.neal, Cfg.gauge_x_m);
R.Storage_Fine_m3(i) = compute_storage(S.fine, G.fine);
R.Storage_Coarse_m3(i) = compute_storage(S.coarse, G.coarse);
R.Storage_Neal_m3(i) = compute_storage(S.neal, G.neal);
R.Depth_Fine_m(i) = channel_depth(S.fine, G.fine, Cfg.gauge_x_m);
R.Depth_Coarse_m(i) = channel_depth(S.coarse, G.coarse, Cfg.gauge_x_m);
R.Depth_Neal_m(i) = channel_depth(S.neal, G.neal, Cfg.gauge_x_m);
R.MidDepth_Fine_m(i) = channel_depth(S.fine, G.fine, Cfg.midpoint_x_m);
R.MidDepth_Coarse_m(i) = channel_depth(S.coarse, G.coarse, Cfg.midpoint_x_m);
R.MidDepth_Neal_m(i) = channel_depth(S.neal, G.neal, Cfg.midpoint_x_m);
end

function q = outlet_discharge(state, grid)
q = sum(state.outlet_flow(:), 'omitnan') * grid.area / 1000 / 3600;
end

function q = section_discharge(state, grid, x_gauge)
[~, col] = min(abs(grid.x - x_gauge));
col = min(col, grid.nx - 1);
q = sum(state.outflow(:, col, 1), 'omitnan') * grid.area / 1000 / 3600;
end

function storage = compute_storage(state, grid)
h = max(state.d_tot ./ 1000, 0);
if grid.mode == "neal"
    V = hp2d_neal_cell_volume(h, grid.River_Width, grid.River_Depth, grid.dx);
else
    V = h .* grid.area;
end
storage = sum(V(:), 'omitnan');
end

function depth = channel_depth(state, grid, x_gauge)
[~, col] = min(abs(grid.x - x_gauge));
rows = grid.channel_rows(:);
if max(state.d_tot(rows, col), [], 'omitnan') <= 1e-6
    depth = 0;
    return;
end
wse = model_wse(state, grid);
physical_bed = grid.z_reference(grid.channel_row, col) - ...
    grid.physical_channel_depth_m;
depth = max(wse(rows, col), [], 'omitnan') - physical_bed;
end

function wse = model_wse(state, grid)
h = max(state.d_tot ./ 1000, 0);
wse = model_bed(grid) + h;
end

function bed = model_bed(grid)
if grid.mode == "neal"
    bed = grid.z - grid.River_Depth;
else
    bed = grid.z;
end
end

function Metrics = score_models(Cfg, T, Mass)
q_ref = T.Outlet_Fine_m3s;
g_ref = T.Gauge_Fine_m3s;
s_ref = T.Storage_Fine_m3;
models = ["FineExplicit10m"; "CoarseOrdinary30m"; "NealSubgrid30m"];
q = {T.Outlet_Fine_m3s, T.Outlet_Coarse_m3s, T.Outlet_Neal_m3s};
g = {T.Gauge_Fine_m3s, T.Gauge_Coarse_m3s, T.Gauge_Neal_m3s};
s = {T.Storage_Fine_m3, T.Storage_Coarse_m3, T.Storage_Neal_m3};

out_rmse = zeros(3,1);
out_nse = ones(3,1);
gauge_rmse = zeros(3,1);
arrival_error = zeros(3,1);
final_q_error = zeros(3,1);
peak_q_error = zeros(3,1);
peak_timing_error = zeros(3,1);
storage_error = zeros(3,1);
storage_rmse = zeros(3,1);
arrival_threshold = 0.5 * max(q_ref);
[~, ref_peak_idx] = max(q_ref);
for i = 1:3
    out_rmse(i) = rmse(q_ref, q{i});
    out_nse(i) = nse(q_ref, q{i});
    gauge_rmse(i) = rmse(g_ref, g{i});
    arrival_error(i) = arrival_time(T.Time_min, q{i}, arrival_threshold) - ...
        arrival_time(T.Time_min, q_ref, arrival_threshold);
    final_q_error(i) = 100 * (q{i}(end) - q_ref(end)) / max(abs(q_ref(end)), 1);
    [peak_q, peak_idx] = max(q{i});
    peak_q_error(i) = 100 * (peak_q - max(q_ref)) / max(max(q_ref), 1);
    peak_timing_error(i) = T.Time_min(peak_idx) - T.Time_min(ref_peak_idx);
    storage_error(i) = 100 * (s{i}(end) - s_ref(end)) / max(abs(s_ref(end)), 1);
    storage_rmse(i) = rmse(s_ref, s{i});
end

Metrics = table(repmat(Cfg.case_id, 3, 1), models, out_rmse, out_nse, ...
    gauge_rmse, arrival_error, peak_timing_error, peak_q_error, final_q_error, ...
    storage_rmse, storage_error, ...
    Mass.MassResidual_pct, ...
    'VariableNames', {'CaseID','Model','OutletRMSE_m3s','OutletNSE', ...
    'InternalGaugeRMSE_m3s','ArrivalTimeError_min','PeakTimingError_min', ...
    'PeakDischargeError_pct','FinalDischargeError_pct','StorageRMSE_m3', ...
    'FinalStorageError_pct','MassResidual_pct'});
end

function Pass = make_pass_table(Cfg, Metrics)
coarse = Metrics(Metrics.Model == "CoarseOrdinary30m", :);
neal = Metrics(Metrics.Model == "NealSubgrid30m", :);
passed = abs(neal.MassResidual_pct) < 0.1 && ...
    neal.OutletRMSE_m3s <= coarse.OutletRMSE_m3s && ...
    neal.InternalGaugeRMSE_m3s <= coarse.InternalGaugeRMSE_m3s && ...
    neal.StorageRMSE_m3 <= coarse.StorageRMSE_m3;
notes = sprintf(['Neal versus ordinary 30 m: outlet RMSE %.3f versus %.3f m3/s; ', ...
    'internal RMSE %.3f versus %.3f m3/s; storage RMSE %.1f versus %.1f m3.'], ...
    neal.OutletRMSE_m3s, coarse.OutletRMSE_m3s, ...
    neal.InternalGaugeRMSE_m3s, coarse.InternalGaugeRMSE_m3s, ...
    neal.StorageRMSE_m3, coarse.StorageRMSE_m3);
Pass = table(Cfg.case_id, logical(passed), string(notes), ...
    'VariableNames', {'CaseID','Pass','Notes'});
end

function plot_hydrographs(Cfg, T, fig_dir)
f = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(T.Time_min, T.Inflow_m3s, 'k--', 'LineWidth', 1.5); hold on;
plot(T.Time_min, T.Outlet_Fine_m3s, 'Color', [0.00 0.35 0.70], 'LineWidth', 2.0);
plot(T.Time_min, T.Outlet_Coarse_m3s, 'Color', [0.85 0.40 0.05], 'LineWidth', 1.8);
plot(T.Time_min, T.Outlet_Neal_m3s, 'Color', [0.10 0.50 0.25], 'LineWidth', 2.0);
grid on; xlim([0 Cfg.duration_min]);
xlabel('Time (min)'); ylabel('Discharge (m^3 s^{-1})');
title('Downstream hydrograph');
legend({'Inflow','Fine explicit, 10 m','Ordinary, 30 m','Neal, 30 m'}, ...
    'Location','southeast');

nexttile;
plot(T.Time_min, T.Gauge_Fine_m3s, 'Color', [0.00 0.35 0.70], 'LineWidth', 2.0); hold on;
plot(T.Time_min, T.Gauge_Coarse_m3s, 'Color', [0.85 0.40 0.05], 'LineWidth', 1.8);
plot(T.Time_min, T.Gauge_Neal_m3s, 'Color', [0.10 0.50 0.25], 'LineWidth', 2.0);
grid on; xlim([0 Cfg.duration_min]);
xlabel('Time (min)'); ylabel('Discharge (m^3 s^{-1})');
title(sprintf('Internal section at x = %.0f m', Cfg.gauge_x_m));
legend({'Fine explicit, 10 m','Ordinary, 30 m','Neal, 30 m'}, ...
    'Location','southeast');
sgtitle(Cfg.forcing_name);
exportgraphics(f, fullfile(fig_dir, 'CompositeChannel_Hydrographs.png'), 'Resolution', 220);
close(f);
end

function q = nash_hydrograph(time_min, peak_m3s, time_to_peak_min, shape)
q = zeros(size(time_min));
positive = time_min > 0;
power = shape - 1;
scaled_time = time_min(positive) ./ time_to_peak_min;
q(positive) = peak_m3s .* scaled_time .^ power .* ...
    exp(power .* (1 - scaled_time));
end

function plot_storage_stage(Cfg, T, fig_dir)
f = figure('Visible','off','Color','w','Position',[100 100 1050 720]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(T.Time_min, T.Storage_Fine_m3, 'Color', [0.00 0.35 0.70], 'LineWidth', 2.0); hold on;
plot(T.Time_min, T.Storage_Coarse_m3, 'Color', [0.85 0.40 0.05], 'LineWidth', 1.8);
plot(T.Time_min, T.Storage_Neal_m3, 'Color', [0.10 0.50 0.25], 'LineWidth', 2.0);
grid on; xlim([0 Cfg.duration_min]);
xlabel('Time (min)'); ylabel('Stored water (m^3)'); title('Domain storage');
legend({'Fine explicit, 10 m','Ordinary, 30 m','Neal, 30 m'}, 'Location','best');

nexttile;
plot(T.Time_min, T.Depth_Fine_m, 'Color', [0.00 0.35 0.70], 'LineWidth', 2.0); hold on;
plot(T.Time_min, T.Depth_Coarse_m, 'Color', [0.85 0.40 0.05], 'LineWidth', 1.8);
plot(T.Time_min, T.Depth_Neal_m, 'Color', [0.10 0.50 0.25], 'LineWidth', 2.0);
yline(Cfg.channel_depth_m, 'k:', 'Bank elevation', 'LineWidth', 1.2);
grid on; xlim([0 Cfg.duration_min]);
xlabel('Time (min)'); ylabel('Depth above channel bed (m)');
title(sprintf('Channel depth at x = %.0f m', Cfg.gauge_x_m));
legend({'Fine explicit, 10 m','Ordinary, 30 m','Neal, 30 m','Bank elevation'}, ...
    'Location','best');
exportgraphics(f, fullfile(fig_dir, 'CompositeChannel_StorageStage.png'), 'Resolution', 220);
close(f);
end

function plot_midpoint_depth(Cfg, T, fig_dir)
f = figure('Visible','off','Color','w','Position',[100 100 1000 520]);
plot(T.Time_min, T.MidDepth_Fine_m, 'Color', [0.00 0.35 0.70], ...
    'LineWidth', 2.0); hold on;
plot(T.Time_min, T.MidDepth_Coarse_m, 'Color', [0.85 0.40 0.05], ...
    'LineWidth', 1.8);
plot(T.Time_min, T.MidDepth_Neal_m, 'Color', [0.10 0.50 0.25], ...
    'LineWidth', 2.0);
yline(Cfg.channel_depth_m, 'k:', 'Bank elevation', 'LineWidth', 1.2);
grid on;
xlim([0 Cfg.duration_min]);
xlabel('Time (min)');
ylabel('Depth above channel bed (m)');
title(sprintf('Center-cell water depth at x = %.0f m', Cfg.midpoint_x_m));
legend({'Fine explicit, 10 m','Ordinary, 30 m','Neal, 30 m', ...
    'Bank elevation'}, 'Location','southeast');
exportgraphics(f, fullfile(fig_dir, 'CompositeChannel_MidCellDepth.png'), ...
    'Resolution', 220);
close(f);
end

function plot_cross_section(Cfg, G, Final, fig_dir)
[~, cf] = min(abs(G.fine.x - Cfg.gauge_x_m));
[~, cc] = min(abs(G.coarse.x - Cfg.gauge_x_m));
[~, cn] = min(abs(G.neal.x - Cfg.gauge_x_m));

bed_f = model_bed(G.fine); bed_c = model_bed(G.coarse);
wse_f = model_wse(Final.fine, G.fine);
wse_c = model_wse(Final.coarse, G.coarse);
wse_n = model_wse(Final.neal, G.neal);
datum = min(bed_f(:,cf));

h_f = max(Final.fine.d_tot(:,cf) ./ 1000, 0);
h_c = max(Final.coarse.d_tot(:,cc) ./ 1000, 0);
h_n = max(Final.neal.d_tot(:,cn) ./ 1000, 0);
wse_f(h_f <= 1e-6,cf) = NaN;
wse_c(h_c <= 1e-6,cc) = NaN;
wse_n(h_n <= 1e-6,cn) = NaN;

[yf,zf] = step_profile(G.fine.y, bed_f(:,cf) - datum, G.fine.dx);
[yc,zc] = step_profile(G.coarse.y, bed_c(:,cc) - datum, G.coarse.dx);
[ywf,swf] = step_profile(G.fine.y, wse_f(:,cf) - datum, G.fine.dx);
[ywc,swc] = step_profile(G.coarse.y, wse_c(:,cc) - datum, G.coarse.dx);
[ywn,swn] = step_profile(G.neal.y, wse_n(:,cn) - datum, G.neal.dx);

f = figure('Visible','off','Color','w','Position',[100 100 1050 550]);
plot(yf, zf, 'k-', 'LineWidth', 1.8); hold on;
plot(yc, zc, 'Color', [0.85 0.40 0.05], 'LineWidth', 1.5);
plot(ywf, swf, 'Color', [0.00 0.35 0.70], 'LineWidth', 2.0);
plot(ywc, swc, 'Color', [0.85 0.40 0.05], 'LineWidth', 1.8);
plot(ywn, swn, 'Color', [0.10 0.50 0.25], 'LineWidth', 2.0);
grid on;
xlabel('Cross-channel distance (m)'); ylabel('Elevation above channel bed (m)');
title(sprintf('Final cross section at x = %.0f m', Cfg.gauge_x_m));
legend({'Physical terrain resolved at 10 m','Ordinary 30 m terrain', ...
    'Fine water surface','Ordinary water surface','Neal water surface'}, ...
    'Location','northoutside','NumColumns',3);
exportgraphics(f, fullfile(fig_dir, 'CompositeChannel_FinalCrossSection.png'), 'Resolution', 220);
close(f);
end

function [xstep,ystep] = step_profile(centers, values, dx)
centers = centers(:);
values = values(:);
left = centers - dx / 2;
right = centers + dx / 2;
xstep = reshape([left'; right'], [], 1);
ystep = repelem(values, 2);
end

function out = block_reduce(A, ratio, reducer)
[ny,nx] = size(A);
if mod(ny, ratio) ~= 0 || mod(nx, ratio) ~= 0
    error('Grid dimensions must be divisible by the aggregation ratio.');
end
out = zeros(ny / ratio, nx / ratio);
for i = 1:size(out,1)
    rows = (i - 1) * ratio + (1:ratio);
    for j = 1:size(out,2)
        cols = (j - 1) * ratio + (1:ratio);
        block = A(rows, cols);
        out(i,j) = reducer(block(:));
    end
end
end

function value = rmse(a,b)
a = double(a(:)); b = double(b(:));
value = sqrt(mean((a-b).^2, 'omitnan'));
end

function value = nse(obs,sim)
obs = double(obs(:)); sim = double(sim(:));
den = sum((obs - mean(obs,'omitnan')).^2, 'omitnan');
if den <= eps
    value = double(all(abs(obs-sim) < 1e-12));
else
    value = 1 - sum((obs-sim).^2, 'omitnan') / den;
end
end

function t = arrival_time(time,q,threshold)
idx = find(q >= threshold, 1, 'first');
if isempty(idx); t = NaN; else; t = time(idx); end
end
