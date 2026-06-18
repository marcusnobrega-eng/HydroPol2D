%RUN_COMPOUND_CHANNEL_SUBGRID_BENCHMARK 10 m full-momentum reference vs 30 m subgrid run.
%
% Geometry:
%   - longitudinal channel slope: 0.5 %
%   - inbank bottom width: 10 m
%   - bank height: 1 m
%   - two overbanks: 40 m each, 1 % lateral slope away from channel
%
% This is a numerical-reference benchmark, not an analytical validation
% case. The 10 m full-momentum run is treated as the baseline reference.

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
geom.slope_long = 0.005;
geom.channel_width = 10.0;
geom.bank_height = 1.0;
geom.overbank_width_each = 40.0;
geom.overbank_lateral_slope = 0.01;
geom.total_width = geom.channel_width + 2 * geom.overbank_width_each;
geom.length = 300.0;
geom.n = 0.035;

fine_dx = 10.0;
coarse_dx = 30.0;
ratio = coarse_dx / fine_dx;
dt_min = 0.02;
dt_s = dt_min * 60;
t_end_min = 30.0;
t_min = (0:dt_min:t_end_min)';
nt = numel(t_min);

Q_peak = 35.0;   % [m3/s], intentionally high enough to activate overbanks
Qin = triangular_hydrograph(t_min, 3.0, 12.0, 24.0, Q_peak);

fine = initialize_model_grid(fine_dx, geom);
coarse = initialize_model_grid(coarse_dx, geom);
SubgridTables = build_lookup_tables_from_fine_dem(fine.z, fine_dx, ratio, geom.n, 0.02, 3.0);
hp2d_validate_subgrid_tables(SubgridTables);

fine_state = initialize_state(fine, dt_min, false, []);
coarse_state = initialize_state(coarse, dt_min, true, SubgridTables);

fine_records = initialize_records(nt);
coarse_records = initialize_records(nt);

for it = 1:nt
    qin_step = Qin(it);

    fine_state = apply_upstream_inflow(fine_state, fine, qin_step, dt_s, false, []);
    coarse_state = apply_upstream_inflow(coarse_state, coarse, qin_step, dt_s, true, SubgridTables);

    fine_state = step_full_momentum(fine_state, fine, dt_min, false, []);
    coarse_state = step_full_momentum(coarse_state, coarse, dt_min, true, SubgridTables);

    fine_records = record_step(fine_records, it, fine_state, fine, false, []);
    coarse_records = record_step(coarse_records, it, coarse_state, coarse, true, SubgridTables);
end

fine_records.cumulative_inflow_m3 = cumtrapz(t_min * 60, Qin);
coarse_records.cumulative_inflow_m3 = fine_records.cumulative_inflow_m3;
fine_records.cumulative_outlet_m3 = cumtrapz(t_min * 60, fine_records.outlet_q_m3s);
coarse_records.cumulative_outlet_m3 = cumtrapz(t_min * 60, coarse_records.outlet_q_m3s);

q_rmse = rmse(coarse_records.outlet_q_m3s, fine_records.outlet_q_m3s);
q_mae = mean(abs(coarse_records.outlet_q_m3s - fine_records.outlet_q_m3s));
q_peak_ref = max(fine_records.outlet_q_m3s);
q_peak_sub = max(coarse_records.outlet_q_m3s);
q_peak_error_pct = rel_pct(q_peak_sub - q_peak_ref, max(q_peak_ref, eps));
[~, i_ref_peak] = max(fine_records.outlet_q_m3s);
[~, i_sub_peak] = max(coarse_records.outlet_q_m3s);
peak_timing_error_min = t_min(i_sub_peak) - t_min(i_ref_peak);
q_nse = nse(fine_records.outlet_q_m3s, coarse_records.outlet_q_m3s);
volume_error_pct = rel_pct(coarse_records.cumulative_outlet_m3(end) - fine_records.cumulative_outlet_m3(end), ...
    max(fine_records.cumulative_outlet_m3(end), eps));
storage_error_pct = rel_pct(coarse_records.storage_m3(end) - fine_records.storage_m3(end), ...
    max(fine_records.storage_m3(end), eps));

profile_time_min = 18.0;
[~, ip] = min(abs(t_min - profile_time_min));
fine_profile = cross_section_profile(fine_state_snapshots(fine, fine_records, ip), fine, false, []);
coarse_profile = cross_section_profile(fine_state_snapshots(coarse, coarse_records, ip), coarse, true, SubgridTables);

metrics = table( ...
    ["outlet_q_rmse_m3s"; "outlet_q_mae_m3s"; "outlet_q_nse"; ...
     "peak_q_error_pct"; "peak_timing_error_min"; "cumulative_outlet_volume_error_pct"; ...
     "final_storage_error_pct"], ...
    [q_rmse; q_mae; q_nse; q_peak_error_pct; peak_timing_error_min; volume_error_pct; storage_error_pct], ...
    [2.0; 1.0; 0.90; 15.0; 3.0; 10.0; 10.0], ...
    [q_rmse < 2.0; q_mae < 1.0; q_nse > 0.90; abs(q_peak_error_pct) < 15.0; ...
     abs(peak_timing_error_min) < 3.0; abs(volume_error_pct) < 10.0; abs(storage_error_pct) < 10.0], ...
    'VariableNames', {'Metric','Value','Threshold','Pass'});

mass_balance = table( ...
    ["fine_10m_full_momentum"; "coarse_30m_lookup_subgrid"], ...
    [fine_records.cumulative_inflow_m3(end); coarse_records.cumulative_inflow_m3(end)], ...
    [fine_records.cumulative_outlet_m3(end); coarse_records.cumulative_outlet_m3(end)], ...
    [fine_records.storage_m3(end); coarse_records.storage_m3(end)], ...
    'VariableNames', {'Run','CumulativeInflow_m3','CumulativeOutlet_m3','FinalStorage_m3'});

pass_fail = table( ...
    "P1-SUBGRID-002", ...
    all(metrics.Pass), ...
    string(sprintf('NSE=%.4f; RMSE=%.3f m3/s; peak_error=%.2f%%; volume_error=%.2f%%', ...
        q_nse, q_rmse, q_peak_error_pct, volume_error_pct)), ...
    'VariableNames', {'CaseID','Pass','Notes'});

time_series = table(t_min, Qin, fine_records.outlet_q_m3s, coarse_records.outlet_q_m3s, ...
    fine_records.storage_m3, coarse_records.storage_m3, ...
    fine_records.cumulative_outlet_m3, coarse_records.cumulative_outlet_m3, ...
    'VariableNames', {'Time_min','Inflow_m3s','Outlet_10m_m3s','Outlet_30mSubgrid_m3s', ...
    'Storage_10m_m3','Storage_30mSubgrid_m3','CumulativeOutlet_10m_m3','CumulativeOutlet_30mSubgrid_m3'});

writetable(metrics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(mass_balance, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(pass_fail, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(time_series, fullfile(out_dir, 'Hydrograph_Storage_TimeSeries.csv'));
writetable(fine_profile, fullfile(out_dir, 'CrossSection_Profile_10m.csv'));
writetable(coarse_profile, fullfile(out_dir, 'CrossSection_Profile_30mSubgrid.csv'));

plot_hydrograph(t_min, Qin, fine_records, coarse_records, fig_dir);
plot_storage(t_min, fine_records, coarse_records, fig_dir);
plot_cross_section(geom, fine_profile, coarse_profile, fig_dir);

fprintf('Compound-channel subgrid benchmark complete. Pass = %d\n', pass_fail.Pass);
fprintf('NSE = %.4f, RMSE = %.3f m3/s, peak error = %.2f%%, volume error = %.2f%%\n', ...
    q_nse, q_rmse, q_peak_error_pct, volume_error_pct);

%% Helpers
function Q = triangular_hydrograph(t_min, t_start, t_peak, t_end, Q_peak)
Q = zeros(size(t_min));
rising = t_min >= t_start & t_min <= t_peak;
falling = t_min > t_peak & t_min <= t_end;
Q(rising) = Q_peak .* (t_min(rising) - t_start) ./ (t_peak - t_start);
Q(falling) = Q_peak .* (1 - (t_min(falling) - t_peak) ./ (t_end - t_peak));
Q(Q < 0) = 0;
end

function grid = initialize_model_grid(dx, geom)
ny = round(geom.total_width / dx);
nx = round(geom.length / dx);
y = ((1:ny) - (ny + 1) / 2) .* dx;
x = ((1:nx) - 0.5) .* dx;
[X,Y] = meshgrid(x, y);
lateral = compound_lateral_elevation(Y, geom);
grid.z = lateral - geom.slope_long .* X;
grid.dx = dx;
grid.ny = ny;
grid.nx = nx;
grid.area = dx^2;
grid.roughness = geom.n .* ones(ny, nx);
grid.outlet = make_outlet_properties(ny, nx);
[~, grid.channel_row] = min(abs(y));
grid.upstream_col = 1;
grid.y = y(:);
grid.x = x(:);
end

function zlat = compound_lateral_elevation(y, geom)
half_channel = geom.channel_width / 2;
dist = abs(y);
zlat = zeros(size(y));
overbank = dist > half_channel;
zlat(overbank) = geom.bank_height + geom.overbank_lateral_slope .* (dist(overbank) - half_channel);
end

function Outlet = make_outlet_properties(ny, nx)
Outlet.face_right = false(ny, nx);
Outlet.face_right(:, nx) = true;
Outlet.face_left = false(ny, nx);
Outlet.face_up = false(ny, nx);
Outlet.face_down = false(ny, nx);
Outlet.flag_HR_full_momentum = 0;
end

function state = initialize_state(grid, dt_min, use_subgrid, SubgridTables)
state.d_tot = zeros(grid.ny, grid.nx);
state.d_p = state.d_tot;
state.outflow = zeros(grid.ny, grid.nx, 5);
state.outlet_flow = zeros(grid.ny, grid.nx);
state.C_a = grid.area .* ones(grid.ny, grid.nx);
state.dt_min = dt_min;
if use_subgrid
    state.C_a = grid.area .* ones(grid.ny, grid.nx);
    state.SubgridTables = SubgridTables;
end
end

function state = apply_upstream_inflow(state, grid, qin_m3s, dt_s, use_subgrid, SubgridTables)
if qin_m3s <= 0
    return;
end
row = grid.channel_row;
col = grid.upstream_col;
vol = qin_m3s * dt_s;
if use_subgrid
    drep = max(state.d_tot ./ 1000, 0);
    V = hp2d_subgrid_lookup_depth(SubgridTables.volume_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
    V(row, col) = V(row, col) + vol;
    drep_new = hp2d_subgrid_inverse_volume(SubgridTables.volume_cell, V, SubgridTables.dz, SubgridTables.maxDepth, grid.dx);
    state.d_tot = 1000 .* drep_new;
else
    state.d_tot(row, col) = state.d_tot(row, col) + 1000 * vol / grid.area;
end
state.d_p = state.d_tot;
end

function state = step_full_momentum(state, grid, dt_min, use_subgrid, SubgridTables)
idx_nan = false(grid.ny, grid.nx);
outlet_index = false(grid.ny, grid.nx);
row_outlet = [];
col_outlet = [];
flag_subgrid = double(use_subgrid);
if ~use_subgrid
    SubgridTables = [];
end
[~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = Full_Momentum_Model_D4( ...
    1, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, grid.z, state.d_tot, state.d_p, grid.roughness, grid.roughness.^2, grid.area, ...
    dt_min, grid.dx, outlet_index, 1, 0.005, row_outlet, col_outlet, 1e-6, ...
    state.outflow, idx_nan, 0, flag_subgrid, [], [], [], [], [], [], [], [], [], ...
    [], 0, 0, SubgridTables, grid.outlet);
state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function rec = initialize_records(nt)
rec.outlet_q_m3s = zeros(nt, 1);
rec.storage_m3 = zeros(nt, 1);
rec.mean_depth_m = zeros(nt, 1);
rec.d_tot_store = [];
end

function rec = record_step(rec, it, state, grid, use_subgrid, SubgridTables)
drep = max(state.d_tot ./ 1000, 0);
if use_subgrid
    V = hp2d_subgrid_lookup_depth(SubgridTables.volume_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
    A = hp2d_subgrid_lookup_depth(SubgridTables.area_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
    C_a = A;
    C_a(~isfinite(C_a) | C_a <= 0) = grid.area;
else
    V = drep .* grid.area;
    C_a = grid.area .* ones(size(drep));
end
rec.storage_m3(it) = sum(V(:), 'omitnan');
rec.mean_depth_m(it) = mean(drep(:), 'omitnan');
rec.outlet_q_m3s(it) = sum(state.outlet_flow(:) .* C_a(:), 'omitnan') / 1000 / 3600;
if it == 1
    rec.d_tot_store = zeros([size(state.d_tot), 3]);
end
snap_idx = 0;
if it == 1; snap_idx = 1; end
if it == round(numel(rec.outlet_q_m3s) * 0.60); snap_idx = 2; end
if it == numel(rec.outlet_q_m3s); snap_idx = 3; end
if snap_idx > 0
    rec.d_tot_store(:,:,snap_idx) = state.d_tot;
end
end

function state = fine_state_snapshots(grid, rec, ip)
state = struct();
if ip <= round(numel(rec.outlet_q_m3s) * 0.60)
    state.d_tot = rec.d_tot_store(:,:,2);
else
    state.d_tot = rec.d_tot_store(:,:,3);
end
if isempty(state.d_tot) || all(state.d_tot(:) == 0)
    state.d_tot = rec.d_tot_store(:,:,end);
end
state.C_a = grid.area .* ones(grid.ny, grid.nx);
end

function profile = cross_section_profile(state, grid, use_subgrid, SubgridTables)
mid_col = round(grid.nx / 2);
drep = max(state.d_tot(:, mid_col) ./ 1000, 0);
if use_subgrid
    invert = SubgridTables.invert_el(:, mid_col);
    volume = hp2d_subgrid_lookup_depth(SubgridTables.volume_cell(:,mid_col,:), drep, SubgridTables.dz, SubgridTables.maxDepth);
    wet_area = hp2d_subgrid_lookup_depth(SubgridTables.area_cell(:,mid_col,:), drep, SubgridTables.dz, SubgridTables.maxDepth);
else
    invert = grid.z(:, mid_col);
    volume = drep .* grid.area;
    wet_area = grid.area .* double(drep > 0);
end
profile = table(grid.y, invert, drep, invert + drep, volume, wet_area, ...
    'VariableNames', {'Y_m','InvertElevation_m','RepresentativeDepth_m','WaterSurface_m','Volume_m3','WettedPlanArea_m2'});
end

function T = build_lookup_tables_from_fine_dem(fine_z, fine_dx, ratio, n, dz, maxDepth)
[ny_f, nx_f] = size(fine_z);
ny = ny_f / ratio;
nx = nx_f / ratio;
if mod(ny_f, ratio) ~= 0 || mod(nx_f, ratio) ~= 0
    error('Fine-grid cross-section must be divisible by the coarse ratio.');
end
% Aggregate 3 by 3 fine cells into each 30 m coarse cell.
depth_axis = 0:dz:maxDepth;
nz = numel(depth_axis);
T.volume_cell = zeros(ny, nx, nz);
T.area_cell = zeros(ny, nx, nz);
T.invert_el = zeros(ny, nx);
for i = 1:ny
    rows = (i-1)*ratio + (1:ratio);
    for j = 1:nx
        cols = (j-1)*ratio + (1:ratio);
        block = fine_z(rows, cols);
        inv = min(block(:));
        T.invert_el(i,j) = inv;
        for k = 1:nz
            eta = inv + depth_axis(k);
            h = max(eta - block, 0);
            T.volume_cell(i,j,k) = sum(h(:)) * fine_dx^2;
            T.area_cell(i,j,k) = nnz(h(:) > 0) * fine_dx^2;
        end
    end
end
T.area_x = zeros(ny, nx-1, nz);
T.invert_x = zeros(ny, nx-1);
for i = 1:ny
    rows = (i-1)*ratio + (1:ratio);
    for j = 1:nx-1
        col_left = j * ratio;
        col_right = col_left + 1;
        face_z = max(fine_z(rows, col_left), fine_z(rows, col_right));
        inv = min(face_z(:));
        T.invert_x(i,j) = inv;
        for k = 1:nz
            eta = inv + depth_axis(k);
            T.area_x(i,j,k) = sum(max(eta - face_z, 0)) * fine_dx;
        end
    end
end
T.area_y = zeros(ny-1, nx, nz);
T.invert_y = zeros(ny-1, nx);
for i = 1:ny-1
    row_up = i * ratio;
    row_down = row_up + 1;
    for j = 1:nx
        cols = (j-1)*ratio + (1:ratio);
        face_z = max(fine_z(row_up, cols), fine_z(row_down, cols));
        inv = min(face_z(:));
        T.invert_y(i,j) = inv;
        for k = 1:nz
            eta = inv + depth_axis(k);
            T.area_y(i,j,k) = sum(max(eta - face_z, 0)) * fine_dx;
        end
    end
end
T.n_x = n .* ones(size(T.area_x));
T.n_y = n .* ones(size(T.area_y));
T.dz = dz;
T.maxDepth = maxDepth;
end

function e = rmse(sim, obs)
d = sim(:) - obs(:);
d = d(isfinite(d));
e = sqrt(mean(d.^2));
end

function v = rel_pct(err, denom)
v = 100 .* err ./ max(abs(denom), eps);
end

function val = nse(obs, sim)
obs = obs(:); sim = sim(:);
ok = isfinite(obs) & isfinite(sim);
obs = obs(ok); sim = sim(ok);
den = sum((obs - mean(obs)).^2);
if den <= eps
    val = double(max(abs(obs - sim)) < 1e-12);
else
    val = 1 - sum((obs - sim).^2) ./ den;
end
end

function plot_hydrograph(t_min, Qin, fine, coarse, fig_dir)
fig = figure('Visible','off');
plot(t_min, Qin, 'k:', 'LineWidth', 1.2); hold on;
plot(t_min, fine.outlet_q_m3s, 'k-', 'LineWidth', 1.6);
plot(t_min, coarse.outlet_q_m3s, 'r--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Discharge (m3/s)');
legend('Inflow', '10 m full momentum', '30 m lookup subgrid', 'Location', 'northwest');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_compound_hydrograph.png'));
close(fig);
end

function plot_storage(t_min, fine, coarse, fig_dir)
fig = figure('Visible','off');
plot(t_min, fine.storage_m3, 'k-', 'LineWidth', 1.6); hold on;
plot(t_min, coarse.storage_m3, 'r--', 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Domain storage (m3)');
legend('10 m full momentum', '30 m lookup subgrid', 'Location', 'best');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_compound_storage.png'));
close(fig);
end

function plot_cross_section(geom, fine_profile, coarse_profile, fig_dir)
yy = linspace(-geom.total_width/2, geom.total_width/2, 361)';
zz = compound_lateral_elevation(yy, geom);
fig = figure('Visible','off');
plot(yy, zz, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4); hold on;
plot(fine_profile.Y_m, fine_profile.WaterSurface_m, 'ko-', 'LineWidth', 1.2);
plot(coarse_profile.Y_m, coarse_profile.WaterSurface_m, 'rs--', 'LineWidth', 1.2);
xlabel('Cross-stream coordinate (m)');
ylabel('Elevation or water surface (m)');
legend('Compound-channel bed', '10 m water surface', '30 m subgrid representative surface', 'Location', 'northwest');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_compound_cross_section.png'));
close(fig);
end
