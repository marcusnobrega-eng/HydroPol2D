%RUN_COMPOUND_CHANNEL_STEADY_LADDER_BENCHMARK
% P1-SUBGRID-003: 10 m local-inertial reference vs 30 m lookup-subgrid
% local-inertial steady-flow ladder for a compound channel.

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
geom.length = 900.0;
geom.n = 0.035;

fine_dx = 10.0;
coarse_dx = 30.0;
ratio = coarse_dx / fine_dx;
dt_min = 0.05;          % 3 s
dt_s = dt_min * 60;
max_time_min = 180.0;
min_time_min = 20.0;
steady_window_min = 10.0;
steady_tol = 1e-3;      % 0.1 percent storage change over the window
window_steps = round(steady_window_min / dt_min);
max_steps = round(max_time_min / dt_min);

Q_ladder = [5; 10; 15; 18; 22; 30; 35];
primary_x = 600.0;
secondary_x = 750.0;

fine = initialize_model_grid(fine_dx, geom);
coarse = initialize_model_grid(coarse_dx, geom);
SubgridTables = build_lookup_tables_from_fine_dem(fine.z, fine_dx, ratio, geom.n, 0.02, 5.0);
hp2d_validate_subgrid_tables(SubgridTables);

fine_state = initialize_state(fine);
coarse_state = initialize_state(coarse);

time_rows = table();
metrics = table();
mass_rows = table();
profile_rows = table();

for iq = 1:numel(Q_ladder)
    Qin = Q_ladder(iq);
    regime = classify_regime(Qin);

    d0 = estimate_normal_depth(Qin, geom);
    fine_state = set_uniform_stage_state(fine_state, fine, d0, false, []);
    coarse_state = set_uniform_stage_state(coarse_state, coarse, d0, true, SubgridTables);

    start_storage_fine = compute_storage(fine_state, fine, false, []);
    start_storage_coarse = compute_storage(coarse_state, coarse, true, SubgridTables);

    storage_fine = zeros(max_steps, 1);
    storage_coarse = zeros(max_steps, 1);
    qout_fine = zeros(max_steps, 1);
    qout_coarse = zeros(max_steps, 1);
    q600_fine = zeros(max_steps, 1);
    q600_coarse = zeros(max_steps, 1);
    q750_fine = zeros(max_steps, 1);
    q750_coarse = zeros(max_steps, 1);

    converged = false;
    nsteps = max_steps;

    for it = 1:max_steps
        fine_state = apply_upstream_inflow(fine_state, fine, Qin, dt_s, false, []);
        coarse_state = apply_upstream_inflow(coarse_state, coarse, Qin, dt_s, true, SubgridTables);

        fine_state = step_local_inertial(fine_state, fine, dt_min, false, []);
        coarse_state = step_subgrid_local_inertial(coarse_state, coarse, dt_min, SubgridTables);

        storage_fine(it) = compute_storage(fine_state, fine, false, []);
        storage_coarse(it) = compute_storage(coarse_state, coarse, true, SubgridTables);
        qout_fine(it) = outlet_discharge(fine_state, fine, false, []);
        qout_coarse(it) = outlet_discharge(coarse_state, coarse, true, SubgridTables);
        q600_fine(it) = internal_discharge(fine_state, fine, primary_x);
        q600_coarse(it) = internal_discharge(coarse_state, coarse, primary_x);
        q750_fine(it) = internal_discharge(fine_state, fine, secondary_x);
        q750_coarse(it) = internal_discharge(coarse_state, coarse, secondary_x);

        if it > window_steps && it * dt_min >= min_time_min
            rel_fine = abs(storage_fine(it) - storage_fine(it-window_steps)) / max(storage_fine(it), 1.0);
            rel_coarse = abs(storage_coarse(it) - storage_coarse(it-window_steps)) / max(storage_coarse(it), 1.0);
            if rel_fine < steady_tol && rel_coarse < steady_tol
                converged = true;
                nsteps = it;
                break;
            end
        end
    end

    idx = (1:nsteps)';
    t_local = idx .* dt_min;
    time_rows = [time_rows; table( ...
        repmat(Qin, nsteps, 1), repmat(string(regime), nsteps, 1), t_local, ...
        storage_fine(idx), storage_coarse(idx), qout_fine(idx), qout_coarse(idx), ...
        q600_fine(idx), q600_coarse(idx), q750_fine(idx), q750_coarse(idx), ...
        'VariableNames', {'Qin_m3s','Regime','Time_min','Storage10m_m3','Storage30mSubgrid_m3', ...
        'Outlet10m_m3s','Outlet30mSubgrid_m3s','Q600_10m_m3s','Q600_30mSubgrid_m3s', ...
        'Q750_10m_m3s','Q750_30mSubgrid_m3s'})]; %#ok<AGROW>

    final_fine = summarize_state(fine_state, fine, primary_x, secondary_x, false, []);
    final_coarse = summarize_state(coarse_state, coarse, primary_x, secondary_x, true, SubgridTables);

    [prof_fine, prof_coarse] = longitudinal_profiles(fine_state, fine, coarse_state, coarse, SubgridTables);
    profile_rmse = rmse(prof_coarse.WaterSurface_m, prof_fine.WaterSurface_m);

    depth_error_primary = final_coarse.Stage600_m - final_fine.Stage600_m;
    depth_error_secondary = final_coarse.Stage750_m - final_fine.Stage750_m;
    depth_rmse = sqrt(mean([depth_error_primary, depth_error_secondary].^2));
    storage_error_pct = rel_pct(final_coarse.Storage_m3 - final_fine.Storage_m3, final_fine.Storage_m3);
    q600_error_pct = rel_pct(final_coarse.Q600_m3s - final_fine.Q600_m3s, max(abs(final_fine.Q600_m3s), eps));
    q750_error_pct = rel_pct(final_coarse.Q750_m3s - final_fine.Q750_m3s, max(abs(final_fine.Q750_m3s), eps));
    internal_q_error_pct = max(abs([q600_error_pct, q750_error_pct]));

    input_volume = Qin * nsteps * dt_s;
    out_volume_fine = sum(qout_fine(idx)) * dt_s;
    out_volume_coarse = sum(qout_coarse(idx)) * dt_s;
    residual_fine_pct = rel_pct(input_volume - out_volume_fine - (final_fine.Storage_m3 - start_storage_fine), input_volume);
    residual_coarse_pct = rel_pct(input_volume - out_volume_coarse - (final_coarse.Storage_m3 - start_storage_coarse), input_volume);
    mass_residual_pct = max(abs([residual_fine_pct, residual_coarse_pct]));

    if Qin <= 15
        depth_threshold = 0.05;
        storage_threshold = 5.0;
        q_threshold = 5.0;
    else
        depth_threshold = 0.10;
        storage_threshold = 10.0;
        q_threshold = 10.0;
    end

    pass_flow = depth_rmse < depth_threshold && profile_rmse < depth_threshold && ...
        abs(storage_error_pct) < storage_threshold && internal_q_error_pct < q_threshold && ...
        mass_residual_pct < 0.1 && converged;

    metrics = [metrics; table( ...
        Qin, string(regime), converged, nsteps * dt_min, depth_rmse, profile_rmse, ...
        storage_error_pct, internal_q_error_pct, mass_residual_pct, logical(pass_flow), ...
        'VariableNames', {'Qin_m3s','Regime','Converged','RunTime_min','GaugeDepthRMSE_m', ...
        'ProfileRMSE_m','StorageError_pct','InternalDischargeError_pct','MassResidual_pct','Pass'})]; %#ok<AGROW>

    mass_rows = [mass_rows; table( ...
        Qin, string(regime), input_volume, out_volume_fine, out_volume_coarse, ...
        final_fine.Storage_m3, final_coarse.Storage_m3, residual_fine_pct, residual_coarse_pct, ...
        'VariableNames', {'Qin_m3s','Regime','InputVolume_m3','OutletVolume10m_m3', ...
        'OutletVolume30mSubgrid_m3','FinalStorage10m_m3','FinalStorage30mSubgrid_m3', ...
        'MassResidual10m_pct','MassResidual30mSubgrid_pct'})]; %#ok<AGROW>

    profile_rows = [profile_rows; add_profile_rows(Qin, regime, prof_fine, prof_coarse)]; %#ok<AGROW>

    fprintf('Q = %.1f m3/s (%s): pass=%d, depth RMSE=%.3f m, profile RMSE=%.3f m, storage err=%.2f%%, q err=%.2f%%, t=%.1f min\n', ...
        Qin, regime, pass_flow, depth_rmse, profile_rmse, storage_error_pct, internal_q_error_pct, nsteps * dt_min);
end

overall_pass = all(metrics.Pass);
pass_fail = table("P1-SUBGRID-003", overall_pass, ...
    string(sprintf('passed_flows=%d/%d; max_depth_rmse=%.3f m; max_storage_error=%.2f%%; max_q_error=%.2f%%', ...
    nnz(metrics.Pass), height(metrics), max(metrics.GaugeDepthRMSE_m), ...
    max(abs(metrics.StorageError_pct)), max(abs(metrics.InternalDischargeError_pct)))), ...
    'VariableNames', {'CaseID','Pass','Notes'});

writetable(metrics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(mass_rows, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(pass_fail, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(time_rows, fullfile(out_dir, 'Steady_Ladder_TimeSeries.csv'));
writetable(profile_rows, fullfile(out_dir, 'Steady_Ladder_FinalProfiles.csv'));

plot_q_stage(metrics, profile_rows, primary_x, secondary_x, fig_dir);
plot_profile_set(profile_rows, [5 18 35], fig_dir);
plot_storage_ladder(metrics, mass_rows, fig_dir);
plot_cross_sections(geom, fine, coarse, SubgridTables, fig_dir);

fprintf('Steady ladder complete. Overall pass = %d\n', overall_pass);

%% Local functions
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
grid.x = x(:);
grid.y = y(:);
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

function state = initialize_state(grid)
state.d_tot = zeros(grid.ny, grid.nx);
state.d_p = state.d_tot;
state.outflow = zeros(grid.ny, grid.nx, 5);
state.outlet_flow = zeros(grid.ny, grid.nx);
state.C_a = grid.area .* ones(grid.ny, grid.nx);
end

function state = set_uniform_stage_state(state, grid, d_thalweg, use_subgrid, SubgridTables)
if use_subgrid
    thalweg = min(SubgridTables.invert_el, [], 1);
    eta = repmat(thalweg + d_thalweg, grid.ny, 1);
    drep = max(eta - SubgridTables.invert_el, 0);
    state.d_tot = 1000 .* drep;
else
    thalweg = min(grid.z, [], 1);
    eta = repmat(thalweg + d_thalweg, grid.ny, 1);
    state.d_tot = 1000 .* max(eta - grid.z, 0);
end
state.d_p = state.d_tot;
state.outflow(:) = 0;
state.outlet_flow(:) = 0;
end

function state = apply_upstream_inflow(state, grid, qin_m3s, dt_s, use_subgrid, SubgridTables)
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
flag_subgrid = double(use_subgrid);
if ~use_subgrid
    SubgridTables = [];
end

[~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = Full_Momentum_Model_D4( ...
    1, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, grid.z, state.d_tot, state.d_p, grid.roughness, grid.roughness.^2, grid.area, ...
    dt_min, grid.dx, outlet_index, 1, 0.005, [], [], 1e-6, ...
    state.outflow, idx_nan, 0, flag_subgrid, [], [], [], [], [], [], [], [], [], ...
    [], 0, 0, SubgridTables, grid.outlet);
state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function state = step_local_inertial(state, grid, dt_min, use_subgrid, SubgridTables)
idx_nan = false(grid.ny, grid.nx);
outlet_index = false(grid.ny, grid.nx);
[row_outlet, col_outlet] = find(grid.outlet.face_right);
roughness = grid.roughness;
legacy_channel_width = zeros(grid.ny, grid.nx);
legacy_channel_depth = zeros(grid.ny, grid.nx);
flag_subgrid = double(use_subgrid);
if ~use_subgrid
    SubgridTables = [];
end
[~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = Local_Inertial_Model_D4( ...
    2, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, grid.z, state.d_tot, state.d_p, roughness, roughness.^2, grid.area, ...
    dt_min, grid.dx, outlet_index, 1, 0.005, row_outlet, col_outlet, 1e-6, ...
    state.outflow, idx_nan, 0, flag_subgrid, roughness, [], legacy_channel_width, legacy_channel_depth, [], [], [], [], ...
    state.C_a, [], 0, 0, SubgridTables);
state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function S = compute_storage(state, grid, use_subgrid, SubgridTables)
drep = max(state.d_tot ./ 1000, 0);
if use_subgrid
    V = hp2d_subgrid_lookup_depth(SubgridTables.volume_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
else
    V = drep .* grid.area;
end
S = sum(V(:), 'omitnan');
end

function q = outlet_discharge(state, grid, use_subgrid, SubgridTables)
drep = max(state.d_tot ./ 1000, 0);
if use_subgrid
    A = hp2d_subgrid_lookup_depth(SubgridTables.area_cell, drep, SubgridTables.dz, SubgridTables.maxDepth);
    C_a = A;
    C_a(~isfinite(C_a) | C_a <= 0) = grid.area;
else
    C_a = grid.area .* ones(size(drep));
end
q = sum(state.outlet_flow(:) .* C_a(:), 'omitnan') / 1000 / 3600;
end

function q = internal_discharge(state, grid, x_gauge)
face_col = max(1, min(grid.nx - 1, round(x_gauge / grid.dx)));
q = sum(state.outflow(:, face_col, 1), 'omitnan') * grid.area / 1000 / 3600;
end

function summary = summarize_state(state, grid, x1, x2, use_subgrid, SubgridTables)
summary = struct();
summary.Storage_m3 = compute_storage(state, grid, use_subgrid, SubgridTables);
summary.Q600_m3s = internal_discharge(state, grid, x1);
summary.Q750_m3s = internal_discharge(state, grid, x2);
summary.Stage600_m = gauge_stage(state, grid, x1, use_subgrid, SubgridTables);
summary.Stage750_m = gauge_stage(state, grid, x2, use_subgrid, SubgridTables);
end

function stage = gauge_stage(state, grid, x_gauge, use_subgrid, SubgridTables)
col = nearest_col(grid, x_gauge);
row = grid.channel_row;
drep = max(state.d_tot(row, col) ./ 1000, 0);
if use_subgrid
    stage = SubgridTables.invert_el(row, col) + drep;
else
    stage = grid.z(row, col) + drep;
end
end

function col = nearest_col(grid, x_gauge)
[~, col] = min(abs(grid.x - x_gauge));
end

function [prof_fine_on_coarse, prof_coarse] = longitudinal_profiles(fine_state, fine, coarse_state, coarse, SubgridTables)
row_f = fine.channel_row;
row_c = coarse.channel_row;
fine_col = arrayfun(@(x) nearest_col(fine, x), coarse.x);
d_f = max(fine_state.d_tot(row_f, fine_col) ./ 1000, 0)';
ws_f = fine.z(row_f, fine_col)' + d_f;
d_c = max(coarse_state.d_tot(row_c, :) ./ 1000, 0)';
ws_c = SubgridTables.invert_el(row_c, :)' + d_c;
prof_fine_on_coarse = table(coarse.x, ws_f, d_f, 'VariableNames', {'X_m','WaterSurface_m','Depth_m'});
prof_coarse = table(coarse.x, ws_c, d_c, 'VariableNames', {'X_m','WaterSurface_m','Depth_m'});
end

function rows = add_profile_rows(Qin, regime, fine_profile, coarse_profile)
n = height(fine_profile);
rows = table(repmat(Qin,n,1), repmat(string(regime),n,1), fine_profile.X_m, ...
    fine_profile.WaterSurface_m, coarse_profile.WaterSurface_m, ...
    fine_profile.Depth_m, coarse_profile.Depth_m, ...
    'VariableNames', {'Qin_m3s','Regime','X_m','WaterSurface10m_m', ...
    'WaterSurface30mSubgrid_m','Depth10m_m','Depth30mSubgrid_m'});
end

function state = step_subgrid_local_inertial(state, grid, dt_min, SubgridTables)
idx_nan = false(grid.ny, grid.nx);
outlet_index = false(grid.ny, grid.nx);
[row_outlet, col_outlet] = find(grid.outlet.face_right);
roughness = grid.roughness;
legacy_channel_width = zeros(grid.ny, grid.nx);
legacy_channel_depth = zeros(grid.ny, grid.nx);
[~,~,~,~,outlet_flow,d_t,~,outflow,~,~,~,~,~,C_a] = Local_Inertial_Model_D4( ...
    2, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, grid.z, state.d_tot, state.d_p, roughness, roughness.^2, grid.area, ...
    dt_min, grid.dx, outlet_index, 1, 0.005, row_outlet, col_outlet, 1e-6, ...
    state.outflow, idx_nan, 0, 1, roughness, [], legacy_channel_width, legacy_channel_depth, [], [], [], [], ...
    state.C_a, [], 0, 0, SubgridTables);
state.d_p = state.d_tot;
state.d_tot = d_t;
state.outflow = outflow;
state.outlet_flow = outlet_flow;
state.C_a = C_a;
end

function T = build_lookup_tables_from_fine_dem(fine_z, fine_dx, ratio, n, dz, maxDepth)
[ny_f, nx_f] = size(fine_z);
ny = ny_f / ratio;
nx = nx_f / ratio;
if mod(ny_f, ratio) ~= 0 || mod(nx_f, ratio) ~= 0
    error('Fine grid dimensions must be divisible by the coarse ratio.');
end
depth_axis = 0:dz:maxDepth;
nz = numel(depth_axis);
coarse_face_width = fine_dx * ratio;
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
T.width_x = zeros(ny, nx-1, nz);
T.wetfrac_x = zeros(ny, nx-1, nz);
T.perimeter_x = zeros(ny, nx-1, nz);
T.Rh_x = zeros(ny, nx-1, nz);
T.K_x = zeros(ny, nx-1, nz);
T.n_x = NaN(ny, nx-1, nz);
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
            [A,W,phi,P,Rh,K,n_eff] = raster_face_hydraulics(face_z(:), eta, fine_dx, n, coarse_face_width);
            T.area_x(i,j,k) = A;
            T.width_x(i,j,k) = W;
            T.wetfrac_x(i,j,k) = phi;
            T.perimeter_x(i,j,k) = P;
            T.Rh_x(i,j,k) = Rh;
            T.K_x(i,j,k) = K;
            T.n_x(i,j,k) = n_eff;
        end
    end
end
T.area_y = zeros(ny-1, nx, nz);
T.width_y = zeros(ny-1, nx, nz);
T.wetfrac_y = zeros(ny-1, nx, nz);
T.perimeter_y = zeros(ny-1, nx, nz);
T.Rh_y = zeros(ny-1, nx, nz);
T.K_y = zeros(ny-1, nx, nz);
T.n_y = NaN(ny-1, nx, nz);
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
            [A,W,phi,P,Rh,K,n_eff] = raster_face_hydraulics(face_z(:), eta, fine_dx, n, coarse_face_width);
            T.area_y(i,j,k) = A;
            T.width_y(i,j,k) = W;
            T.wetfrac_y(i,j,k) = phi;
            T.perimeter_y(i,j,k) = P;
            T.Rh_y(i,j,k) = Rh;
            T.K_y(i,j,k) = K;
            T.n_y(i,j,k) = n_eff;
        end
    end
end
T.dz = dz;
T.maxDepth = maxDepth;
end

function [A,W,phi,P,Rh,K,n_eff] = raster_face_hydraulics(z_face, eta, dx, n, coarse_face_width)
d = eta - z_face(:);
d(d < 0) = 0;
wet = d > 0;
A = sum(d) * dx;
W = sum(wet) * dx;
phi = min(max(W / max(coarse_face_width, eps), 0), 1);
if A <= 0 || ~any(wet)
    P = 0;
    Rh = 0;
    K = 0;
    n_eff = NaN;
    return;
end
P = W;
if wet(1)
    P = P + d(1);
end
if wet(end)
    P = P + d(end);
end
for ii = 1:(numel(z_face)-1)
    lower_z = min(z_face(ii), z_face(ii+1));
    upper_z = max(z_face(ii), z_face(ii+1));
    if eta > lower_z
        P = P + max(min(eta, upper_z) - lower_z, 0);
    end
end
Rh = A / max(P, eps);
K = (1/n) * A * Rh^(2/3);
H_G = A / max(coarse_face_width, eps);
denom = mean(d.^(5/3) ./ n);
if H_G > 0 && denom > 0 && isfinite(denom)
    n_eff = H_G^(5/3) / denom;
else
    n_eff = NaN;
end
end

function d = estimate_normal_depth(Q, geom)
lo = 1e-4;
hi = 5.0;
for k = 1:80
    mid = 0.5 * (lo + hi);
    if compound_manning_q(mid, geom) < Q
        lo = mid;
    else
        hi = mid;
    end
end
d = 0.5 * (lo + hi);
end

function Q = compound_manning_q(d, geom)
[A, P] = compound_area_perimeter(d, geom);
R = A / max(P, eps);
Q = (1 / geom.n) * A * R^(2/3) * sqrt(geom.slope_long);
end

function [A, P] = compound_area_perimeter(d, geom)
B = geom.channel_width;
H = geom.bank_height;
S = geom.overbank_lateral_slope;
W = geom.overbank_width_each;
if d <= H
    A = B * d;
    P = B + 2*d;
    return;
end
df = d - H;
L = min(W, df / S);
A_one = df * L - 0.5 * S * L^2;
if df > S * W
    A_one = A_one + (df - S * W) * W;
end
A = B * d + 2 * A_one;
P = B + 2 * H + 2 * L * sqrt(1 + S^2);
if df > S * W
    P = P + 2 * (df - S * W);
end
end

function regime = classify_regime(Q)
if Q <= 15
    regime = "inbank";
elseif Q <= 18
    regime = "near_bankfull";
else
    regime = "overbank";
end
end

function e = rmse(sim, obs)
d = sim(:) - obs(:);
d = d(isfinite(d));
e = sqrt(mean(d.^2));
end

function v = rel_pct(err, denom)
v = 100 .* err ./ max(abs(denom), eps);
end

function plot_q_stage(metrics, profile_rows, primary_x, secondary_x, fig_dir)
fig = figure('Visible','off');
for k = 1:2
    if k == 1
        xg = primary_x;
        ttl = 'Primary gauge';
    else
        xg = secondary_x;
        ttl = 'Secondary gauge';
    end
    subplot(1,2,k);
    Qs = unique(profile_rows.Qin_m3s);
    s10 = zeros(size(Qs));
    s30 = zeros(size(Qs));
    for i = 1:numel(Qs)
        rows = profile_rows(profile_rows.Qin_m3s == Qs(i), :);
        [~, ix] = min(abs(rows.X_m - xg));
        s10(i) = rows.WaterSurface10m_m(ix);
        s30(i) = rows.WaterSurface30mSubgrid_m(ix);
    end
    plot(Qs, s10, 'ko-', 'LineWidth', 1.3); hold on;
    plot(Qs, s30, 'rs--', 'LineWidth', 1.3);
    xlabel('Inflow (m3/s)');
    ylabel('Water surface elevation (m)');
    title(ttl);
    grid on;
end
legend('10 m reference', '30 m lookup subgrid', 'Location', 'best');
saveas(fig, fullfile(fig_dir, 'p1_subgrid_steady_ladder_q_stage.png'));
close(fig);
end

function plot_profile_set(profile_rows, Qset, fig_dir)
fig = figure('Visible','off');
tiledlayout(numel(Qset),1);
for i = 1:numel(Qset)
    rows = profile_rows(profile_rows.Qin_m3s == Qset(i), :);
    nexttile;
    plot(rows.X_m, rows.WaterSurface10m_m, 'k-', 'LineWidth', 1.3); hold on;
    plot(rows.X_m, rows.WaterSurface30mSubgrid_m, 'r--', 'LineWidth', 1.3);
    ylabel('WSE (m)');
    title(sprintf('Q = %.0f m3/s', Qset(i)));
    grid on;
end
xlabel('Distance downstream (m)');
legend('10 m reference', '30 m lookup subgrid', 'Location', 'best');
saveas(fig, fullfile(fig_dir, 'p1_subgrid_steady_ladder_profiles.png'));
close(fig);
end

function plot_storage_ladder(metrics, mass_rows, fig_dir)
fig = figure('Visible','off');
yyaxis left;
plot(mass_rows.Qin_m3s, mass_rows.FinalStorage10m_m3, 'ko-', 'LineWidth', 1.3); hold on;
plot(mass_rows.Qin_m3s, mass_rows.FinalStorage30mSubgrid_m3, 'rs--', 'LineWidth', 1.3);
ylabel('Final storage (m3)');
yyaxis right;
plot(metrics.Qin_m3s, metrics.StorageError_pct, 'b^-', 'LineWidth', 1.2);
ylabel('Storage error (%)');
xlabel('Inflow (m3/s)');
legend('10 m storage', '30 m subgrid storage', 'Storage error', 'Location', 'northwest');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_steady_ladder_storage.png'));
close(fig);
end

function plot_cross_sections(geom, fine, coarse, SubgridTables, fig_dir)
yy = linspace(-geom.total_width/2, geom.total_width/2, 361)';
zz = compound_lateral_elevation(yy, geom);
fig = figure('Visible','off');
plot(yy, zz, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4); hold on;
plot(fine.y, compound_lateral_elevation(fine.y, geom), 'ko', 'MarkerSize', 4);
plot(coarse.y, SubgridTables.invert_el(:, round(coarse.nx/2)), 'rs', 'MarkerSize', 6);
xlabel('Cross-stream coordinate (m)');
ylabel('Bed or invert elevation (m)');
legend('Analytical bed', '10 m cell centers', '30 m lookup cell inverts', 'Location', 'northwest');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_steady_ladder_cross_section.png'));
close(fig);
end
