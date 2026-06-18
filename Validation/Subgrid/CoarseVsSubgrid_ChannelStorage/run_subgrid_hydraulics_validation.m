%RUN_SUBGRID_HYDRAULICS_VALIDATION Phase 1 lookup-table subgrid validation.
%
% This driver validates the shared-face lookup-table subgrid hydraulics with
% controlled synthetic tables. The legacy River_Width/River_Depth overbank
% pathway is intentionally excluded from this Phase 1 case.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(fileparts(fileparts(case_dir)));
func_dir = fullfile(model_dir, 'HydroPol2D_Functions');
addpath(func_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
fig_dir = fullfile(case_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

Resolution = 20.0;        % [m]
channel_width = 5.0;      % [m]
bankfull_depth = 1.0;     % [m]
manning_n = 0.035;        % [s m^(-1/3)]
dz = 0.01;                % [m]
maxDepth = 2.0;           % [m]
ny = 3;
nx = 4;

T_incised = make_incised_tables(ny, nx, Resolution, channel_width, ...
    bankfull_depth, manning_n, dz, maxDepth);
T_flat = make_flat_tables(ny, nx, Resolution, manning_n, dz, maxDepth);

hp2d_validate_subgrid_tables(T_incised);
hp2d_validate_subgrid_tables(T_flat);

metrics = table();
mass_rows = table();
pass_rows = table();

%% P1-SUBGRID-001A: storage curve and inverse storage.
depth_samples = reshape([0 0.05 0.25 0.80 1.00 1.35 1.90 2.00 2.20 0.63 1.12 1.75], ny, nx);
depth_samples_clamped = min(depth_samples, maxDepth);
V_model = hp2d_subgrid_lookup_depth(T_incised.volume_cell, depth_samples, dz, maxDepth);
A_model = hp2d_subgrid_lookup_depth(T_incised.area_cell, depth_samples, dz, maxDepth);
V_truth = rectangular_incised_volume(depth_samples_clamped, Resolution, channel_width, bankfull_depth);
A_truth = rectangular_incised_plan_area(depth_samples_clamped, Resolution, channel_width, bankfull_depth);
d_inverse = hp2d_subgrid_inverse_volume(T_incised.volume_cell, V_truth, dz, maxDepth, Resolution);

vol_rmse = rmse(V_model, V_truth);
area_rmse = rmse(A_model, A_truth);
inverse_rmse = rmse(d_inverse, depth_samples_clamped);
max_vol_error = max(abs(V_model(:) - V_truth(:)));
mono_violations = count_negative_diff(T_incised.volume_cell);

metrics = append_metric(metrics, 'P1-SUBGRID-001A', 'volume_rmse_m3', vol_rmse, 1e-9, vol_rmse < 1e-9);
metrics = append_metric(metrics, 'P1-SUBGRID-001A', 'wetted_area_rmse_m2', area_rmse, 1e-9, area_rmse < 1e-9);
metrics = append_metric(metrics, 'P1-SUBGRID-001A', 'inverse_depth_rmse_m', inverse_rmse, 1e-9, inverse_rmse < 1e-9);
metrics = append_metric(metrics, 'P1-SUBGRID-001A', 'max_volume_error_m3', max_vol_error, 1e-8, max_vol_error < 1e-8);
metrics = append_metric(metrics, 'P1-SUBGRID-001A', 'monotonicity_violations', mono_violations, 0, mono_violations == 0);

%% P1-SUBGRID-001B: shared-face conveyance.
face_depth = 0.75 * ones(ny, nx-1);
face_eta = face_depth;
face_area = hp2d_subgrid_lookup_eta(T_incised.area_x, face_eta, T_incised.invert_x, dz, maxDepth);
face_width = hp2d_subgrid_lookup_eta(T_incised.width_x, face_eta, T_incised.invert_x, dz, maxDepth);
face_perimeter = hp2d_subgrid_lookup_eta(T_incised.perimeter_x, face_eta, T_incised.invert_x, dz, maxDepth);
face_Rh = hp2d_subgrid_lookup_eta(T_incised.Rh_x, face_eta, T_incised.invert_x, dz, maxDepth);
face_n = hp2d_subgrid_lookup_eta(T_incised.n_x, face_eta, T_incised.invert_x, dz, maxDepth);
K_model = hp2d_subgrid_lookup_eta(T_incised.K_x, face_eta, T_incised.invert_x, dz, maxDepth);
slope = 0.001;
Q_model = K_model .* sqrt(slope);
face_area_truth = channel_width .* face_depth;
face_width_truth = channel_width .* ones(size(face_depth));
face_perimeter_truth = channel_width + 2 .* face_depth;
face_Rh_truth = face_area_truth ./ face_perimeter_truth;
K_truth = (1 / manning_n) .* face_area_truth .* face_Rh_truth.^(2/3);
Q_truth = K_truth .* sqrt(slope);

face_area_error_pct = rel_pct(max(abs(face_area(:) - face_area_truth(:))), max(face_area_truth(:)));
face_width_error_pct = rel_pct(max(abs(face_width(:) - face_width_truth(:))), max(face_width_truth(:)));
face_perimeter_error_pct = rel_pct(max(abs(face_perimeter(:) - face_perimeter_truth(:))), max(face_perimeter_truth(:)));
face_Rh_error_pct = rel_pct(max(abs(face_Rh(:) - face_Rh_truth(:))), max(face_Rh_truth(:)));
conveyance_error_pct = rel_pct(max(abs(K_model(:) - K_truth(:))), max(K_truth(:)));
discharge_error_pct = rel_pct(max(abs(Q_model(:) - Q_truth(:))), max(Q_truth(:)));

metrics = append_metric(metrics, 'P1-SUBGRID-001B', 'face_area_error_pct', face_area_error_pct, 0.1, face_area_error_pct < 0.1);
metrics = append_metric(metrics, 'P1-SUBGRID-001B', 'face_width_error_pct', face_width_error_pct, 0.1, face_width_error_pct < 0.1);
metrics = append_metric(metrics, 'P1-SUBGRID-001B', 'face_perimeter_error_pct', face_perimeter_error_pct, 0.1, face_perimeter_error_pct < 0.1);
metrics = append_metric(metrics, 'P1-SUBGRID-001B', 'face_Rh_error_pct', face_Rh_error_pct, 0.1, face_Rh_error_pct < 0.1);
metrics = append_metric(metrics, 'P1-SUBGRID-001B', 'conveyance_error_pct', conveyance_error_pct, 0.1, conveyance_error_pct < 0.1);
metrics = append_metric(metrics, 'P1-SUBGRID-001B', 'discharge_error_pct', discharge_error_pct, 0.1, discharge_error_pct < 0.1);

%% P1-SUBGRID-001C: areal source/sink volume conversion.
d_old = reshape([0.20 0.40 0.75 1.05 1.30 1.70 0.15 0.95 1.50 0.55 1.10 1.85], ny, nx);
delta_mm = reshape([20 -5 12 -10 8 -7 15 -8 5 -3 11 -4], ny, nx);
[d_new, V_old, V_new, ~, residual] = hp2d_subgrid_apply_volume_change( ...
    d_old, delta_mm, T_incised, Resolution^2);
expected_delta = delta_mm ./ 1000 .* Resolution^2;
volume_residual = max(abs((V_new(:) - V_old(:)) - expected_delta(:)));
d_truth = hp2d_subgrid_inverse_volume(T_incised.volume_cell, V_old + expected_delta, dz, maxDepth, Resolution);
depth_error = max(abs(d_new(:) - d_truth(:)));
negative_storage_count = nnz(V_new(:) < -1e-12);

metrics = append_metric(metrics, 'P1-SUBGRID-001C', 'source_sink_volume_residual_m3', volume_residual, 1e-9, volume_residual < 1e-9);
metrics = append_metric(metrics, 'P1-SUBGRID-001C', 'representative_depth_error_m', depth_error, 1e-9, depth_error < 1e-9);
metrics = append_metric(metrics, 'P1-SUBGRID-001C', 'negative_storage_count', negative_storage_count, 0, negative_storage_count == 0);

%% P1-SUBGRID-001D: local-inertial shared-face coupling.
eta_line = repmat([1.00 0.96 0.92 0.88], ny, 1);
q_prev = zeros(ny, nx, 2);
idx_nan = false(ny, nx);
dt_s = 0.25;
h_min = 1e-6;
[q_face, Hf_x, ~, Wf_x, ~] = subgrid_topography_model( ...
    1, eta_line, q_prev, manning_n, Resolution, dt_s, h_min, ...
    idx_nan, [], T_incised);
Sx = (eta_line(:,2:end) - eta_line(:,1:end-1)) ./ Resolution;
face_A_for_q = hp2d_subgrid_lookup_eta(T_incised.area_x, max(eta_line(:,1:end-1), eta_line(:,2:end)), T_incised.invert_x, dz, maxDepth);
expected_q = -9.81 .* (face_A_for_q ./ Wf_x(:,1:end-1)) .* dt_s .* Sx;
local_q_error = max(abs(q_face(:,1:end-1,1) - expected_q), [], 'all');
Vol_Flux = dt_s .* ( ...
    [zeros(ny,1), q_face(:,1:end-1,1) .* Wf_x(:,1:end-1)] - ...
    [q_face(:,1:end-1,1) .* Wf_x(:,1:end-1), zeros(ny,1)]);
closed_volume_residual = abs(sum(Vol_Flux(:)));

metrics = append_metric(metrics, 'P1-SUBGRID-001D', 'local_inertial_q_error_m2_s', local_q_error, 1e-12, local_q_error < 1e-12);
metrics = append_metric(metrics, 'P1-SUBGRID-001D', 'closed_domain_volume_residual_m3', closed_volume_residual, 1e-9, closed_volume_residual < 1e-9);
mass_rows = append_mass(mass_rows, 'P1-SUBGRID-001D', 0, 0, closed_volume_residual, 0);

%% P1-SUBGRID-001E: full-momentum flat lookup equivalence.
d0_mm = 1000 .* repmat([0.30 0.25 0.20 0.15], ny, 1);
[d_coarse, out_coarse] = run_full_momentum_once(false, T_flat, d0_mm, Resolution, manning_n);
[d_subgrid_flat, out_subgrid_flat] = run_full_momentum_once(true, T_flat, d0_mm, Resolution, manning_n);

flat_depth_rmse = rmse(d_subgrid_flat ./ 1000, d_coarse ./ 1000);
flat_flux_rmse = rmse(out_subgrid_flat(:,:,1:2), out_coarse(:,:,1:2));
flat_nse = nse(out_coarse(:,:,1), out_subgrid_flat(:,:,1));

metrics = append_metric(metrics, 'P1-SUBGRID-001E', 'full_momentum_flat_depth_rmse_m', flat_depth_rmse, 1e-10, flat_depth_rmse < 1e-10);
metrics = append_metric(metrics, 'P1-SUBGRID-001E', 'full_momentum_flat_hydrograph_rmse_mmh', flat_flux_rmse, 1e-8, flat_flux_rmse < 1e-8);
metrics = append_metric(metrics, 'P1-SUBGRID-001E', 'full_momentum_flat_flux_nse', flat_nse, 0.999, flat_nse > 0.999);

%% P1-SUBGRID-001F: full-momentum incised-channel volume diagnostic.
d_incised0_mm = 1000 .* repmat([0.75 0.65 0.55 0.45], ny, 1);
[d_incised1, ~] = run_full_momentum_once(true, T_incised, d_incised0_mm, Resolution, manning_n);
V0 = hp2d_subgrid_lookup_depth(T_incised.volume_cell, d_incised0_mm ./ 1000, dz, maxDepth);
V1 = hp2d_subgrid_lookup_depth(T_incised.volume_cell, d_incised1 ./ 1000, dz, maxDepth);
incised_mass_residual = abs(sum(V1(:)) - sum(V0(:)));
incised_mass_error_pct = rel_pct(incised_mass_residual, sum(V0(:)));

metrics = append_metric(metrics, 'P1-SUBGRID-001F', 'full_momentum_incised_mass_residual_m3', incised_mass_residual, 1e-6, incised_mass_residual < 1e-6);
metrics = append_metric(metrics, 'P1-SUBGRID-001F', 'full_momentum_incised_mass_error_pct', incised_mass_error_pct, 0.1, incised_mass_error_pct < 0.1);
mass_rows = append_mass(mass_rows, 'P1-SUBGRID-001F', sum(V0(:)), sum(V1(:)), incised_mass_residual, incised_mass_error_pct);

%% P1-SUBGRID-001G: low-friction local-inertial stability diagnostic.
T_low_n = make_flat_tables(ny, nx, Resolution, 0.015, dz, maxDepth);
hp2d_validate_subgrid_tables(T_low_n);
eta_low = repmat(linspace(0.35, 0.20, nx), ny, 1);
q_prev_low = zeros(ny, nx, 2);
[q_low, ~, ~, W_low_x, W_low_y] = subgrid_topography_model( ...
    2, eta_low, q_prev_low, 0.015, Resolution, 1.0, h_min, idx_nan, [], T_low_n);
Q_low = zeros(size(q_low));
Q_low(:,:,1) = q_low(:,:,1) .* W_low_x;
Q_low(:,:,2) = q_low(:,:,2) .* W_low_y;
Vol_low = 1.0 .* ( ...
    [zeros(ny,1), Q_low(:,1:(nx-1),1)] - Q_low(:,:,1) ...
    - Q_low(:,:,2) + [Q_low(2:end,:,2); zeros(1,nx)] );
low_finite = all(isfinite(q_low(:)));
low_closed_residual = abs(sum(Vol_low(:), 'omitnan'));
low_max_q = max(abs(q_low(:)), [], 'omitnan');

metrics = append_metric(metrics, 'P1-SUBGRID-001G', 'low_friction_finite_flag', double(low_finite), 1, low_finite);
metrics = append_metric(metrics, 'P1-SUBGRID-001G', 'low_friction_closed_volume_residual_m3', low_closed_residual, 1e-9, low_closed_residual < 1e-9);
metrics = append_metric(metrics, 'P1-SUBGRID-001G', 'low_friction_max_unit_q_m2_s', low_max_q, 10, low_max_q < 10);
mass_rows = append_mass(mass_rows, 'P1-SUBGRID-001G', 0, 0, low_closed_residual, 0);

%% Write outputs and figures.
case_ids = unique(metrics.CaseID, 'stable');
for k = 1:numel(case_ids)
    idx = strcmp(metrics.CaseID, case_ids{k});
    case_pass = all(metrics.Pass(idx));
    pass_rows = append_pass(pass_rows, case_ids{k}, case_pass, ...
        join(string(metrics.Metric(idx)) + "=" + string(metrics.Value(idx)), "; "));
end

nan_count = sum(any(ismissing(metrics), 2)) + nnz(~isfinite(metrics.Value));
pass_rows = append_pass(pass_rows, 'P1-SUBGRID-001', nan_count == 0 && all(pass_rows.Pass), ...
    sprintf('nan_count=%d; all_subcases_pass=%d', nan_count, all(pass_rows.Pass)));

writetable(metrics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(mass_rows, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(pass_rows, fullfile(out_dir, 'Pass_Fail.csv'));

plot_storage_curve(T_incised, Resolution, channel_width, bankfull_depth, fig_dir);
plot_conveyance(face_depth, Q_model, Q_truth, fig_dir);
plot_full_momentum_equivalence(d_coarse, d_subgrid_flat, d_incised0_mm, d_incised1, fig_dir);

fprintf('Subgrid validation complete. Overall pass = %d\n', all(pass_rows.Pass));

%% Local helpers
function T = make_incised_tables(ny, nx, R, W, B, n, dz, maxDepth)
depth_axis = 0:dz:maxDepth;
nz = numel(depth_axis);
T.volume_cell = zeros(ny, nx, nz);
T.area_cell = zeros(ny, nx, nz);
T.area_x = zeros(ny, nx-1, nz);
T.area_y = zeros(ny-1, nx, nz);
T.n_x = NaN(ny, nx-1, nz);
T.n_y = NaN(ny-1, nx, nz);
T.width_x = zeros(ny, nx-1, nz);
T.width_y = zeros(ny-1, nx, nz);
T.wetfrac_x = zeros(ny, nx-1, nz);
T.wetfrac_y = zeros(ny-1, nx, nz);
T.perimeter_x = zeros(ny, nx-1, nz);
T.perimeter_y = zeros(ny-1, nx, nz);
T.Rh_x = zeros(ny, nx-1, nz);
T.Rh_y = zeros(ny-1, nx, nz);
T.K_x = zeros(ny, nx-1, nz);
T.K_y = zeros(ny-1, nx, nz);
for kk = 1:nz
    d = depth_axis(kk);
    T.volume_cell(:,:,kk) = rectangular_incised_volume(d, R, W, B);
    T.area_cell(:,:,kk) = rectangular_incised_plan_area(d, R, W, B);
    [Af,Wf,phif,Pf,Rhf,Kf,n_eff] = rectangular_incised_face_hydraulics(d, R, W, B, n);
    T.area_x(:,:,kk) = Af;
    T.area_y(:,:,kk) = Af;
    T.width_x(:,:,kk) = Wf;
    T.width_y(:,:,kk) = Wf;
    T.wetfrac_x(:,:,kk) = phif;
    T.wetfrac_y(:,:,kk) = phif;
    T.perimeter_x(:,:,kk) = Pf;
    T.perimeter_y(:,:,kk) = Pf;
    T.Rh_x(:,:,kk) = Rhf;
    T.Rh_y(:,:,kk) = Rhf;
    T.K_x(:,:,kk) = Kf;
    T.K_y(:,:,kk) = Kf;
    T.n_x(:,:,kk) = n_eff;
    T.n_y(:,:,kk) = n_eff;
end
T.invert_el = zeros(ny, nx);
T.invert_x = zeros(ny, nx-1);
T.invert_y = zeros(ny-1, nx);
T.dz = dz;
T.maxDepth = maxDepth;
end

function T = make_flat_tables(ny, nx, R, n, dz, maxDepth)
depth_axis = 0:dz:maxDepth;
nz = numel(depth_axis);
T.volume_cell = zeros(ny, nx, nz);
T.area_cell = zeros(ny, nx, nz);
T.area_x = zeros(ny, nx-1, nz);
T.area_y = zeros(ny-1, nx, nz);
T.n_x = n .* ones(ny, nx-1, nz);
T.n_y = n .* ones(ny-1, nx, nz);
T.width_x = zeros(ny, nx-1, nz);
T.width_y = zeros(ny-1, nx, nz);
T.wetfrac_x = zeros(ny, nx-1, nz);
T.wetfrac_y = zeros(ny-1, nx, nz);
T.perimeter_x = zeros(ny, nx-1, nz);
T.perimeter_y = zeros(ny-1, nx, nz);
T.Rh_x = zeros(ny, nx-1, nz);
T.Rh_y = zeros(ny-1, nx, nz);
T.K_x = zeros(ny, nx-1, nz);
T.K_y = zeros(ny-1, nx, nz);
for kk = 1:nz
    d = depth_axis(kk);
    T.volume_cell(:,:,kk) = R^2 * d;
    T.area_cell(:,:,kk) = R^2 * double(d > 0);
    T.area_x(:,:,kk) = R * d;
    T.area_y(:,:,kk) = R * d;
    width = R * double(d > 0);
    perimeter = (R + 2*d) * double(d > 0);
    Rh = (R*d) / max(perimeter, eps);
    K = (1/n) * (R*d) * Rh^(2/3) * double(d > 0);
    T.width_x(:,:,kk) = width;
    T.width_y(:,:,kk) = width;
    T.wetfrac_x(:,:,kk) = double(d > 0);
    T.wetfrac_y(:,:,kk) = double(d > 0);
    T.perimeter_x(:,:,kk) = perimeter;
    T.perimeter_y(:,:,kk) = perimeter;
    T.Rh_x(:,:,kk) = Rh;
    T.Rh_y(:,:,kk) = Rh;
    T.K_x(:,:,kk) = K;
    T.K_y(:,:,kk) = K;
end
T.invert_el = zeros(ny, nx);
T.invert_x = zeros(ny, nx-1);
T.invert_y = zeros(ny-1, nx);
T.dz = dz;
T.maxDepth = maxDepth;
end

function V = rectangular_incised_volume(d, R, W, B)
d = max(d, 0);
V = R .* W .* min(d, B) + R.^2 .* max(d - B, 0);
end

function A = rectangular_incised_plan_area(d, R, W, B)
A = zeros(size(d));
A(d > 0 & d <= B) = R .* W;
A(d > B) = R.^2;
end

function A = rectangular_incised_face_area(d, R, W, B)
d = max(d, 0);
A = W .* min(d, B) + R .* max(d - B, 0);
end

function [A,Wf,phi,P,Rh,K,n_eff] = rectangular_incised_face_hydraulics(d, R, W, B, n)
A = rectangular_incised_face_area(d, R, W, B);
if d <= 0
    Wf = 0;
    phi = 0;
    P = 0;
    Rh = 0;
    K = 0;
    n_eff = NaN;
elseif d <= B
    Wf = W;
    phi = W / R;
    P = W + 2*d;
    Rh = A / P;
    K = (1/n) * A * Rh^(2/3);
    H_G = A / R;
    denom = phi * d^(5/3) / n;
    n_eff = H_G^(5/3) / denom;
else
    Wf = R;
    phi = 1;
    P = R + 2*d;
    Rh = A / P;
    K = (1/n) * A * Rh^(2/3);
    d_channel = d;
    d_overbank = d - B;
    mean_d53_over_n = (W/R) * d_channel^(5/3) / n + ((R-W)/R) * d_overbank^(5/3) / n;
    H_G = A / R;
    n_eff = H_G^(5/3) / mean_d53_over_n;
end
end

function [d_t, outflow] = run_full_momentum_once(use_subgrid, T, d0_mm, R, n)
ny = size(d0_mm, 1);
nx = size(d0_mm, 2);
z = zeros(ny, nx);
d_p = d0_mm;
roughness = n .* ones(ny, nx);
roughness_squared = roughness.^2;
cell_area = R^2;
time_step_min = 0.002;
outlet_index = false(ny, nx);
outlet_type = 1;
slope_outlet = 0.001;
row_outlet = [];
col_outlet = [];
d_tolerance = 1e-6;
outflow0 = zeros(ny, nx, 5);
idx_nan = false(ny, nx);
flag_subgrid = double(use_subgrid);
if ~use_subgrid
    T = [];
end
[~,~,~,~,~,d_t,~,outflow] = Full_Momentum_Model_D4( ...
    1, [], [], [], [], [], [], [], [], [], [], [], [], ...
    0, z, d0_mm, d_p, roughness, roughness_squared, cell_area, ...
    time_step_min, R, outlet_index, outlet_type, slope_outlet, ...
    row_outlet, col_outlet, d_tolerance, outflow0, idx_nan, 0, ...
    flag_subgrid, [], [], [], [], [], [], [], [], [], ...
    [], 0, 0, T, struct());
end

function out = append_metric(in, case_id, metric, value, threshold, pass)
row = table(string(case_id), string(metric), value, threshold, logical(pass), ...
    'VariableNames', {'CaseID','Metric','Value','Threshold','Pass'});
out = [in; row];
end

function out = append_mass(in, case_id, initial_volume, final_volume, residual, residual_pct)
row = table(string(case_id), initial_volume, final_volume, residual, residual_pct, ...
    'VariableNames', {'CaseID','InitialVolume_m3','FinalVolume_m3','Residual_m3','Residual_pct'});
out = [in; row];
end

function out = append_pass(in, case_id, pass, notes)
row = table(string(case_id), logical(pass), string(notes), ...
    'VariableNames', {'CaseID','Pass','Notes'});
out = [in; row];
end

function e = rmse(a, b)
diffv = a(:) - b(:);
diffv = diffv(isfinite(diffv));
e = sqrt(mean(diffv.^2));
end

function v = rel_pct(err, denom)
v = 100 .* err ./ max(abs(denom), eps);
end

function val = nse(obs, sim)
obs = obs(:);
sim = sim(:);
ok = isfinite(obs) & isfinite(sim);
obs = obs(ok);
sim = sim(ok);
den = sum((obs - mean(obs)).^2);
if den <= eps
    val = double(max(abs(obs - sim)) < 1e-12);
else
    val = 1 - sum((obs - sim).^2) ./ den;
end
end

function nbad = count_negative_diff(V)
dV = diff(V, 1, 3);
nbad = nnz(dV(:) < -1e-12);
end

function plot_storage_curve(T, R, W, B, fig_dir)
depth = reshape(0:T.dz:T.maxDepth, 1, 1, []);
V = squeeze(T.volume_cell(1,1,:));
Vtruth = squeeze(rectangular_incised_volume(depth, R, W, B));
fig = figure('Visible','off');
plot(squeeze(depth), Vtruth, 'k-', 'LineWidth', 1.8); hold on;
plot(squeeze(depth), V, 'ro', 'MarkerSize', 3);
xlabel('Representative depth above invert (m)');
ylabel('Stored volume (m3)');
legend('Analytical', 'Lookup table', 'Location', 'northwest');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_storage_curve.png'));
close(fig);
end

function plot_conveyance(face_depth, Q_model, Q_truth, fig_dir)
fig = figure('Visible','off');
plot(face_depth(:), Q_truth(:), 'k.', 'MarkerSize', 16); hold on;
plot(face_depth(:), Q_model(:), 'rx', 'MarkerSize', 8);
xlabel('Face representative depth (m)');
ylabel('Discharge for S = 0.001 (m3/s)');
legend('Analytical', 'Lookup conveyance', 'Location', 'northwest');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_conveyance.png'));
close(fig);
end

function plot_full_momentum_equivalence(d_coarse, d_subgrid_flat, d0_incised, d1_incised, fig_dir)
fig = figure('Visible','off');
tiledlayout(1,2);
nexttile;
plot(d_coarse(2,:) ./ 1000, 'k-o', 'LineWidth', 1.4); hold on;
plot(d_subgrid_flat(2,:) ./ 1000, 'r--x', 'LineWidth', 1.4);
xlabel('Column');
ylabel('Depth (m)');
title('Flat-table equivalence');
legend('Coarse', 'Lookup subgrid', 'Location', 'best');
grid on;
nexttile;
plot(d0_incised(2,:) ./ 1000, 'k-o', 'LineWidth', 1.4); hold on;
plot(d1_incised(2,:) ./ 1000, 'b--x', 'LineWidth', 1.4);
xlabel('Column');
ylabel('Representative depth (m)');
title('Incised-channel diagnostic');
legend('Initial', 'After one step', 'Location', 'best');
grid on;
saveas(fig, fullfile(fig_dir, 'p1_subgrid_full_momentum.png'));
close(fig);
end
