%% ============================================================
% Evaluate HydroPol2D non-breaking wave benchmark
% Standalone script
%
% Requires existing workspace variables:
%   Paths
%   running_control
%   saver_memory_maps
%   DEM_raster
%
% Reads:
%   Paths.Temp/save_map_hydro_*.mat
%
% Uses:
%   Maps.Hydro.d [mm]
%
% Diagnostics included:
%   1. Analytical vs HydroPol2D centerline profiles
%   2. RMSE / MAE / max error / relative L2 error
%   3. Wetting-front error
%   4. Full saved-map stored-volume series
%   5. Analytical expected stored volume over active domain
%   6. Model - analytical storage residual
%   7. Optional outlet hydrograph integration diagnostic
%
% Notes:
%   For this benchmark, the most useful mass/volume diagnostic is:
%
%       V_model_stored(t) - V_analytical_stored(t)
%
%   because the analytical non-breaking wave profile defines the expected
%   water volume in the computational domain at each time.
%
% ============================================================

close all; clc;

%% ================= BENCHMARK PARAMETERS =================

% Grid / hydraulic parameters
dx_user = 25;          % [m], fallback if DEM_raster.cellsize is unavailable
n_manning = 0.02;    % [s/m^(1/3)]
u_wave = 1;       % [m/s]

% Times to compare [s]
compare_times_sec = [12 36 48 60]*60;
compare_times_min = compare_times_sec / 60;

% Depth threshold used for wetting-front detection [m]
depth_threshold = 1e-4;

% Coordinate convention.
% true  = finite-volume cell centers: x = (j - 0.5) dx
% false = original script convention: x = (j - 1) dx
use_cell_centers = true;

% Optional hydrograph diagnostic.
% This is only a diagnostic; the main benchmark compares stored volume
% against the analytical stored volume.
include_outlet_hydrograph_diagnostic = true;

%% ================= OUTPUT FOLDERS =================

outDir = fullfile(Paths.Results, 'NonBreakingWave_Validation');
figDir = fullfile(outDir, 'Figures');
tabDir = fullfile(outDir, 'Tables');

if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(tabDir), mkdir(tabDir); end

tempDir = Paths.Temp;

%% ================= MODEL TIME =================

t_records = gather(running_control.time_records(:));

if isdatetime(t_records)

    t_model_min = minutes(t_records - t_records(1));

elseif isduration(t_records)

    t_model_min = minutes(t_records);

else

    % HydroPol2D commonly stores elapsed time in minutes
    t_model_min = double(t_records(:));

end

t_model_min = double(t_model_min(:));
t_model_sec = t_model_min * 60;

%% ================= DOMAIN =================

ny = size(DEM_raster.Z,1);
nx = size(DEM_raster.Z,2);

if isfield(DEM_raster, 'cellsize') && ~isempty(DEM_raster.cellsize)
    dx = DEM_raster.cellsize;
else
    dx = dx_user;
end

A_cell = dx^2;  % [m2]

if use_cell_centers
    x_profile = ((0:nx-1)' + 0.5) * dx;  % [m]
else
    x_profile = ((0:nx-1)' * dx);        % [m]
end

row_mid = round(ny/2);

% Active computational cells
if exist('idx_nan','var') && ~isempty(idx_nan)

    active_mask = ~idx_nan;

else

    active_mask = isfinite(DEM_raster.Z);

end

active_mask = logical(active_mask);

% Number of active cells per x-column.
% This lets us compute analytical volume over the actual active domain,
% even if the mask is not perfectly rectangular.
active_cells_per_col = sum(active_mask, 1)';  % [nx x 1]

A_active_total = nnz(active_mask) * A_cell;

fprintf('\n============= NON-BREAKING WAVE DOMAIN =============\n');
fprintf('dx                         = %.6f m\n', dx);
fprintf('Cell area                  = %.6f m2\n', A_cell);
fprintf('Grid size                  = %d rows x %d cols\n', ny, nx);
fprintf('Active cells               = %d\n', nnz(active_mask));
fprintf('Active domain area         = %.6f m2\n', A_active_total);
fprintf('Coordinate convention      = %s\n', string_if(use_cell_centers, 'cell centers', 'left-edge/original'));
fprintf('====================================================\n\n');

%% ================= OPTIONAL OUTLET HYDROGRAPH DIAGNOSTIC =================

has_outlet_hydrograph = false;

if include_outlet_hydrograph_diagnostic && ...
        exist('outlet_states','var') && ...
        isfield(outlet_states,'outlet_hydrograph') && ...
        ~isempty(outlet_states.outlet_hydrograph) && ...
        exist('running_control','var') && ...
        isfield(running_control,'time_hydrograph') && ...
        ~isempty(running_control.time_hydrograph)

    q_out_model = gather(outlet_states.outlet_hydrograph(:));  % [m3/s]
    t_q = gather(running_control.time_hydrograph(:));

    if isdatetime(t_q)

        t_q_min = minutes(t_q - t_q(1));

    elseif isduration(t_q)

        t_q_min = minutes(t_q);

    else

        t_q_min = double(t_q(:));  % usually elapsed minutes

    end

    t_q_min = double(t_q_min(:));
    t_q_sec = t_q_min * 60;
    q_out_model = double(q_out_model(:));

    valid_q = isfinite(t_q_sec) & isfinite(q_out_model);

    t_q_sec = t_q_sec(valid_q);
    t_q_min = t_q_min(valid_q);
    q_out_model = q_out_model(valid_q);

    [t_q_sec, order_q] = sort(t_q_sec);
    t_q_min = t_q_min(order_q);
    q_out_model = q_out_model(order_q);

    [t_q_sec, ia_q] = unique(t_q_sec, 'stable');
    t_q_min = t_q_min(ia_q);
    q_out_model = q_out_model(ia_q);

    q_out_model = max(q_out_model, 0);

    if numel(t_q_sec) >= 2

        has_outlet_hydrograph = true;

        % Two interpretations:
        % 1. trapz: q is an instantaneous sample at time stamp.
        % 2. rect: q(i) represents the interval ending at t(i).
        V_out_cum_trapz = cumtrapz(t_q_sec, q_out_model);

        dt_q_sec = [0; diff(t_q_sec)];
        V_out_cum_rect = cumsum(q_out_model .* dt_q_sec);

    else

        has_outlet_hydrograph = false;
        V_out_cum_trapz = [];
        V_out_cum_rect = [];

    end

else

    t_q_min = [];
    t_q_sec = [];
    q_out_model = [];
    V_out_cum_trapz = [];
    V_out_cum_rect = [];

end

fprintf('\n============= OPTIONAL OUTLET HYDROGRAPH CHECK =============\n');
fprintf('has_outlet_hydrograph = %d\n', has_outlet_hydrograph);

if has_outlet_hydrograph

    fprintf('Hydrograph points          = %d\n', numel(t_q_sec));
    fprintf('t_q min/max                = %.6f / %.6f min\n', min(t_q_min), max(t_q_min));
    fprintf('q_out min/max/mean         = %.6e / %.6e / %.6e m3/s\n', ...
        min(q_out_model), max(q_out_model), mean(q_out_model));
    fprintf('Outlet volume trapz final  = %.12e m3\n', V_out_cum_trapz(end));
    fprintf('Outlet volume rect final   = %.12e m3\n', V_out_cum_rect(end));
    fprintf('rect - trapz final         = %.12e m3\n', V_out_cum_rect(end) - V_out_cum_trapz(end));

else

    fprintf('Outlet hydrograph unavailable or not used.\n');

end

fprintf('=============================================================\n\n');

%% ================= FULL SAVED-MAP VOLUME SERIES =================

nSaved = numel(t_model_min);

V_model_stored_series = NaN(nSaved,1);
V_analytical_stored_series = NaN(nSaved,1);
V_residual_series = NaN(nSaved,1);
V_error_pct_series = NaN(nSaved,1);
V_error_pct_initial_series = NaN(nSaved,1);

FrontAnalytical_series_m = NaN(nSaved,1);
FrontModel_series_m = NaN(nSaved,1);
FrontError_series_m = NaN(nSaved,1);

MaxDepthModel_series_m = NaN(nSaved,1);
MaxDepthAnalytical_series_m = NaN(nSaved,1);

% Initial analytical/model volumes for reference
h_ana_t0 = analytical_nonbreaking_depth_script(x_profile, 0, n_manning, u_wave);
V_analytical_initial = sum(h_ana_t0(:) .* active_cells_per_col(:), 'omitnan') * A_cell;

[D0_mm, ~, ~] = load_depth_map_from_temp_script(1, tempDir, saver_memory_maps);
D0_m = double(gather(D0_mm)) / 1000;
D0_m(~isfinite(D0_m)) = 0;
D0_m = max(D0_m, 0);

V_model_initial = sum(D0_m(active_mask), 'omitnan') * A_cell;

for jj = 1:nSaved

    t_sec = t_model_sec(jj);

    %% ---- Load model saved map ----
    [Dj_mm, ~, ~] = load_depth_map_from_temp_script(jj, tempDir, saver_memory_maps);

    Dj_m = double(gather(Dj_mm)) / 1000;
    Dj_m(~isfinite(Dj_m)) = 0;
    Dj_m = max(Dj_m, 0);

    %% ---- Model stored volume ----
    V_model_stored_series(jj) = sum(Dj_m(active_mask), 'omitnan') * A_cell;

    %% ---- Analytical full-domain stored volume ----
    h_ana_j = analytical_nonbreaking_depth_script(x_profile, t_sec, n_manning, u_wave);

    V_analytical_stored_series(jj) = ...
        sum(h_ana_j(:) .* active_cells_per_col(:), 'omitnan') * A_cell;

    %% ---- Residual: model minus analytical ----
    V_residual_series(jj) = V_model_stored_series(jj) - V_analytical_stored_series(jj);

    V_error_pct_series(jj) = ...
        100 * V_residual_series(jj) / max(V_analytical_stored_series(jj), eps);

    V_error_pct_initial_series(jj) = ...
        100 * V_residual_series(jj) / max(V_analytical_initial, eps);

    %% ---- Front diagnostics from centerline ----
    h_model_mid = Dj_m(row_mid,:)';

    threshold = depth_threshold;

    front_ana_idx = find(h_ana_j > threshold, 1, 'last');
    front_mod_idx = find(h_model_mid > threshold, 1, 'last');

    if isempty(front_ana_idx)
        front_ana = NaN;
    else
        front_ana = x_profile(front_ana_idx);
    end

    if isempty(front_mod_idx)
        front_mod = NaN;
    else
        front_mod = x_profile(front_mod_idx);
    end

    FrontAnalytical_series_m(jj) = front_ana;
    FrontModel_series_m(jj) = front_mod;
    FrontError_series_m(jj) = front_mod - front_ana;

    MaxDepthModel_series_m(jj) = max(Dj_m(:), [], 'omitnan');
    MaxDepthAnalytical_series_m(jj) = max(h_ana_j(:), [], 'omitnan');

end

FullVolumeSeries = table( ...
    t_model_min, ...
    t_model_sec, ...
    V_model_stored_series, ...
    V_analytical_stored_series, ...
    V_residual_series, ...
    V_error_pct_series, ...
    V_error_pct_initial_series, ...
    FrontAnalytical_series_m, ...
    FrontModel_series_m, ...
    FrontError_series_m, ...
    MaxDepthModel_series_m, ...
    MaxDepthAnalytical_series_m, ...
    'VariableNames', { ...
        'Time_min', ...
        'Time_s', ...
        'ModelStoredVolume_m3', ...
        'AnalyticalStoredVolume_m3', ...
        'Residual_ModelMinusAnalytical_m3', ...
        'Residual_pct_of_AnalyticalVolume', ...
        'Residual_pct_of_InitialAnalyticalVolume', ...
        'FrontAnalytical_m', ...
        'FrontModel_m', ...
        'FrontError_m', ...
        'MaxDepthModel_m', ...
        'MaxDepthAnalytical_m'});

writetable(FullVolumeSeries, fullfile(tabDir, 'NonBreakingWave_Full_VolumeSeries.csv'));

fprintf('\n============= FULL SAVED-MAP VOLUME SUMMARY =============\n');
fprintf('Analytical initial volume     = %.12e m3\n', V_analytical_initial);
fprintf('Model initial volume          = %.12e m3\n', V_model_initial);
fprintf('Initial model - analytical    = %.12e m3\n', V_model_initial - V_analytical_initial);
fprintf('Initial bias                  = %.12e %% of analytical initial\n', ...
    100 * (V_model_initial - V_analytical_initial) / max(V_analytical_initial, eps));
fprintf('Max abs volume residual       = %.12e m3\n', max(abs(V_residual_series), [], 'omitnan'));
fprintf('Max abs residual %% analytical= %.12e %%\n', max(abs(V_error_pct_series), [], 'omitnan'));
fprintf('Final volume residual         = %.12e m3\n', V_residual_series(end));
fprintf('Final residual %% analytical  = %.12e %%\n', V_error_pct_series(end));
fprintf('==========================================================\n\n');

%% ================= STORAGE FOR PROFILE COMPARISON =================

nCompare = numel(compare_times_min);

Profiles = table();
Profiles.x_m = x_profile;

Metrics = table('Size',[nCompare 13], ...
    'VariableTypes', { ...
        'double','double','double','double','double','double','double', ...
        'double','double','double','double','double','double'}, ...
    'VariableNames', { ...
        'TargetTime_min', ...
        'SavedTime_min', ...
        'RMSE_m', ...
        'MAE_m', ...
        'MaxAbsError_m', ...
        'L2Relative', ...
        'FrontAnalytical_m', ...
        'FrontModel_m', ...
        'FrontError_m', ...
        'ModelStoredVolume_m3', ...
        'AnalyticalStoredVolume_m3', ...
        'VolumeResidual_m3', ...
        'VolumeResidual_pct'});

%% ================= MAIN PROFILE FIGURE =================

figure('Color','w','Position',[100 100 1200 650]);
hold on;

colors = lines(nCompare);

for k = 1:nCompare

    target_min = compare_times_min(k);
    target_sec = compare_times_sec(k);

    [~, idx_time] = min(abs(t_model_min - target_min));
    saved_min = t_model_min(idx_time);

    %% ---- Load modeled depth map ----
    [D_mm, local_i, store] = load_depth_map_from_temp_script( ...
        idx_time, tempDir, saver_memory_maps);

    D_m = double(gather(D_mm)) / 1000;
    D_m(~isfinite(D_m)) = 0;
    D_m = max(D_m, 0);

    h_model = D_m(row_mid,:)';  % [m]
    h_model(~isfinite(h_model)) = NaN;

    %% ---- Analytical profile ----
    h_ana = analytical_nonbreaking_depth_script( ...
        x_profile, target_sec, n_manning, u_wave);

    %% ---- Metrics ----
    valid = isfinite(h_model) & isfinite(h_ana);

    err = h_model(valid) - h_ana(valid);

    RMSE = sqrt(mean(err.^2));
    MAE = mean(abs(err));
    MaxAbsError = max(abs(err));
    L2Relative = sqrt(sum(err.^2)) / max(sqrt(sum(h_ana(valid).^2)), eps);

    threshold = depth_threshold;

    front_ana_idx = find(h_ana > threshold, 1, 'last');
    front_mod_idx = find(h_model > threshold, 1, 'last');

    if isempty(front_ana_idx)
        front_ana = NaN;
    else
        front_ana = x_profile(front_ana_idx);
    end

    if isempty(front_mod_idx)
        front_mod = NaN;
    else
        front_mod = x_profile(front_mod_idx);
    end

    front_error = front_mod - front_ana;

    %% ---- Volume diagnostics at comparison time ----
    V_model_compare = sum(D_m(active_mask), 'omitnan') * A_cell;

    V_ana_compare = sum(h_ana(:) .* active_cells_per_col(:), 'omitnan') * A_cell;

    V_res_compare = V_model_compare - V_ana_compare;
    V_res_pct_compare = 100 * V_res_compare / max(V_ana_compare, eps);

    %% ---- Store results ----
    Metrics.TargetTime_min(k) = target_min;
    Metrics.SavedTime_min(k) = saved_min;
    Metrics.RMSE_m(k) = RMSE;
    Metrics.MAE_m(k) = MAE;
    Metrics.MaxAbsError_m(k) = MaxAbsError;
    Metrics.L2Relative(k) = L2Relative;
    Metrics.FrontAnalytical_m(k) = front_ana;
    Metrics.FrontModel_m(k) = front_mod;
    Metrics.FrontError_m(k) = front_error;
    Metrics.ModelStoredVolume_m3(k) = V_model_compare;
    Metrics.AnalyticalStoredVolume_m3(k) = V_ana_compare;
    Metrics.VolumeResidual_m3(k) = V_res_compare;
    Metrics.VolumeResidual_pct(k) = V_res_pct_compare;

    Profiles.(sprintf('Analytical_h_t%05ds_m', target_sec)) = h_ana;
    Profiles.(sprintf('HydroPol2D_h_t%05ds_m', target_sec)) = h_model;

    %% ---- Plot ----
    plot(x_profile, h_ana, '-', ...
        'Color', colors(k,:), ...
        'LineWidth', 2.2);

    plot(x_profile, h_model, '--', ...
        'Color', colors(k,:), ...
        'LineWidth', 2.2);

    fprintf(['t = %.1f min | saved = %.1f min | store = %d | local_i = %d | ', ...
             'RMSE = %.5f m | front error = %.2f m | ', ...
             'Vmodel - Vana = %.6e m3 (%.6e %%)\n'], ...
        target_min, saved_min, store, local_i, RMSE, front_error, ...
        V_res_compare, V_res_pct_compare);

end

%% ================= FORMAT PROFILE FIGURE =================

grid on; box on;

xlabel('x [m]');
ylabel('Water depth h [m]');
title('Non-Breaking Wave: Analytical vs HydroPol2D Centerline Profiles');

legend_entries = strings(2*nCompare,1);
ii = 1;

for k = 1:nCompare
    legend_entries(ii) = sprintf('Analytical t = %.0f min', compare_times_min(k)); ii = ii + 1;
    legend_entries(ii) = sprintf('HydroPol2D t = %.0f min', compare_times_min(k)); ii = ii + 1;
end

legend(legend_entries, 'Location','northeastoutside');

set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5,'TickDir','out');

exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Profile_Comparison.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Profile_Comparison.png'), 'Resolution',300);

%% ================= ERROR METRICS FIGURE =================

figure('Color','w','Position',[100 100 1200 450]);

subplot(1,4,1)
plot(Metrics.TargetTime_min, Metrics.RMSE_m, 'o-', 'LineWidth',2);
grid on; box on;
xlabel('Time [min]');
ylabel('RMSE [m]');
title('RMSE');

subplot(1,4,2)
plot(Metrics.TargetTime_min, Metrics.L2Relative, 's-', 'LineWidth',2);
grid on; box on;
xlabel('Time [min]');
ylabel('Relative L2 error [-]');
title('Relative L2 error');

subplot(1,4,3)
plot(Metrics.TargetTime_min, Metrics.FrontError_m, 'd-', 'LineWidth',2);
grid on; box on;
xlabel('Time [min]');
ylabel('Front error [m]');
title('Wetting-front error');

subplot(1,4,4)
plot(Metrics.TargetTime_min, Metrics.VolumeResidual_pct, '^-', 'LineWidth',2);
grid on; box on;
xlabel('Time [min]');
ylabel('Volume residual [%]');
title('Storage residual');
ytickformat('%.2e');

exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Error_Metrics.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Error_Metrics.png'), 'Resolution',300);

%% ================= FULL VOLUME DIAGNOSTIC FIGURES =================

figure('Color','w','Position',[150 150 1200 500]);

subplot(1,2,1)

plot(t_model_min, V_analytical_stored_series, 'k--', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Analytical stored volume');
hold on;

plot(t_model_min, V_model_stored_series, 'b-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 4, ...
    'DisplayName', 'Model stored volume');

grid on; box on;
xlabel('Time [min]');
ylabel('Stored volume [m^3]');
title('Non-breaking wave stored volume');
legend('Location','best');

subplot(1,2,2)

plot(t_model_min, V_residual_series, 'r-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 4);

grid on; box on;
xlabel('Time [min]');
ylabel('Model - analytical volume [m^3]');
title('Stored-volume residual');

exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Volume_Diagnostic.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Volume_Diagnostic.png'), 'Resolution',300);

figure('Color','w','Position',[200 200 1200 500]);

subplot(1,2,1)

plot(t_model_min, V_error_pct_series, 'm-o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 4);

grid on; box on;
xlabel('Time [min]');
ylabel('Residual [% of analytical volume]');
title(sprintf('Storage residual, max = %.3e %%', ...
    max(abs(V_error_pct_series), [], 'omitnan')));
ytickformat('%.2e');

subplot(1,2,2)

plot(t_model_min, FrontError_series_m, 'c-s', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 4);

grid on; box on;
xlabel('Time [min]');
ylabel('Front error [m]');
title('Wetting-front error over all saved maps');

exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Residual_and_Front_Series.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_Residual_and_Front_Series.png'), 'Resolution',300);

%% ================= OPTIONAL OUTLET HYDROGRAPH FIGURES =================

if has_outlet_hydrograph

    figure('Color','w','Position',[250 250 1000 500]);

    plot(t_q_min, V_out_cum_trapz, 'r--', ...
        'LineWidth', 2, ...
        'DisplayName', 'Outlet volume trapz');
    hold on;

    plot(t_q_min, V_out_cum_rect, 'b-', ...
        'LineWidth', 2, ...
        'DisplayName', 'Outlet volume rectangular');

    plot(t_q_min, V_out_cum_rect - V_out_cum_trapz, 'k:', ...
        'LineWidth', 2, ...
        'DisplayName', 'Rect - trapz');

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Cumulative outlet volume [m^3]');
    title('Outlet Volume Integration Method Check');
    legend('Location','best');

    exportgraphics(gcf, fullfile(figDir,'NonBreakingWave_OutletVolume_IntegrationCheck.png'), 'Resolution',300);

    OutletHydrographDiagnostic = table( ...
        t_q_min, ...
        t_q_sec, ...
        q_out_model, ...
        V_out_cum_trapz, ...
        V_out_cum_rect, ...
        V_out_cum_rect - V_out_cum_trapz, ...
        'VariableNames', { ...
            'Time_min', ...
            'Time_s', ...
            'OutletDischarge_m3s', ...
            'OutletVolume_trapz_m3', ...
            'OutletVolume_rect_m3', ...
            'RectMinusTrapz_m3'});

    writetable(OutletHydrographDiagnostic, ...
        fullfile(tabDir,'NonBreakingWave_OutletHydrograph_Diagnostic.csv'));

end

%% ================= EXPORT TABLES =================

writetable(Profiles, fullfile(tabDir,'NonBreakingWave_Profile_Comparison.csv'));
writetable(Metrics,  fullfile(tabDir,'NonBreakingWave_Error_Metrics.csv'));
writetable(FullVolumeSeries, fullfile(tabDir,'NonBreakingWave_Full_VolumeSeries.csv'));

fprintf('\nNon-breaking wave validation exported to:\n%s\n', outDir);

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function h = analytical_nonbreaking_depth_script(x, t, n, u)

    arg = -(7/3) * n^2 * u^2 .* (x - u*t);

    h = zeros(size(x));

    wet = arg > 0;
    h(wet) = arg(wet).^(3/7);

    h(~isfinite(h)) = 0;
    h(h < 0) = 0;

end

function [D_mm, local_i, store] = load_depth_map_from_temp_script(global_i, tempDir, saver_memory_maps)

    store = ceil(global_i / saver_memory_maps);
    local_i = global_i - (store-1)*saver_memory_maps;

    f1 = fullfile(tempDir, sprintf('save_map_hydro_%d.mat', store));
    f2 = fullfile(tempDir, sprintf('save_map_hydro_%d', store));

    if isfile(f1)

        mapFile = f1;

    elseif isfile(f2)

        mapFile = f2;

    else

        error('Could not find temporary map file:\n%s\nor\n%s', f1, f2);

    end

    S = load(mapFile, 'Maps');

    if ~isfield(S,'Maps') || ~isfield(S.Maps,'Hydro') || ~isfield(S.Maps.Hydro,'d')
        error('File does not contain Maps.Hydro.d:\n%s', mapFile);
    end

    D_mm = S.Maps.Hydro.d(:,:,local_i);

end

function out = string_if(cond, a, b)

    if cond
        out = a;
    else
        out = b;
    end

end