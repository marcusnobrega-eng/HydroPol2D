%% ============================================================
% RITTER DAM-BREAK: Analytical vs HydroPol2D depth profiles
% + full-domain mass-balance diagnostic
%
% Reads HydroPol2D temporary map files:
%   save_map_hydro_1, save_map_hydro_2, ...
%
% Uses:
%   Maps.Hydro.d [mm]
%
% Mass balance:
%   For Ritter dam-break there is no rainfall input.
%
%   If the domain is closed:
%       V_initial = V_stored(t)
%
%   If the domain has outlet/open-boundary flow:
%       V_initial = V_stored(t) + V_outlet(t)
%
% This script reports:
%   1) Storage-only error
%   2) Stored + explicit outlet error
%   3) Full saved-map volume series
%
% IMPORTANT:
%   V_initial is taken from the first saved HydroPol2D map, so the mass
%   conservation check is based on the actual model initial water volume.
%
%   The analysis is restricted to analysis_end_sec to avoid including
%   post-run/reset/empty maps, such as the problematic 60 s map.
% ============================================================

close all; clc;

%% ================= USER SETTINGS =================

% Case/static folder where Ritter_Analytical_Profiles.csv was saved
staticDir = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/RitterDam/Static';

% Analytical parameters
h0    = 1.0;      % [m]
x_dam = 250;      % [m]
g     = 9.81;     % [m/s2]

% Times to compare [s]
compare_times_sec = [10 20 30 40 50];

% Analyze only saved maps and comparison times up to this time [s]
analysis_end_sec = 50;

% If your run saved maps every 10 s, this will match exactly.
% If not, the script finds the nearest saved map.
use_nearest_saved_time = true;

% Depth threshold used to locate wet/dry front [m]
depth_threshold = 1e-4;

%% ================= PATHS =================

tempDir = Paths.Temp;
outDir  = fullfile(Paths.Results, 'Ritter_Validation');

if ~isfolder(outDir)
    mkdir(outDir);
end

figDir = fullfile(outDir, 'Figures');
tabDir = fullfile(outDir, 'Tables');

if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(tabDir), mkdir(tabDir); end

%% ================= MODEL TIMES FOR SAVED MAPS =================

t_records = gather(running_control.time_records(:));

% Convert model saved-map times to seconds
if isdatetime(t_records)

    t_model_sec = seconds(t_records - t_records(1));

elseif isduration(t_records)

    t_model_sec = seconds(t_records);

else

    % HydroPol2D usually stores elapsed time in minutes when flag_elapsed_time = 1
    if exist('flags','var') && isfield(flags,'flag_elapsed_time') && flags.flag_elapsed_time == 1
        t_model_sec = double(t_records) * 60;  % min -> s
    else
        % Conservative fallback: assume minutes
        t_model_sec = double(t_records) * 60;
    end

end

t_model_sec = double(t_model_sec(:));

%% ================= ANALYSIS TIME WINDOW =================

valid_saved_window = t_model_sec <= analysis_end_sec;

if ~any(valid_saved_window)
    error('No saved maps found at or before %.3f s.', analysis_end_sec);
end

% Restrict comparison times to the analysis window
compare_times_sec = compare_times_sec(compare_times_sec <= analysis_end_sec);

if isempty(compare_times_sec)
    error('No compare_times_sec values are <= analysis_end_sec.');
end

fprintf('\n============= ANALYSIS WINDOW =============\n');
fprintf('Analysis end time = %.3f s\n', analysis_end_sec);
fprintf('Saved maps used   = %d / %d\n', nnz(valid_saved_window), numel(t_model_sec));
fprintf('Compare times [s] = ');
fprintf('%.3f ', compare_times_sec);
fprintf('\n===========================================\n\n');

%% ================= OUTLET HYDROGRAPH FOR MASS BALANCE =================

has_outlet_hydrograph = exist('outlet_states','var') && ...
    isfield(outlet_states,'outlet_hydrograph') && ...
    ~isempty(outlet_states.outlet_hydrograph);

if has_outlet_hydrograph

    q_out_model = gather(outlet_states.outlet_hydrograph(:)); % [m3/s]

    % Use running_control.time_hydrograph if available
    if isfield(running_control,'time_hydrograph') && ~isempty(running_control.time_hydrograph)

        t_q = gather(running_control.time_hydrograph(:));

        if isdatetime(t_q)

            t_q_sec = seconds(t_q - t_q(1));

        elseif isduration(t_q)

            t_q_sec = seconds(t_q);

        else

            % HydroPol2D usually stores elapsed time in minutes
            t_q_sec = double(t_q) * 60;

        end

    else

        warning('running_control.time_hydrograph not found. Outlet volume will not be computed.');
        has_outlet_hydrograph = false;

    end

    if has_outlet_hydrograph

        t_q_sec = double(t_q_sec(:));
        q_out_model = double(q_out_model(:));

        valid_q = isfinite(t_q_sec) & isfinite(q_out_model);

        t_q_sec = t_q_sec(valid_q);
        q_out_model = q_out_model(valid_q);

        [t_q_sec, order_q] = sort(t_q_sec);
        q_out_model = q_out_model(order_q);

        [t_q_sec, ia_q] = unique(t_q_sec, 'stable');
        q_out_model = q_out_model(ia_q);

        q_out_model = max(q_out_model, 0);

        if numel(t_q_sec) >= 2
            cumV_out_model = cumtrapz(t_q_sec, q_out_model); % [m3]
        else
            warning('Outlet hydrograph has fewer than 2 valid points. Outlet volume assumed zero.');
            has_outlet_hydrograph = false;
            cumV_out_model = [];
        end

    end

else

    warning('outlet_states.outlet_hydrograph not found. Outlet volume will be assumed zero.');
    t_q_sec = [];
    q_out_model = [];
    cumV_out_model = [];

end

%% ================= OUTLET HYDROGRAPH DEBUG =================

fprintf('\n============= OUTLET HYDROGRAPH DEBUG =============\n');
fprintf('has_outlet_hydrograph = %d\n', has_outlet_hydrograph);

if has_outlet_hydrograph

    fprintf('Number hydrograph points = %d\n', numel(t_q_sec));
    fprintf('t_q_sec min/max          = %.6f / %.6f s\n', min(t_q_sec), max(t_q_sec));
    fprintf('q_out min/max/mean       = %.6e / %.6e / %.6e m3/s\n', ...
        min(q_out_model), max(q_out_model), mean(q_out_model));
    fprintf('cum outlet volume final  = %.12f m3\n', cumV_out_model(end));

else

    fprintf('Outlet hydrograph not available or not usable.\n');

end

fprintf('====================================================\n\n');

%% ================= DOMAIN GEOMETRY =================

% Model grid
ny = size(DEM_raster.Z, 1);
nx = size(DEM_raster.Z, 2);

dx = DEM_raster.cellsize;        % [m]
A_cell = dx^2;                   % [m2]

% Finite-volume cell-center coordinates
x_profile = ((0:nx-1)' + 0.5) * dx;    % [m]

% Centerline row
row_mid = round(ny/2);

%% ================= MASS BALANCE SETUP =================

% Active computational domain
if exist('idx_nan','var') && ~isempty(idx_nan)

    active_mask = ~idx_nan;

else

    active_mask = isfinite(DEM_raster.Z);

end

domain_area = nnz(active_mask) * A_cell;

% Geometric initial reservoir region.
% This assumes dam is at a cell face x = x_dam and cells are represented by centers.
initial_wet_cols = x_profile < x_dam;

initial_wet_mask = false(ny,nx);
initial_wet_mask(:,initial_wet_cols) = true;
initial_wet_mask = initial_wet_mask & active_mask;

% Geometric initial reservoir volume [m3]
V_initial_geom = h0 * nnz(initial_wet_mask) * A_cell;

% Model-based initial volume from first saved HydroPol2D map.
% This is the reference used for numerical mass conservation.
[D0_mm, ~, ~] = load_hydropol_depth_map(1, tempDir, saver_memory_maps);

D0_m = double(D0_mm) / 1000;
D0_m(~isfinite(D0_m)) = 0;
D0_m = max(D0_m, 0);

V_initial_model = sum(D0_m(active_mask), 'omitnan') * A_cell;

% Use model initial volume for conservation error.
V_initial = V_initial_model;

setup_volume_bias_pct = 100 * (V_initial_model - V_initial_geom) / max(V_initial_geom, eps);

fprintf('\n============= RITTER MASS SETUP =============\n');
fprintf('Cell size dx                  = %.6f m\n', dx);
fprintf('Cell area                     = %.6f m2\n', A_cell);
fprintf('Domain area                   = %.6f m2\n', domain_area);
fprintf('Initial wet column count geom = %d\n', nnz(initial_wet_cols));
fprintf('Grid wet length geom          = %.6f m\n', nnz(initial_wet_cols) * dx);
fprintf('Geometric initial volume      = %.12f m3\n', V_initial_geom);
fprintf('Model-map initial volume      = %.12f m3\n', V_initial_model);
fprintf('Setup volume bias             = %.6f %%\n', setup_volume_bias_pct);

if has_outlet_hydrograph
    fprintf('Outlet hydrograph volume      = included\n');
else
    fprintf('Outlet hydrograph volume      = not included / assumed zero\n');
end

fprintf('=============================================\n\n');

%% ================= DOMAIN / FRONT CHECK =================

domain_xmax = max(x_profile);

fprintf('\n============= RITTER FRONT CHECK =============\n');
fprintf('Domain xmax = %.3f m\n', domain_xmax);

for ii_check = 1:numel(compare_times_sec)
    t_check = compare_times_sec(ii_check);
    front_check = x_dam + 2*sqrt(g*h0)*t_check;
    fprintf('Ritter front at %.2f s = %.3f m | distance to domain end = %.3f m\n', ...
        t_check, front_check, domain_xmax - front_check);
end

fprintf('==============================================\n\n');

%% ================= LOAD ANALYTICAL CSV IF AVAILABLE =================

analyticalFile = fullfile(staticDir, 'Ritter_Analytical_Profiles.csv');

if isfile(analyticalFile)

    Analytical = readtable(analyticalFile);

else

    warning('Analytical CSV not found. Analytical profiles will be computed internally.');
    Analytical = table();
    Analytical.x_m = x_profile;

end

%% ================= RITTER ANALYTICAL FUNCTION =================

c0 = sqrt(g*h0);
x_rel = x_profile - x_dam;

compute_ritter_h = @(t) local_ritter_depth(x_rel, t, h0, c0, g);

%% ================= FULL SAVED-MAP MASS SERIES =================

saved_ids = find(valid_saved_window);
t_saved_analysis = t_model_sec(saved_ids);

nSaved = numel(saved_ids);

V_stored_series = NaN(nSaved,1);
V_outlet_series = NaN(nSaved,1);
V_total_series  = NaN(nSaved,1);
StorageOnlyError_series_pct = NaN(nSaved,1);
TotalMassError_series_pct = NaN(nSaved,1);

for jj_local = 1:nSaved

    jj = saved_ids(jj_local);

    [Dj_mm, ~, ~] = load_hydropol_depth_map(jj, tempDir, saver_memory_maps);

    Dj_m = double(Dj_mm) / 1000;
    Dj_m(~isfinite(Dj_m)) = 0;
    Dj_m = max(Dj_m, 0);

    V_stored_series(jj_local) = sum(Dj_m(active_mask), 'omitnan') * A_cell;

    if has_outlet_hydrograph && numel(t_q_sec) >= 2

        if t_model_sec(jj) <= max(t_q_sec)
            V_outlet_series(jj_local) = interp1(t_q_sec, cumV_out_model, t_model_sec(jj), 'linear', 'extrap');
        else
            V_outlet_series(jj_local) = cumV_out_model(end);
        end

        V_outlet_series(jj_local) = max(V_outlet_series(jj_local), 0);

    else

        V_outlet_series(jj_local) = 0;

    end

    V_total_series(jj_local) = V_stored_series(jj_local) + V_outlet_series(jj_local);

    StorageOnlyError_series_pct(jj_local) = ...
        100 * (V_stored_series(jj_local) - V_initial) / max(V_initial, eps);

    TotalMassError_series_pct(jj_local) = ...
        100 * (V_total_series(jj_local) - V_initial) / max(V_initial, eps);

    fprintf('Saved map %d | t = %.3f s | max depth = %.6e m | stored volume = %.6f m3 | total mass error = %.6f %%\n', ...
        jj, t_model_sec(jj), max(Dj_m(:),[],'omitnan'), V_stored_series(jj_local), TotalMassError_series_pct(jj_local));

end

FullMassSeries = table( ...
    t_saved_analysis, ...
    V_stored_series, ...
    V_outlet_series, ...
    V_total_series, ...
    StorageOnlyError_series_pct, ...
    TotalMassError_series_pct, ...
    'VariableNames', { ...
        'Time_s', ...
        'StoredVolume_m3', ...
        'OutletVolume_m3', ...
        'StoredPlusOutletVolume_m3', ...
        'StorageOnlyError_pct', ...
        'TotalMassError_pct'});

writetable(FullMassSeries, fullfile(tabDir, 'Ritter_Full_SavedMap_MassSeries.csv'));

%% ================= EXTRACT MODEL PROFILES =================

nCompare = numel(compare_times_sec);

Profiles = table();
Profiles.x_m = x_profile;

Metrics = table('Size',[nCompare 15], ...
    'VariableTypes', { ...
        'double','double','double','double','double','double', ...
        'double','double','double','double','double','double', ...
        'double','double','double'}, ...
    'VariableNames', { ...
        'TargetTime_s', ...
        'SavedTime_s', ...
        'RMSE_m', ...
        'MAE_m', ...
        'MaxAbsError_m', ...
        'L2Relative', ...
        'FrontAnalytical_m', ...
        'FrontModel_m', ...
        'InitialVolume_m3', ...
        'ModelStoredVolume_m3', ...
        'OutletVolume_m3', ...
        'TotalAccountedVolume_m3', ...
        'StorageOnlyError_pct', ...
        'TotalMassError_m3', ...
        'TotalMassError_pct'});

figure('Color','w','Position',[100 100 1200 750]);
hold on;

colors = lines(nCompare);

for kk = 1:nCompare

    target_t = compare_times_sec(kk);

    if use_nearest_saved_time

        valid_candidate_ids = find(t_model_sec <= analysis_end_sec);
        [~, local_idx_time] = min(abs(t_model_sec(valid_candidate_ids) - target_t));
        idx_time = valid_candidate_ids(local_idx_time);

    else

        idx_time = find(abs(t_model_sec - target_t) < 1e-9 & t_model_sec <= analysis_end_sec, 1, 'first');

        if isempty(idx_time)
            error('No exact saved map found for t = %.3f s within analysis window.', target_t);
        end

    end

    saved_t = t_model_sec(idx_time);

    %% ---- Load HydroPol2D map from temporary chunks ----
    [D_mm, local_i, store] = load_hydropol_depth_map(idx_time, tempDir, saver_memory_maps);

    % Full-domain depth map [m]
    D_full_m = double(D_mm) / 1000;
    D_full_m(~isfinite(D_full_m)) = 0;
    D_full_m = max(D_full_m, 0);

    %% ---- Full-domain model stored volume ----
    V_model_stored = sum(D_full_m(active_mask), 'omitnan') * A_cell;  % [m3]

    %% ---- Outlet volume up to saved_t ----
    if has_outlet_hydrograph && numel(t_q_sec) >= 2

        if saved_t <= max(t_q_sec)
            V_outlet_to_saved = interp1(t_q_sec, cumV_out_model, saved_t, 'linear', 'extrap');
        else
            V_outlet_to_saved = cumV_out_model(end);
        end

        V_outlet_to_saved = max(V_outlet_to_saved, 0);

    else

        V_outlet_to_saved = 0;

    end

    %% ---- Mass balance diagnostics ----
    StorageOnlyError_pct = 100 * (V_model_stored - V_initial) / max(V_initial, eps);

    V_total_accounted = V_model_stored + V_outlet_to_saved;

    TotalMassError_m3 = V_total_accounted - V_initial;
    TotalMassError_pct = 100 * TotalMassError_m3 / max(V_initial, eps);

    %% ---- Centerline model profile ----
    h_model = D_full_m(row_mid,:)';  % [m]
    h_model(~isfinite(h_model)) = NaN;

    %% ---- Analytical profile ----
    colName = sprintf('h_t%03ds_m', round(target_t));

    if ismember(colName, Analytical.Properties.VariableNames)

        h_analytical = Analytical.(colName);
        h_analytical = h_analytical(:);

    else

        h_analytical = compute_ritter_h(target_t);

    end

    %% ---- Ensure same length ----
    nMin = min(numel(h_model), numel(h_analytical));

    x_use = x_profile(1:nMin);
    h_model = h_model(1:nMin);
    h_analytical = h_analytical(1:nMin);

    %% ---- Metrics ----
    valid = isfinite(h_model) & isfinite(h_analytical);

    err = h_model(valid) - h_analytical(valid);

    RMSE = sqrt(mean(err.^2));
    MAE  = mean(abs(err));
    MaxAbsError = max(abs(err));
    L2Relative = sqrt(sum(err.^2)) / max(sqrt(sum(h_analytical(valid).^2)), eps);

    % Analytical front location
    front_analytical = x_dam + 2*c0*target_t;

    % Model front location: last x where depth exceeds threshold
    wet_idx = find(h_model > depth_threshold);

    if isempty(wet_idx)
        front_model = NaN;
    else
        front_model = x_use(max(wet_idx));
    end

    %% ---- Store metrics ----
    Metrics.TargetTime_s(kk)            = target_t;
    Metrics.SavedTime_s(kk)             = saved_t;
    Metrics.RMSE_m(kk)                  = RMSE;
    Metrics.MAE_m(kk)                   = MAE;
    Metrics.MaxAbsError_m(kk)           = MaxAbsError;
    Metrics.L2Relative(kk)              = L2Relative;
    Metrics.FrontAnalytical_m(kk)       = front_analytical;
    Metrics.FrontModel_m(kk)            = front_model;
    Metrics.InitialVolume_m3(kk)        = V_initial;
    Metrics.ModelStoredVolume_m3(kk)    = V_model_stored;
    Metrics.OutletVolume_m3(kk)         = V_outlet_to_saved;
    Metrics.TotalAccountedVolume_m3(kk) = V_total_accounted;
    Metrics.StorageOnlyError_pct(kk)    = StorageOnlyError_pct;
    Metrics.TotalMassError_m3(kk)       = TotalMassError_m3;
    Metrics.TotalMassError_pct(kk)      = TotalMassError_pct;

    %% ---- Store profile table ----
    Profiles.(sprintf('Analytical_h_t%03ds_m', round(target_t))) = h_analytical;
    Profiles.(sprintf('Model_h_t%03ds_m', round(target_t)))      = h_model;

    %% ---- Plot ----
    plot(x_use, h_analytical, '-', ...
        'Color', colors(kk,:), ...
        'LineWidth', 2);

    plot(x_use, h_model, '--', ...
        'Color', colors(kk,:), ...
        'LineWidth', 2);

    fprintf(['t target = %6.2f s | saved = %6.2f s | store = %d | local_i = %d | ', ...
             'RMSE = %.5f m | front error = %.2f m | ', ...
             'Vstored = %.6f m3 | Vout = %.6f m3 | total mass error = %.4f %%\n'], ...
        target_t, saved_t, store, local_i, RMSE, front_model - front_analytical, ...
        V_model_stored, V_outlet_to_saved, TotalMassError_pct);

end

%% ================= FIGURE FORMAT =================

grid on; box on;
xlabel('x [m]');
ylabel('Water depth h [m]');
title(sprintf('Ritter Dam-Break: Analytical vs HydroPol2D Centerline Profiles, t \\leq %.0f s', analysis_end_sec));

legend_entries = strings(2*nCompare,1);

ii = 1;

for kk = 1:nCompare

    legend_entries(ii) = sprintf('Analytical t = %d s', compare_times_sec(kk));
    ii = ii + 1;

    legend_entries(ii) = sprintf('HydroPol2D t = %d s', compare_times_sec(kk));
    ii = ii + 1;

end

legend(legend_entries, 'Location','northeastoutside');

set(gca,'FontName','Helvetica','FontSize',12,'LineWidth',1.5,'TickDir','out');

exportgraphics(gcf, fullfile(figDir, 'Ritter_Analytical_vs_HydroPol2D_Profiles.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir, 'Ritter_Analytical_vs_HydroPol2D_Profiles.png'), 'Resolution',300);

%% ================= ERROR PLOTS =================

figure('Color','w','Position',[100 100 1000 500]);

subplot(1,2,1)
plot(Metrics.TargetTime_s, Metrics.RMSE_m, 'o-', 'LineWidth',2);
grid on; box on;
xlabel('Time [s]');
ylabel('RMSE [m]');
title('Depth-profile RMSE');

subplot(1,2,2)
plot(Metrics.TargetTime_s, Metrics.FrontModel_m - Metrics.FrontAnalytical_m, 's-', 'LineWidth',2);
grid on; box on;
xlabel('Time [s]');
ylabel('Front position error [m]');
title('Wave-front error');

exportgraphics(gcf, fullfile(figDir, 'Ritter_Error_Metrics.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir, 'Ritter_Error_Metrics.png'), 'Resolution',300);

%% ================= MASS BALANCE PLOTS AT COMPARE TIMES =================

figure('Color','w','Position',[100 100 1200 500]);

subplot(1,2,1)

plot(Metrics.TargetTime_s, Metrics.ModelStoredVolume_m3, 'o-', ...
    'LineWidth',2, ...
    'DisplayName','Stored volume');
hold on;

plot(Metrics.TargetTime_s, Metrics.OutletVolume_m3, 's-', ...
    'LineWidth',2, ...
    'DisplayName','Cumulative outlet volume');

plot(Metrics.TargetTime_s, Metrics.TotalAccountedVolume_m3, 'd-', ...
    'LineWidth',2, ...
    'DisplayName','Stored + outlet');

yline(V_initial, 'k--', ...
    'LineWidth',1.8, ...
    'DisplayName','Initial volume');

grid on; box on;
xlabel('Time [s]');
ylabel('Volume [m^3]');
title(sprintf('Ritter volume accounting at compare times, t \\leq %.0f s', analysis_end_sec));
legend('Location','best');

subplot(1,2,2)

plot(Metrics.TargetTime_s, Metrics.StorageOnlyError_pct, 'o--', ...
    'LineWidth',1.5, ...
    'DisplayName','Storage-only error');
hold on;

plot(Metrics.TargetTime_s, Metrics.TotalMassError_pct, 's-', ...
    'LineWidth',2, ...
    'DisplayName','Stored + outlet error');

grid on; box on;
xlabel('Time [s]');
ylabel('Mass error [%]');
title(sprintf('Mass conservation error at compare times, t \\leq %.0f s', analysis_end_sec));
legend('Location','best');

exportgraphics(gcf, fullfile(figDir, 'Ritter_Mass_Balance.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir, 'Ritter_Mass_Balance.png'), 'Resolution',300);

%% ================= FULL SAVED-MAP MASS SERIES PLOTS =================

figure('Color','w','Position',[100 100 1200 500]);

subplot(1,2,1)

plot(t_saved_analysis, V_stored_series, 'b-o', ...
    'LineWidth',1.5, ...
    'DisplayName','Stored volume');
hold on;

plot(t_saved_analysis, V_outlet_series, 'r-s', ...
    'LineWidth',1.5, ...
    'DisplayName','Cumulative outlet volume');

plot(t_saved_analysis, V_total_series, 'k-d', ...
    'LineWidth',1.5, ...
    'DisplayName','Stored + outlet');

yline(V_initial, 'k--', ...
    'LineWidth',1.5, ...
    'DisplayName','Initial volume');

grid on; box on;
xlabel('Time [s]');
ylabel('Volume [m^3]');
title(sprintf('Saved-map volume accounting, t \\leq %.0f s', analysis_end_sec));
legend('Location','best');

subplot(1,2,2)

plot(t_saved_analysis, StorageOnlyError_series_pct, 'b-o', ...
    'LineWidth',1.3, ...
    'DisplayName','Storage-only error');
hold on;

plot(t_saved_analysis, TotalMassError_series_pct, 'm-s', ...
    'LineWidth',1.5, ...
    'DisplayName','Stored + outlet error');

grid on; box on;
xlabel('Time [s]');
ylabel('Mass error [%]');
title(sprintf('Mass error across saved maps, t \\leq %.0f s', analysis_end_sec));
legend('Location','best');

exportgraphics(gcf, fullfile(figDir, 'Ritter_Full_SavedMap_MassSeries.pdf'), 'ContentType','vector');
exportgraphics(gcf, fullfile(figDir, 'Ritter_Full_SavedMap_MassSeries.png'), 'Resolution',300);

%% ================= EXPORT TABLES =================

writetable(Profiles, fullfile(tabDir, 'Ritter_Depth_Profiles_Comparison.csv'));
writetable(Metrics,  fullfile(tabDir, 'Ritter_Profile_Error_Metrics.csv'));

MassBalance = Metrics(:, { ...
    'TargetTime_s', ...
    'SavedTime_s', ...
    'InitialVolume_m3', ...
    'ModelStoredVolume_m3', ...
    'OutletVolume_m3', ...
    'TotalAccountedVolume_m3', ...
    'StorageOnlyError_pct', ...
    'TotalMassError_m3', ...
    'TotalMassError_pct'});

writetable(MassBalance, fullfile(tabDir, 'Ritter_Mass_Balance.csv'));

fprintf('\nRitter validation exported to:\n%s\n', outDir);

%% ================= ZOOMED VOLUME ACCOUNTING =================

figure('Color','w','Position',[100 100 1200 500]);

subplot(1,2,1)

plot(t_saved_analysis, V_stored_series, 'b-o', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5, ...
    'DisplayName', 'Stored volume');
hold on;

plot(t_saved_analysis, V_total_series, 'k-d', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5, ...
    'DisplayName', 'Stored + outlet');

yline(V_initial, 'k--', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'Initial volume');

grid on; box on;
xlabel('Time [s]');
ylabel('Volume [m^3]');
title('Zoomed volume accounting');
legend('Location','best');

V_all = [V_stored_series(:); V_total_series(:); V_initial];
Vmin = min(V_all, [], 'omitnan');
Vmax = max(V_all, [], 'omitnan');
ypad = 0.05 * max(Vmax - Vmin, 1);
ylim([Vmin - ypad, Vmax + ypad]);

subplot(1,2,2)

V_deficit = V_initial - V_total_series;
V_deficit_pct = 100 * V_deficit ./ max(V_initial, eps);

plot(t_saved_analysis, V_deficit, 'r-o', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 5);

grid on; box on;
xlabel('Time [s]');
ylabel('Volume deficit [m^3]');
title('Missing volume');

yyaxis right
plot(t_saved_analysis, V_deficit_pct, 'm-s', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5);
ylabel('Volume deficit [%]');

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function h = local_ritter_depth(x_rel, t, h0, c0, g)

    h = zeros(size(x_rel));

    if t <= 0
        h(x_rel <= 0) = h0;
        return;
    end

    left_region = x_rel < -c0*t;
    fan_region  = x_rel >= -c0*t & x_rel <= 2*c0*t;

    h(left_region) = h0;

    h(fan_region) = (4/(9*g)) .* ...
        (c0 - x_rel(fan_region)./(2*t)).^2;

    h(h < 0) = 0;

end

function [D_mm, local_i, store] = load_hydropol_depth_map(global_i, tempDir, saver_memory_maps)

    store = ceil(global_i / saver_memory_maps);
    local_i = global_i - (store-1)*saver_memory_maps;

    % Try both naming conventions used in HydroPol2D scripts
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

    if ~isfield(S, 'Maps') || ~isfield(S.Maps, 'Hydro') || ~isfield(S.Maps.Hydro, 'd')
        error('Temporary map file does not contain Maps.Hydro.d: %s', mapFile);
    end

    D_mm = S.Maps.Hydro.d(:,:,local_i);

end