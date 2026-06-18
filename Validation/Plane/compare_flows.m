%% ============================================================
% Benchmark vs modeled outlet hydrograph + NSE
% Tilted plane analytical kinematic-wave comparison
% + comprehensive mass-balance diagnostic
% ============================================================

%% ================= USER-EDITABLE BENCHMARK PARAMETERS =================

% Plane / hydraulic parameters
S0 = 0.01;                % bed slope [-]
n_manning = 0.015;        % Manning roughness [s/m^(1/3)]

% Rainfall excess intensity
% For the classical tilted-plane benchmark that gives Qss = 55.5556 m3/s
% over a 2000 m x 1000 m plane, this is usually 100 mm/h.
rain_excess_mmhr = 100;   % [mm/h]

% Rainfall duration.
% For the analytical hydrograph below, rainfall is assumed continuous over
% the benchmark window unless you specify otherwise.
rain_duration_min = 180;  % [min]

% Geometry
% Default: infer from model grid.
% If your analytical plane length/width are different, overwrite these.
dx = Wshed_Properties.Resolution;   % [m]

flow_length_m = ny * dx;            % [m] downslope length
plane_width_m = nx * dx;            % [m] lateral width

% Benchmark time settings
benchmark_dt_min = 10;              % [min]
benchmark_end_min = 180;            % [min]

% Kinematic-wave exponent for wide rectangular overland flow with Manning
m_kw = 5/3;

% Optional: include saved-map storage series if available
use_saved_maps_for_mass = true;

%% ================= OUTPUT PATHS =================

if exist('Dirs','var') && isfield(Dirs,'Tables') && ~isempty(Dirs.Tables)
    tabDir = Dirs.Tables;
else
    tabDir = pwd;
end

if exist('Paths','var') && isfield(Paths,'Results') && ~isempty(Paths.Results)
    figDir = fullfile(Paths.Results, 'TiltedPlane_Validation', 'Figures');
else
    figDir = fullfile(pwd, 'TiltedPlane_Validation', 'Figures');
end

if ~isfolder(figDir)
    mkdir(figDir);
end

%% ================= DOMAIN / ACTIVE AREA =================

% Active cells
if exist('idx_nan','var') && ~isempty(idx_nan)

    active_mask = ~idx_nan;

elseif exist('DEM_raster','var') && isfield(DEM_raster,'Z')

    active_mask = isfinite(DEM_raster.Z);

elseif exist('z','var')

    active_mask = isfinite(z);

else

    active_mask = true(ny,nx);

end

% Cell area
if exist('cell_area','var') && isscalar(cell_area)
    A_cell = cell_area;
else
    A_cell = dx^2;
end

A_catchment = nnz(active_mask) * A_cell;  % [m2]

fprintf('\n========================================\n');
fprintf('Tilted plane domain / area\n');
fprintf('dx                       = %.6f m\n', dx);
fprintf('Cell area                = %.6f m2\n', A_cell);
fprintf('Active cells             = %d\n', nnz(active_mask));
fprintf('Active catchment area    = %.6f m2\n', A_catchment);
fprintf('Flow length              = %.3f m\n', flow_length_m);
fprintf('Plane width              = %.3f m\n', plane_width_m);
fprintf('========================================\n\n');

%% ================= ANALYTICAL KINEMATIC-WAVE HYDROGRAPH =================

% Convert rainfall excess to [m/s]
rain_excess_mps = rain_excess_mmhr / 1000 / 3600;

% Manning kinematic-wave coefficient for unit discharge:
% q = alpha * h^m
% q [m2/s], h [m]
alpha_kw = (1 / n_manning) * sqrt(S0);

% Steady unit discharge at outlet [m2/s]
q_unit_steady = rain_excess_mps * flow_length_m;

% Steady total outlet discharge [m3/s]
Q_steady = q_unit_steady * plane_width_m;

% Time to equilibrium / time of concentration [s]
% tc = (L / (alpha * r^(m-1)))^(1/m)
tc_sec = (flow_length_m / (alpha_kw * rain_excess_mps^(m_kw - 1)))^(1 / m_kw);
tc_min = tc_sec / 60;

% Benchmark time vector [min]
t_bench = (0:benchmark_dt_min:benchmark_end_min)';

% Convert time to seconds
t_bench_sec = t_bench * 60;

% Analytical unit discharge [m2/s]
q_unit_bench = alpha_kw .* (rain_excess_mps .* t_bench_sec).^m_kw;

% Do not exceed steady discharge
q_unit_bench = min(q_unit_bench, q_unit_steady);

% Total outlet discharge [m3/s]
q_bench = plane_width_m .* q_unit_bench;

% Ensure exact zero at t = 0
q_bench(t_bench_sec == 0) = 0;

fprintf('\n========================================\n');
fprintf('Tilted plane analytical benchmark\n');
fprintf('Slope S0                  = %.6f [-]\n', S0);
fprintf('Manning n                 = %.6f s/m^(1/3)\n', n_manning);
fprintf('Rainfall excess           = %.6f mm/h\n', rain_excess_mmhr);
fprintf('Rainfall duration         = %.6f min\n', rain_duration_min);
fprintf('Steady discharge Qss      = %.6f m3/s\n', Q_steady);
fprintf('Time of concentration tc  = %.3f min\n', tc_min);
fprintf('========================================\n\n');

%% ================= MODELED HYDROGRAPH =================

t_model = gather(running_control.time_hydrograph(:));   % [min or datetime]
q_model = gather(outlet_states.outlet_hydrograph(:));   % [m3/s]

% Convert datetime/duration to minutes if needed
if isdatetime(t_model)
    t_model = minutes(t_model - t_model(1));
elseif isduration(t_model)
    t_model = minutes(t_model);
end

t_model = double(t_model(:));
q_model = double(q_model(:));

% Remove invalid values
valid = isfinite(t_model) & isfinite(q_model);
t_model = t_model(valid);
q_model = q_model(valid);

% Sort model data in case time vector is not strictly ordered
[t_model, sort_id] = sort(t_model);
q_model = q_model(sort_id);

% Remove duplicate times if present
[t_model, unique_id] = unique(t_model, 'stable');
q_model = q_model(unique_id);

% Outlet discharge cannot be negative
q_model = max(q_model, 0);

% Convert model hydrograph time to seconds
t_model_sec = t_model * 60;

%% Cumulative outlet volume from modeled hydrograph [m3]
% Two possible interpretations:
% 1) trapz: q_model is instantaneous at timestamps.
% 2) backward rectangular: q_model(i) is the mean/end-step discharge
%    over interval [t(i-1), t(i)].
%
% HydroPol2D hydrographs often behave more like interval/output-step values,
% so compare both.

if numel(t_model_sec) >= 2

    % Instantaneous-point interpretation
    V_out_model_cum_trapz = cumtrapz(t_model_sec, q_model);

    % End-of-interval / backward-rectangle interpretation
    dt_hydro_sec = [0; diff(t_model_sec)];
    V_out_model_cum_rect = cumsum(q_model .* dt_hydro_sec);

    % Use rectangular integration for mass balance.
    V_out_model_cum = V_out_model_cum_rect;

else

    V_out_model_cum_trapz = zeros(size(t_model_sec));
    V_out_model_cum_rect  = zeros(size(t_model_sec));
    V_out_model_cum       = zeros(size(t_model_sec));

end

fprintf('\n============= OUTLET VOLUME INTEGRATION CHECK =============\n');
fprintf('Outlet volume, trapz final = %.12e m3\n', V_out_model_cum_trapz(end));
fprintf('Outlet volume, rect final  = %.12e m3\n', V_out_model_cum_rect(end));
fprintf('Difference rect - trapz    = %.12e m3\n', ...
    V_out_model_cum_rect(end) - V_out_model_cum_trapz(end));
fprintf('===========================================================\n\n');

%% ================= INTERPOLATION =================

% Align model to benchmark times.
% Use NaN outside model range instead of extrapolating.
q_model_i = interp1(t_model, q_model, t_bench, 'linear', NaN);

% Preserve initial model value if available
if ~isempty(q_model_i) && ~isempty(q_model)
    q_model_i(1) = q_model(1);
end

% Use only times where both analytical and modeled values exist
valid_eval = isfinite(q_model_i) & isfinite(q_bench);

t_eval = t_bench(valid_eval);
q_bench_eval = q_bench(valid_eval);
q_model_eval = q_model_i(valid_eval);

%% ================= NSE =================

den_nse = sum((q_bench_eval - mean(q_bench_eval)).^2);

if den_nse > 0
    NSE = 1 - sum((q_model_eval - q_bench_eval).^2) / den_nse;
else
    NSE = NaN;
end

fprintf('\n========================================\n');
fprintf('Tilted plane benchmark NSE = %.4f\n', NSE);
fprintf('========================================\n\n');

%% ================= ADDITIONAL ERROR DIAGNOSTICS =================

error_abs = q_model_i - q_bench;
abs_error = abs(error_abs);

rel_error = abs_error ./ max(abs(q_bench), 1e-6);

RMSE = sqrt(mean((q_model_eval - q_bench_eval).^2, 'omitnan'));
MAE  = mean(abs(q_model_eval - q_bench_eval), 'omitnan');
BIAS = mean(q_model_eval - q_bench_eval, 'omitnan');

fprintf('RMSE = %.6f m3/s\n', RMSE);
fprintf('MAE  = %.6f m3/s\n', MAE);
fprintf('BIAS = %.6f m3/s\n\n', BIAS);

%% ================= MASS BALANCE SETUP =================

% Initial stored volume from first saved map if available; otherwise from
% current/initial depth variables if available; otherwise assume dry.
V_initial = NaN;
initial_source = "unavailable";

if use_saved_maps_for_mass && ...
        exist('Paths','var') && isfield(Paths,'Temp') && ...
        exist('saver_memory_maps','var') && ...
        exist('running_control','var') && isfield(running_control,'time_records')

    try
        [D0_mm, ~, ~] = load_hydropol_depth_map(1, Paths.Temp, saver_memory_maps);
        D0_m = double(gather(D0_mm)) / 1000;
        D0_m(~isfinite(D0_m)) = 0;
        D0_m = max(D0_m, 0);

        V_initial = sum(D0_m(active_mask), 'omitnan') * A_cell;
        initial_source = "first saved map";

    catch ME
        warning('Could not read first saved map for initial volume: %s', ME.message);
    end

end

if ~isfinite(V_initial)

    if exist('d_initial','var')
        d0_mm = gather(d_initial);
        d0_m = max(double(d0_mm),0) / 1000;
        V_initial = sum(d0_m(active_mask), 'omitnan') * A_cell;
        initial_source = "d_initial";

    elseif exist('depths','var') && isfield(depths,'d_0')
        d0_mm = gather(depths.d_0);
        d0_m = max(double(d0_mm),0) / 1000;
        V_initial = sum(d0_m(active_mask), 'omitnan') * A_cell;
        initial_source = "depths.d_0";

    else
        V_initial = 0;
        initial_source = "assumed dry";
    end

end

fprintf('\n============= MASS BALANCE SETUP =============\n');
fprintf('Initial volume source     = %s\n', initial_source);
fprintf('Initial stored volume     = %.12e m3\n', V_initial);
fprintf('Catchment area            = %.12e m2\n', A_catchment);
fprintf('Rainfall intensity        = %.6f mm/h\n', rain_excess_mmhr);
fprintf('Rainfall duration         = %.6f min\n', rain_duration_min);
fprintf('==============================================\n\n');

%% ================= MASS BALANCE AT HYDROGRAPH TIMES =================

% Cumulative rainfall volume at model hydrograph times
rain_time_model_min = min(max(t_model, 0), rain_duration_min);
V_rain_model_cum = rain_excess_mps .* (rain_time_model_min * 60) .* A_catchment;

% Stored volume at hydrograph times is only available if a final/current
% depth field is available. Saved-map series is handled separately below.
V_stored_final = NaN;
stored_final_source = "unavailable";

if exist('depths','var') && isfield(depths,'d_t')

    d_final_mm = gather(depths.d_t);
    d_final_m = max(double(d_final_mm), 0) / 1000;

    V_stored_final = sum(d_final_m(active_mask), 'omitnan') * A_cell;
    stored_final_source = "depths.d_t";

elseif exist('d_t','var')

    d_final_mm = gather(d_t);
    d_final_m = max(double(d_final_mm), 0) / 1000;

    V_stored_final = sum(d_final_m(active_mask), 'omitnan') * A_cell;
    stored_final_source = "d_t";

end

V_rain_total = rain_excess_mps * (rain_duration_min * 60) * A_catchment;
V_out_model_total = 0;

if ~isempty(V_out_model_cum)
    V_out_model_total = V_out_model_cum(end);
end

if isfinite(V_stored_final)

    V_accounted_final = V_stored_final + V_out_model_total;
    V_expected_final = V_initial + V_rain_total;

    V_residual_final = V_accounted_final - V_expected_final;
    mass_error_final_pct = 100 * V_residual_final / max(V_expected_final, eps);

else

    V_accounted_final = NaN;
    V_expected_final = V_initial + V_rain_total;
    V_residual_final = NaN;
    mass_error_final_pct = NaN;

end

fprintf('============= FINAL EVENT MASS BALANCE =============\n');
fprintf('Rainfall volume total       = %.12e m3\n', V_rain_total);
fprintf('Outlet volume total         = %.12e m3\n', V_out_model_total);
fprintf('Final stored volume source  = %s\n', stored_final_source);
fprintf('Final stored volume         = %.12e m3\n', V_stored_final);
fprintf('Expected final volume basis = initial + rain = %.12e m3\n', V_expected_final);
fprintf('Accounted final volume      = stored + outlet = %.12e m3\n', V_accounted_final);
fprintf('Final residual              = %.12e m3\n', V_residual_final);
fprintf('Final mass error            = %.12e %%\n', mass_error_final_pct);
fprintf('====================================================\n\n');

%% ================= SAVED-MAP MASS BALANCE SERIES =================

has_saved_mass_series = false;

if use_saved_maps_for_mass && ...
        exist('Paths','var') && isfield(Paths,'Temp') && ...
        exist('saver_memory_maps','var') && ...
        exist('running_control','var') && isfield(running_control,'time_records')

    try

        t_records = gather(running_control.time_records(:));

        if isdatetime(t_records)
            t_saved_min = minutes(t_records - t_records(1));
        elseif isduration(t_records)
            t_saved_min = minutes(t_records);
        else
            % HydroPol2D commonly stores elapsed time in minutes
            t_saved_min = double(t_records);
        end

        t_saved_min = double(t_saved_min(:));
        t_saved_sec = t_saved_min * 60;

        nSaved = numel(t_saved_min);

        V_stored_series = NaN(nSaved,1);
        V_rain_series   = NaN(nSaved,1);
        V_out_series    = NaN(nSaved,1);
        V_expected_series = NaN(nSaved,1);
        V_accounted_series = NaN(nSaved,1);
        V_residual_series = NaN(nSaved,1);
        MassError_series_pct = NaN(nSaved,1);

        for jj = 1:nSaved

            [Dj_mm, ~, ~] = load_hydropol_depth_map(jj, Paths.Temp, saver_memory_maps);

            Dj_m = double(gather(Dj_mm)) / 1000;
            Dj_m(~isfinite(Dj_m)) = 0;
            Dj_m = max(Dj_m, 0);

            V_stored_series(jj) = sum(Dj_m(active_mask), 'omitnan') * A_cell;

            rain_t_min = min(max(t_saved_min(jj), 0), rain_duration_min);
            V_rain_series(jj) = rain_excess_mps * (rain_t_min * 60) * A_catchment;

            if numel(t_model_sec) >= 2
                if t_saved_sec(jj) <= max(t_model_sec)
                    V_out_series(jj) = interp1(t_model_sec, V_out_model_cum, t_saved_sec(jj), 'linear', 'extrap');
                else
                    V_out_series(jj) = V_out_model_cum(end);
                end
                V_out_series(jj) = max(V_out_series(jj), 0);
            else
                V_out_series(jj) = 0;
            end

            V_expected_series(jj) = V_initial + V_rain_series(jj);
            V_accounted_series(jj) = V_stored_series(jj) + V_out_series(jj);

            V_residual_series(jj) = V_accounted_series(jj) - V_expected_series(jj);

            MassError_series_pct(jj) = ...
                100 * V_residual_series(jj) / max(V_expected_series(jj), eps);

        end

        has_saved_mass_series = true;

        FullMassSeries = table( ...
            t_saved_min, ...
            V_rain_series, ...
            V_stored_series, ...
            V_out_series, ...
            V_expected_series, ...
            V_accounted_series, ...
            V_residual_series, ...
            MassError_series_pct, ...
            'VariableNames', { ...
                'Time_min', ...
                'RainVolume_m3', ...
                'StoredVolume_m3', ...
                'OutletVolume_m3', ...
                'ExpectedVolume_initial_plus_rain_m3', ...
                'AccountedVolume_stored_plus_outlet_m3', ...
                'Residual_m3', ...
                'MassError_pct'});

        full_mass_file = fullfile(tabDir, 'TiltedPlane_Full_SavedMap_MassSeries.csv');
        writetable(FullMassSeries, full_mass_file);

        fprintf('Saved-map mass series exported to:\n%s\n', full_mass_file);

        fprintf('\n============= SAVED-MAP MASS SERIES SUMMARY =============\n');
        fprintf('Max abs mass error = %.12e %%\n', max(abs(MassError_series_pct), [], 'omitnan'));
        fprintf('Final saved-map mass error = %.12e %%\n', MassError_series_pct(end));
        fprintf('==========================================================\n\n');

    catch ME

        warning('Saved-map mass balance series could not be computed: %s', ME.message);
        has_saved_mass_series = false;

    end

end

%% ================= HYDROGRAPH PLOT =================

figure('Color','w','Position',[100 100 900 500]);

plot(t_bench, q_bench, 'k-o', ...
    'LineWidth', 2, ...
    'MarkerSize', 5, ...
    'DisplayName', 'Analytical (Kinematic wave)'); 
hold on;

plot(t_bench, q_model_i, 'r-s', ...
    'LineWidth', 2, ...
    'MarkerSize', 5, ...
    'DisplayName', 'HydroPol2D');

plot(t_model, q_model, 'r:', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'HydroPol2D raw');

grid on; box on;

xlabel('Time [min]');
ylabel('Outlet discharge [m^3/s]');

title(sprintf('Tilted Plane Hydrograph Comparison (NSE = %.3f)', NSE));

legend('Location', 'best');

exportgraphics(gcf, fullfile(figDir, 'TiltedPlane_Hydrograph_Comparison.png'), 'Resolution',300);

%% ================= ERROR PLOTS =================

figure('Color','w','Position',[150 150 900 450]);

plot(t_bench, rel_error, 'b-', 'LineWidth', 2);
grid on; box on;

xlabel('Time [min]');
ylabel('Relative error [-]');
title('Relative Error vs Time');

exportgraphics(gcf, fullfile(figDir, 'TiltedPlane_Relative_Error.png'), 'Resolution',300);

figure('Color','w','Position',[200 200 900 450]);

plot(t_bench, error_abs, 'm-', 'LineWidth', 2);
grid on; box on;

xlabel('Time [min]');
ylabel('Model - Analytical [m^3/s]');
title('Signed Error vs Time');

exportgraphics(gcf, fullfile(figDir, 'TiltedPlane_Signed_Error.png'), 'Resolution',300);

%% ================= CUMULATIVE OUTLET / RAINFALL PLOT =================

figure('Color','w','Position',[250 250 1000 500]);

plot(t_model, V_out_model_cum, 'r-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Cumulative outlet volume');
hold on;

plot(t_model, V_rain_model_cum, 'b--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Cumulative rainfall volume');

grid on; box on;

xlabel('Time [min]');
ylabel('Volume [m^3]');
title('Cumulative Rainfall and Outlet Volume');
legend('Location','best');

exportgraphics(gcf, fullfile(figDir, 'TiltedPlane_Cumulative_Rainfall_Outlet.png'), 'Resolution',300);

%% ================= SAVED-MAP MASS BALANCE PLOTS =================

if has_saved_mass_series

    figure('Color','w','Position',[300 300 1200 500]);

    subplot(1,2,1)

    plot(t_saved_min, V_rain_series, 'b--', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Rainfall volume');
    hold on;

    plot(t_saved_min, V_stored_series, 'g-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5, ...
        'DisplayName', 'Stored volume');

    plot(t_saved_min, V_out_series, 'r-s', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5, ...
        'DisplayName', 'Outlet volume');

    plot(t_saved_min, V_accounted_series, 'k-d', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5, ...
        'DisplayName', 'Stored + outlet');

    plot(t_saved_min, V_expected_series, 'k--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Initial + rain');

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Volume [m^3]');
    title('Tilted Plane Volume Accounting');
    legend('Location','best');

    subplot(1,2,2)

    plot(t_saved_min, MassError_series_pct, 'm-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5);

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Mass error [%]');
    title(sprintf('Mass Error, max = %.3e %%', ...
        max(abs(MassError_series_pct), [], 'omitnan')));
    ytickformat('%.2e');

    exportgraphics(gcf, fullfile(figDir, 'TiltedPlane_Mass_Balance_SavedMaps.png'), 'Resolution',300);

    %% Zoomed residual / deficit plot

    figure('Color','w','Position',[350 350 1200 500]);

    subplot(1,2,1)

    plot(t_saved_min, V_residual_series, 'r-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5);

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Residual volume [m^3]');
    title('Mass residual: stored + outlet - initial - rain');

    subplot(1,2,2)

    plot(t_saved_min, MassError_series_pct, 'm-s', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 5);

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Mass error [%]');
    title('Mass error');
    ytickformat('%.2e');

    exportgraphics(gcf, fullfile(figDir, 'TiltedPlane_Mass_Residual_Zoom.png'), 'Resolution',300);

end



%% ================= EXPORT HYDROGRAPH TABLE =================

Hydrograph_Comparison = table( ...
    t_bench, ...
    q_bench, ...
    q_model_i, ...
    error_abs, ...
    abs_error, ...
    rel_error, ...
    'VariableNames', { ...
        'Time_min', ...
        'Analytical_m3s', ...
        'Model_m3s', ...
        'Error_m3s', ...
        'AbsError_m3s', ...
        'RelError'});

export_file = fullfile(tabDir, 'TiltedPlane_Benchmark_vs_Model.csv');
writetable(Hydrograph_Comparison, export_file);

fprintf('Hydrograph comparison table exported to:\n%s\n', export_file);

%% ================= EXPORT EVENT MASS BALANCE TABLE =================

EventMassBalance = table( ...
    S0, ...
    n_manning, ...
    rain_excess_mmhr, ...
    rain_duration_min, ...
    A_catchment, ...
    V_initial, ...
    V_rain_total, ...
    V_out_model_total, ...
    V_stored_final, ...
    V_expected_final, ...
    V_accounted_final, ...
    V_residual_final, ...
    mass_error_final_pct, ...
    'VariableNames', { ...
        'Slope', ...
        'Manning_n', ...
        'Rain_excess_mmhr', ...
        'Rain_duration_min', ...
        'Catchment_area_m2', ...
        'Initial_volume_m3', ...
        'Rain_volume_total_m3', ...
        'Outlet_volume_total_m3', ...
        'Final_stored_volume_m3', ...
        'Expected_final_volume_initial_plus_rain_m3', ...
        'Accounted_final_volume_stored_plus_outlet_m3', ...
        'Final_residual_m3', ...
        'Final_mass_error_pct'});

mass_file = fullfile(tabDir, 'TiltedPlane_Event_Mass_Balance.csv');
writetable(EventMassBalance, mass_file);

fprintf('Event mass balance table exported to:\n%s\n', mass_file);

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

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