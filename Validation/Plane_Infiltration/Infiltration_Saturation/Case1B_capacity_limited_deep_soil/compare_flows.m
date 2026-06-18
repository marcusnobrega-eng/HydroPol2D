%% ============================================================
% Tilted Plane Flow Comparison + Infiltration Evolution
% ============================================================
%
% Purpose:
%   1) Keep the outlet hydrograph comparison against the no-infiltration
%      analytical kinematic-wave benchmark.
%
%   2) Add infiltration diagnostics from saved maps:
%        Maps.Hydro.I_t   [mm]    soil / bucket storage
%        Maps.Hydro.f     [mm/h]  actual infiltration flux
%        Maps.Hydro.C     [mm/h]  infiltration capacity
%
%   3) Test the hypothesis:
%        "The subsurface fills, then infiltration decreases/stops."
%
% Required after HydroPol2D run:
%   running_control.time_hydrograph
%   outlet_states.outlet_hydrograph
%   running_control.time_records
%   saved map files in Paths.Temp, or current Maps in memory
%
% Notes:
%   - This script compares the current HydroPol2D run against a
%     no-infiltration analytical benchmark.
%   - If your HydroPol2D run includes infiltration, the model hydrograph
%     should be lower/delayed relative to the no-infiltration benchmark.
%   - The infiltration diagnostics are computed from saved map stacks.
%
% ============================================================

%% ================= USER PARAMETERS =================

% Rainfall input used in the analytical no-infiltration benchmark
rain_excess_mmhr = 20;      % [mm/h]

% Rainfall duration / benchmark duration
rain_duration_min = 12 * 60;  % [min]
benchmark_end_min = 12 * 60;  % [min]
benchmark_dt_min  = 10;      % [min]

% Manning roughness
n_manning = 0.015;           % [s/m^(1/3)]

% If the DEM is available, slope is inferred from the DEM.
% This fallback matches the uploaded synthetic catchment script.
S0_fallback = 0.01;          % [-]
% If you truly generated a 0.01 percent slope, this value would be 0.0001.
% However, the uploaded generator code uses S0 = 0.01.

% Kinematic-wave exponent for Manning wide-sheet flow
m_kw = 5/3;

% Saved maps
use_saved_maps = true;

% If saver_memory_maps is not in the workspace, this matches your save logic.
if ~exist('saver_memory_maps','var') || isempty(saver_memory_maps)
    saver_memory_maps = 12;
end

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

dx = Wshed_Properties.Resolution;   % [m]

if exist('ny','var') && exist('nx','var')
    ny0 = ny;
    nx0 = nx;
elseif exist('DEM_raster','var') && isfield(DEM_raster,'Z')
    [ny0, nx0] = size(DEM_raster.Z);
elseif exist('elevation','var')
    [ny0, nx0] = size(elevation);
else
    error('Cannot infer nx/ny. Please provide ny/nx, DEM_raster.Z, or elevation.');
end

% Active cells
if exist('idx_nan','var') && ~isempty(idx_nan)
    active_mask = ~idx_nan;
elseif exist('DEM_raster','var') && isfield(DEM_raster,'Z')
    active_mask = isfinite(double(gather(DEM_raster.Z)));
elseif exist('elevation','var')
    active_mask = isfinite(double(gather(elevation)));
else
    active_mask = true(ny0, nx0);
end

% Cell area
if exist('cell_area','var') && isscalar(cell_area)
    A_cell = cell_area;
else
    A_cell = dx^2;
end

A_catchment = nnz(active_mask) * A_cell;  % [m2]

flow_length_m = ny0 * dx;       % [m]
plane_width_m = nx0 * dx;       % [m]

%% ================= INFER SLOPE FROM DEM IF POSSIBLE =================

S0 = S0_fallback;

DEM_for_slope = [];

if exist('DEM_raster','var') && isfield(DEM_raster,'Z')
    DEM_for_slope = double(gather(DEM_raster.Z));
elseif exist('elevation','var')
    DEM_for_slope = double(gather(elevation));
end

if ~isempty(DEM_for_slope)
    row_mean_z = mean(DEM_for_slope, 2, 'omitnan');

    if numel(row_mean_z) >= 2
        dz = row_mean_z(end) - row_mean_z(1);
        dy = max((numel(row_mean_z) - 1) * dx, dx);

        S0_inferred = abs(dz) / dy;

        if isfinite(S0_inferred) && S0_inferred > 0
            S0 = S0_inferred;
        end
    end
end

fprintf('\n========================================\n');
fprintf('Tilted plane domain / benchmark setup\n');
fprintf('dx                       = %.6f m\n', dx);
fprintf('Active cells             = %d\n', nnz(active_mask));
fprintf('Active catchment area    = %.6f m2\n', A_catchment);
fprintf('Flow length              = %.3f m\n', flow_length_m);
fprintf('Plane width              = %.3f m\n', plane_width_m);
fprintf('Slope used S0            = %.8f [-]\n', S0);
fprintf('Manning n                = %.6f\n', n_manning);
fprintf('Rainfall intensity       = %.6f mm/h\n', rain_excess_mmhr);
fprintf('Rainfall duration        = %.6f min\n', rain_duration_min);
fprintf('========================================\n\n');

%% ================= ANALYTICAL NO-INFILTRATION HYDROGRAPH =================

% Convert rainfall excess to [m/s]
rain_excess_mps = rain_excess_mmhr / 1000 / 3600;

% Manning kinematic-wave coefficient for unit discharge:
% q = alpha h^m
alpha_kw = (1 / n_manning) * sqrt(S0);

% Steady unit discharge at outlet [m2/s]
q_unit_steady = rain_excess_mps * flow_length_m;

% Steady total outlet discharge [m3/s]
Q_steady_noinf = q_unit_steady * plane_width_m;

% Time to equilibrium / time of concentration [s]
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
fprintf('No-infiltration analytical benchmark\n');
fprintf('Steady no-infiltration Qss = %.6f m3/s\n', Q_steady_noinf);
fprintf('Time of concentration tc   = %.3f min\n', tc_min);
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

% Sort model data
[t_model, sort_id] = sort(t_model);
q_model = q_model(sort_id);

% Remove duplicate times
[t_model, unique_id] = unique(t_model, 'stable');
q_model = q_model(unique_id);

% Outlet discharge cannot be negative
q_model = max(q_model, 0);

% Convert model hydrograph time to seconds
t_model_sec = t_model * 60;

% Interpolate model to benchmark times
q_model_i = interp1(t_model, q_model, t_bench, 'linear', NaN);

if ~isempty(q_model_i) && ~isempty(q_model)
    q_model_i(1) = q_model(1);
end

% NSE against no-infiltration analytical benchmark
valid_eval = isfinite(q_model_i) & isfinite(q_bench);

if nnz(valid_eval) >= 2
    q_bench_eval = q_bench(valid_eval);
    q_model_eval = q_model_i(valid_eval);

    den_nse = sum((q_bench_eval - mean(q_bench_eval)).^2);

    if den_nse > 0
        NSE = 1 - sum((q_model_eval - q_bench_eval).^2) / den_nse;
    else
        NSE = NaN;
    end
else
    NSE = NaN;
end

fprintf('\n========================================\n');
fprintf('NSE against no-infiltration benchmark = %.4f\n', NSE);
fprintf('========================================\n\n');

%% ================= OUTLET VOLUME =================

if numel(t_model_sec) >= 2
    dt_hydro_sec = [0; diff(t_model_sec)];
    V_out_model_cum = cumsum(q_model .* dt_hydro_sec);
else
    V_out_model_cum = zeros(size(t_model_sec));
end

V_rain_model_cum = rain_excess_mps .* ...
    (min(max(t_model,0), rain_duration_min) * 60) .* A_catchment;

%% ================= READ SAVED INFILTRATION MAPS =================

has_infiltration_series = false;

t_saved_min = [];
mean_It_mm = [];
median_It_mm = [];
max_It_mm = [];
mean_f_mmhr = [];
max_f_mmhr = [];
mean_C_mmhr = [];
max_C_mmhr = [];
mean_depth_mm = [];
soil_storage_volume_m3 = [];
surface_storage_volume_m3 = [];
remaining_capacity_volume_m3 = [];
mean_storage_fill_frac = [];

if use_saved_maps && exist('running_control','var') && isfield(running_control,'time_records')

    t_records_raw = gather(running_control.time_records(:));

    if isdatetime(t_records_raw)
        t_all_records_min = minutes(t_records_raw - t_records_raw(1));
    elseif isduration(t_records_raw)
        t_all_records_min = minutes(t_records_raw);
    else
        t_all_records_min = double(t_records_raw);
    end

    t_all_records_min = double(t_all_records_min(:));

    % Only use records within the modeled simulation window.
    if ~isempty(t_model)
        tmax_for_maps = max(t_model);
    else
        tmax_for_maps = benchmark_end_min;
    end

    valid_record = isfinite(t_all_records_min) & ...
                   t_all_records_min >= 0 & ...
                   t_all_records_min <= tmax_for_maps + 1e-9;

    global_record_ids = find(valid_record);
    t_saved_min = t_all_records_min(valid_record);

    nSaved = numel(t_saved_min);

    mean_It_mm = NaN(nSaved,1);
    median_It_mm = NaN(nSaved,1);
    max_It_mm = NaN(nSaved,1);

    mean_f_mmhr = NaN(nSaved,1);
    max_f_mmhr = NaN(nSaved,1);

    mean_C_mmhr = NaN(nSaved,1);
    max_C_mmhr = NaN(nSaved,1);

    mean_depth_mm = NaN(nSaved,1);

    soil_storage_volume_m3 = NaN(nSaved,1);
    surface_storage_volume_m3 = NaN(nSaved,1);
    remaining_capacity_volume_m3 = NaN(nSaved,1);
    mean_storage_fill_frac = NaN(nSaved,1);

    % Estimate unsaturated-zone storage capacity [mm]
    UZ_capacity_mm = infer_uz_capacity_mm( ...
        [ny0 nx0], ...
        active_mask, ...
        Soil_Properties, ...
        DEM_for_slope, ...
        exist_var('elevation'), ...
        exist_var('BC_States'), ...
        exist_var('GW_table'));

    UZ_capacity_mm = double(gather(UZ_capacity_mm));
    UZ_capacity_mm(~isfinite(UZ_capacity_mm)) = NaN;

    if exist('Paths','var') && isfield(Paths,'Temp') && ~isempty(Paths.Temp)
        tempDir = Paths.Temp;
    else
        tempDir = '';
    end

    if exist('Maps','var')
        Maps_fallback = Maps;
    else
        Maps_fallback = [];
    end

    last_It_map = [];
    last_f_map = [];
    last_C_map = [];

    for jj = 1:nSaved

        global_i = global_record_ids(jj);

        try
            [HydroMaps, local_i, ~] = load_hydropol_hydro_stack( ...
                global_i, tempDir, saver_memory_maps, Maps_fallback);

            % Depth map [mm]
            [D_mm, ok_d] = read_hydro_field(HydroMaps, 'd', local_i);

            if ok_d
                D_mm = to_double_cpu(D_mm);
                D_mm(~isfinite(D_mm)) = NaN;
                D_mm = max(D_mm, 0);

                mean_depth_mm(jj) = mean(D_mm(active_mask), 'omitnan');
                surface_storage_volume_m3(jj) = ...
                    sum((D_mm(active_mask) / 1000) * A_cell, 'omitnan');
            end

            % Soil storage I_t [mm]
            [It_mm, ok_it] = read_hydro_field(HydroMaps, 'I_t', local_i);

            if ok_it
                It_mm = to_double_cpu(It_mm);
                It_mm(~isfinite(It_mm)) = NaN;
                It_mm = max(It_mm, 0);

                mean_It_mm(jj) = mean(It_mm(active_mask), 'omitnan');
                median_It_mm(jj) = median(It_mm(active_mask), 'omitnan');
                max_It_mm(jj) = max(It_mm(active_mask), [], 'omitnan');

                soil_storage_volume_m3(jj) = ...
                    sum((It_mm(active_mask) / 1000) * A_cell, 'omitnan');

                rem_mm = max(UZ_capacity_mm - It_mm, 0);

                remaining_capacity_volume_m3(jj) = ...
                    sum((rem_mm(active_mask) / 1000) * A_cell, 'omitnan');

                fill_frac = It_mm ./ max(UZ_capacity_mm, 1e-12);
                fill_frac(~isfinite(fill_frac)) = NaN;
                fill_frac = min(max(fill_frac, 0), 1);

                mean_storage_fill_frac(jj) = mean(fill_frac(active_mask), 'omitnan');

                last_It_map = It_mm;
            end

            % Actual infiltration flux f [mm/h]
            [F_mmhr, ok_f] = read_hydro_field(HydroMaps, 'f', local_i);

            if ok_f
                F_mmhr = to_double_cpu(F_mmhr);
                F_mmhr(~isfinite(F_mmhr)) = NaN;
                F_mmhr = max(F_mmhr, 0);

                mean_f_mmhr(jj) = mean(F_mmhr(active_mask), 'omitnan');
                max_f_mmhr(jj) = max(F_mmhr(active_mask), [], 'omitnan');

                last_f_map = F_mmhr;
            end

            % Infiltration capacity C [mm/h]
            [C_mmhr, ok_c] = read_hydro_field(HydroMaps, 'C', local_i);

            if ok_c
                C_mmhr = to_double_cpu(C_mmhr);
                C_mmhr(~isfinite(C_mmhr)) = NaN;
                C_mmhr = max(C_mmhr, 0);

                mean_C_mmhr(jj) = mean(C_mmhr(active_mask), 'omitnan');
                max_C_mmhr(jj) = max(C_mmhr(active_mask), [], 'omitnan');

                last_C_map = C_mmhr;
            end

        catch ME
            warning('Could not read saved infiltration maps for record %d: %s', ...
                global_i, ME.message);
        end
    end

    has_infiltration_series = any(isfinite(mean_It_mm)) || ...
                              any(isfinite(mean_f_mmhr)) || ...
                              any(isfinite(mean_C_mmhr));

end

%% ================= INFILTRATION VOLUME / EQUIVALENT FLOW =================

if has_infiltration_series

    t_saved_sec = t_saved_min * 60;

    if numel(t_saved_sec) >= 2
        dt_saved_sec = [0; diff(t_saved_sec)];
    else
        dt_saved_sec = zeros(size(t_saved_sec));
    end

    % Area-equivalent infiltration rate [m3/s]
    Q_infiltration_equiv_m3s = ...
        (mean_f_mmhr / 1000 / 3600) * A_catchment;

    % Cumulative infiltrated volume estimated from mean f [m3]
    V_infiltration_cum_from_f = cumsum(Q_infiltration_equiv_m3s .* dt_saved_sec);

    % Instantaneous rainfall volume rate [m3/s]
    Q_rain_m3s = rain_excess_mps * A_catchment;

    % Instantaneous rainfall-excess-equivalent rate after infiltration.
    % This is not routed; it is a diagnostic water-balance rate.
    Q_rain_minus_infiltration_m3s = ...
        max(Q_rain_m3s - Q_infiltration_equiv_m3s, 0);

    % Interpolate outlet hydrograph to saved-map times
    if numel(t_model) >= 2
        q_model_saved = interp1(t_model, q_model, t_saved_min, 'linear', NaN);
        V_out_saved = interp1(t_model, V_out_model_cum, t_saved_min, 'linear', NaN);
    else
        q_model_saved = NaN(size(t_saved_min));
        V_out_saved = NaN(size(t_saved_min));
    end

    % Cumulative rainfall volume at saved-map times
    V_rain_saved = rain_excess_mps .* ...
        (min(max(t_saved_min,0), rain_duration_min) * 60) .* A_catchment;

    fprintf('\n================ INFILTRATION SERIES SUMMARY ================\n');
    fprintf('Final mean I_t                  = %.6f mm\n', mean_It_mm(end));
    fprintf('Final mean f                    = %.6f mm/h\n', mean_f_mmhr(end));
    fprintf('Final mean C                    = %.6f mm/h\n', mean_C_mmhr(end));
    fprintf('Final mean storage fill fraction= %.6f [-]\n', mean_storage_fill_frac(end));
    fprintf('Final cumulative infiltration   = %.6e m3\n', V_infiltration_cum_from_f(end));
    fprintf('Final soil storage volume       = %.6e m3\n', soil_storage_volume_m3(end));
    fprintf('Final remaining UZ capacity     = %.6e m3\n', remaining_capacity_volume_m3(end));
    fprintf('=============================================================\n\n');

end

%% ============================================================
% PLOT 1: HYDROGRAPH COMPARISON ONLY
% ============================================================

figure('Color','w','Position',[100 100 1000 520]);

plot(t_bench, q_bench, 'k-', ...
    'LineWidth', 2.2, ...
    'DisplayName', 'Analytical no-infiltration benchmark');
hold on;

plot(t_bench, q_model_i, 'r-s', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 4, ...
    'DisplayName', 'HydroPol2D interpolated');

plot(t_model, q_model, 'r:', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'HydroPol2D raw');

yline(Q_steady_noinf, 'k--', ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf('No-infiltration Q_{ss}=%.3f m^3/s', Q_steady_noinf));

grid on; box on;

xlabel('Time [min]');
ylabel('Outlet discharge [m^3/s]');

title(sprintf('Tilted Plane Hydrograph: Model vs No-Infiltration Benchmark, NSE = %.3f', NSE));

legend('Location','best');

exportgraphics(gcf, ...
    fullfile(figDir, 'TiltedPlane_Hydrograph_Comparison_WithInfiltration.png'), ...
    'Resolution', 300);

%% ============================================================
% PLOT 2: INFILTRATION EVOLUTION
% ============================================================

if has_infiltration_series

    figure('Color','w','Position',[150 150 1300 850]);

    % --------------------------------------------------------
    % Mean f and C
    % --------------------------------------------------------
    subplot(2,2,1)

    plot(t_saved_min, mean_f_mmhr, 'b-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 4, ...
        'DisplayName', 'Mean actual infiltration f');
    hold on;

    plot(t_saved_min, mean_C_mmhr, 'm-s', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 4, ...
        'DisplayName', 'Mean capacity C');

    yline(rain_excess_mmhr, 'k--', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'Rainfall');

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Rate [mm/h]');
    title('Infiltration Rate and Capacity');
    legend('Location','best');

    % --------------------------------------------------------
    % Storage filling
    % --------------------------------------------------------
    subplot(2,2,2)

    yyaxis left
    plot(t_saved_min, mean_It_mm, 'g-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 4, ...
        'DisplayName', 'Mean I_t');
    ylabel('Mean soil storage I_t [mm]');

    yyaxis right
    plot(t_saved_min, 100 * mean_storage_fill_frac, 'k-s', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 4, ...
        'DisplayName', 'Mean storage fill');
    ylabel('Mean storage fill [%]');

    grid on; box on;
    xlabel('Time [min]');
    title('Subsurface Storage Filling');

    % --------------------------------------------------------
    % Volumes
    % --------------------------------------------------------
    subplot(2,2,3)

    plot(t_saved_min, V_rain_saved, 'b--', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Cumulative rainfall volume');
    hold on;

    plot(t_saved_min, V_infiltration_cum_from_f, 'g-', ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Cumulative infiltration from f');

    plot(t_saved_min, soil_storage_volume_m3, 'c-', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Soil storage volume I_t');

    plot(t_saved_min, remaining_capacity_volume_m3, 'k:', ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Remaining UZ capacity');

    if any(isfinite(V_out_saved))
        plot(t_saved_min, V_out_saved, 'r-', ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Cumulative outlet volume');
    end

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Volume [m^3]');
    title('Rainfall, Infiltration, Storage, and Outlet Volumes');
    legend('Location','best');

    % --------------------------------------------------------
    % Outlet discharge vs rainfall-minus-infiltration diagnostic
    % --------------------------------------------------------
    subplot(2,2,4)

    plot(t_model, q_model, 'r-', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Outlet discharge');
    hold on;

    plot(t_saved_min, Q_rain_minus_infiltration_m3s, 'b--', ...
        'LineWidth', 1.8, ...
        'DisplayName', 'Area rainfall - mean infiltration');

    plot(t_saved_min, Q_infiltration_equiv_m3s, 'g:', ...
        'LineWidth', 2.0, ...
        'DisplayName', 'Mean infiltration as Q-equivalent');

    yline(Q_steady_noinf, 'k--', ...
        'LineWidth', 1.0, ...
        'DisplayName', 'No-infiltration Q_{ss}');

    grid on; box on;
    xlabel('Time [min]');
    ylabel('Flow equivalent [m^3/s]');
    title('Outlet Flow and Infiltration-Adjusted Water Supply');
    legend('Location','best');

    sgtitle('Infiltration Evolution Diagnostics');

    exportgraphics(gcf, ...
        fullfile(figDir, 'TiltedPlane_Infiltration_Evolution_Diagnostics.png'), ...
        'Resolution', 300);

else

    warning('No infiltration map series found. Check that Maps.Hydro.I_t, Maps.Hydro.f, and Maps.Hydro.C are being saved.');

end

%% ============================================================
% PLOT 3: FINAL INFILTRATION MAPS
% ============================================================

if has_infiltration_series && ...
        exist('last_It_map','var') && ~isempty(last_It_map)

    figure('Color','w','Position',[200 200 1500 480]);

    subplot(1,3,1)
    imagesc_masked(last_It_map, active_mask);
    axis image off;
    colorbar;
    title('Final I_t [mm]');

    subplot(1,3,2)
    if exist('last_f_map','var') && ~isempty(last_f_map)
        imagesc_masked(last_f_map, active_mask);
        title('Final f [mm/h]');
    else
        imagesc(nan(size(last_It_map)));
        title('Final f unavailable');
    end
    axis image off;
    colorbar;

    subplot(1,3,3)
    if exist('last_C_map','var') && ~isempty(last_C_map)
        imagesc_masked(last_C_map, active_mask);
        title('Final C [mm/h]');
    else
        imagesc(nan(size(last_It_map)));
        title('Final C unavailable');
    end
    axis image off;
    colorbar;

    sgtitle('Final Infiltration State Maps');

    exportgraphics(gcf, ...
        fullfile(figDir, 'TiltedPlane_Final_Infiltration_Maps.png'), ...
        'Resolution', 300);

end

%% ================= EXPORT TABLES =================

Hydrograph_Comparison = table( ...
    t_bench, ...
    q_bench, ...
    q_model_i, ...
    q_model_i - q_bench, ...
    abs(q_model_i - q_bench), ...
    abs(q_model_i - q_bench) ./ max(abs(q_bench), 1e-6), ...
    'VariableNames', { ...
        'Time_min', ...
        'Analytical_NoInfiltration_m3s', ...
        'HydroPol2D_m3s', ...
        'Error_m3s', ...
        'AbsError_m3s', ...
        'RelError'});

hydrograph_file = fullfile(tabDir, 'TiltedPlane_Hydrograph_Comparison_WithInfiltration.csv');
writetable(Hydrograph_Comparison, hydrograph_file);

fprintf('Hydrograph comparison table exported to:\n%s\n', hydrograph_file);

if has_infiltration_series

    Infiltration_Evolution = table( ...
        t_saved_min, ...
        mean_depth_mm, ...
        mean_It_mm, ...
        median_It_mm, ...
        max_It_mm, ...
        mean_f_mmhr, ...
        max_f_mmhr, ...
        mean_C_mmhr, ...
        max_C_mmhr, ...
        100 * mean_storage_fill_frac, ...
        soil_storage_volume_m3, ...
        surface_storage_volume_m3, ...
        remaining_capacity_volume_m3, ...
        V_rain_saved, ...
        V_infiltration_cum_from_f, ...
        V_out_saved, ...
        Q_infiltration_equiv_m3s, ...
        Q_rain_minus_infiltration_m3s, ...
        q_model_saved, ...
        'VariableNames', { ...
            'Time_min', ...
            'Mean_SurfaceDepth_mm', ...
            'Mean_It_mm', ...
            'Median_It_mm', ...
            'Max_It_mm', ...
            'Mean_f_mmhr', ...
            'Max_f_mmhr', ...
            'Mean_C_mmhr', ...
            'Max_C_mmhr', ...
            'Mean_StorageFill_pct', ...
            'SoilStorageVolume_m3', ...
            'SurfaceStorageVolume_m3', ...
            'RemainingUZCapacity_m3', ...
            'RainVolume_m3', ...
            'CumulativeInfiltrationFromMeanF_m3', ...
            'OutletVolume_m3', ...
            'InfiltrationEquivalentFlow_m3s', ...
            'RainMinusInfiltrationEquivalentFlow_m3s', ...
            'OutletFlow_m3s'});

    infil_file = fullfile(tabDir, 'TiltedPlane_Infiltration_Evolution.csv');
    writetable(Infiltration_Evolution, infil_file);

    fprintf('Infiltration evolution table exported to:\n%s\n', infil_file);

end

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function x = to_double_cpu(x)
    x = double(gather(x));
end

function tf = exist_var(varname)
    tf = evalin('caller', sprintf('exist(''%s'',''var'')', varname)) == 1;
end

function Aout = expand_to_grid(Ain, targetSize)

    Ain = double(gather(Ain));

    if isscalar(Ain)
        Aout = Ain * ones(targetSize);
    elseif isequal(size(Ain), targetSize)
        Aout = Ain;
    else
        error('Cannot expand variable of size [%s] to target size [%s].', ...
            num2str(size(Ain)), num2str(targetSize));
    end

end

function UZ_capacity_mm = infer_uz_capacity_mm( ...
    targetSize, active_mask, Soil_Properties, DEM_for_slope, ...
    has_elevation, has_BC_States, has_GW_table)

    theta_r = expand_to_grid(Soil_Properties.theta_r, targetSize);
    theta_s = expand_to_grid(Soil_Properties.theta_sat, targetSize);

    if isfield(Soil_Properties, 'Soil_Depth') && ~isempty(Soil_Properties.Soil_Depth)
        soil_depth = expand_to_grid(Soil_Properties.Soil_Depth, targetSize); % [m]
    else
        warning('Soil_Properties.Soil_Depth not found. Assuming 1 m.');
        soil_depth = ones(targetSize);
    end

    elevation_grid = [];

    if has_elevation
        elevation_grid = evalin('caller', 'double(gather(elevation))');
    elseif ~isempty(DEM_for_slope)
        elevation_grid = DEM_for_slope;
    elseif evalin('caller', 'exist(''DEM_raster'',''var'') && isfield(DEM_raster,''Z'')')
        elevation_grid = evalin('caller', 'double(gather(DEM_raster.Z))');
    end

    h_t = [];

    if has_BC_States
        if evalin('caller', 'isfield(BC_States,''h_t'') && ~isempty(BC_States.h_t)')
            h_t = evalin('caller', 'double(gather(BC_States.h_t))');
        end
    end

    if isempty(h_t) && has_GW_table
        h_t = evalin('caller', 'double(gather(GW_table))');
    end

    if ~isempty(h_t) && ~isempty(elevation_grid)

        h_t = expand_to_grid(h_t, targetSize);
        elevation_grid = expand_to_grid(elevation_grid, targetSize);

        % Same logic as infiltration module:
        % GW_Depth = h_t - (elevation - Soil_Depth)
        % zwt      = Soil_Depth - GW_Depth
        GW_saturated_thickness = h_t - (elevation_grid - soil_depth);
        zwt = soil_depth - GW_saturated_thickness;

    else

        % Fallback for synthetic plane with water table at bedrock:
        % unsaturated thickness equals soil depth.
        zwt = soil_depth;

    end

    zwt = max(zwt, 0);

    UZ_capacity_mm = zwt .* max(theta_s - theta_r, 0) * 1000;

    UZ_capacity_mm(~active_mask) = NaN;

end

function [HydroMaps, local_i, store] = load_hydropol_hydro_stack( ...
    global_i, tempDir, saver_memory_maps, Maps_fallback)

    store = ceil(global_i / saver_memory_maps);
    local_i = global_i - (store - 1) * saver_memory_maps;

    HydroMaps = [];

    if ~isempty(tempDir)

        f1 = fullfile(tempDir, sprintf('save_map_hydro_%d.mat', store));
        f2 = fullfile(tempDir, sprintf('save_map_hydro_%d', store));

        if isfile(f1)
            S = load(f1, 'Maps');
            if isfield(S, 'Maps') && isfield(S.Maps, 'Hydro')
                HydroMaps = S.Maps.Hydro;
            end
        elseif isfile(f2)
            S = load(f2, 'Maps');
            if isfield(S, 'Maps') && isfield(S.Maps, 'Hydro')
                HydroMaps = S.Maps.Hydro;
            end
        end

    end

    % Fallback to current in-memory Maps, useful for the latest unflushed stack.
    if isempty(HydroMaps)
        if ~isempty(Maps_fallback) && isfield(Maps_fallback, 'Hydro')
            HydroMaps = Maps_fallback.Hydro;
        else
            error('Could not find Hydro maps for global record %d.', global_i);
        end
    end

end

function [A, ok] = read_hydro_field(HydroMaps, fieldName, local_i)

    A = [];
    ok = false;

    if ~isfield(HydroMaps, fieldName)
        return
    end

    stack = HydroMaps.(fieldName);

    if ndims(stack) >= 3
        if size(stack,3) >= local_i
            A = stack(:,:,local_i);
            ok = true;
        end
    elseif ismatrix(stack)
        if local_i == 1
            A = stack;
            ok = true;
        end
    end

end

function imagesc_masked(A, active_mask)

    A = double(gather(A));
    A(~active_mask) = NaN;
    imagesc(A);
    set(gca, 'YDir', 'normal');

end