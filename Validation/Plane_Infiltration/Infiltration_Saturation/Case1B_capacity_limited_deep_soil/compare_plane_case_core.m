function Results = compare_plane_case_core(Cfg)
%% ============================================================
% compare_plane_case_core
%
% Shared diagnostic core for plane infiltration cases.
%
% Updated for the new HydroPol2D infiltration closure:
%
%   Ksoil = wK * Ks + (1 - wK) * K(theta_top)
%
% Therefore the expected infiltration is generally a RANGE:
%
%   C_min = wK * Ks * G
%   C_max = Ks * G
%
% where:
%
%   G = dh_max / Ltop + 1
%
% For Case 1A, rainfall is less than C_min, so expected f = rainfall.
% For Case 1B, f should lie between C_min and C_max.
%
% Reads from current HydroPol2D workspace:
%   running_control.time_hydrograph
%   outlet_states.outlet_hydrograph
%   Maps.Hydro.I_t, Maps.Hydro.f, Maps.Hydro.C, Maps.Hydro.d
%
% Also reads:
%   CaseName/Static/Analytical_KinematicWave_Hydrograph.csv
%
% Outputs:
%   - hydrograph comparison figure
%   - infiltration evolution figure
%   - CSV diagnostic tables
% ============================================================

Cfg = fill_defaults(Cfg);

caseRootDir = fullfile(Cfg.rootDir, Cfg.caseFolder);
staticDir   = fullfile(caseRootDir, 'Static');

if evalin('base', 'exist(''Paths'',''var'')')
    Paths = evalin('base', 'Paths');
else
    Paths = struct();
end

if isfield(Paths, 'Results') && ~isempty(Paths.Results)
    outRoot = fullfile(Paths.Results, 'Plane_Infiltration_Diagnostics', Cfg.caseFolder);
else
    outRoot = fullfile(caseRootDir, 'Diagnostics');
end

figDir = fullfile(outRoot, 'Figures');
tabDir = fullfile(outRoot, 'Tables');

if ~isfolder(figDir); mkdir(figDir); end
if ~isfolder(tabDir); mkdir(tabDir); end

%% ================= DOMAIN CONSTANTS =================

A_domain = Cfg.Lx_m * Cfg.Ly_m;

Qss_noinf = Cfg.rain_mm_h / 1000 / 3600 * A_domain;

%% ================= EXPECTED INFILTRATION RANGE =================
%
% With the new closure:
%
%   Ksoil = wK*Ks + (1-wK)*Ktop
%
% Since 0 <= Ktop <= Ks:
%
%   wK*Ks <= Ksoil <= Ks
%
% Therefore:
%
%   C_low  = wK*Ks*G
%   C_high = Ks*G
%
% Then:
%
%   f_low  = min(P, C_low)
%   f_high = min(P, C_high)
%
% and discharge lies between:
%
%   Q_low  = flow using f_high
%   Q_high = flow using f_low
%
% because larger infiltration means smaller discharge.

if isfinite(Cfg.expected_f_low_mm_h) && isfinite(Cfg.expected_f_high_mm_h)

    expected_f_low_mm_h  = Cfg.expected_f_low_mm_h;
    expected_f_high_mm_h = Cfg.expected_f_high_mm_h;

elseif isfinite(Cfg.ksat_mm_h) && isfinite(Cfg.Ksurf_weight) && ...
       isfinite(Cfg.Ltop_m) && isfinite(Cfg.dh_max_m)

    G_expected = Cfg.dh_max_m / Cfg.Ltop_m + 1;

    C_low_mm_h  = Cfg.Ksurf_weight * Cfg.ksat_mm_h * G_expected;
    C_high_mm_h = Cfg.ksat_mm_h * G_expected;

    expected_f_low_mm_h  = min(Cfg.rain_mm_h, C_low_mm_h);
    expected_f_high_mm_h = min(Cfg.rain_mm_h, C_high_mm_h);

else

    expected_f_low_mm_h  = NaN;
    expected_f_high_mm_h = NaN;

end

% Make sure low <= high
if isfinite(expected_f_low_mm_h) && isfinite(expected_f_high_mm_h)
    ftmp1 = min(expected_f_low_mm_h, expected_f_high_mm_h);
    ftmp2 = max(expected_f_low_mm_h, expected_f_high_mm_h);
    expected_f_low_mm_h  = ftmp1;
    expected_f_high_mm_h = ftmp2;
end

if isfinite(expected_f_low_mm_h) && isfinite(expected_f_high_mm_h)

    expected_excess_high_mm_h = max(Cfg.rain_mm_h - expected_f_low_mm_h, 0);
    expected_excess_low_mm_h  = max(Cfg.rain_mm_h - expected_f_high_mm_h, 0);

    Q_expected_low_m3s = expected_excess_low_mm_h  / 1000 / 3600 * A_domain;
    Q_expected_high_m3s = expected_excess_high_mm_h / 1000 / 3600 * A_domain;

else

    Q_expected_low_m3s  = NaN;
    Q_expected_high_m3s = NaN;

end

%% ================= LOAD MODEL HYDROGRAPH =================

[t_model, q_model] = get_model_hydrograph();

if isempty(t_model)
    error('Could not read model hydrograph from running_control/outlet_states.');
end

%% ================= LOAD / BUILD ANALYTICAL BENCHMARK =================

benchFile = fullfile(staticDir, 'Analytical_KinematicWave_Hydrograph.csv');

if isfile(benchFile)

    Bench = readtable(benchFile);

    if ismember('time_min', Bench.Properties.VariableNames)
        t_bench = Bench.time_min;
    else
        error('Benchmark file exists but does not contain time_min.');
    end

    qVar = find(contains(lower(Bench.Properties.VariableNames), 'q'), 1, 'first');

    if isempty(qVar)
        error('Benchmark file exists but no Q column was found.');
    end

    q_bench = Bench.(Bench.Properties.VariableNames{qVar});

else

    warning('Benchmark CSV not found. Building analytical benchmark internally.');
    [t_bench, q_bench] = build_noinf_benchmark(Cfg);

end

t_bench = double(t_bench(:));
q_bench = double(q_bench(:));

%% ================= INTERPOLATE MODEL TO BENCHMARK =================

if numel(t_model) >= 2
    q_model_i = interp1(t_model, q_model, t_bench, 'linear', NaN);
else
    q_model_i = NaN(size(t_bench));
end

valid_eval = isfinite(q_model_i) & isfinite(q_bench);

if nnz(valid_eval) >= 2
    NSE = compute_nse(q_bench(valid_eval), q_model_i(valid_eval));
    RMSE = sqrt(mean((q_model_i(valid_eval) - q_bench(valid_eval)).^2, 'omitnan'));
    BIAS = mean(q_model_i(valid_eval) - q_bench(valid_eval), 'omitnan');
else
    NSE = NaN;
    RMSE = NaN;
    BIAS = NaN;
end

%% ================= HYDROGRAPH PLOT =================

fig = figure('Color','w','Position',[100 100 1050 560]);

plot(t_bench, q_bench, 'k-', ...
    'LineWidth', 2.2, ...
    'DisplayName', 'Analytical no-infiltration');
hold on;

plot(t_model, q_model, 'r-', ...
    'LineWidth', 1.8, ...
    'DisplayName', 'HydroPol2D');

plot(t_bench, q_model_i, 'ro', ...
    'MarkerSize', 4, ...
    'DisplayName', 'HydroPol2D interpolated');

yline(Qss_noinf, 'k--', ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf('No-inf Q_{ss}=%.3f m^3/s', Qss_noinf));

% Expected range using new Ksurf interface closure
t_end_plot = max([t_bench(:); t_model(:)], [], 'omitnan');

if isfinite(Q_expected_low_m3s) && isfinite(Q_expected_high_m3s)

    if abs(Q_expected_high_m3s - Q_expected_low_m3s) < 1e-9

        yline(Q_expected_low_m3s, 'b--', ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('Expected Q=%.3f m^3/s', Q_expected_low_m3s));

    else

        patch([0 t_end_plot t_end_plot 0], ...
              [Q_expected_low_m3s Q_expected_low_m3s Q_expected_high_m3s Q_expected_high_m3s], ...
              [0.80 0.88 1.00], ...
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.35, ...
              'DisplayName', sprintf('Expected Q range %.3f–%.3f m^3/s', ...
                  Q_expected_low_m3s, Q_expected_high_m3s));

        yline(Q_expected_low_m3s, 'b--', ...
            'LineWidth', 1.1, ...
            'HandleVisibility', 'off');

        yline(Q_expected_high_m3s, 'b--', ...
            'LineWidth', 1.1, ...
            'HandleVisibility', 'off');

    end

end

if isfinite(Cfg.expected_fill_time_h)
    xline(Cfg.expected_fill_time_h * 60, 'm:', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Approx. fill time %.2f h', Cfg.expected_fill_time_h));
end

grid on; box on;

xlabel('Time [min]');
ylabel('Outlet discharge [m^3/s]');

title(sprintf('%s: Outlet Hydrograph, NSE vs no-inf = %.3f', ...
    strrep(Cfg.caseFolder, '_', '\_'), NSE));

legend('Location','best');

exportgraphics(fig, fullfile(figDir, '01_hydrograph_comparison.png'), 'Resolution', 300);

%% ================= HYDROGRAPH TABLE =================

HydrographTable = table( ...
    t_bench, ...
    q_bench, ...
    q_model_i, ...
    q_model_i - q_bench, ...
    abs(q_model_i - q_bench), ...
    'VariableNames', { ...
        'Time_min', ...
        'Analytical_NoInfiltration_m3s', ...
        'HydroPol2D_Interpolated_m3s', ...
        'SignedError_m3s', ...
        'AbsError_m3s'});

writetable(HydrographTable, fullfile(tabDir, 'hydrograph_comparison.csv'));

%% ================= OUTLET VOLUME =================

t_model_sec = t_model * 60;

if numel(t_model_sec) >= 2
    dt_sec = [0; diff(t_model_sec)];
    V_out_cum = cumsum(q_model .* dt_sec);
else
    V_out_cum = zeros(size(t_model_sec));
end

V_rain_cum = Cfg.rain_mm_h / 1000 / 3600 .* ...
    (min(max(t_model, 0), Cfg.rain_duration_min) * 60) .* A_domain;

%% ================= INFILTRATION DIAGNOSTICS =================

hasInfiltrationSeries = false;
InfilTable = table();

if Cfg.flag_infiltration == 1

    It_stack = collect_hydro_field('I_t');
    f_stack  = collect_hydro_field('f');
    C_stack  = collect_hydro_field('C');
    d_stack  = collect_hydro_field('d');

    if ~isempty(It_stack) || ~isempty(f_stack) || ~isempty(C_stack)

        [activeMask, ~] = get_active_mask_from_workspace(It_stack, f_stack, C_stack, Cfg);

        nSeriesRaw = max([stack_len(It_stack), stack_len(f_stack), stack_len(C_stack), stack_len(d_stack)]);

        [t_saved_min, nSeries] = get_saved_times(nSeriesRaw, Cfg);

        mean_It_mm = NaN(nSeries,1);
        max_It_mm  = NaN(nSeries,1);

        mean_f_mm_h = NaN(nSeries,1);
        max_f_mm_h  = NaN(nSeries,1);

        mean_C_mm_h = NaN(nSeries,1);
        max_C_mm_h  = NaN(nSeries,1);

        mean_d_mm = NaN(nSeries,1);

        for ii = 1:nSeries

            if stack_len(It_stack) >= ii
                A = to_double(It_stack(:,:,ii));
                A(~isfinite(A)) = NaN;
                A = max(A, 0);
                mean_It_mm(ii) = mean(A(activeMask), 'omitnan');
                max_It_mm(ii)  = max(A(activeMask), [], 'omitnan');
            end

            if stack_len(f_stack) >= ii
                A = to_double(f_stack(:,:,ii));
                A(~isfinite(A)) = NaN;
                A = max(A, 0);
                mean_f_mm_h(ii) = mean(A(activeMask), 'omitnan');
                max_f_mm_h(ii)  = max(A(activeMask), [], 'omitnan');
            end

            if stack_len(C_stack) >= ii
                A = to_double(C_stack(:,:,ii));
                A(~isfinite(A)) = NaN;
                A = max(A, 0);
                mean_C_mm_h(ii) = mean(A(activeMask), 'omitnan');
                max_C_mm_h(ii)  = max(A(activeMask), [], 'omitnan');
            end

            if stack_len(d_stack) >= ii
                A = to_double(d_stack(:,:,ii));
                A(~isfinite(A)) = NaN;
                A = max(A, 0);
                mean_d_mm(ii) = mean(A(activeMask), 'omitnan');
            end

        end

        %% ----- Storage convention / event-fill diagnostic -----
        %
        % HydroPol2D versions may store I_t as either:
        %   1) storage above residual:
        %        I_t = (theta - theta_r) * zwt * 1000
        %
        % or
        %   2) total water storage:
        %        I_t = theta * zwt * 1000
        %
        % The infiltration module currently treats I_t as storage above
        % residual, but the diagnostic auto-detects the plotted convention
        % using the first saved I_t value.

        first_valid_it = mean_It_mm(find(isfinite(mean_It_mm), 1, 'first'));

        expected_initial_effective_mm = ...
            max(Cfg.theta_i - Cfg.theta_r, 0) * Cfg.zwt0_m * 1000;

        expected_initial_total_mm = ...
            Cfg.theta_i * Cfg.zwt0_m * 1000;

        if isempty(first_valid_it) || ~isfinite(first_valid_it)

            It_convention = "unknown";
            max_storage_mm = (Cfg.theta_sat - Cfg.theta_r) * Cfg.zwt0_m * 1000;

        else

            diff_effective = abs(first_valid_it - expected_initial_effective_mm);
            diff_total     = abs(first_valid_it - expected_initial_total_mm);

            if diff_total < diff_effective
                It_convention = "total_water_storage";
                max_storage_mm = Cfg.theta_sat * Cfg.zwt0_m * 1000;
            else
                It_convention = "above_residual_storage";
                max_storage_mm = (Cfg.theta_sat - Cfg.theta_r) * Cfg.zwt0_m * 1000;
            end

        end

        % Use the first valid mapped I_t as the event baseline.
        It0_obs_mm = first_valid_it;

        if isempty(It0_obs_mm) || ~isfinite(It0_obs_mm)
            It0_obs_mm = NaN;
        end

        event_storage_fill_frac = ...
            (mean_It_mm - It0_obs_mm) ./ max(max_storage_mm - It0_obs_mm, eps);

        event_storage_fill_frac = min(max(event_storage_fill_frac, 0), 1);

        remaining_storage_mm = max(max_storage_mm - mean_It_mm, 0);

        %% ----- Cumulative infiltration from mean f -----

        t_saved_sec = t_saved_min * 60;

        if numel(t_saved_sec) >= 2
            dt_saved_sec = [0; diff(t_saved_sec)];
        else
            dt_saved_sec = zeros(size(t_saved_sec));
        end

        Q_infiltration_equiv_m3s = mean_f_mm_h / 1000 / 3600 * A_domain;

        V_infiltration_cum = cumsum(Q_infiltration_equiv_m3s .* dt_saved_sec);

        Q_rain_m3s = Cfg.rain_mm_h / 1000 / 3600 * A_domain;

        Q_rain_minus_infiltration_m3s = ...
            max(Q_rain_m3s - Q_infiltration_equiv_m3s, 0);

        InfilTable = table( ...
            t_saved_min, ...
            mean_d_mm, ...
            mean_It_mm, ...
            max_It_mm, ...
            mean_f_mm_h, ...
            max_f_mm_h, ...
            mean_C_mm_h, ...
            max_C_mm_h, ...
            100 * event_storage_fill_frac, ...
            remaining_storage_mm, ...
            Q_infiltration_equiv_m3s, ...
            Q_rain_minus_infiltration_m3s, ...
            V_infiltration_cum, ...
            repmat(It_convention, nSeries, 1), ...
            'VariableNames', { ...
                'Time_min', ...
                'Mean_SurfaceDepth_mm', ...
                'Mean_It_mm', ...
                'Max_It_mm', ...
                'Mean_f_mm_h', ...
                'Max_f_mm_h', ...
                'Mean_C_mm_h', ...
                'Max_C_mm_h', ...
                'Event_StorageFill_pct', ...
                'Mean_RemainingStorage_mm', ...
                'InfiltrationEquivalentFlow_m3s', ...
                'RainMinusInfiltrationEquivalentFlow_m3s', ...
                'CumulativeInfiltration_m3', ...
                'Detected_It_Convention'});

        writetable(InfilTable, fullfile(tabDir, 'infiltration_evolution.csv'));

        hasInfiltrationSeries = true;

        %% ---------- Infiltration plot ----------

        fig2 = figure('Color','w','Position',[150 150 1300 850]);

        subplot(2,2,1)

        plot(t_saved_min, mean_f_mm_h, 'b-o', ...
            'LineWidth', 1.7, ...
            'MarkerSize', 4, ...
            'DisplayName', 'Mean f');
        hold on;

        plot(t_saved_min, mean_C_mm_h, 'm-s', ...
            'LineWidth', 1.7, ...
            'MarkerSize', 4, ...
            'DisplayName', 'Mean C');

        yline(Cfg.rain_mm_h, 'k--', ...
            'LineWidth', 1.2, ...
            'DisplayName', 'Rainfall');

        if isfinite(expected_f_low_mm_h) && isfinite(expected_f_high_mm_h)

            if abs(expected_f_high_mm_h - expected_f_low_mm_h) < 1e-9

                yline(expected_f_low_mm_h, 'b--', ...
                    'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Expected f=%.2f mm/h', expected_f_low_mm_h));

            else

                yline(expected_f_low_mm_h, 'b:', ...
                    'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Expected f low=%.2f', expected_f_low_mm_h));

                yline(expected_f_high_mm_h, 'b--', ...
                    'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Expected f high=%.2f', expected_f_high_mm_h));

            end

        end

        grid on; box on;
        xlabel('Time [min]');
        ylabel('Rate [mm/h]');
        title('Infiltration Rate and Capacity');
        legend('Location','best');

        subplot(2,2,2)

        yyaxis left
        plot(t_saved_min, mean_It_mm, 'g-o', ...
            'LineWidth', 1.7, ...
            'MarkerSize', 4);
        ylabel('Mean I_t [mm]');

        yyaxis right
        plot(t_saved_min, 100 * event_storage_fill_frac, 'k-s', ...
            'LineWidth', 1.7, ...
            'MarkerSize', 4);
        ylabel('Event storage fill [%]');

        if isfinite(Cfg.expected_fill_time_h)
            xline(Cfg.expected_fill_time_h * 60, 'm:', ...
                'LineWidth', 1.4);
        end

        grid on; box on;
        xlabel('Time [min]');
        title(sprintf('Subsurface Storage Evolution (%s)', It_convention));

        subplot(2,2,3)

        plot(t_saved_min, V_infiltration_cum, 'g-', ...
            'LineWidth', 2.0, ...
            'DisplayName', 'Cumulative infiltration');
        hold on;

        plot(t_model, V_rain_cum, 'b--', ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Cumulative rainfall');

        plot(t_model, V_out_cum, 'r-', ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Cumulative outlet');

        grid on; box on;
        xlabel('Time [min]');
        ylabel('Volume [m^3]');
        title('Cumulative Volumes');
        legend('Location','best');

        subplot(2,2,4)

        plot(t_model, q_model, 'r-', ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Outlet Q');
        hold on;

        plot(t_saved_min, Q_rain_minus_infiltration_m3s, 'b--', ...
            'LineWidth', 1.8, ...
            'DisplayName', 'Rain - mean infiltration');

        plot(t_saved_min, Q_infiltration_equiv_m3s, 'g:', ...
            'LineWidth', 2.0, ...
            'DisplayName', 'Infiltration equivalent Q');

        yline(Qss_noinf, 'k--', ...
            'LineWidth', 1.1, ...
            'DisplayName', 'No-inf Qss');

        if isfinite(Q_expected_low_m3s) && isfinite(Q_expected_high_m3s)

            if abs(Q_expected_high_m3s - Q_expected_low_m3s) < 1e-9

                yline(Q_expected_low_m3s, 'b:', ...
                    'LineWidth', 1.4, ...
                    'DisplayName', 'Expected Q');

            else

                yline(Q_expected_low_m3s, 'b:', ...
                    'LineWidth', 1.2, ...
                    'DisplayName', 'Expected Q low');

                yline(Q_expected_high_m3s, 'b--', ...
                    'LineWidth', 1.2, ...
                    'DisplayName', 'Expected Q high');

            end

        end

        grid on; box on;
        xlabel('Time [min]');
        ylabel('Flow [m^3/s]');
        title('Flow Equivalent Diagnostics');
        legend('Location','best');

        sgtitle(strrep(Cfg.caseFolder, '_', '\_'));

        exportgraphics(fig2, fullfile(figDir, '02_infiltration_evolution.png'), 'Resolution', 300);

    else

        warning('No infiltration map stacks found for %s.', Cfg.caseFolder);

    end

end

%% ================= SUMMARY =================

Summary = table();

Summary.CaseFolder = string(Cfg.caseFolder);
Summary.Rain_mm_h = Cfg.rain_mm_h;
Summary.FlagInfiltration = Cfg.flag_infiltration;

Summary.Expected_f_low_mm_h = expected_f_low_mm_h;
Summary.Expected_f_high_mm_h = expected_f_high_mm_h;
Summary.Expected_Q_low_m3s = Q_expected_low_m3s;
Summary.Expected_Q_high_m3s = Q_expected_high_m3s;

Summary.NoInfiltration_Qss_m3s = Qss_noinf;
Summary.NSE_vs_NoInfiltration = NSE;
Summary.RMSE_vs_NoInfiltration_m3s = RMSE;
Summary.BIAS_vs_NoInfiltration_m3s = BIAS;
Summary.Final_Model_Q_m3s = q_model(end);

if hasInfiltrationSeries
    Summary.Final_Mean_f_mm_h = InfilTable.Mean_f_mm_h(end);
    Summary.Final_Mean_C_mm_h = InfilTable.Mean_C_mm_h(end);
    Summary.Final_Mean_It_mm = InfilTable.Mean_It_mm(end);
    Summary.Final_Event_StorageFill_pct = InfilTable.Event_StorageFill_pct(end);
    Summary.Detected_It_Convention = string(InfilTable.Detected_It_Convention(end));
else
    Summary.Final_Mean_f_mm_h = NaN;
    Summary.Final_Mean_C_mm_h = NaN;
    Summary.Final_Mean_It_mm = NaN;
    Summary.Final_Event_StorageFill_pct = NaN;
    Summary.Detected_It_Convention = "unknown";
end

writetable(Summary, fullfile(tabDir, 'case_summary.csv'));

Results = struct();
Results.Cfg = Cfg;
Results.Summary = Summary;
Results.HydrographTable = HydrographTable;
Results.InfilTable = InfilTable;
Results.figDir = figDir;
Results.tabDir = tabDir;

fprintf('\n============================================================\n');
fprintf('Diagnostics complete for %s\n', Cfg.caseFolder);
fprintf('Figures written to:\n%s\n', figDir);
fprintf('Tables written to:\n%s\n', tabDir);
fprintf('============================================================\n');

end

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function Cfg = fill_defaults(Cfg)

if ~isfield(Cfg, 'rootDir')
    Cfg.rootDir = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Plane_Infiltration/Infiltration_Saturation';
end

if ~isfield(Cfg, 'Lx_m'); Cfg.Lx_m = 2000; end
if ~isfield(Cfg, 'Ly_m'); Cfg.Ly_m = 1000; end
if ~isfield(Cfg, 'dx_m'); Cfg.dx_m = 20; end
if ~isfield(Cfg, 'S0'); Cfg.S0 = 0.01; end
if ~isfield(Cfg, 'n_manning'); Cfg.n_manning = 0.015; end

if ~isfield(Cfg, 'theta_sat'); Cfg.theta_sat = 0.385; end
if ~isfield(Cfg, 'theta_r'); Cfg.theta_r = 0.068; end
if ~isfield(Cfg, 'theta_i'); Cfg.theta_i = 0.0997; end

if ~isfield(Cfg, 'ksat_mm_h'); Cfg.ksat_mm_h = NaN; end
if ~isfield(Cfg, 'Ltop_m'); Cfg.Ltop_m = 0.05; end
if ~isfield(Cfg, 'dh_max_m'); Cfg.dh_max_m = 1.0; end
if ~isfield(Cfg, 'Ksurf_weight'); Cfg.Ksurf_weight = 0.5; end

if ~isfield(Cfg, 'expected_f_low_mm_h'); Cfg.expected_f_low_mm_h = NaN; end
if ~isfield(Cfg, 'expected_f_high_mm_h'); Cfg.expected_f_high_mm_h = NaN; end
if ~isfield(Cfg, 'expected_fill_time_h'); Cfg.expected_fill_time_h = NaN; end

end

function [t_model, q_model] = get_model_hydrograph()

if ~evalin('base', 'exist(''running_control'',''var'')') || ...
   ~evalin('base', 'exist(''outlet_states'',''var'')')
    t_model = [];
    q_model = [];
    return
end

running_control = evalin('base', 'running_control');
outlet_states   = evalin('base', 'outlet_states');

t_raw = running_control.time_hydrograph(:);

if isdatetime(t_raw)
    t_model = minutes(t_raw - t_raw(1));
elseif isduration(t_raw)
    t_model = minutes(t_raw);
else
    t_model = to_double(t_raw);
end

q_model = to_double(outlet_states.outlet_hydrograph(:));

valid = isfinite(t_model) & isfinite(q_model);

t_model = t_model(valid);
q_model = q_model(valid);

[t_model, sort_id] = sort(t_model);
q_model = q_model(sort_id);

[t_model, unique_id] = unique(t_model, 'stable');
q_model = q_model(unique_id);

q_model = max(q_model, 0);

end

function [t_bench, q_bench] = build_noinf_benchmark(Cfg)

r = Cfg.rain_mm_h / 1000 / 3600;

alpha = sqrt(Cfg.S0) / Cfg.n_manning;
m_kw = 5/3;

Lx = Cfg.Lx_m;
Ly = Cfg.Ly_m;

q_unit_steady = r * Lx;

dt_min = 10;
t_bench = (0:dt_min:Cfg.sim_duration_min)';

t_sec = t_bench * 60;

q_unit = alpha .* (r .* t_sec).^m_kw;
q_unit = min(q_unit, q_unit_steady);
q_unit(t_sec == 0) = 0;

q_bench = Ly .* q_unit;

end

function NSE = compute_nse(obs, sim)

den = sum((obs - mean(obs)).^2);

if den > 0
    NSE = 1 - sum((sim - obs).^2) / den;
else
    NSE = NaN;
end

end

function A = collect_hydro_field(fieldName)

A = [];

% Prefer saved files if available.
if evalin('base', 'exist(''Paths'',''var'')')
    Paths = evalin('base', 'Paths');

    if isfield(Paths, 'Temp') && ~isempty(Paths.Temp) && isfolder(Paths.Temp)

        files = dir(fullfile(Paths.Temp, 'save_map_hydro_*'));

        if ~isempty(files)

            fileNums = NaN(numel(files),1);

            for i = 1:numel(files)
                tok = regexp(files(i).name, 'save_map_hydro_(\d+)', 'tokens', 'once');
                if ~isempty(tok)
                    fileNums(i) = str2double(tok{1});
                end
            end

            [~, order] = sort(fileNums);
            files = files(order);

            for i = 1:numel(files)

                fpath = fullfile(files(i).folder, files(i).name);

                try
                    S = load(fpath);

                    if isfield(S, 'Maps') && ...
                       isfield(S.Maps, 'Hydro') && ...
                       isfield(S.Maps.Hydro, fieldName)

                        B = S.Maps.Hydro.(fieldName);

                        if ndims(B) == 2
                            B = reshape(B, size(B,1), size(B,2), 1);
                        end

                        A = cat(3, A, B);

                    end

                catch
                    % Ignore unreadable files.
                end

            end

        end

    end
end

% Fallback to in-memory Maps.
if isempty(A) && evalin('base', 'exist(''Maps'',''var'')')

    Maps = evalin('base', 'Maps');

    if isfield(Maps, 'Hydro') && isfield(Maps.Hydro, fieldName)

        A = Maps.Hydro.(fieldName);

        if ndims(A) == 2
            A = reshape(A, size(A,1), size(A,2), 1);
        end

    end

end

end

function n = stack_len(A)

if isempty(A)
    n = 0;
elseif ndims(A) == 2
    n = 1;
else
    n = size(A,3);
end

end

function [t_saved_min, nSeries] = get_saved_times(nSeriesRaw, Cfg)

t_saved_min = [];

if evalin('base', 'exist(''running_control'',''var'')')

    running_control = evalin('base', 'running_control');

    if isfield(running_control, 'time_records') && ~isempty(running_control.time_records)

        t_raw = running_control.time_records(:);

        if isdatetime(t_raw)
            t_saved_min = minutes(t_raw - t_raw(1));
        elseif isduration(t_raw)
            t_saved_min = minutes(t_raw);
        else
            t_saved_min = to_double(t_raw);
        end

        t_saved_min = t_saved_min(isfinite(t_saved_min));
        t_saved_min = t_saved_min(t_saved_min >= 0 & t_saved_min <= Cfg.sim_duration_min + 1e-9);

    end

end

if isempty(t_saved_min)
    nSeries = nSeriesRaw;
    t_saved_min = linspace(0, Cfg.sim_duration_min, nSeries)';
else
    nSeries = min(nSeriesRaw, numel(t_saved_min));
    t_saved_min = t_saved_min(1:nSeries);
end

end

function [activeMask, A_cell] = get_active_mask_from_workspace(It_stack, f_stack, C_stack, Cfg)

targetSize = [];

if ~isempty(It_stack)
    targetSize = size(It_stack(:,:,1));
elseif ~isempty(f_stack)
    targetSize = size(f_stack(:,:,1));
elseif ~isempty(C_stack)
    targetSize = size(C_stack(:,:,1));
else
    targetSize = [round(Cfg.Lx_m / Cfg.dx_m), round(Cfg.Ly_m / Cfg.dx_m)];
end

if evalin('base', 'exist(''idx_nan'',''var'')')

    idx_nan = evalin('base', 'idx_nan');
    activeMask = ~logical(idx_nan);

elseif evalin('base', 'exist(''DEM_raster'',''var'')')

    DEM_raster = evalin('base', 'DEM_raster');

    if isfield(DEM_raster, 'Z')
        activeMask = isfinite(to_double(DEM_raster.Z));
    else
        activeMask = true(targetSize);
    end

elseif evalin('base', 'exist(''elevation'',''var'')')

    elevation = evalin('base', 'elevation');
    activeMask = isfinite(to_double(elevation));

else

    activeMask = true(targetSize);

end

if ~isequal(size(activeMask), targetSize)
    activeMask = true(targetSize);
end

if evalin('base', 'exist(''Wshed_Properties'',''var'')')

    Wshed_Properties = evalin('base', 'Wshed_Properties');

    if isfield(Wshed_Properties, 'Resolution')
        dx = Wshed_Properties.Resolution;
    else
        dx = Cfg.dx_m;
    end

else

    dx = Cfg.dx_m;

end

A_cell = dx^2;

end

function A = to_double(A)

try
    A = gather(A);
catch
end

A = double(A);

end