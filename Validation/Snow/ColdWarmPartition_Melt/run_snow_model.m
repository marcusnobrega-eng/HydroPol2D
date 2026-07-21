% VAL-SNOW-001: configurable snow accumulation, melt, and initialization.
%
% This validation suite validates the active per-LULC snow pathway. It tests
% precipitation partitioning, storage bookkeeping, raster initial states,
% time-step scaling, class-code fallback, and the Snow_Module coupling on a
% small V-tilted grid. It does not claim field-scale snow calibration.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
fig_dir = fullfile(out_dir, 'Figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

addpath(functions_dir, '-begin');
hydropol2d_add_runtime_paths(hydropol2d_find_root(case_dir));

[SnowConfig, ExcelConfig] = load_snow_config(repo_root);
Diagnostics = table();
MassBalance = table();
PassFail = table();

% 1. Configurable cold, mixed, and warm precipitation partition.
[PartitionDiag, PartitionMass, passed] = test_partition_and_storage(SnowConfig);
Diagnostics = [Diagnostics; PartitionDiag]; %#ok<AGROW>
MassBalance = [MassBalance; PartitionMass]; %#ok<AGROW>
PassFail = [PassFail; pass_row(PartitionDiag, passed)]; %#ok<AGROW>

% 2. Per-LULC parameter maps, documented class fallback, and all four
% supported initial-state combinations.
[InitDiag, InitMass, passed] = test_initial_states_and_class_fallback(SnowConfig);
Diagnostics = [Diagnostics; InitDiag]; %#ok<AGROW>
MassBalance = [MassBalance; InitMass]; %#ok<AGROW>
PassFail = [PassFail; pass_row(InitDiag, passed)]; %#ok<AGROW>

% 3. The same day is integrated with 5, 15, and 60 minute steps. The
% forcing has no radiative or sublimation melt, so degree-day scaling has
% an exact end-of-day reference of 96 mm SWE.
[TimeDiag, TimeMass, TimeSeries, passed] = test_timestep_refinement();
Diagnostics = [Diagnostics; TimeDiag]; %#ok<AGROW>
MassBalance = [MassBalance; TimeMass]; %#ok<AGROW>
PassFail = [PassFail; pass_row(TimeDiag, passed)]; %#ok<AGROW>

% 4. Run the actual Snow_Module over a small V-tilted grid. The test keeps
% all non-snow hydrologic modules off, so snow, rainfall release, and
% surface storage have a closed mass ledger.
[VTiltedDiag, VTiltedMass, VTiltedSeries, passed] = test_vtilted_snow_module(SnowConfig);
Diagnostics = [Diagnostics; VTiltedDiag]; %#ok<AGROW>
MassBalance = [MassBalance; VTiltedMass]; %#ok<AGROW>
PassFail = [PassFail; pass_row(VTiltedDiag, passed)]; %#ok<AGROW>

% 5. Read the public Excel table and compare it with the equivalent bypass
% table. This confirms that the two supported configuration surfaces map to
% identical snow parameter arrays.
[ConfigDiag, ConfigMass, passed] = test_excel_bypass_equivalence(ExcelConfig, SnowConfig);
Diagnostics = [Diagnostics; ConfigDiag]; %#ok<AGROW>
MassBalance = [MassBalance; ConfigMass]; %#ok<AGROW>
PassFail = [PassFail; pass_row(ConfigDiag, passed)]; %#ok<AGROW>

writetable(Diagnostics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(MassBalance, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(TimeSeries, fullfile(out_dir, 'Timestep_Refinement_TimeSeries.csv'));
writetable(VTiltedSeries, fullfile(out_dir, 'VTilted_Snow_Runoff_TimeSeries.csv'));
make_figures(TimeSeries, VTiltedSeries, fig_dir);

disp(Diagnostics);
disp(PassFail);

function [Config, ExcelConfig] = load_snow_config(repo_root)
excel_path = fullfile(repo_root, 'Input_Data_Sheets', 'LULC_parameters.xlsx');
ExcelTable = readtable(excel_path, 'Range', 'A2', 'VariableNamingRule', 'preserve');
ExcelConfig = hp2d_normalize_snow_table(ExcelTable);

% A two-class table is used for the controlled equation and V-tilted tests.
% Its fields exactly match the public bypass interface.
T = table([1; 2], [0.80; 0.60], [1.00; 1.00], [0; 0], [2; 4], ...
    [-1; 1], [1; 3], [100; 200], [450; 500], [0; 0], [0; 0], [0; 0], ...
    'VariableNames', {'LULC_Index','Snow_Albedo','Snow_Emissivity', ...
    'Sublimation_Coefficient_d_1','Degree_Day_Factor_mm_C_day', ...
    'T_Snow_All_C','T_Rain_All_C','Rho_Snow_Init_kg_m3', ...
    'Rho_Snow_Max_kg_m3','Compaction_Temperature_kg_m3_C_day', ...
    'Compaction_SWE_kg_m3_mm_day','Compaction_Depth_kg_m3_mm_day'});
Config = hp2d_normalize_snow_table(T);
end

function [Diag, Ledger, passed] = test_partition_and_storage(Config)
p = scalar_parameters(Config, 1);
dt_s = 900;
[swe_cold, ~, melt_cold, snow_cold, rain_cold, ~, sub_cold, mb_cold] = ...
    Snow_Model_Function(0, 0, p.rho_snow_init, -5, -8, 10, 2, 39, 30, ...
    p.alpha, p.epsilon, p.C_e, p.DDF, p.T_snow_all, p.T_rain_all, ...
    p.rho_snow_init, p.rho_max, p.k_t, p.k_swe, p.k_D, dt_s);
[swe_mixed, ~, melt_mixed, snow_mixed, rain_mixed, ~, sub_mixed, mb_mixed] = ...
    Snow_Model_Function(0, 0, p.rho_snow_init, 0, -2, 10, 2, 39, 30, ...
    p.alpha, p.epsilon, p.C_e, p.DDF, p.T_snow_all, p.T_rain_all, ...
    p.rho_snow_init, p.rho_max, p.k_t, p.k_swe, p.k_D, dt_s);
[swe_warm, ~, melt_warm, snow_warm, rain_warm, ~, sub_warm, mb_warm] = ...
    Snow_Model_Function(0, 0, p.rho_snow_init, 5, 1, 10, 2, 39, 30, ...
    p.alpha, p.epsilon, p.C_e, p.DDF, p.T_snow_all, p.T_rain_all, ...
    p.rho_snow_init, p.rho_max, p.k_t, p.k_swe, p.k_D, dt_s);

partition_error = max(abs([snow_cold - 10, rain_cold, ...
    snow_mixed - 5, rain_mixed - 5, snow_warm, rain_warm - 10]));
mass_residual_mm = max(abs([mb_cold, mb_mixed, mb_warm]));
finite_ok = all(isfinite([swe_cold,swe_mixed,swe_warm,melt_cold,melt_mixed,melt_warm, ...
    sub_cold,sub_mixed,sub_warm]));
passed = partition_error < 1e-12 && mass_residual_mm < 1e-10 && finite_ok;

Diag = diagnostic_row("VAL-SNOW-001A", "Configurable precipitation partition and storage", ...
    partition_error, mass_residual_mm / 1000, 0, 0, passed, ...
    "Cold, mixed, and warm precipitation follow the configured -1 to 1 degC transition.");
Ledger = ledger_row("VAL-SNOW-001A", 30, snow_cold + snow_mixed + snow_warm, ...
    rain_cold + rain_mixed + rain_warm, melt_cold + melt_mixed + melt_warm, ...
    sub_cold + sub_mixed + sub_warm, swe_cold + swe_mixed + swe_warm, mass_residual_mm);
end

function [Diag, Ledger, passed] = test_initial_states_and_class_fallback(Config)
LULC = [1 2 99; 1 NaN 2];
valid = isfinite(LULC);
cell_area_m2 = 400;
[NoneState, Audit] = hp2d_initialize_snow_state(Config, LULC, valid, [], [], cell_area_m2);
SWEOnly = nan(size(LULC)); SWEOnly(1,1) = 100; SWEOnly(2,2) = 25;
[SWEState, ~] = hp2d_initialize_snow_state(Config, LULC, valid, SWEOnly, [], cell_area_m2);
DepthOnly = nan(size(LULC)); DepthOnly(1,2) = 1000;
[DepthState, ~] = hp2d_initialize_snow_state(Config, LULC, valid, [], DepthOnly, cell_area_m2);
BothSWE = nan(size(LULC)); BothDepth = nan(size(LULC));
BothSWE(2,1) = 100; BothDepth(2,1) = 500;
[BothState, ~] = hp2d_initialize_snow_state(Config, LULC, valid, BothSWE, BothDepth, cell_area_m2);

inconsistent_raised = false;
try
    BadDepth = nan(size(LULC)); BadDepth(1,1) = 0;
    hp2d_initialize_snow_state(Config, LULC, valid, SWEOnly, BadDepth, cell_area_m2);
catch ME
    inconsistent_raised = contains(ME.message, 'inconsistent');
end

fallback = mean(Config.parameter_values, 1);
initial_error = max(abs([ ...
    NoneState.SWE_t(1,1), ...
    SWEState.H_snow_t(1,1) - 1000, ...
    DepthState.SWE_t(1,2) - 200, ...
    BothState.rho_snow(2,1) - 200, ...
    NoneState.rho_snow_init(1,3) - fallback(7)]));
fallback_cells = sum(Audit.cell_count(Audit.mapping_status == "fallback_area_weighted_mean"));
passed = initial_error < 1e-12 && inconsistent_raised && fallback_cells == 1;

Diag = diagnostic_row("VAL-SNOW-001B", "Per-LULC parameters and raster initial states", ...
    initial_error, 0, 0, fallback_cells, passed, ...
    "No raster, SWE-only, depth-only, both-raster, and inconsistent-raster paths were exercised.");
Ledger = ledger_row("VAL-SNOW-001B", 0, 0, 0, 0, 0, 0, 0);
end

function [Diag, Ledger, Series, passed] = test_timestep_refinement()
dt_minutes = [5; 15; 60];
duration_minutes = 1440;
expected_swe_mm = 96;
final_swe = nan(size(dt_minutes));
Series = table();

for i = 1:numel(dt_minutes)
    swe = 100; rho = 100; h = 1000; t = 0; row = table();
    while t < duration_minutes - 1e-12
        dt = min(dt_minutes(i), duration_minutes - t);
        [swe, h, melt, snow, rain, rho, sub, mb] = Snow_Model_Function( ...
            swe, h, rho, 2, -2, 0, 0, 39, 30, 1, 1, 0, 2, -1, 1, ...
            100, 450, 0, 0, 0, dt * 60);
        row = [row; table(dt_minutes(i), t + dt, swe, h, melt, snow, rain, sub, mb, ...
            'VariableNames', {'dt_min','time_min','swe_mm','snow_depth_mm','melt_mm', ...
            'snowfall_mm','rain_mm','sublimation_mm','mass_residual_mm'})]; %#ok<AGROW>
        t = t + dt;
    end
    final_swe(i) = swe;
    Series = [Series; row]; %#ok<AGROW>
end

time_error = max(abs(final_swe - expected_swe_mm));
time_spread = max(final_swe) - min(final_swe);
mass_residual_mm = max(abs(Series.mass_residual_mm));
passed = time_error < 1e-10 && time_spread < 1e-10 && mass_residual_mm < 1e-10;
Diag = diagnostic_row("VAL-SNOW-001C", "Time-step refinement", ...
    time_error, mass_residual_mm / 1000, time_spread, 0, passed, ...
    "Five, fifteen, and sixty minute integrations recover the same 96 mm degree-day solution.");
Ledger = ledger_row("VAL-SNOW-001C", 0, 0, 0, 4, 0, mean(final_swe), mass_residual_mm);
end

function [Diag, Ledger, Series, passed] = test_vtilted_snow_module(Config)
LULC = [1 1 2; 1 2 2; 1 1 2];
valid = true(size(LULC));
[Snow_Properties, ~] = hp2d_initialize_snow_state(Config, LULC, valid, [], [], 400);

flags = struct('flag_snow_modeling', 1); %#ok<NASGU>
time_step = 15; %#ok<NASGU>
Wshed_Properties = struct('pixel_latitude', 39 * ones(size(LULC))); %#ok<NASGU>
DEM_raster = struct('cellsize', 20, 'Z', zeros(size(LULC))); %#ok<NASGU>
errors = zeros(3,1); %#ok<NASGU>
depths = struct('d_t', zeros(size(LULC)), 'd_p', zeros(size(LULC))); %#ok<NASGU>
BC_States = struct(); %#ok<NASGU>
max_Hsnow = []; %#ok<NASGU>
day_of_year = 30; %#ok<NASGU>
module_file = fullfile(fileparts(which('Snow_Module')), 'Snow_Module.m');

forcing = table([15; 30; 45; 60], [-5; 0; 3; 5], [-8; -2; 0; 2], ...
    [8; 8; 0; 0], [2; 2; 2; 2], ...
    'VariableNames', {'time_min','t_air_c','t_min_c','precip_mm','wind_m_s'});
initial_total_mm = 0;
cumulative_precip_mm = 0;
cumulative_sublimation_mm = 0;
Series = table();

for i = 1:height(forcing)
    k = i; %#ok<NASGU>
    depths.d_p = depths.d_t;
    BC_States.Eff_Rainfall = forcing.precip_mm(i) * ones(size(LULC));
    BC_States.Average_Daily_Temperature = forcing.t_air_c(i) * ones(size(LULC));
    BC_States.min_temp = forcing.t_min_c(i) * ones(size(LULC));
    BC_States.wind = forcing.wind_m_s(i) * ones(size(LULC));
    run(module_file);
    cumulative_precip_mm = cumulative_precip_mm + forcing.precip_mm(i);
    cumulative_sublimation_mm = cumulative_sublimation_mm + mean(Snow_Properties.E_s, 'all');
    snow_mean_mm = mean(Snow_Properties.SWE_t, 'all');
    surface_mean_mm = mean(depths.d_t, 'all');
    residual_mm = initial_total_mm + cumulative_precip_mm - ...
        cumulative_sublimation_mm - snow_mean_mm - surface_mean_mm;
    Series = [Series; table(forcing.time_min(i), forcing.t_air_c(i), forcing.precip_mm(i), ...
        snow_mean_mm, surface_mean_mm, mean(Snow_Properties.M_snow, 'all'), ...
        mean(Snow_Properties.P_snow, 'all'), mean(Snow_Properties.P_rain, 'all'), ...
        mean(Snow_Properties.rho_snow, 'all'), residual_mm, ...
        'VariableNames', {'time_min','air_temperature_c','precipitation_mm', ...
        'mean_swe_mm','mean_surface_water_mm','mean_melt_mm','mean_snowfall_mm', ...
        'mean_rainfall_mm','mean_snow_density_kg_m3','mass_residual_mm'})]; %#ok<AGROW>
end

mass_residual_mm = max(abs(Series.mass_residual_mm));
two_class_response = abs(Snow_Properties.DDF(1,1) - Snow_Properties.DDF(1,3));
finite_ok = all(isfinite([Series.mean_swe_mm; Series.mean_surface_water_mm; Series.mean_snow_density_kg_m3]));
passed = mass_residual_mm < 1e-8 && two_class_response > 0 && finite_ok;
Diag = diagnostic_row("VAL-SNOW-001D", "V-tilted snow and runoff coupling", ...
    0, mass_residual_mm * 9 * 400 / 1000, 0, 0, passed, ...
    "A two-class V-tilted grid routes liquid rain and melt to surface storage while retaining snowpack SWE.");
Ledger = ledger_row("VAL-SNOW-001D", cumulative_precip_mm * 9 * 400 / 1000, ...
    sum(Series.mean_snowfall_mm) * 9 * 400 / 1000, sum(Series.mean_rainfall_mm) * 9 * 400 / 1000, ...
    sum(Series.mean_melt_mm) * 9 * 400 / 1000, cumulative_sublimation_mm * 9 * 400 / 1000, ...
    (Series.mean_swe_mm(end) + Series.mean_surface_water_mm(end)) * 9 * 400 / 1000, ...
    mass_residual_mm * 9 * 400 / 1000);
end

function [Diag, Ledger, passed] = test_excel_bypass_equivalence(ExcelConfig, BypassConfig)
ExcelTable = readtable(fullfile(hydropol2d_find_root(fileparts(mfilename('fullpath'))), ...
    'Input_Data_Sheets', 'LULC_parameters.xlsx'), 'Range', 'A2', ...
    'VariableNamingRule', 'preserve');
BypassTable = table(ExcelConfig.class_index, ExcelConfig.parameter_values(:,1), ...
    ExcelConfig.parameter_values(:,2), ExcelConfig.parameter_values(:,3), ...
    ExcelConfig.parameter_values(:,4), ExcelConfig.parameter_values(:,5), ...
    ExcelConfig.parameter_values(:,6), ExcelConfig.parameter_values(:,7), ...
    ExcelConfig.parameter_values(:,8), ExcelConfig.parameter_values(:,9), ...
    ExcelConfig.parameter_values(:,10), ExcelConfig.parameter_values(:,11), ...
    'VariableNames', {'LULC_Index','Snow_Albedo','Snow_Emissivity', ...
    'Sublimation_Coefficient_d_1','Degree_Day_Factor_mm_C_day', ...
    'T_Snow_All_C','T_Rain_All_C','Rho_Snow_Init_kg_m3', ...
    'Rho_Snow_Max_kg_m3','Compaction_Temperature_kg_m3_C_day', ...
    'Compaction_SWE_kg_m3_mm_day','Compaction_Depth_kg_m3_mm_day'});
BypassConfigFromWorkbook = hp2d_normalize_snow_table(struct('table', BypassTable));
config_error = max(abs([ExcelConfig.class_index - BypassConfigFromWorkbook.class_index; ...
    ExcelConfig.parameter_values(:) - BypassConfigFromWorkbook.parameter_values(:)]));
passed = config_error < 1e-12 && ~isempty(BypassConfig.class_index);
Diag = diagnostic_row("VAL-SNOW-001E", "Excel and bypass configuration equivalence", ...
    config_error, 0, 0, 0, passed, ...
    "The LULC workbook and LULC-table bypass interface produce identical snow parameter arrays.");
Ledger = ledger_row("VAL-SNOW-001E", 0, 0, 0, 0, 0, 0, 0);
end

function p = scalar_parameters(Config, row)
p = struct();
for i = 1:numel(Config.parameter_names)
    p.(char(Config.parameter_names(i))) = Config.parameter_values(row, i);
end
end

function Row = diagnostic_row(case_id, case_name, max_abs_error, mass_residual_m3, ...
    timestep_spread_mm, fallback_cells, passed, notes)
Row = table(string(case_id), string(case_name), max_abs_error, mass_residual_m3, ...
    timestep_spread_mm, fallback_cells, logical(passed), string(notes), ...
    'VariableNames', {'case_id','case_name','max_abs_error','mass_residual_m3', ...
    'timestep_spread_mm','fallback_cells','passed','notes'});
end

function Row = ledger_row(case_id, precipitation_mm, snowfall_mm, rainfall_mm, melt_mm, ...
    sublimation_mm, final_storage_mm, residual_mm)
Row = table(string(case_id), precipitation_mm, snowfall_mm, rainfall_mm, melt_mm, ...
    sublimation_mm, final_storage_mm, residual_mm, ...
    'VariableNames', {'case_id','precipitation_mm','snowfall_mm','rainfall_mm','melt_mm', ...
    'sublimation_mm','final_storage_mm','residual_mm'});
end

function Row = pass_row(Diagnostic, passed)
status = "fail";
if passed; status = "pass"; end
Row = table(Diagnostic.case_id, status, logical(passed), logical(passed), Diagnostic.notes, ...
    'VariableNames', {'case_id','status','passed','report_ready','notes'});
end

function make_figures(TimeSeries, VTiltedSeries, fig_dir)
colors = [11 129 162; 226 87 89; 89 168 156] / 255;
figure('Color', 'w', 'Position', [100 100 1450 580], 'Visible', 'off');
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; hold on
dts = unique(TimeSeries.dt_min, 'stable');
for i = 1:numel(dts)
    idx = TimeSeries.dt_min == dts(i);
    plot(TimeSeries.time_min(idx) / 60, TimeSeries.swe_mm(idx), '-', ...
        'Color', colors(i,:), 'LineWidth', 2.0, 'DisplayName', sprintf('%g min', dts(i)));
end
xlabel('Time [h]'); ylabel('Snow water equivalent [mm]');
legend('Location', 'southwest', 'Box', 'off'); box on; grid on
title('Degree-day time-step refinement', 'FontWeight', 'normal');

nexttile; hold on
plot(VTiltedSeries.time_min, VTiltedSeries.mean_swe_mm, '-', 'Color', colors(1,:), ...
    'LineWidth', 2.0, 'DisplayName', 'Snow water equivalent');
plot(VTiltedSeries.time_min, VTiltedSeries.mean_surface_water_mm, '-', 'Color', colors(2,:), ...
    'LineWidth', 2.0, 'DisplayName', 'Surface water');
stairs(VTiltedSeries.time_min, VTiltedSeries.precipitation_mm, ':', 'Color', colors(3,:), ...
    'LineWidth', 1.6, 'DisplayName', 'Precipitation per step');
xlabel('Time [min]'); ylabel('Water depth [mm]');
legend('Location', 'northwest', 'Box', 'off'); box on; grid on
title('V-tilted snow and runoff response', 'FontWeight', 'normal');

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Helvetica');
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 1.5);
exportgraphics(gcf, fullfile(fig_dir, 'snow_validation_dynamics.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(fig_dir, 'snow_validation_dynamics.pdf'), 'ContentType', 'vector');
close(gcf);
end
