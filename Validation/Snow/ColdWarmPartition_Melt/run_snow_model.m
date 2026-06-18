% P1-SNOW-001: run HydroPol2D Snow_Model_Function against the Phase 1 reference.
%
% Run from the repository root or from this case folder in MATLAB:
%   run('HydroPol2D_Model/Validation/Snow/ColdWarmPartition_Melt/run_snow_model.m')

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
reference_file = fullfile(repo_root, 'HydroPol2D_Model', 'Validation', ...
    'Reference_Outputs', 'Phase1', 'snow_degree_day', ...
    'P1-SNOW-001_reference.csv');
output_dir = fullfile(case_dir, 'Outputs', 'Validation');

if ~exist(functions_dir, 'dir')
    error('HydroPol2D functions directory not found: %s', functions_dir);
end
if ~exist(reference_file, 'file')
    error('Reference file not found. Generate it with phase1_reference_solutions.py: %s', reference_file);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

addpath(functions_dir);

ref = readtable(reference_file, 'TextType', 'string');
model_rows = table();

scenarios = unique(ref.scenario, 'stable');
for s = 1:numel(scenarios)
    scenario = scenarios(s);
    idx = ref.scenario == scenario;
    ref_s = ref(idx, :);

    swe_prev = 0;
    h_snow_prev = 0;

    for i = 1:height(ref_s)
        T_air = ref_s.temperature_c(i);
        T_min = ref_s.t_min_c(i);
        P = ref_s.precipitation_mm(i);
        wind = ref_s.wind_m_s(i);
        lat = ref_s.lat_deg(i);
        DOY = ref_s.doy(i);
        H_snow_t0 = ref_s.h_snow_t0_mm(i);
        alpha = ref_s.alpha(i);
        epsilon = ref_s.epsilon(i);
        C_e = ref_s.c_e(i);
        DDF = ref_s.ddf_mm_c_day(i);
        T_thresh = ref_s.t_thresh_c(i);
        rho_snow_init = ref_s.rho_snow_init_kg_m3(i);
        rho_max = ref_s.rho_max_kg_m3(i);
        k_t = ref_s.k_t(i);
        k_swe = ref_s.k_swe(i);
        k_D = ref_s.k_d(i);

        [SWE_t, H_snow_t, M_snow, P_snow, P_rain, rho_snow, E_s, mb_error] = ...
            Snow_Model_Function(swe_prev, h_snow_prev, T_air, T_min, P, ...
            wind, lat, DOY, H_snow_t0, alpha, epsilon, C_e, DDF, T_thresh, ...
            rho_snow_init, rho_max, k_t, k_swe, k_D);

        next_row = table( ...
            scenario, ...
            ref_s.zone_id(i), ...
            ref_s.zone_cells(i), ...
            ref_s.time_day(i), ...
            ref_s.cell_area_m2(i), ...
            T_air, ...
            T_min, ...
            P, ...
            wind, ...
            lat, ...
            DOY, ...
            P_rain, ...
            P_snow, ...
            M_snow, ...
            E_s, ...
            SWE_t, ...
            H_snow_t, ...
            rho_snow, ...
            mb_error, ...
            'VariableNames', {'scenario','zone_id','zone_cells','time_day', ...
            'cell_area_m2','temperature_c','t_min_c','precipitation_mm', ...
            'wind_m_s','lat_deg','doy','rain_mm','snowfall_mm','melt_mm', ...
            'sublimation_mm','swe_mm','snow_depth_mm','rho_snow_kg_m3', ...
            'mass_residual_mm'});

        model_rows = [model_rows; next_row]; %#ok<AGROW>

        swe_prev = SWE_t;
        h_snow_prev = H_snow_t;
    end
end

writetable(model_rows, fullfile(output_dir, 'HydroPol2D_Snow_Model.csv'));

disp('Wrote HydroPol2D snow model output to:');
disp(fullfile(output_dir, 'HydroPol2D_Snow_Model.csv'));
disp('Run compare_snow_model.py to generate metrics and pass/fail diagnostics.');
