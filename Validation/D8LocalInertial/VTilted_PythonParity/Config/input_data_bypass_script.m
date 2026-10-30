% VAL-D8LI-001 full-model configuration (D8 local-inertial routing).
% run_d8_local_inertial_check.m supplies InputPaths and ValidationCase.mode
% ('hybrid' or 'full'). The catchment, rain and parameters are those of
% make_python_reference.py, part 2 (the Python run_case runs).

if ~exist('InputPaths', 'var') || ~isstruct(InputPaths)
    error('VAL-D8LI-001 requires InputPaths from run_d8_local_inertial_check.m.');
end
if ~exist('ValidationCase', 'var') || ~isfield(ValidationCase, 'mode')
    error('VAL-D8LI-001 requires ValidationCase.mode from run_d8_local_inertial_check.m.');
end

InputData_Bypass = struct();

G = struct();
G.time_step_model = 1;                  % min
G.min_time_step = 0.01;                 % s (Python routing.min_dt_s)
G.max_time_step = 2;                    % s (Python routing.max_dt_s)
G.time_step_increments = 0.001;
G.time_step_change = 0.0001;            % refresh the D4 step every step, as Python does
G.alfa_min = 0.50;                      % Python routing.courant
G.alfa_max = 0.50;
G.date_begin = datetime(2026, 10, 30, 0, 0, 0);
G.date_end = datetime(2026, 10, 30, 0, 30, 0);
G.routing_time = 30;                    % min
G.etp_input_spreadsheet = '';
G.slope_outlet = 0.004;
G.n_outlets_data = 1;                   % the lowest perimeter cell (bottom centre)
G.record_time_maps = 2;
G.record_time_hydrographs = 1;
G.record_time_spatial_rainfall = 2;
G.time_save_ETP = 60;
G.record_time_spatial_ETP = 60;
G.Pol_min = 0.001;
G.depth_wse = 0.01;
G.flag_wse = 0;
G.Krs_ETP = 0.13;
G.albedo = 0.25;
G.volume_error = 1e6;
G.factor_reduction = 2;
G.CA_depth_tolerance = 1e-4;
G.alfa_1 = 0.15;                        % bank depth 0.15 A^0.2 (A in km2)
G.alfa_2 = 0.2;
G.beta_1 = 5;                           % width 5 A^0.2
G.beta_2 = 0.2;
G.Manning = 0.04;                       % channel n
G.River_K_coeff = 1;
G.ADD = 0;
G.min_Bt = 0;
G.Bmin = 0;
G.Bmax = 0;
G.min_area = 0.0012;                    % km2: rivers drain 3 cells or more
G.tau = 0;
G.K_value = 0;
G.sl = 0;
G.slope_DTM = 0;
G.RP = 1;
G.Rainfall_Duration = 10;
G.K = 1;
G.a = 1;
G.b = 1;
G.c = 1;
G.dt_design = 10;
G.d8_min_depth_m = 1e-3;
G.d8_hydraulic_table_path = '';
G.d8_max_stage_above_bed_m = 100;
InputData_Bypass.general = G;

flag_names = { ...
    'flag_rainfall','flag_spatial_rainfall','flag_ETP','flag_input_rainfall_map', ...
    'flag_rainfall_multiple_runs','flag_data_source','flag_inflow','flag_satellite_rainfall', ...
    'flag_alternated_blocks','flag_huff','flag_stage_hydrograph','flag_input_ETP_map', ...
    'flag_timestep','flag_infiltration','flag_critical','flag_D8','flag_CA','flag_inertial', ...
    'flag_full_momentum','flag_waterbalance','flag_waterquality','flag_reservoir','flag_wq_model', ...
    'flag_groundwater_modeling','flag_real_time_satellite_rainfall','flag_dam_break', ...
    'flag_human_instability','flag_boundary','flag_numerical_scheme','flag_outlet_type', ...
    'flag_adaptive_timestepping','flag_neglect_infiltration_river','flag_subgrid', ...
    'flag_spatial_albedo','flag_river_rasters','flag_baseflow','flag_kinematic', ...
    'flag_diffusive','flag_DTM','flag_abstraction','flag_overbanks','flag_snow_modeling', ...
    'flag_WQ_Rasters','flag_GPU','flag_single','flag_warmup','flag_initial_buildup', ...
    'flag_resample','flag_smoothening','flag_trunk','flag_fill_DEM','flag_smooth_cells', ...
    'flag_reduce_DEM','flag_export_maps','flag_river_heigth_compensation', ...
    'flag_dashboard','flag_elapsed_time','flag_obs_gauges'};
InputData_Bypass.flags = cell2struct(num2cell(zeros(numel(flag_names), 1)), flag_names, 1);
InputData_Bypass.flags.flag_rainfall = 1;
InputData_Bypass.flags.flag_timestep = 1;
InputData_Bypass.flags.flag_inertial = 1;
InputData_Bypass.flags.flag_numerical_scheme = 1;
InputData_Bypass.flags.flag_outlet_type = 1;
InputData_Bypass.flags.flag_adaptive_timestepping = 1;
InputData_Bypass.flags.flag_d8_local_inertial = 1;
InputData_Bypass.flags.flag_d8_full_domain = double(strcmp(ValidationCase.mode, 'full'));

LULC = table({'V-tilted'}, 1, 0.04, 0, 0, 0, 0, 0, 0, NaN, 0.3, 0.8, 0.98, 0, 2, -1, 1, ...
    100, 450, 0, 0, 0, ...
    'VariableNames', {'LULC','Index','roughness','h_0_mm','d_0_mm','C1','C2','C3','C4', ...
    'index_impervious','root_depth_m','Snow_Albedo','Snow_Emissivity', ...
    'Sublimation_Coefficient_d_1','Degree_Day_Factor_mm_C_day','T_Snow_All_C', ...
    'T_Rain_All_C','Rho_Snow_Init_kg_m3','Rho_Snow_Max_kg_m3', ...
    'Compaction_Temperature_kg_m3_C_day','Compaction_SWE_kg_m3_mm_day', ...
    'Compaction_Depth_kg_m3_mm_day'});
InputData_Bypass.LULC = struct('table', LULC);

SOIL = table({'V-tilted soil'}, 1, 0.001, 1.6, 2.5, 0.45, 0.05, 0.20, 0.20, 0.001, ...
    2.0, 1, 1, 1, 0.10, 0.25, 0.50, ...
    'VariableNames', {'Soil_type','Index','ksat_mm_h','n_vg','alpha_vg_1_m','theta_sat', ...
    'theta_r','theta_i','Sy','ksat_gw_mm_h','Soil_Depth_m', ...
    'Ks_multiplier_near_surface','Ks_multiplier_root_zone', ...
    'Ks_multiplier_transmission','Ltop_m','dh_max_m','l_vg'});
InputData_Bypass.SOIL = struct('table', SOIL);

% Each time starts an interval: 47 mm/h from 0 to 10 min, then none.
InputData_Bypass.Rainfall_Parameters = struct( ...
    'time_rainfall', [0; 10; 20; 30], ...
    'intensity_rainfall', [47; 0; 0; 0], ...
    'time_step_rainfall', 10, ...
    'rainfall_duration', 30, ...
    'n_obs_rainfall', 4);
