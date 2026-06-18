% V-tilted infiltration validation bypass inputs.
%
% The parent runner must define ValidationCase from
% Config/Infiltration_Case_Registry.csv before HydroPol2D_preprocessing calls
% this script.

if ~exist('ValidationCase','var') || ~istable(ValidationCase)
    error('ValidationCase table row is required before running this bypass script.');
end

Cfg = ValidationCase;
InputData_Bypass = struct();

InputData_Bypass.general = struct();
InputData_Bypass.general.time_step_model      = Cfg.dt_min;  % [min]
InputData_Bypass.general.min_time_step        = 1; % [s]
InputData_Bypass.general.max_time_step        = max(1, Cfg.dt_min * 60); % [s]
InputData_Bypass.general.time_step_increments = 0.001;
InputData_Bypass.general.time_step_change     = 0.0001;
InputData_Bypass.general.alfa_min             = 0.50;
InputData_Bypass.general.alfa_max             = 0.50;
InputData_Bypass.general.date_begin           = datetime(2025,5,1,0,0,0);
InputData_Bypass.general.date_end             = InputData_Bypass.general.date_begin + hours(Cfg.duration_h);
InputData_Bypass.general.routing_time         = Cfg.duration_h * 60;
InputData_Bypass.general.dt_rainfall_maps_min = Cfg.dt_min;
InputData_Bypass.general.dt_transpiration_maps_min = 1440;
InputData_Bypass.general.dt_evaporation_maps_min = 1440;
InputData_Bypass.general.rainfall_filename_example = 'rain_2025_05_01_00_00.tif';
InputData_Bypass.general.slope_outlet = 0.02;
InputData_Bypass.general.n_outlets_data = 1;
InputData_Bypass.general.record_time_maps = max(Cfg.dt_min, 10);
InputData_Bypass.general.record_time_hydrographs = max(Cfg.dt_min, 10);
InputData_Bypass.general.Pol_min = 0.001;
InputData_Bypass.general.depth_wse = 0.01;
InputData_Bypass.general.flag_wse = 0;
InputData_Bypass.general.record_time_spatial_rainfall = max(Cfg.dt_min, 10);
InputData_Bypass.general.time_save_ETP = 1440;
InputData_Bypass.general.record_time_spatial_ETP = 1440;
InputData_Bypass.general.Krs_ETP = 0.13;
InputData_Bypass.general.albedo = 0.25;
InputData_Bypass.general.volume_error = 1e9;
InputData_Bypass.general.factor_reduction = 2;
InputData_Bypass.general.CA_depth_tolerance = 1e-4;
InputData_Bypass.general.alfa_1 = 1;
InputData_Bypass.general.alfa_2 = 1;
InputData_Bypass.general.beta_1 = 1;
InputData_Bypass.general.beta_2 = 1;
InputData_Bypass.general.Manning = 0.04;
InputData_Bypass.general.River_K_coeff = 1;
InputData_Bypass.general.ADD = 0;
InputData_Bypass.general.min_Bt = 0;
InputData_Bypass.general.Bmin = 0;
InputData_Bypass.general.Bmax = 0;
InputData_Bypass.general.min_area = 0.0001;
InputData_Bypass.general.tau = 0;
InputData_Bypass.general.K_value = 0;
InputData_Bypass.general.sl = 0;
InputData_Bypass.general.slope_DTM = 0;
InputData_Bypass.general.topo_path = InputPaths.topo_path;
InputData_Bypass.general.RP = 1;
InputData_Bypass.general.Rainfall_Duration = Cfg.duration_h * 60;
InputData_Bypass.general.K = 1;
InputData_Bypass.general.a = 1;
InputData_Bypass.general.b = 1;
InputData_Bypass.general.c = 1;
InputData_Bypass.general.dt_design = Cfg.dt_min;

InputData_Bypass.flags = struct();
InputData_Bypass.flags.flag_rainfall = double(Cfg.rainfall_mm_h > 0);
InputData_Bypass.flags.flag_spatial_rainfall = 0;
InputData_Bypass.flags.flag_ETP = 0;
InputData_Bypass.flags.flag_input_rainfall_map = 0;
InputData_Bypass.flags.flag_rainfall_multiple_runs = 0;
InputData_Bypass.flags.flag_data_source = 0;
InputData_Bypass.flags.flag_inflow = 0;
InputData_Bypass.flags.flag_satellite_rainfall = 0;
InputData_Bypass.flags.flag_alternated_blocks = 0;
InputData_Bypass.flags.flag_huff = 0;
InputData_Bypass.flags.flag_stage_hydrograph = 0;
InputData_Bypass.flags.flag_input_ETP_map = 0;
InputData_Bypass.flags.flag_timestep = 1;
InputData_Bypass.flags.flag_infiltration = double(Cfg.infiltration_enabled);
InputData_Bypass.flags.flag_critical = 0;
InputData_Bypass.flags.flag_D8 = 0;
InputData_Bypass.flags.flag_CA = 0;
InputData_Bypass.flags.flag_inertial = 1;
InputData_Bypass.flags.flag_full_momentum = 0;
InputData_Bypass.flags.flag_waterbalance = 0;
InputData_Bypass.flags.flag_waterquality = 0;
InputData_Bypass.flags.flag_reservoir = 0;
InputData_Bypass.flags.flag_wq_model = 0;
InputData_Bypass.flags.flag_groundwater_modeling = double(Cfg.enable_recharge);
InputData_Bypass.flags.flag_real_time_satellite_rainfall = 0;
InputData_Bypass.flags.flag_dam_break = 0;
InputData_Bypass.flags.flag_human_instability = 0;
InputData_Bypass.flags.flag_boundary = 0;
InputData_Bypass.flags.flag_numerical_scheme = 1;
InputData_Bypass.flags.flag_outlet_type = 1;
InputData_Bypass.flags.flag_adaptive_timestepping = 1;
InputData_Bypass.flags.flag_neglect_infiltration_river = 0;
InputData_Bypass.flags.flag_subgrid = 0;
InputData_Bypass.flags.flag_spatial_albedo = 0;
InputData_Bypass.flags.flag_river_rasters = 0;
InputData_Bypass.flags.flag_baseflow = 0;
InputData_Bypass.flags.flag_kinematic = 0;
InputData_Bypass.flags.flag_diffusive = 0;
InputData_Bypass.flags.flag_DTM = 0;
InputData_Bypass.flags.flag_abstraction = 0;
InputData_Bypass.flags.flag_overbanks = 0;
InputData_Bypass.flags.flag_snow_modeling = 0;
InputData_Bypass.flags.flag_WQ_Rasters = 0;
InputData_Bypass.flags.flag_GPU = 0;
InputData_Bypass.flags.flag_single = 0;
InputData_Bypass.flags.flag_warmup = 1;
InputData_Bypass.flags.flag_initial_buildup = 0;
InputData_Bypass.flags.flag_resample = 0;
InputData_Bypass.flags.flag_smoothening = 0;
InputData_Bypass.flags.flag_trunk = 0;
InputData_Bypass.flags.flag_fill_DEM = 0;
InputData_Bypass.flags.flag_smooth_cells = 0;
InputData_Bypass.flags.flag_reduce_DEM = 0;
InputData_Bypass.flags.flag_export_maps = 1;
InputData_Bypass.flags.flag_river_heigth_compensation = 0;
InputData_Bypass.flags.flag_dashboard = 0;
InputData_Bypass.flags.flag_elapsed_time = 0;
InputData_Bypass.flags.flag_obs_gauges = 0;

LULC_table = table();
LULC_table.LC = {'Left hillslope'; 'Channel strip'; 'Right hillslope'};
LULC_table.Index = [1; 2; 3];
LULC_table.roughness = [0.060; 0.035; 0.060];
LULC_table.h_0_mm = zeros(3,1);
LULC_table.d_0_mm = zeros(3,1);
LULC_table.C1 = zeros(3,1);
LULC_table.C2 = zeros(3,1);
LULC_table.C3 = zeros(3,1);
LULC_table.C4 = zeros(3,1);
LULC_table.index_impervious = [999; NaN; NaN];
LULC_table.root_depth_m = Cfg.root_depth_m * ones(3,1);
InputData_Bypass.LULC.table = LULC_table;

ks1 = Cfg.ksat_mm_h;
ks2 = Cfg.ksat_mm_h;
ks3 = Cfg.ksat_mm_h;
if Cfg.regime == "spatial_soil_contrast"
    ks1 = Cfg.ksat_mm_h;
    ks2 = Cfg.ksat_alt_mm_h;
    ks3 = Cfg.ksat_alt_mm_h;
end

SOIL_table = table();
SOIL_table.Soil_type = {'Validation soil 1'; 'Validation soil 2'; 'Validation soil 3'};
SOIL_table.Index = [1; 2; 3];
SOIL_table.ksat_mm_h = [ks1; ks2; ks3];
SOIL_table.n_vg = Cfg.n_vg * ones(3,1);
SOIL_table.alpha_vg_1_m = Cfg.alpha_vg * ones(3,1);
SOIL_table.theta_sat = Cfg.theta_sat * ones(3,1);
SOIL_table.theta_r = Cfg.theta_r * ones(3,1);
SOIL_table.theta_i = Cfg.theta_i * ones(3,1);
SOIL_table.Sy = max(Cfg.theta_sat - Cfg.theta_r, 1e-6) * ones(3,1);
SOIL_table.ksat_gw_mm_h = [ks1; ks2; ks3];
SOIL_table.Soil_Depth_m = Cfg.soil_depth_m * ones(3,1);
SOIL_table.Ks_multiplier_near_surface = Cfg.ks_mult_near * ones(3,1);
SOIL_table.Ks_multiplier_root_zone = Cfg.ks_mult_root * ones(3,1);
SOIL_table.Ks_multiplier_transmission = Cfg.ks_mult_trans * ones(3,1);
InputData_Bypass.SOIL.table = SOIL_table;

Rainfall_Parameters = struct();
Rainfall_Parameters.time_rainfall = [0; Cfg.duration_h * 60];
Rainfall_Parameters.intensity_rainfall = [Cfg.rainfall_mm_h; Cfg.rainfall_mm_h];
Rainfall_Parameters.time_step_rainfall = Cfg.duration_h * 60;
Rainfall_Parameters.rainfall_duration = Cfg.duration_h * 60;
Rainfall_Parameters.n_obs_rainfall = 2;
InputData_Bypass.Rainfall_Parameters = Rainfall_Parameters;

InputData_Bypass.Snow_Properties = struct();
