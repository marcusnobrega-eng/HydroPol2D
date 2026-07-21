function Cal = hydropol2d_calibration_config_template()
%HYDROPOL2D_CALIBRATION_CONFIG_TEMPLATE Template for new catchments.
%
% Copy this file into a catchment-specific calibration folder and rename it.
% Then edit the sections marked USER INPUT. The generic GA engine does not
% know about any particular catchment; all model-specific behavior enters
% through this configuration struct.

%% USER INPUT: paths and identity
model_root = '/path/to/HydroPol2D';
case_root = fullfile(model_root, 'Examples', 'YourCatchment');
addpath(fullfile(model_root, 'Calibration', 'GenericGA'));
addpath(case_root);

Cal = struct();
Cal.name = 'YourCatchment_Outlet_GA';
Cal.random_seed = 42;
Cal.output_root = fullfile(case_root, 'Outputs', 'Calibration', ...
    char(datetime("now", "Format", "yyyyMMdd_HHmmss")));

%% USER INPUT: GA controls
Cal.ga = struct();
Cal.ga.n_generations = 5;
Cal.ga.n_individuals = 10;
Cal.ga.elite_fraction = 0.20;
Cal.ga.tournament_size = 3;
Cal.ga.crossover_fraction = 0.85;
Cal.ga.mutation_probability = 0.25;
Cal.ga.mutation_scale_fraction = 0.15;

%% USER INPUT: single-objective weights
Cal.objective = struct();
Cal.objective.weights = struct( ...
    'rmse', 1.0, ...
    'volume', 1.0, ...
    'peak', 0.75, ...
    'timing', 0.25, ...
    'bias', 0.25);
Cal.objective.timing_scale_min = 120;

%% USER INPUT: baseline class parameters
% These are the starting values the GA perturbs. Include every LULC/SOIL
% class active in the rasters if you want global multipliers to affect all
% classes. Add/remove classes freely; the generic catalog expands to the
% number of classes listed here.
Cal.base_overrides = struct();

Cal.base_overrides.LULC = struct();
Cal.base_overrides.LULC.Index = [10; 20; 30];           % raster class codes
Cal.base_overrides.LULC.roughness = [0.08; 0.05; 0.035];
Cal.base_overrides.LULC.root_depth_m = [1.0; 0.6; 0.25];

Cal.base_overrides.SOIL = struct();
Cal.base_overrides.SOIL.Index = [4; 8; 9; 11; 12];      % raster class codes
Cal.base_overrides.SOIL.ksat_mm_h = [1; 3.4; 10.9; 29.9; 117.8];
Cal.base_overrides.SOIL.n_vg = [1.31; 1.56; 1.89; 1.75; 2.68];
Cal.base_overrides.SOIL.alpha_vg_1_m = [1.9; 3.6; 6.0; 11.0; 14.5];
Cal.base_overrides.SOIL.theta_sat = [0.309; 0.399; 0.387; 0.390; 0.430];
Cal.base_overrides.SOIL.theta_r = [0.095; 0.078; 0.100; 0.049; 0.045];
Cal.base_overrides.SOIL.theta_i = [0.260; 0.335; 0.330; 0.322; 0.353];
Cal.base_overrides.SOIL.Sy = [0.08; 0.18; 0.22; 0.25; 0.30];
Cal.base_overrides.SOIL.ksat_gw_mm_h = [20; 68; 218; 598; 2356];
Cal.base_overrides.SOIL.Soil_Depth_m = ones(5,1);
Cal.base_overrides.SOIL.Ks_multiplier_near_surface = ones(5,1);
Cal.base_overrides.SOIL.Ks_multiplier_root_zone = ones(5,1);
Cal.base_overrides.SOIL.Ks_multiplier_transmission = ones(5,1);

%% USER INPUT: optional direct InputData_Bypass overrides
% Example:
% Cal.input_data_overrides.general.slope_outlet = 0.02;
% Cal.input_data_overrides.flags.flag_infiltration = 1;
Cal.input_data_overrides = struct();

%% USER INPUT: one or more rainfall/event cases
case1 = struct();
case1.name = 'Event_001';
case1.resolution = 100;
case1.simulation_minutes = 720;
case1.record_maps_minutes = 720;
case1.record_hydrographs_minutes = 5;
case1.static_folder = fullfile(case_root, 'Static');
case1.rainfall_timeseries_path = fullfile(case_root, 'Forcing', 'Rainfall', 'rain.csv');
case1.use_spatial_rainfall_gauges = false;
case1.event_start = datetime(2021, 1, 1, 0, 0, 0);
case1.observed_runoff_csv = fullfile(case_root, 'Forcing', 'Observations', 'outlet.csv');
Cal.cases = case1;

%% USER INPUT: catchment-specific runner
% This function must accept:
%   Eval = runner(Cal, caseDef, candidate, info)
% and return Eval.metrics fields used by the objective function.
Cal.model_runner = @your_catchment_calibration_candidate_runner;

%% GENERIC PARAMETER CATALOG
Catalog = hydropol2d_default_parameter_catalog(Cal.base_overrides);

% Optional: edit Catalog before building specs. Examples:
% Catalog.SOIL.Fields(contains({Catalog.SOIL.Fields.field}, 'theta_i')).lower = 0.05;
% Catalog.InputData(end+1) = struct('name','my_parameter', 'path','general.x', ...);

Cal.parameter_specs = hydropol2d_build_parameter_specs(Cal.base_overrides, Catalog);

% Enable exact names or wildcard patterns. Examples:
%   "soil_*_theta_i" enables class-specific initial moisture for all soils.
%   "lulc_*_roughness" enables class-specific roughness for all LULC classes.
enableNames = [
    "lulc_roughness_multiplier_all"
    "soil_theta_i_multiplier_all"
    "soil_ks_near_surface_multiplier_all"
    "soil_ks_root_zone_multiplier_all"
    "soil_ks_transmission_multiplier_all"
    "soil_depth_multiplier_all"];
Cal.parameter_specs = hydropol2d_enable_parameter_specs(Cal.parameter_specs, enableNames);
Cal.parameter_specs = hydropol2d_finalize_parameter_specs(Cal.parameter_specs);
end
