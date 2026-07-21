function PF = pf_config_vtilted(scenario)
%PF_CONFIG_VTILTED Configuration for the V-tilted augmented particle filter.
%
% This file intentionally keeps the DA prototype configuration in MATLAB so
% parameters, priors, and truth values can be edited without touching the
% production HydroPol2D input spreadsheets or model source.
if nargin < 1 || strlength(string(scenario)) == 0
    scenario = "constant_rain";
end

case_dir = fileparts(mfilename('fullpath'));
model_root = fullfile(case_dir, '..', '..');

PF = struct();
PF.case_dir = case_dir;
PF.model_root = model_root;
PF.static_source_dir = fullfile(model_root, 'Validation', 'Phase1_VTilted_Catchment', 'Static');
PF.runs_dir = fullfile(case_dir, 'Runs');
PF.outputs_dir = fullfile(case_dir, 'Outputs');
PF.figures_dir = fullfile(case_dir, 'Figures');

PF.random_seed = 42;
PF.n_particles = 100;
PF.assimilation_window_min = 30;
PF.simulation_duration_min = 360;
PF.dt_min = 1;
PF.resample_threshold_fraction = 0.5;
PF.forecast_backend = 'compact_vtilted_augmented';
PF.state_update_method = 'identity';
PF.spinup_min = 0;

PF.resampling = struct();
PF.resampling.resample_every_window = false;
PF.resampling.inflation_lambda = 1.5;

PF.weighting = struct();
PF.weighting.method = 'paper_beta_normal_pdf';
PF.weighting.group_weights = struct( ...
    'discharge', 1/3, ...
    'groundwater_depth', 1/3, ...
    'soil_moisture', 1/3);
PF.weighting.sigma = struct( ...
    'discharge_mm_h', 0.05, ...
    'groundwater_depth_m', 0.025, ...
    'soil_moisture_m3m3', 0.025);
PF.observation_noise_scale = 1;

PF.forcing = struct();
PF.forcing.rainfall_intensity_mm_h = 100;
PF.forcing.rainfall_duration_min = inf;

PF.groundwater = struct();
PF.groundwater.target_dt_min = 1440;

PF.routing = struct();
PF.routing.min_slope = 1e-5;
PF.routing.max_drain_fraction = 0.85;
PF.routing.surface_depth_epsilon_m = 1e-5;

PF.base = struct();
PF.base.lulc(1) = struct('name', 'open_pervious', 'manning_n', 0.035, 'root_depth_m', 0.45);
PF.base.lulc(2) = struct('name', 'rough_compacted', 'manning_n', 0.070, 'root_depth_m', 0.25);

PF.base.soil(1) = struct( ...
    'name', 'sandy_loam', ...
    'ksat_surface_mm_h', 18, ...
    'ksat_rootzone_mm_h', 8, ...
    'storage_capacity_mm', 120, ...
    'theta_s', 0.43, ...
    'theta_r', 0.07, ...
    'initial_sm_fraction', 0.45);
PF.base.soil(2) = struct( ...
    'name', 'clay_loam', ...
    'ksat_surface_mm_h', 4, ...
    'ksat_rootzone_mm_h', 1.5, ...
    'storage_capacity_mm', 180, ...
    'theta_s', 0.47, ...
    'theta_r', 0.10, ...
    'initial_sm_fraction', 0.55);

PF.base.gw(1) = struct( ...
    'name', 'single_aquifer', ...
    'ksat_mm_h', 0.8, ...
    'specific_yield', 0.15, ...
    'initial_wtd_m', 0.8);

% The source V-tilted raster has left hillslope, channel, and right hillslope
% classes. The two hillslopes share parameter class 1; the channel is class 2.
PF.class_map = struct('lulc', [1 2 1], 'soil', [1 2 1]);

PF.observations = [
    pf_obs('Q_outlet', 'discharge', nan, nan, 5, 0.12, 'm3/s', true)
    pf_obs('SM_left', 'soil_moisture', 0.55, 0.25, 15, 0.025, 'm3/m3', true)
    pf_obs('SM_right', 'soil_moisture', 0.55, 0.75, 15, 0.025, 'm3/m3', true)
    pf_obs('GW_channel', 'groundwater_depth', 0.70, 0.50, 30, 0.08, 'm', false)
    pf_obs('GW_left', 'groundwater_depth', 0.55, 0.25, 60, 0.06, 'm', false)
    pf_obs('GW_right', 'groundwater_depth', 0.55, 0.75, 60, 0.06, 'm', false)
    pf_obs('GW_divide', 'groundwater_depth', 0.25, 0.50, 60, 0.08, 'm', false)
];

PF.parameters = [
    pf_param('lulc1_manning_mult', 'LULC', 1, 'manning_mult', 1.0, 1.25, 0.35, 2.50, 0.50, 0.05, 'paper_multiplicative', true)
    pf_param('lulc1_root_depth_mult', 'LULC', 1, 'root_depth_mult', 1.0, 1.00, 0.35, 1.80, 0.35, 0.04, 'paper_multiplicative', false)
    pf_param('lulc2_manning_mult', 'LULC', 2, 'manning_mult', 1.0, 0.80, 0.35, 2.50, 0.50, 0.05, 'paper_multiplicative', true)
    pf_param('lulc2_root_depth_mult', 'LULC', 2, 'root_depth_mult', 1.0, 1.00, 0.35, 1.80, 0.35, 0.04, 'paper_multiplicative', false)
    pf_param('soil1_ksat_surface_mult', 'SOIL', 1, 'ksat_surface_mult', 1.0, 1.35, 0.10, 6.00, 0.90, 0.10, 'paper_multiplicative', true)
    pf_param('soil1_ksat_rootzone_mult', 'SOIL', 1, 'ksat_rootzone_mult', 1.0, 1.00, 0.10, 6.00, 0.90, 0.10, 'paper_multiplicative', false)
    pf_param('soil1_storage_mult', 'SOIL', 1, 'storage_mult', 1.0, 1.00, 0.50, 1.60, 0.35, 0.04, 'paper_multiplicative', false)
    pf_param('soil1_initial_sm_mult', 'SOIL', 1, 'initial_sm_mult', 1.0, 1.00, 0.30, 1.70, 0.50, 0.05, 'paper_multiplicative', false)
    pf_param('soil2_ksat_surface_mult', 'SOIL', 2, 'ksat_surface_mult', 1.0, 0.70, 0.10, 6.00, 0.90, 0.10, 'paper_multiplicative', true)
    pf_param('soil2_ksat_rootzone_mult', 'SOIL', 2, 'ksat_rootzone_mult', 1.0, 1.00, 0.10, 6.00, 0.90, 0.10, 'paper_multiplicative', false)
    pf_param('soil2_storage_mult', 'SOIL', 2, 'storage_mult', 1.0, 1.00, 0.50, 1.60, 0.35, 0.04, 'paper_multiplicative', false)
    pf_param('soil2_initial_sm_mult', 'SOIL', 2, 'initial_sm_mult', 1.0, 1.00, 0.30, 1.70, 0.50, 0.05, 'paper_multiplicative', false)
    pf_param('gw1_ksat_mult', 'GW', 1, 'ksat_mult', 1.0, 1.00, 0.10, 6.00, 0.90, 0.10, 'paper_multiplicative', false)
    pf_param('gw1_specific_yield_mult', 'GW', 1, 'specific_yield_mult', 1.0, 1.00, 0.35, 1.80, 0.50, 0.05, 'paper_multiplicative', false)
    pf_param('ic1_initial_wtd_mult', 'IC', 1, 'initial_wtd_mult', 1.0, 1.00, 0.05, 1.50, 0.80, 0.00, 'paper_multiplicative', false)
];

PF = apply_scenario(PF, string(scenario));
end

function PF = apply_scenario(PF, scenario)
PF.scenario = scenario;
switch lower(scenario)
    case "gw_recession"
        PF.n_particles = 20;
        PF.assimilation_window_min = 60;
        PF.simulation_duration_min = 480;
        PF.forcing.rainfall_intensity_mm_h = 20;
        PF.forcing.rainfall_duration_min = 120;
        PF.groundwater.target_dt_min = 60;
        PF.weighting.sigma.discharge_mm_h = 0.03;
        PF.weighting.sigma.groundwater_depth_m = 0.03;
        PF.weighting.sigma.soil_moisture_m3m3 = 0.025;

        PF = set_obs_enabled(PF, "GW_channel", true);

        PF = set_param(PF, 'gw1_ksat_mult', 1.8, true);
        PF = set_param(PF, 'gw1_specific_yield_mult', 0.70, true);
        PF = set_param(PF, 'ic1_initial_wtd_mult', 0.55, true);
    case {"gw_no_rain_recession", "gw_no_rain_recession_exact_obs"}
        PF.n_particles = 20;
        PF.assimilation_window_min = 60;
        PF.simulation_duration_min = 720;
        PF.forcing.rainfall_intensity_mm_h = 0;
        PF.forcing.rainfall_duration_min = 0;
        PF.groundwater.target_dt_min = 30;
        PF.resample_threshold_fraction = 0.70;
        PF.weighting.group_weights.discharge = 0.45;
        PF.weighting.group_weights.groundwater_depth = 0.55;
        PF.weighting.group_weights.soil_moisture = 0;
        PF.weighting.sigma.discharge_mm_h = 0.01;
        PF.weighting.sigma.groundwater_depth_m = 0.04;
        if lower(scenario) == "gw_no_rain_recession_exact_obs"
            PF.observation_noise_scale = 0;
            PF.weighting.sigma.discharge_mm_h = 0.001;
            PF.weighting.sigma.groundwater_depth_m = 0.005;
        end

        PF = set_obs_type_enabled(PF, "soil_moisture", false);
        PF = set_obs_type_enabled(PF, "groundwater_depth", true);

        PF = set_param(PF, 'lulc1_manning_mult', 1.0, false);
        PF = set_param(PF, 'lulc2_manning_mult', 1.0, false);
        PF = set_param(PF, 'soil1_ksat_surface_mult', 1.0, false);
        PF = set_param(PF, 'soil2_ksat_surface_mult', 1.0, false);
        PF = set_param(PF, 'gw1_ksat_mult', 3.0, true);
        PF = set_param(PF, 'gw1_specific_yield_mult', 1.25, true);
        PF = set_param_prior(PF, 'ic1_initial_wtd_mult', 0.20, -0.05, -0.20, 0.80, 2.00, 0.00, true);
    case "constant_rain"
    otherwise
        error('pf_config_vtilted:scenario', 'Unknown PF scenario: %s', scenario);
end
end

function PF = set_obs_enabled(PF, obs_id, enabled)
idx = [PF.observations.obs_id] == string(obs_id);
for i = find(idx)
    PF.observations(i).enabled = enabled;
end
end

function PF = set_obs_type_enabled(PF, type, enabled)
idx = [PF.observations.type] == string(type);
for i = find(idx)
    PF.observations(i).enabled = enabled;
end
end

function PF = set_param(PF, name, truth, enabled)
idx = [PF.parameters.name] == string(name);
PF.parameters(idx).truth = truth;
PF.parameters(idx).enabled = enabled;
end

function PF = set_param_prior(PF, name, baseline, truth, lower, upper, initial_std, perturb_std, enabled)
idx = [PF.parameters.name] == string(name);
PF.parameters(idx).baseline = baseline;
PF.parameters(idx).truth = truth;
PF.parameters(idx).lower = lower;
PF.parameters(idx).upper = upper;
PF.parameters(idx).initial_std = initial_std;
PF.parameters(idx).perturb_std = perturb_std;
PF.parameters(idx).enabled = enabled;
end

function S = pf_param(name, target_type, class_id, property, baseline, truth, lower, upper, initial_std, perturb_std, transform, enabled)
S = struct( ...
    'name', string(name), ...
    'target_type', string(target_type), ...
    'class_id', class_id, ...
    'property', string(property), ...
    'baseline', baseline, ...
    'truth', truth, ...
    'lower', lower, ...
    'upper', upper, ...
    'initial_std', initial_std, ...
    'perturb_std', perturb_std, ...
    'transform', string(transform), ...
    'enabled', logical(enabled));
end

function S = pf_obs(obs_id, type, row_fraction, col_fraction, sample_interval_min, sigma, units, enabled)
S = struct( ...
    'obs_id', string(obs_id), ...
    'type', string(type), ...
    'row_fraction', row_fraction, ...
    'col_fraction', col_fraction, ...
    'sample_interval_min', sample_interval_min, ...
    'sigma', sigma, ...
    'units', string(units), ...
    'enabled', logical(enabled));
end
