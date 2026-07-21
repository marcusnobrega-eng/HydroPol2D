function Results = run_vtilted_augmented_pf_full_hydropol2d(varargin)
%RUN_VTILTED_AUGMENTED_PF_FULL_HYDROPOL2D Particle filter using HydroPol2D.
%
% This is the DA runner that uses the real HydroPol2D preprocessing,
% Hydrological_Model, Groundwater_Module, and selected routing solver. It does
% not use the compact surrogate backend in run_vtilted_augmented_pf.m.

p = inputParser;
addParameter(p, 'SmokeTest', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'NParticles', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(p, 'NWindows', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(p, 'Routing', 'full_momentum', @(x) ischar(x) || isstring(x));
addParameter(p, 'Scenario', 'constant_rain', @(x) ischar(x) || isstring(x));
addParameter(p, 'TruthOnly', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'TruthSeriesIntervalMin', 5, @(x) isscalar(x) && x > 0);
parse(p, varargin{:});

PF = pf_config_vtilted(p.Results.Scenario);
PF.forecast_backend = 'hydropol2d_full';
PF.full_model = struct();
PF.full_model.routing = canonical_routing(p.Results.Routing);
PF.full_model.output_dir = fullfile(PF.runs_dir, 'FullHydroPol2D');

if p.Results.SmokeTest
    PF.n_particles = 2;
    PF.simulation_duration_min = PF.assimilation_window_min;
end
if ~isempty(p.Results.NParticles)
    PF.n_particles = round(p.Results.NParticles);
end
if ~isempty(p.Results.NWindows)
    PF.simulation_duration_min = round(p.Results.NWindows) * PF.assimilation_window_min;
end
if p.Results.TruthOnly
    PF.assimilation_window_min = p.Results.TruthSeriesIntervalMin;
end

rng(PF.random_seed);
if ~exist(PF.outputs_dir, 'dir'); mkdir(PF.outputs_dir); end

BaseState = hp2d_da_prepare_base_state(PF);
TruthTheta = local_make_theta(PF.parameters, 'truth');
TruthState = hp2d_da_apply_theta(BaseState, TruthTheta, PF);
TruthState = hp2d_da_force_initial_conditions(TruthState, TruthTheta);

if p.Results.TruthOnly
    Results = run_truth_series_only(PF, TruthState);
    return;
end

Particles = repmat(struct('theta', [], 'state', [], 'weight', [], ...
    'objective', [], 'log_likelihood', []), PF.n_particles, 1);
for i = 1:PF.n_particles
    Particles(i).theta = local_make_theta(PF.parameters, 'sample');
    Particles(i).state = hp2d_da_apply_theta(BaseState, Particles(i).theta, PF);
    Particles(i).state = hp2d_da_force_initial_conditions(Particles(i).state, Particles(i).theta);
    Particles(i).weight = 1 / PF.n_particles;
end

ObsRows = struct([]);
FitRows = struct([]);
EssRows = struct([]);
ParamRows = struct([]);
n_windows = round(PF.simulation_duration_min / PF.assimilation_window_min);
ParamRows = append_param_rows(ParamRows, PF, Particles, TruthTheta, 0, 0, [Particles.weight]');

for w = 1:n_windows
    t_end_min = w * PF.assimilation_window_min;
    TruthState = hp2d_da_forecast_window(TruthState, t_end_min);
    ObsTemplate = hp2d_da_observation_template(PF, TruthState, w, t_end_min);
    y_truth = hp2d_da_extract_observations(TruthState, ObsTemplate);
    y_obs = y_truth + PF.observation_noise_scale .* ObsTemplate.sigma .* randn(height(ObsTemplate), 1);
    y_obs = hp2d_da_bound_observations(y_obs, ObsTemplate);
    ObsTemplate.value = y_obs;

    sim = nan(PF.n_particles, height(ObsTemplate));
    for i = 1:PF.n_particles
        Particles(i).state = hp2d_da_forecast_window(Particles(i).state, t_end_min);
        sim(i, :) = hp2d_da_extract_observations(Particles(i).state, ObsTemplate)';
    end

    prior_weights = [Particles.weight]';
    [~, objective, logL, WeightInfo] = pf_paper_weights(PF, hp2d_da_domain_for_weights(TruthState), ObsTemplate, sim);
    weights = normalize_particle_weights(prior_weights .* WeightInfo.total_likelihood(:));
    Neff = 1 / sum(weights .^ 2);
    did_resample = Neff < PF.resample_threshold_fraction * PF.n_particles;

    for i = 1:PF.n_particles
        Particles(i).weight = weights(i);
        Particles(i).objective = objective(i);
        Particles(i).log_likelihood = logL(i);
    end

    ObsRows = append_obs_rows(ObsRows, ObsTemplate, y_truth, w);
    FitRows = append_fit_rows(FitRows, ObsTemplate, sim, weights, w);
    ParamRows = append_param_rows(ParamRows, PF, Particles, TruthTheta, w, t_end_min, weights);
    EssRows(end+1).window = w; %#ok<AGROW>
    EssRows(end).time_min = t_end_min;
    EssRows(end).effective_sample_size = Neff;
    EssRows(end).resampled = did_resample;
    EssRows(end).weighted_objective = sum(weights .* objective);
    write_pf_csvs(PF, ObsRows, FitRows, ParamRows, EssRows);

    if did_resample
        sigma_by_name = inflated_parameter_sigma_local(PF, Particles, weights);
        idx = pf_systematic_resample(weights);
        source_idx = idx;
        Particles = Particles(idx);
        copies_seen = zeros(PF.n_particles, 1);
        for i = 1:PF.n_particles
            copies_seen(source_idx(i)) = copies_seen(source_idx(i)) + 1;
            if copies_seen(source_idx(i)) > 1
                Particles(i).theta = pf_perturb_theta_inflated(Particles(i).theta, PF.parameters, sigma_by_name);
                Particles(i).state = hp2d_da_apply_theta(Particles(i).state, Particles(i).theta, PF);
            end
            Particles(i).weight = 1 / PF.n_particles;
        end
    end

    fprintf('Full HydroPol2D PF window %02d/%02d: Neff=%.2f, resampled=%d\n', ...
        w, n_windows, Neff, did_resample);
end

Results = struct();
Results.PF = PF;
Results.Observations = struct2table_safe_local(ObsRows);
Results.Observation_Fit = struct2table_safe_local(FitRows);
Results.Parameter_Evolution = struct2table_safe_local(ParamRows);
Results.Effective_Sample_Size = struct2table_safe_local(EssRows);
Results.Particles = Particles;

write_pf_csvs(PF, ObsRows, FitRows, ParamRows, EssRows);
save(fullfile(PF.outputs_dir, 'FullHydroPol2D_PF_Results.mat'), 'Results', '-v7.3');
end

function write_pf_csvs(PF, ObsRows, FitRows, ParamRows, EssRows)
% ponytail: cheap window checkpoint; add MAT checkpoints only if runs need restart.
writetable(struct2table_safe_local(ObsRows), fullfile(PF.outputs_dir, 'FullHydroPol2D_Synthetic_Observations.csv'));
writetable(struct2table_safe_local(FitRows), fullfile(PF.outputs_dir, 'FullHydroPol2D_Observation_Fit.csv'));
writetable(struct2table_safe_local(ParamRows), fullfile(PF.outputs_dir, 'FullHydroPol2D_Parameter_Evolution.csv'));
writetable(struct2table_safe_local(EssRows), fullfile(PF.outputs_dir, 'FullHydroPol2D_Effective_Sample_Size.csv'));
end

function Results = run_truth_series_only(PF, TruthState)
TruthRows = struct([]);
n_windows = ceil(PF.simulation_duration_min / PF.assimilation_window_min);
for w = 1:n_windows
    t_end_min = min(w * PF.assimilation_window_min, PF.simulation_duration_min);
    TruthState = hp2d_da_forecast_window(TruthState, t_end_min);
    ObsTemplate = hp2d_da_observation_template(PF, TruthState, w, t_end_min);
    y_truth = hp2d_da_extract_observations(TruthState, ObsTemplate);
    for i = 1:height(ObsTemplate)
        TruthRows(end+1).window = w; %#ok<AGROW>
        TruthRows(end).time_min = t_end_min;
        TruthRows(end).obs_id = ObsTemplate.obs_id(i);
        TruthRows(end).type = ObsTemplate.type(i);
        TruthRows(end).truth_value = y_truth(i);
        TruthRows(end).units = ObsTemplate.units(i);
    end
    fprintf('Full HydroPol2D truth window %02d/%02d at %.1f min\n', w, n_windows, t_end_min);
end
Results = struct();
Results.PF = PF;
Results.Truth = struct2table_safe_local(TruthRows);
Results.TruthState = TruthState;
writetable(Results.Truth, fullfile(PF.outputs_dir, 'FullHydroPol2D_Synthetic_Truth_Series.csv'));
save(fullfile(PF.outputs_dir, 'FullHydroPol2D_Synthetic_Truth_Series.mat'), 'Results', '-v7.3');
plot_truth_series(Results.Truth, PF.figures_dir);
end

function plot_truth_series(T, fig_dir)
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
obs_ids = unique(string(T.obs_id), 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 max(300, 260*numel(obs_ids))]);
tiledlayout(numel(obs_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(obs_ids)
    nexttile;
    idx = string(T.obs_id) == obs_ids(i);
    plot(T.time_min(idx), T.truth_value(idx), 'k-', 'LineWidth', 1.4);
    ylabel(sprintf('%s [%s]', char(obs_ids(i)), char(T.units(find(idx,1)))), 'Interpreter', 'none');
    grid on; box on;
end
xlabel('Time [min]');
exportgraphics(fig, fullfile(fig_dir, 'FullHydroPol2D_synthetic_truth_series.png'), 'Resolution', 200);
close(fig);
end

function State = hp2d_da_prepare_base_state(PF)
case_dir = fileparts(mfilename('fullpath'));
model_root = hydropol2d_find_root(case_dir);
functions_dir = fullfile(model_root, 'HydroPol2D_Functions');
static_dir = fullfile(model_root, 'Validation', 'Phase1_VTilted_Catchment', 'Static');
da_routing = PF.full_model.routing;
da_output_dir = PF.full_model.output_dir;
addpath(functions_dir, '-begin');
hydropol2d_add_runtime_paths(model_root);
addpath(case_dir);
setenv('HYDROPOL2D_VTILT_ROUTING', char(da_routing));
setenv('HYDROPOL2D_DA_TOTAL_DURATION_MIN', sprintf('%.15g', PF.simulation_duration_min));
setenv('HYDROPOL2D_DA_RAIN_MM_H', sprintf('%.15g', PF.forcing.rainfall_intensity_mm_h));
setenv('HYDROPOL2D_DA_RAIN_DURATION_MIN', sprintf('%.15g', PF.forcing.rainfall_duration_min));
setenv('HYDROPOL2D_DA_GW_TARGET_DT_MIN', sprintf('%.15g', PF.groundwater.target_dt_min));

Paths = make_da_paths(fullfile(da_output_dir, 'Base'), true);
InputPaths = make_vtilted_input_paths(static_dir);
input_data_bypass_script_path = fullfile(case_dir, 'pf_input_data_bypass_vtilted.m');
use_inputpaths_bypass = 1; %#ok<NASGU>
use_inputdata_bypass = 1; %#ok<NASGU>
clean_output_folder = true; %#ok<NASGU>
run_postprocessing = false; %#ok<NASGU>
enable_logging = false; %#ok<NASGU>
model_folder = ''; %#ok<NASGU>
GD = []; %#ok<NASGU>
resultsDir = Paths.Results; %#ok<NASGU>
export_root_dir = Paths.Root; %#ok<NASGU>

HydroPol2D_preprocessing;

% Production preprocessing does not retain its soil-class masks. Rebuild
% them from the DA source raster, then carry them inside the isolated state.
Soil_Properties.idx_soil = hp2d_da_source_soil_masks( ...
    InputPaths.SOIL_path, size(LULC_Properties.idx_lulc, 3), idx_nan);

da_routing_after_preprocessing = string(getenv('HYDROPOL2D_VTILT_ROUTING'));
running_control = hp2d_da_force_duration(running_control);
Rainfall_Parameters = hp2d_da_force_rainfall(Rainfall_Parameters, running_control);
flags = hp2d_da_force_process_flags(flags);
flags.flag_subgrid = 0;
flags.flag_full_momentum = double(da_routing_after_preprocessing == "full_momentum");
flags.flag_inertial = double(da_routing_after_preprocessing == "local_inertial");
flags.flag_CA = 0;
flags.flag_diffusive = 0;
flags.flag_kinematic = 0;
flags.flag_groundwater_modeling = 1;
flags.flag_groundwater_async = 1;
flags.flag_capillary_rise = 1;
flags.flag_baseflow = 1;
flags.flag_dashboard = 0;
gw_target_dt_min = str2double(getenv('HYDROPOL2D_DA_GW_TARGET_DT_MIN'));
if isfinite(gw_target_dt_min) && gw_target_dt_min > 0
    flags.groundwater_target_dt_min = gw_target_dt_min;
end

State = pack_caller_workspace();
State.DA_Base = hp2d_da_base_parameters(State, PF);
end

function running_control = hp2d_da_force_duration(running_control)
duration_min = str2double(getenv('HYDROPOL2D_DA_TOTAL_DURATION_MIN'));
if ~isfinite(duration_min) || duration_min <= 0
    return;
end
running_control.routing_time = duration_min;
if isfield(running_control, 'time_step_model') && running_control.time_step_model > 0
    running_control.steps = ceil(duration_min / running_control.time_step_model);
end
if isfield(running_control, 'record_time_maps') && running_control.record_time_maps > 0
    running_control.time_records = 0:running_control.record_time_maps:duration_min;
    if running_control.time_records(end) < duration_min
        running_control.time_records(end+1) = duration_min;
    end
    running_control.number_of_records = max(numel(running_control.time_records) - 1, 1);
end
if isfield(running_control, 'record_time_hydrographs') && running_control.record_time_hydrographs > 0
    running_control.time_record_hydrograph = 0:running_control.record_time_hydrographs:duration_min;
    if running_control.time_record_hydrograph(end) < duration_min
        running_control.time_record_hydrograph(end+1) = duration_min;
    end
    running_control.time_hydrograph = running_control.time_record_hydrograph(:);
end
end

function Rainfall_Parameters = hp2d_da_force_rainfall(Rainfall_Parameters, running_control)
rain_mm_h = str2double(getenv('HYDROPOL2D_DA_RAIN_MM_H'));
rain_duration_min = str2double(getenv('HYDROPOL2D_DA_RAIN_DURATION_MIN'));
if ~isfinite(rain_mm_h) || rain_mm_h < 0
    return;
end
if ~isfinite(rain_duration_min) || rain_duration_min <= 0
    rain_duration_min = running_control.routing_time;
end
if rain_duration_min >= running_control.routing_time
    Rainfall_Parameters.time_rainfall = [0; running_control.routing_time];
    Rainfall_Parameters.intensity_rainfall = [rain_mm_h; rain_mm_h];
else
    Rainfall_Parameters.time_rainfall = [0; rain_duration_min; running_control.routing_time];
    Rainfall_Parameters.intensity_rainfall = [rain_mm_h; rain_mm_h; 0];
end
Rainfall_Parameters.time_step_rainfall = min(rain_duration_min, running_control.routing_time);
Rainfall_Parameters.rainfall_duration = rain_duration_min;
Rainfall_Parameters.rainfall_end_time = rain_duration_min;
Rainfall_Parameters.n_obs_rainfall = numel(Rainfall_Parameters.time_rainfall);
Rainfall_Parameters.index_aggregation = 1;
end

function State = hp2d_da_forecast_window(State, t_end_min)
restore_workspace(State);
DA_window_end_h = t_end_min; %#ok<NASGU>
HydroPol2D_Main_While_DA_Window;
State = pack_caller_workspace();
end

function flags = hp2d_da_force_process_flags(flags)
flags.flag_rainfall = 1;
flags.flag_spatial_rainfall = 0;
flags.flag_input_rainfall_map = 0;
flags.flag_satellite_rainfall = 0;
flags.flag_inflow = 0;
flags.flag_stage_hydrograph = 0;
flags.flag_infiltration = 1;
flags.flag_ETP = 0;
flags.flag_input_ETP_map = 0;
flags.flag_waterquality = 0;
flags.flag_reservoir = 0;
flags.flag_boundary = 0;
end

function State = hp2d_da_apply_theta(State, theta, PF)
S = State;
if ~isfield(S, 'DA_Base')
    S.DA_Base = hp2d_da_base_parameters(S, PF);
end
is_initial_state = hp2d_da_is_initial_state(S);

for k = 1:numel(PF.base.lulc)
    mask = hp2d_da_lulc_mask(S, PF, k);
    if isfield(theta.lulc(k), 'manning_mult')
        S.LULC_Properties.roughness(mask) = S.DA_Base.lulc(k).manning_n * theta.lulc(k).manning_mult;
    end
    if isfield(theta.lulc(k), 'root_depth_mult')
        if ~isfield(S.LULC_Properties, 'root_depth_m') || isempty(S.LULC_Properties.root_depth_m)
            S.LULC_Properties.root_depth_m = nan(size(S.Elevation_Properties.elevation_cell));
        end
        S.LULC_Properties.root_depth_m(mask) = S.DA_Base.lulc(k).root_depth_m * theta.lulc(k).root_depth_mult;
    end
end

for k = 1:numel(PF.base.soil)
    mask = hp2d_da_soil_mask(S, PF, k);
    if isfield(theta.soil(k), 'ksat_surface_mult') && isfield(S.Soil_Properties, 'Ks_multiplier_near_surface')
        S.Soil_Properties.Ks_multiplier_near_surface(mask) = theta.soil(k).ksat_surface_mult;
    end
    if isfield(theta.soil(k), 'ksat_rootzone_mult') && isfield(S.Soil_Properties, 'Ks_multiplier_root_zone')
        S.Soil_Properties.Ks_multiplier_root_zone(mask) = theta.soil(k).ksat_rootzone_mult;
    end
    if isfield(theta.soil(k), 'storage_mult') && isfield(S.Soil_Properties, 'theta_sat')
        S.Soil_Properties.theta_sat(mask) = min(0.90, S.DA_Base.soil(k).theta_sat * theta.soil(k).storage_mult);
    end
    if is_initial_state && isfield(theta.soil(k), 'initial_sm_mult') && ...
            ~isfield(S.Soil_Properties, 'Layers')
        S.Soil_Properties.I_t(mask) = S.DA_Base.soil(k).I_t_mm * theta.soil(k).initial_sm_mult;
        S.Soil_Properties.I_p(mask) = S.Soil_Properties.I_t(mask);
    end
end

if isfield(theta.gw(1), 'ksat_mult') && isfield(S.Soil_Properties, 'ksat_gw')
    S.Soil_Properties.ksat_gw = S.DA_Base.gw.ksat_gw .* theta.gw(1).ksat_mult;
end
if isfield(theta.gw(1), 'specific_yield_mult') && isfield(S.Soil_Properties, 'Sy')
    S.Soil_Properties.Sy = S.DA_Base.gw.Sy .* theta.gw(1).specific_yield_mult;
end
if is_initial_state && isfield(S, 'BC_States') && isfield(S.BC_States, 'h_t')
    if isfield(theta, 'ic') && isfield(theta.ic(1), 'initial_wtd_mult')
        wtd = S.DA_Base.gw.wtd_m .* theta.ic(1).initial_wtd_mult;
        S.BC_States.h_t = S.Elevation_Properties.elevation_cell - wtd;
        S.BC_States.h_0 = S.BC_States.h_t;
    elseif isfield(theta.gw(1), 'initial_wtd_mult')
        wtd = S.DA_Base.gw.wtd_m .* theta.gw(1).initial_wtd_mult;
        S.BC_States.h_t = S.Elevation_Properties.elevation_cell - wtd;
        S.BC_States.h_0 = S.BC_States.h_t;
    end
end

if isfield(S.Soil_Properties, 'Layers') && isstruct(S.Soil_Properties.Layers)
    [S.Soil_Properties, S.LULC_Properties] = derive_layered_soil_profile( ...
        S.Soil_Properties, S.LULC_Properties, S.BC_States, ...
        S.Elevation_Properties.elevation_cell, S.idx_nan, struct());

    if is_initial_state
        for k = 1:numel(PF.base.soil)
            mask = hp2d_da_soil_mask(S, PF, k);
            if isfield(theta.soil(k), 'initial_sm_mult') && ...
                    isfield(S.DA_Base.soil(k), 'Layers')
                storage_fields = {'near_surface_storage_mm', 'root_zone_storage_mm', 'transmission_storage_mm'};
                for f = 1:numel(storage_fields)
                    field_name = storage_fields{f};
                    if isfield(S.Soil_Properties.Layers, field_name) && ...
                            isfield(S.DA_Base.soil(k).Layers, field_name)
                        S.Soil_Properties.Layers.(field_name)(mask) = ...
                            S.DA_Base.soil(k).Layers.(field_name)(mask) .* ...
                            theta.soil(k).initial_sm_mult;
                    end
                end
            end
        end
    end

    S.Soil_Properties = sync_layered_soil_storage(S.Soil_Properties, S.idx_nan);
end

S.DA_initial_conditions_applied = true;
State = S;
end

function S = hp2d_da_force_initial_conditions(S, theta)
if isfield(theta, 'ic') && isfield(theta.ic(1), 'initial_wtd_mult') && ...
        isfield(S, 'BC_States') && isfield(S.BC_States, 'h_t')
    wtd = S.DA_Base.gw.wtd_m .* theta.ic(1).initial_wtd_mult;
    S.BC_States.h_t = S.Elevation_Properties.elevation_cell - wtd;
    S.BC_States.h_0 = S.BC_States.h_t;
    if isfield(S.Soil_Properties, 'Layers') && isstruct(S.Soil_Properties.Layers)
        [S.Soil_Properties, S.LULC_Properties] = derive_layered_soil_profile( ...
            S.Soil_Properties, S.LULC_Properties, S.BC_States, ...
            S.Elevation_Properties.elevation_cell, S.idx_nan, struct());
        S.Soil_Properties = sync_layered_soil_storage(S.Soil_Properties, S.idx_nan);
    end
end
S.DA_initial_conditions_applied = true;
end

function tf = hp2d_da_is_initial_state(S)
if isfield(S, 'DA_initial_conditions_applied')
    tf = ~S.DA_initial_conditions_applied;
    return;
end
if ~isfield(S, 't') || isempty(S.t) || ~isfinite(S.t)
    tf = true;
    return;
end
tf = S.t <= 1e-3;
end

function Base = hp2d_da_base_parameters(S, PF)
Base = struct();
for k = 1:numel(PF.base.lulc)
    mask = hp2d_da_lulc_mask(S, PF, k);
    Base.lulc(k).manning_n = median(S.LULC_Properties.roughness(mask), 'omitnan');
    if isfield(S.LULC_Properties, 'root_depth_m')
        Base.lulc(k).root_depth_m = median(S.LULC_Properties.root_depth_m(mask), 'omitnan');
    else
        Base.lulc(k).root_depth_m = 1;
    end
end
for k = 1:numel(PF.base.soil)
    mask = hp2d_da_soil_mask(S, PF, k);
    Base.soil(k).theta_sat = median(S.Soil_Properties.theta_sat(mask), 'omitnan');
    Base.soil(k).I_t_mm = median(S.Soil_Properties.I_t(mask), 'omitnan');
    if isfield(S.Soil_Properties, 'Layers') && isstruct(S.Soil_Properties.Layers)
        storage_fields = {'near_surface_storage_mm', 'root_zone_storage_mm', 'transmission_storage_mm'};
        for f = 1:numel(storage_fields)
            field_name = storage_fields{f};
            if isfield(S.Soil_Properties.Layers, field_name)
                Base.soil(k).Layers.(field_name) = S.Soil_Properties.Layers.(field_name);
            end
        end
    end
end
Base.gw.ksat_gw = S.Soil_Properties.ksat_gw;
Base.gw.Sy = S.Soil_Properties.Sy;
Base.gw.wtd_m = S.Elevation_Properties.elevation_cell - S.BC_States.h_t;
end

function mask = hp2d_da_lulc_mask(S, PF, parameter_class)
source_classes = hp2d_da_source_classes(PF, 'lulc', ...
    size(S.LULC_Properties.idx_lulc, 3), parameter_class);
mask = any(S.LULC_Properties.idx_lulc(:, :, source_classes), 3);
end

function mask = hp2d_da_soil_mask(S, PF, parameter_class)
if isfield(S, 'idx_soil')
    idx_soil = S.idx_soil;
elseif isfield(S.Soil_Properties, 'idx_soil')
    idx_soil = S.Soil_Properties.idx_soil;
else
    error('DA state does not contain class-aware soil masks.');
end
source_classes = hp2d_da_source_classes(PF, 'soil', size(idx_soil, 3), parameter_class);
mask = any(idx_soil(:, :, source_classes), 3);
end

function idx_soil = hp2d_da_source_soil_masks(soil_path, n_source_classes, idx_nan)
[soil_codes, ~] = readgeoraster(soil_path);
soil_codes = double(soil_codes);
idx_soil = false([size(soil_codes), n_source_classes]);
for source_class = 1:n_source_classes
    idx_soil(:, :, source_class) = soil_codes == source_class;
end
idx_soil(repmat(idx_nan, 1, 1, n_source_classes)) = false;
end

function source_classes = hp2d_da_source_classes(PF, property_name, n_source_classes, parameter_class)
if ~isfield(PF, 'class_map') || ~isfield(PF.class_map, property_name)
    class_map = 1:n_source_classes;
else
    class_map = PF.class_map.(property_name)(:).';
end
if numel(class_map) ~= n_source_classes || any(~isfinite(class_map)) || ...
        any(class_map < 1) || any(abs(class_map - round(class_map)) > 1e-9)
    error('PF.class_map.%s must map each source class to a positive integer parameter class.', property_name);
end
source_classes = find(class_map == parameter_class);
if isempty(source_classes)
    error('No source %s classes map to particle parameter class %d.', property_name, parameter_class);
end
end

function T = hp2d_da_observation_template(PF, State, window, time_min)
rows = struct([]);
for i = 1:numel(PF.observations)
    obs = PF.observations(i);
    if ~obs.enabled
        continue;
    end
    rows(end+1).window = window; %#ok<AGROW>
    rows(end).time_min = time_min;
    rows(end).obs_id = obs.obs_id;
    rows(end).type = obs.type;
    rows(end).value = nan;
    rows(end).sigma = obs_sigma(PF, State, obs.type, obs.sigma);
    rows(end).units = obs.units;
    [rows(end).row, rows(end).col] = obs_location(State, obs);
end
T = struct2table_safe_local(rows);
end

function y = hp2d_da_extract_observations(State, Obs)
y = nan(height(Obs), 1);
for i = 1:height(Obs)
    switch char(Obs.type(i))
        case 'discharge'
            y(i) = sum(State.outlet_states.outlet_flow(:), 'omitnan') / 1000 / 3600 * State.Wshed_Properties.cell_area;
        case 'soil_moisture'
            r = Obs.row(i); c = Obs.col(i);
            y(i) = hp2d_da_surface_soil_moisture(State.Soil_Properties, r, c);
        case 'groundwater_depth'
            r = Obs.row(i); c = Obs.col(i);
            y(i) = State.Elevation_Properties.elevation_cell(r, c) - State.BC_States.h_t(r, c);
    end
end
end

function theta = hp2d_da_surface_soil_moisture(Soil_Properties, r, c)
if isfield(Soil_Properties, 'Layers') && isstruct(Soil_Properties.Layers) && ...
        isfield(Soil_Properties.Layers, 'near_surface_storage_mm')
    L = Soil_Properties.Layers;
    thick = max(L.near_surface_thickness_m(r, c), eps);
    theta = L.theta_r_near_surface(r, c) + ...
        L.near_surface_storage_mm(r, c) / 1000 / thick;
    theta = min(max(theta, L.theta_r_near_surface(r, c)), L.theta_sat_near_surface(r, c));
else
    theta = Soil_Properties.theta_r(r, c) + ...
        Soil_Properties.I_t(r, c) / max(Soil_Properties.Soil_Depth(r, c) * 1000, eps);
    theta = min(max(theta, Soil_Properties.theta_r(r, c)), Soil_Properties.theta_sat(r, c));
end
end

function y = hp2d_da_bound_observations(y, Obs)
for i = 1:height(Obs)
    switch char(Obs.type(i))
        case 'discharge'
            y(i) = max(y(i), 0);
        case 'soil_moisture'
            y(i) = min(max(y(i), 0), 1);
        case 'groundwater_depth'
            y(i) = max(y(i), 0);
    end
end
end

function [row, col] = obs_location(State, obs)
if obs.type == "discharge"
    row = nan; col = nan; return
end
[ny, nx] = size(State.Elevation_Properties.elevation_cell);
row = min(max(round(1 + obs.row_fraction * (ny - 1)), 1), ny);
col = min(max(round(1 + obs.col_fraction * (nx - 1)), 1), nx);
end

function sigma = obs_sigma(PF, State, type, fallback)
switch char(type)
    case 'discharge'
        area_m2 = sum(~isnan(State.Elevation_Properties.elevation_cell), 'all') * State.Wshed_Properties.cell_area;
        sigma = PF.weighting.sigma.discharge_mm_h / 1000 / 3600 * area_m2;
    case 'groundwater_depth'
        sigma = PF.weighting.sigma.groundwater_depth_m;
    case 'soil_moisture'
        sigma = PF.weighting.sigma.soil_moisture_m3m3;
    otherwise
        sigma = fallback;
end
end

function Domain = hp2d_da_domain_for_weights(State)
Domain.valid = ~isnan(State.Elevation_Properties.elevation_cell);
Domain.cell_area_m2 = State.Wshed_Properties.cell_area;
end

function S = pack_caller_workspace()
names = evalin('caller', 'who');
skip = ["S","State","names","skip","ans","t_end_min","DA_window_end_h","DA_is_final_window"];
S = struct();
for i = 1:numel(names)
    name = string(names{i});
    if any(name == skip)
        continue;
    end
    S.(names{i}) = evalin('caller', names{i});
end
end

function restore_workspace(S)
fn = fieldnames(S);
for i = 1:numel(fn)
    assignin('caller', fn{i}, S.(fn{i}));
end
end

function theta = local_make_theta(specs, mode)
mode = string(mode);
targets = [specs.target_type];
theta.lulc = repmat(struct(), 1, max([specs(targets == "LULC").class_id]));
theta.soil = repmat(struct(), 1, max([specs(targets == "SOIL").class_id]));
theta.gw = repmat(struct(), 1, max([specs(targets == "GW").class_id]));
theta.ic = repmat(struct(), 1, max([specs(targets == "IC").class_id]));
for i = 1:numel(specs)
    spec = specs(i);
    if mode == "truth"
        value = spec.truth;
    elseif ~spec.enabled
        value = spec.baseline;
    else
        value = min(max(spec.baseline * (1 + sqrt(spec.initial_std) * randn()), spec.lower), spec.upper);
    end
    theta = set_theta(theta, spec, value);
end
end

function theta = set_theta(theta, spec, value)
switch char(lower(spec.target_type))
    case 'lulc'
        theta.lulc(spec.class_id).(char(spec.property)) = value;
    case 'soil'
        theta.soil(spec.class_id).(char(spec.property)) = value;
    case 'gw'
        theta.gw(spec.class_id).(char(spec.property)) = value;
    case 'ic'
        theta.ic(spec.class_id).(char(spec.property)) = value;
end
end

function sigma_by_name = inflated_parameter_sigma_local(PF, Particles, weights)
weights = normalize_particle_weights(weights);
specs = PF.parameters([PF.parameters.enabled]);
[names, first] = theta_values_local(Particles(1).theta, specs);
values = nan(numel(Particles), numel(first));
values(1, :) = first;
for i = 2:numel(Particles)
    [~, values(i, :)] = theta_values_local(Particles(i).theta, specs);
end
sigma_by_name = struct();
for p = 1:numel(names)
    mu = sum(weights .* values(:, p), 'omitnan');
    sd = sqrt(sum(weights .* (values(:, p) - mu).^2, 'omitnan'));
    sigma_by_name.(matlab.lang.makeValidName(char(names(p)))) = PF.resampling.inflation_lambda * sd;
end
end

function [names, values] = theta_values_local(theta, specs)
names = strings(1, numel(specs));
values = nan(1, numel(specs));
for i = 1:numel(specs)
    spec = specs(i);
    names(i) = spec.name;
    switch char(lower(spec.target_type))
        case 'lulc'
            values(i) = theta.lulc(spec.class_id).(char(spec.property));
        case 'soil'
            values(i) = theta.soil(spec.class_id).(char(spec.property));
        case 'gw'
            values(i) = theta.gw(spec.class_id).(char(spec.property));
        case 'ic'
            values(i) = theta.ic(spec.class_id).(char(spec.property));
    end
end
end

function Rows = append_obs_rows(Rows, Obs, truth, window)
for i = 1:height(Obs)
    Rows(end+1).window = window; %#ok<AGROW>
    Rows(end).time_min = Obs.time_min(i);
    Rows(end).obs_id = Obs.obs_id(i);
    Rows(end).type = Obs.type(i);
    Rows(end).truth_value = truth(i);
    Rows(end).observed = Obs.value(i);
    Rows(end).sigma = Obs.sigma(i);
    Rows(end).units = Obs.units(i);
end
end

function Rows = append_fit_rows(Rows, Obs, sim, weights, window)
for j = 1:height(Obs)
    for i = 1:size(sim, 1)
        Rows(end+1).window = window; %#ok<AGROW>
        Rows(end).time_min = Obs.time_min(j);
        Rows(end).particle_id = i;
        Rows(end).obs_id = Obs.obs_id(j);
        Rows(end).type = Obs.type(j);
        Rows(end).observed = Obs.value(j);
        Rows(end).simulated = sim(i, j);
        Rows(end).weight = weights(i);
        Rows(end).units = Obs.units(j);
    end
end
end

function Rows = append_param_rows(Rows, PF, Particles, TruthTheta, window, time_min, weights)
specs = PF.parameters([PF.parameters.enabled]);
[names, truth] = theta_values_local(TruthTheta, specs);
for i = 1:numel(Particles)
    [~, values] = theta_values_local(Particles(i).theta, specs);
    for p = 1:numel(names)
        Rows(end+1).window = window; %#ok<AGROW>
        Rows(end).time_min = time_min;
        Rows(end).particle_id = i;
        Rows(end).parameter = names(p);
        Rows(end).value = values(p);
        Rows(end).truth_value = truth(p);
        Rows(end).weight = weights(i);
    end
end
end

function weights = normalize_particle_weights(weights)
weights = double(weights(:));
weights(~isfinite(weights) | weights < 0) = 0;
if sum(weights) <= 0
    weights(:) = 1 / numel(weights);
else
    weights = weights ./ sum(weights);
end
end

function mode = canonical_routing(mode)
mode = lower(strtrim(string(mode)));
switch mode
    case {"full_momentum", "fullmomentum", "fm"}
        mode = "full_momentum";
    case {"local_inertial", "localinertial", "li"}
        mode = "local_inertial";
    otherwise
        error('Only full_momentum and local_inertial are wired for full HydroPol2D DA currently.');
end
end

function InputPaths = make_vtilted_input_paths(static_dir)
InputPaths = struct();
InputPaths.case_root = fileparts(static_dir);
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_dir, 'DTB.tif');
InputPaths.GW_table_path = fullfile(static_dir, 'GW_table.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_dir, 'Initial_SM.tif');
InputPaths.LAI_path = fullfile(static_dir, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_dir, 'Albedo.tif');
InputPaths.Subgrid_DEM_path = '';
InputPaths.RiverWidths_path = '';
InputPaths.RiverDepths_path = '';
InputPaths.Warmup_Depth_path = '';
InputPaths.Initial_Buildup_path = '';
InputPaths.B1_path = '';
InputPaths.B2_path = '';
InputPaths.W1_path = '';
InputPaths.W2_path = '';
InputPaths.Rainfall_Rasters_Folder = '';
InputPaths.Transpiration_Rasters_Folder = '';
InputPaths.Evaporation_Rasters_Folder = '';
InputPaths.Rainfall_Raster_Files = {};
InputPaths.Transpiration_Raster_Files = {};
InputPaths.Evaporation_Raster_Files = {};
InputPaths.Inflow_Hydrograph_CSV = '';
InputPaths.Stage_Hydrograph_CSV = '';
InputPaths.Observed_Gauges_CSV = '';
InputPaths.ETP_input_spreadsheet = '';
InputPaths.Rainfall_Timeseries_File = '/Users/mngomes/Downloads/Rainfall_Intensity_Data.xlsx';
InputPaths.Outlet_Cells_CSV = '';
end

function Paths = make_da_paths(root_dir, clean_output)
if exist(root_dir, 'dir') && clean_output
    rmdir(root_dir, 's');
end
Paths.Root = root_dir;
Paths.Results = fullfile(root_dir, 'Modeling_Results');
Paths.Temp = fullfile(root_dir, 'Temporary_Files');
Paths.Logs = fullfile(root_dir, 'Logs');
Paths.FigPDF = fullfile(Paths.Results, 'Figures_PDF');
Paths.FigFIG = fullfile(Paths.Results, 'Figures_FIG');
Paths.Tables = fullfile(Paths.Results, 'Tables_CSV');
Paths.RastersWD = fullfile(Paths.Results, 'Rasters_Water_Depths');
Paths.RastersWSE = fullfile(Paths.Results, 'Rasters_WSE');
Paths.RastersStatic = fullfile(Paths.Results, 'Rasters_Static');
Paths.RastersVelocity = fullfile(Paths.Results, 'Rasters_Velocity');
Paths.RastersHazard = fullfile(Paths.Results, 'Rasters_Hazard');
Paths.WQMaps = fullfile(Paths.Results, 'Rasters_WQ');
Paths.HRMaps = fullfile(Paths.Results, 'Rasters_Human_Risk');
Paths.Anim = fullfile(Paths.Results, 'GIFs_MP4');
Paths.Shapes = fullfile(Paths.Results, 'Shapefiles');
fields = fieldnames(Paths);
for i = 1:numel(fields)
    if ~exist(Paths.(fields{i}), 'dir')
        mkdir(Paths.(fields{i}));
    end
end
end

function T = struct2table_safe_local(S)
if isempty(S)
    T = table();
else
    T = struct2table(S);
end
end
