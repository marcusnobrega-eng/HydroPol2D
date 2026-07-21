function Results = run_vtilted_augmented_pf(varargin)
%RUN_VTILTED_AUGMENTED_PF Augmented particle filter prototype for V-tilted.
%
% The prototype is isolated from production HydroPol2D source. Particles carry
% theta + state explicitly; v1 updates theta and resamples state, while the
% analysis update for state is identity.
%
% Optional name-value pairs:
%   'SmokeTest'  true/false
%   'NParticles' positive scalar
%   'NWindows'   positive scalar
%   'MakePlots'  true/false

p = inputParser;
addParameter(p, 'SmokeTest', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'NParticles', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(p, 'NWindows', [], @(x) isempty(x) || (isscalar(x) && x > 0));
addParameter(p, 'MakePlots', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

PF = pf_config_vtilted();
if p.Results.SmokeTest
    PF.n_particles = 4;
    PF.simulation_duration_min = 2 * PF.assimilation_window_min;
end
if ~isempty(p.Results.NParticles)
    PF.n_particles = round(p.Results.NParticles);
end
if ~isempty(p.Results.NWindows)
    PF.simulation_duration_min = round(p.Results.NWindows) * PF.assimilation_window_min;
end
PF.make_plots = logical(p.Results.MakePlots);
validate_pf_timing(PF);

rng(PF.random_seed);
prepare_output_dirs(PF);

Domain = load_vtilted_domain(PF);
write_two_class_static_copy(PF, Domain);

truth_theta = make_theta(PF.parameters, 'truth');
truth_state = make_initial_state(PF, Domain, truth_theta);
[Truth, Observations] = generate_truth_observations(PF, Domain, truth_state, truth_theta);

Particles = initialize_particles(PF, Domain);
n_windows = PF.simulation_duration_min / PF.assimilation_window_min;
if abs(n_windows - round(n_windows)) > 1e-9
    error('Simulation duration must be an integer multiple of assimilation_window_min.');
end
n_windows = round(n_windows);

ParamRows = struct([]);
ObjectiveRows = struct([]);
ObsRows = struct([]);
EssRows = struct([]);
PosteriorRows = struct([]);
StateRows = struct([]);

initial_weights = ones(PF.n_particles, 1) ./ PF.n_particles;
[ParamRows, PosteriorRows] = append_parameter_rows( ...
    ParamRows, PosteriorRows, PF, Particles, truth_theta, 0, 0, initial_weights);

for w = 1:n_windows
    t_start = (w - 1) * PF.assimilation_window_min;
    t_end = w * PF.assimilation_window_min;
    obs_w = Observations(Observations.window == w, :);

    sim_matrix = nan(PF.n_particles, height(obs_w));

    for i = 1:PF.n_particles
        [forecast_state, window_diag, y_sim] = forecast_window(PF, Domain, Particles(i).state, Particles(i).theta, t_start, t_end, obs_w);
        analysis_state = pf_apply_state_update_identity(forecast_state, obs_w, window_diag);
        Particles(i).state = analysis_state;
        sim_matrix(i, :) = y_sim(:)';
    end

    prior_weights = [Particles.weight]';
    [~, J, logL, WeightInfo] = pf_paper_weights(PF, Domain, obs_w, sim_matrix);
    weights = normalize_particle_weights(prior_weights .* WeightInfo.total_likelihood(:));
    for i = 1:PF.n_particles
        Particles(i).weight = weights(i);
        Particles(i).objective = J(i);
        Particles(i).log_likelihood = logL(i);
        Particles(i).history(w).time_min = t_end;
        Particles(i).history(w).objective = J(i);
        Particles(i).history(w).log_likelihood = logL(i);
        Particles(i).history(w).outlet_q_m3s = Particles(i).state.diagnostics.last_outlet_q_m3s;
    end

    Neff = 1 / sum(weights .^ 2);
    did_resample = t_end >= PF.spinup_min && Neff < PF.resample_threshold_fraction * PF.n_particles;

    [ParamRows, PosteriorRows] = append_parameter_rows( ...
        ParamRows, PosteriorRows, PF, Particles, truth_theta, w, t_end, weights);
    ObjectiveRows = append_objective_rows(ObjectiveRows, PF, Particles, w, t_end, weights, J, logL, WeightInfo);
    ObsRows = append_observation_rows(ObsRows, PF, obs_w, sim_matrix, weights, w, t_end);
    EssRows = append_ess_row(EssRows, w, t_end, Neff, did_resample, sum(weights .* J));
    StateRows = append_state_rows(StateRows, PF, Particles, w, t_end, weights);

    if did_resample
        sigma_by_name = inflated_parameter_sigma(PF, Particles, weights);
        idx = pf_systematic_resample(weights);
        source_idx = idx;
        Particles = Particles(idx);
        copies_seen = zeros(PF.n_particles, 1);
        for i = 1:PF.n_particles
            Particles(i).id = i;
            copies_seen(source_idx(i)) = copies_seen(source_idx(i)) + 1;
            if copies_seen(source_idx(i)) > 1
                Particles(i).theta = pf_perturb_theta_inflated(Particles(i).theta, PF.parameters, sigma_by_name);
            end
            Particles(i).weight = 1 / PF.n_particles;
        end
    end

    fprintf('PF window %02d/%02d complete: Neff=%.2f, resampled=%d, weighted J=%.3f\n', ...
        w, n_windows, Neff, did_resample, sum(weights .* J));
end

Results = struct();
Results.PF = PF;
Results.Domain = compact_domain_for_save(Domain);
Results.Truth = Truth;
Results.Observations = Observations;
Results.Particles = Particles;
Results.Parameter_Evolution = struct2table_safe(ParamRows);
Results.Objective_Function = struct2table_safe(ObjectiveRows);
Results.Effective_Sample_Size = struct2table_safe(EssRows);
Results.Observation_Fit = struct2table_safe(ObsRows);
Results.Observation_Ensemble_Summary = observation_ensemble_summary(Results.Observation_Fit);
Results.Posterior_Summary = struct2table_safe(PosteriorRows);
Results.Posterior_Summary_Physical = physical_posterior_table(PF, Results.Posterior_Summary);
Results.State_Summary = struct2table_safe(StateRows);

write_results(PF, Results);
if PF.make_plots
    make_pf_plots(PF, Results);
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

function prepare_output_dirs(PF)
dirs = {PF.runs_dir, PF.outputs_dir, PF.figures_dir};
for i = 1:numel(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i});
    end
end
end

function validate_pf_timing(PF)
if rem(PF.assimilation_window_min, PF.dt_min) ~= 0
    error('assimilation_window_min must be an integer multiple of dt_min.');
end
if rem(PF.simulation_duration_min, PF.assimilation_window_min) ~= 0
    error('simulation_duration_min must be an integer multiple of assimilation_window_min.');
end
for i = 1:numel(PF.observations)
    obs = PF.observations(i);
    if ~obs.enabled
        continue;
    end
    dt_obs = obs.sample_interval_min;
    if rem(dt_obs, PF.dt_min) ~= 0
        error('Observation %s sample_interval_min must be an integer multiple of dt_min.', obs.obs_id);
    end
    a_over_o = PF.assimilation_window_min / dt_obs;
    o_over_a = dt_obs / PF.assimilation_window_min;
    if abs(a_over_o - round(a_over_o)) > 1e-9 && abs(o_over_a - round(o_over_a)) > 1e-9
        error('Observation %s interval must divide or be divided by assimilation_window_min.', obs.obs_id);
    end
end
end

function Domain = load_vtilted_domain(PF)
[DEM, R] = readgeoraster(fullfile(PF.static_source_dir, 'DEM.tif'));
DEM = double(DEM);
[LULC0, ~] = readgeoraster(fullfile(PF.static_source_dir, 'LULC.tif'));
[SOIL0, ~] = readgeoraster(fullfile(PF.static_source_dir, 'SOIL.tif'));
[DTB, ~] = readgeoraster(fullfile(PF.static_source_dir, 'DTB.tif'));
[Zone0, ~] = readgeoraster(fullfile(PF.static_source_dir, 'Zone_ID.tif'));

valid = isfinite(DEM);
LULC = ones(size(DEM));
LULC(Zone0 ~= 1) = 2;
LULC(~valid) = nan;
SOIL = ones(size(DEM));
SOIL(Zone0 ~= 1) = 2;
SOIL(~valid) = nan;

dx = raster_cellsize(R);
[receiver, distance, route_order] = build_receivers(DEM, valid, dx);

Domain = struct();
Domain.DEM = DEM;
Domain.R = R;
Domain.valid = valid;
Domain.LULC_original = double(LULC0);
Domain.SOIL_original = double(SOIL0);
Domain.Zone_ID = double(Zone0);
Domain.LULC = LULC;
Domain.SOIL = SOIL;
Domain.DTB = double(DTB);
Domain.dx = dx;
Domain.cell_area_m2 = dx ^ 2;
Domain.ny = size(DEM, 1);
Domain.nx = size(DEM, 2);
Domain.receiver = receiver;
Domain.receiver_distance_m = distance;
Domain.route_order = route_order;
Domain.obs_locations = observation_locations(PF, Domain);
end

function dx = raster_cellsize(R)
if isprop(R, 'CellExtentInWorldX')
    dx = abs(R.CellExtentInWorldX);
elseif isfield(R, 'CellExtentInWorldX')
    dx = abs(R.CellExtentInWorldX);
else
    dx = 20;
end
end

function [receiver, distance, route_order] = build_receivers(DEM, valid, dx)
[ny, nx] = size(DEM);
receiver = zeros(ny, nx);
distance = nan(ny, nx);
offsets = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1];

for r = 1:ny
    for c = 1:nx
        if ~valid(r, c)
            continue;
        end
        best_slope = -inf;
        best_idx = 0;
        best_dist = dx;
        for k = 1:size(offsets, 1)
            rr = r + offsets(k, 1);
            cc = c + offsets(k, 2);
            if rr < 1 || rr > ny || cc < 1 || cc > nx || ~valid(rr, cc)
                continue;
            end
            dist = dx * hypot(offsets(k, 1), offsets(k, 2));
            slope = (DEM(r, c) - DEM(rr, cc)) / dist;
            if slope > best_slope
                best_slope = slope;
                best_idx = sub2ind([ny, nx], rr, cc);
                best_dist = dist;
            end
        end
        if best_slope > 0
            receiver(r, c) = best_idx;
            distance(r, c) = best_dist;
        elseif r == 1
            receiver(r, c) = 0;
            distance(r, c) = dx;
        else
            receiver(r, c) = 0;
            distance(r, c) = dx;
        end
    end
end

[~, order] = sort(DEM(:), 'descend', 'MissingPlacement', 'last');
route_order = order(valid(order));
end

function locations = observation_locations(PF, Domain)
locations = struct([]);
for i = 1:numel(PF.observations)
    obs = PF.observations(i);
    locations(i).obs_id = obs.obs_id;
    locations(i).type = obs.type;
    if obs.type == "discharge"
        locations(i).row = nan;
        locations(i).col = nan;
        locations(i).idx = nan;
    else
        row = min(max(round(1 + obs.row_fraction * (Domain.ny - 1)), 1), Domain.ny);
        col = min(max(round(1 + obs.col_fraction * (Domain.nx - 1)), 1), Domain.nx);
        locations(i).row = row;
        locations(i).col = col;
        locations(i).idx = sub2ind([Domain.ny, Domain.nx], row, col);
    end
end
end

function write_two_class_static_copy(PF, Domain)
static_dir = fullfile(PF.runs_dir, 'Static_TwoClass');
if ~exist(static_dir, 'dir')
    mkdir(static_dir);
end
write_tif_safe(fullfile(static_dir, 'DEM.tif'), Domain.DEM, Domain.R);
write_tif_safe(fullfile(static_dir, 'LULC_two_class.tif'), Domain.LULC, Domain.R);
write_tif_safe(fullfile(static_dir, 'SOIL_two_class.tif'), Domain.SOIL, Domain.R);
write_tif_safe(fullfile(static_dir, 'Zone_ID.tif'), Domain.Zone_ID, Domain.R);
end

function write_tif_safe(path, A, R)
try
    geotiffwrite(path, single(A), R);
catch
    save(strrep(path, '.tif', '.mat'), 'A', 'R');
end
end

function theta = make_theta(specs, mode)
mode = string(mode);
theta = struct();
theta.lulc = repmat(struct(), 1, max_class(specs, "LULC"));
theta.soil = repmat(struct(), 1, max_class(specs, "SOIL"));
theta.gw = repmat(struct(), 1, max_class(specs, "GW"));
theta.ic = repmat(struct(), 1, max_class(specs, "IC"));

for i = 1:numel(specs)
    spec = specs(i);
    if mode == "truth"
        value = spec.truth;
    elseif mode == "baseline"
        value = spec.baseline;
    elseif mode == "sample"
        value = sample_prior_value(spec);
    else
        error('Unknown theta mode: %s', mode);
    end
    theta = set_theta(theta, spec, value);
end
end

function n = max_class(specs, target)
idx = [specs.target_type] == target;
if any(idx)
    n = max([specs(idx).class_id]);
else
    n = 1;
end
end

function value = sample_prior_value(spec)
if spec.transform == "paper_multiplicative"
    value = spec.baseline * (1 + sqrt(spec.initial_std) * randn());
elseif spec.transform == "log"
    value = exp(log(max(spec.baseline, eps)) + spec.initial_std * randn());
else
    value = spec.baseline + spec.initial_std * randn();
end
value = min(max(value, spec.lower), spec.upper);
end

function theta = set_theta(theta, spec, value)
property = char(spec.property);
switch char(lower(spec.target_type))
    case 'lulc'
        theta.lulc(spec.class_id).(property) = value;
    case 'soil'
        theta.soil(spec.class_id).(property) = value;
    case 'gw'
        theta.gw(spec.class_id).(property) = value;
    case 'ic'
        theta.ic(spec.class_id).(property) = value;
end
end

function Particles = initialize_particles(PF, Domain)
Particles = repmat(empty_particle(), PF.n_particles, 1);
for i = 1:PF.n_particles
    theta = make_theta(PF.parameters, 'sample');
    state = make_initial_state(PF, Domain, theta);
    Particles(i).id = i;
    Particles(i).theta = theta;
    Particles(i).state = state;
    Particles(i).weight = 1 / PF.n_particles;
    Particles(i).log_likelihood = 0;
    Particles(i).objective = nan;
    Particles(i).history = struct([]);
end
end

function P = empty_particle()
P = struct('id', [], 'theta', [], 'state', [], 'weight', [], ...
    'log_likelihood', [], 'objective', [], 'history', struct([]));
end

function state = make_initial_state(PF, Domain, theta)
state = struct();
state.time = 0;
state.depths = struct();
state.depths.d_t = zeros(Domain.ny, Domain.nx);
state.depths.d_tot = state.depths.d_t;
state.depths.d_p = state.depths.d_t;

state.Soil_Properties = struct();
state.Soil_Properties.Soil_Depth = Domain.DTB;
state.Soil_Properties.Layers = struct();
state.Soil_Properties.idx_soil = Domain.SOIL;

state.BC_States = struct();
state.Hydro_States = struct('f', zeros(Domain.ny, Domain.nx));
state.GW_States = struct();
state.Snow_Properties = struct();
state.WQ_States = struct();
state.Reservoir_Data = struct();
state.cumulative_fluxes = struct( ...
    'rain_m3', 0, ...
    'infiltration_m3', 0, ...
    'recharge_m3', 0, ...
    'exfiltration_m3', 0, ...
    'outlet_m3', 0);
state.diagnostics = struct( ...
    'last_outlet_q_m3s', 0, ...
    'last_mean_velocity_m_s', 0, ...
    'last_rainfall_mm_h', 0, ...
    'last_infiltration_mm_h', 0, ...
    'last_recharge_mm_h', 0, ...
    'last_exfiltration_mm_h', 0);

state = apply_theta_to_state(PF, Domain, state, theta, true);
end

function state = apply_theta_to_state(PF, Domain, state, theta, is_initial)
roughness = nan(Domain.ny, Domain.nx);
root_depth = nan(Domain.ny, Domain.nx);
for k = 1:numel(PF.base.lulc)
    idx = Domain.LULC == k;
    roughness(idx) = PF.base.lulc(k).manning_n * theta.lulc(k).manning_mult;
    root_depth(idx) = PF.base.lulc(k).root_depth_m * theta.lulc(k).root_depth_mult;
end

ksat_surface = nan(Domain.ny, Domain.nx);
ksat_rootzone = nan(Domain.ny, Domain.nx);
storage_capacity = nan(Domain.ny, Domain.nx);
theta_s = nan(Domain.ny, Domain.nx);
theta_r = nan(Domain.ny, Domain.nx);
initial_fraction = nan(Domain.ny, Domain.nx);
for k = 1:numel(PF.base.soil)
    idx = Domain.SOIL == k;
    ksat_surface(idx) = PF.base.soil(k).ksat_surface_mm_h * theta.soil(k).ksat_surface_mult;
    ksat_rootzone(idx) = PF.base.soil(k).ksat_rootzone_mm_h * theta.soil(k).ksat_rootzone_mult;
    storage_capacity(idx) = PF.base.soil(k).storage_capacity_mm * theta.soil(k).storage_mult;
    theta_s(idx) = PF.base.soil(k).theta_s;
    theta_r(idx) = PF.base.soil(k).theta_r;
    initial_fraction(idx) = min(max(PF.base.soil(k).initial_sm_fraction * theta.soil(k).initial_sm_mult, 0), 1);
end

near_cap = min(25, 0.20 .* storage_capacity);
root_fraction = min(max(root_depth ./ max(state.Soil_Properties.Soil_Depth, 0.1), 0.20), 0.70);
root_cap = min(max(storage_capacity - near_cap, 0), root_fraction .* storage_capacity);
trans_cap = max(storage_capacity - near_cap - root_cap, 0);

state.Soil_Properties.roughness = roughness;
state.Soil_Properties.root_depth_m = root_depth;
state.Soil_Properties.ksat_surface_mm_h = ksat_surface;
state.Soil_Properties.ksat_rootzone_mm_h = ksat_rootzone;
state.Soil_Properties.storage_capacity_mm = storage_capacity;
state.Soil_Properties.theta_s = theta_s;
state.Soil_Properties.theta_r = theta_r;
state.Soil_Properties.Layers.near_surface_capacity_mm = near_cap;
state.Soil_Properties.Layers.root_zone_capacity_mm = root_cap;
state.Soil_Properties.Layers.transmission_capacity_mm = trans_cap;

if is_initial
    total_initial = initial_fraction .* storage_capacity;
    near = min(total_initial, 0.70 .* near_cap);
    rem = max(total_initial - near, 0);
    root = min(rem, root_cap);
    rem = max(rem - root, 0);
    trans = min(rem, trans_cap);

    state.Soil_Properties.Layers.near_surface_storage_mm = near;
    state.Soil_Properties.Layers.root_zone_storage_mm = root;
    state.Soil_Properties.Layers.transmission_storage_mm = trans;

    wtd = PF.base.gw(1).initial_wtd_m * theta.ic(1).initial_wtd_mult * ones(Domain.ny, Domain.nx);
    wtd(~Domain.valid) = nan;
    state.BC_States.gw_depth_m = wtd;
    state.BC_States.h_t = Domain.DEM - wtd;
    state.BC_States.h_0 = state.BC_States.h_t;
else
    [state, excess_mm] = clamp_soil_layers(state);
    state.depths.d_t = state.depths.d_t + excess_mm / 1000;
end

state.GW_States.ksat_mm_h = PF.base.gw(1).ksat_mm_h * theta.gw(1).ksat_mult * ones(Domain.ny, Domain.nx);
state.GW_States.Sy = PF.base.gw(1).specific_yield * theta.gw(1).specific_yield_mult * ones(Domain.ny, Domain.nx);
state.GW_States.ksat_mm_h(~Domain.valid) = nan;
state.GW_States.Sy(~Domain.valid) = nan;
state.BC_States.h_t = Domain.DEM - state.BC_States.gw_depth_m;
end

function [state, excess_mm] = clamp_soil_layers(state)
L = state.Soil_Properties.Layers;
excess_mm = zeros(size(L.near_surface_storage_mm));

over = max(L.near_surface_storage_mm - L.near_surface_capacity_mm, 0);
L.near_surface_storage_mm = L.near_surface_storage_mm - over;
L.root_zone_storage_mm = L.root_zone_storage_mm + over;

over = max(L.root_zone_storage_mm - L.root_zone_capacity_mm, 0);
L.root_zone_storage_mm = L.root_zone_storage_mm - over;
L.transmission_storage_mm = L.transmission_storage_mm + over;

over = max(L.transmission_storage_mm - L.transmission_capacity_mm, 0);
L.transmission_storage_mm = L.transmission_storage_mm - over;
excess_mm = excess_mm + over;

state.Soil_Properties.Layers = L;
end

function [Truth, Observations] = generate_truth_observations(PF, Domain, state, theta)
n_windows = round(PF.simulation_duration_min / PF.assimilation_window_min);
Schedule = observation_schedule(PF, Domain);
TruthRows = struct([]);
ObsRows = struct([]);

for w = 1:n_windows
    t_start = (w - 1) * PF.assimilation_window_min;
    t_end = w * PF.assimilation_window_min;
    obs_template = Schedule(Schedule.window == w, :);
    [state, window_diag, y_true] = forecast_window(PF, Domain, state, theta, t_start, t_end, obs_template);

    for j = 1:height(obs_template)
        obs_val = y_true(j) + obs_template.sigma(j) * randn();
        TruthRows(end+1).window = w; %#ok<AGROW>
        TruthRows(end).time_min = obs_template.time_min(j);
        TruthRows(end).obs_id = obs_template.obs_id(j);
        TruthRows(end).type = obs_template.type(j);
        TruthRows(end).truth_value = y_true(j);
        TruthRows(end).sample_interval_min = obs_template.sample_interval_min(j);
        TruthRows(end).units = obs_template.units(j);
        TruthRows(end).outlet_q_m3s = window_diag.outlet_q_m3s;

        ObsRows(end+1).window = w; %#ok<AGROW>
        ObsRows(end).time_min = obs_template.time_min(j);
        ObsRows(end).obs_id = obs_template.obs_id(j);
        ObsRows(end).type = obs_template.type(j);
        ObsRows(end).value = obs_val;
        ObsRows(end).sample_interval_min = obs_template.sample_interval_min(j);
        ObsRows(end).sigma = obs_template.sigma(j);
        ObsRows(end).units = obs_template.units(j);
    end
end

Truth = struct2table_safe(TruthRows);
Observations = struct2table_safe(ObsRows);
end

function T = observation_schedule(PF, Domain)
rows = struct([]);
for i = 1:numel(PF.observations)
    obs = PF.observations(i);
    if ~obs.enabled
        continue;
    end
    times = obs.sample_interval_min:obs.sample_interval_min:PF.simulation_duration_min;
    for k = 1:numel(times)
        rows(end+1).window = ceil(times(k) / PF.assimilation_window_min); %#ok<AGROW>
        rows(end).time_min = times(k);
        rows(end).obs_id = obs.obs_id;
        rows(end).type = obs.type;
        rows(end).value = nan;
        rows(end).sample_interval_min = obs.sample_interval_min;
        rows(end).sigma = observation_sigma_in_output_units(PF, Domain, obs.type, obs.sigma);
        rows(end).units = obs.units;
    end
end
T = struct2table_safe(rows);
if ~isempty(T)
    T = sortrows(T, {'window', 'time_min', 'obs_id'});
end
end

function sigma = observation_sigma_in_output_units(PF, Domain, type, fallback_sigma)
switch char(type)
    case 'discharge'
        area_m2 = sum(Domain.valid(:)) * Domain.cell_area_m2;
        sigma = PF.weighting.sigma.discharge_mm_h / 1000 / 3600 * area_m2;
    case 'groundwater_depth'
        sigma = PF.weighting.sigma.groundwater_depth_m;
    case 'soil_moisture'
        sigma = PF.weighting.sigma.soil_moisture_m3m3;
    otherwise
        sigma = fallback_sigma;
end
end

function [state, window_diag, obs_sim] = forecast_window(PF, Domain, state, theta, t_start, t_end, obs_table)
if nargin < 7
    obs_table = table();
end
state = apply_theta_to_state(PF, Domain, state, theta, false);

n_steps = round((t_end - t_start) / PF.dt_min);
q = zeros(n_steps, 1);
mean_velocity = zeros(n_steps, 1);
infiltration_mm = zeros(n_steps, 1);
recharge_mm = zeros(n_steps, 1);
exfiltration_mm = zeros(n_steps, 1);
obs_sim = nan(height(obs_table), 1);

for s = 1:n_steps
    t_now = t_start + (s - 1) * PF.dt_min;
    rain_rate = rainfall_rate(PF, t_now);
    [state, step_diag] = forecast_step(PF, Domain, state, rain_rate);
    q(s) = step_diag.outlet_q_m3s;
    mean_velocity(s) = step_diag.mean_velocity_m_s;
    infiltration_mm(s) = step_diag.mean_infiltration_mm;
    recharge_mm(s) = step_diag.mean_recharge_mm;
    exfiltration_mm(s) = step_diag.mean_exfiltration_mm;
    state.diagnostics.last_outlet_q_m3s = step_diag.outlet_q_m3s;
    state.diagnostics.last_mean_velocity_m_s = step_diag.mean_velocity_m_s;
    state.diagnostics.last_rainfall_mm_h = rain_rate;
    state.diagnostics.last_infiltration_mm_h = step_diag.mean_infiltration_mm * 60 / PF.dt_min;
    state.diagnostics.last_recharge_mm_h = step_diag.mean_recharge_mm * 60 / PF.dt_min;
    state.diagnostics.last_exfiltration_mm_h = step_diag.mean_exfiltration_mm * 60 / PF.dt_min;

    t_after = t_start + s * PF.dt_min;
    if height(obs_table) > 0
        hit = abs(obs_table.time_min - t_after) < 1e-9;
        if any(hit)
            hit_rows = find(hit);
            for jj = 1:numel(hit_rows)
                row = hit_rows(jj);
                if obs_table.type(row) == "discharge"
                    n_avg = max(1, round(obs_table.sample_interval_min(row) / PF.dt_min));
                    obs_sim(row) = mean(q(max(1, s - n_avg + 1):s), 'omitnan');
                else
                    obs_sim(row) = extract_observation_vector(PF, Domain, state, obs_table(row, :));
                end
            end
        end
    end
end

state.time = t_end;
state.diagnostics.last_outlet_q_m3s = q(end);
state.diagnostics.last_mean_velocity_m_s = mean(mean_velocity, 'omitnan');
state.diagnostics.last_rainfall_mm_h = rainfall_rate(PF, max(t_end - PF.dt_min, 0));
state.diagnostics.last_infiltration_mm_h = mean(infiltration_mm, 'omitnan') * 60 / PF.dt_min;
state.diagnostics.last_recharge_mm_h = mean(recharge_mm, 'omitnan') * 60 / PF.dt_min;
state.diagnostics.last_exfiltration_mm_h = mean(exfiltration_mm, 'omitnan') * 60 / PF.dt_min;

window_diag = struct();
window_diag.time_min = t_end;
window_diag.outlet_q_m3s = q(end);
window_diag.mean_outlet_q_m3s = mean(q, 'omitnan');
window_diag.mean_velocity_m_s = mean(mean_velocity, 'omitnan');
end

function rain_rate = rainfall_rate(PF, time_min)
if time_min < PF.forcing.rainfall_duration_min
    rain_rate = PF.forcing.rainfall_intensity_mm_h;
else
    rain_rate = 0;
end
end

function [state, diag] = forecast_step(PF, Domain, state, rain_rate_mm_h)
dt_h = PF.dt_min / 60;
dt_s = PF.dt_min * 60;
area = Domain.cell_area_m2;
valid = Domain.valid;

rain_step = rain_rate_mm_h * dt_h;
depth_mm = state.depths.d_t * 1000;
depth_mm(valid) = depth_mm(valid) + rain_step;

L = state.Soil_Properties.Layers;
near_deficit = max(L.near_surface_capacity_mm - L.near_surface_storage_mm, 0);
near_deficit_frac = near_deficit ./ max(L.near_surface_capacity_mm, eps);
infil_capacity_mm_h = state.Soil_Properties.ksat_surface_mm_h .* (0.15 + 0.85 .* near_deficit_frac);
infil_mm = min(depth_mm, infil_capacity_mm_h * dt_h);
infil_mm(~valid) = 0;
depth_mm = depth_mm - infil_mm;

[L, recharge_mm] = update_soil_layers(L, infil_mm, state.Soil_Properties.ksat_rootzone_mm_h, dt_h);
state.Soil_Properties.Layers = L;

[state, exfiltration_mm] = update_groundwater(PF, Domain, state, recharge_mm, dt_h);
depth_mm = depth_mm + exfiltration_mm;

state.depths.d_t = max(depth_mm / 1000, 0);
[state, outlet_volume_m3, mean_velocity] = route_surface(PF, Domain, state, dt_s);

state.depths.d_tot = state.depths.d_t;
state.depths.d_p = state.depths.d_t;
state.Hydro_States.f = infil_mm / dt_h;
state.cumulative_fluxes.rain_m3 = state.cumulative_fluxes.rain_m3 + sum(rain_step / 1000 * area .* valid, 'all');
state.cumulative_fluxes.infiltration_m3 = state.cumulative_fluxes.infiltration_m3 + sum(infil_mm / 1000 * area, 'all', 'omitnan');
state.cumulative_fluxes.recharge_m3 = state.cumulative_fluxes.recharge_m3 + sum(recharge_mm / 1000 * area, 'all', 'omitnan');
state.cumulative_fluxes.exfiltration_m3 = state.cumulative_fluxes.exfiltration_m3 + sum(exfiltration_mm / 1000 * area, 'all', 'omitnan');
state.cumulative_fluxes.outlet_m3 = state.cumulative_fluxes.outlet_m3 + outlet_volume_m3;

diag = struct();
diag.outlet_q_m3s = outlet_volume_m3 / dt_s;
diag.mean_velocity_m_s = mean_velocity;
diag.mean_infiltration_mm = mean(infil_mm(valid), 'omitnan');
diag.mean_recharge_mm = mean(recharge_mm(valid), 'omitnan');
diag.mean_exfiltration_mm = mean(exfiltration_mm(valid), 'omitnan');
end

function [L, recharge_mm] = update_soil_layers(L, infil_mm, ksat_rootzone_mm_h, dt_h)
near_space = max(L.near_surface_capacity_mm - L.near_surface_storage_mm, 0);
add_near = min(infil_mm, near_space);
L.near_surface_storage_mm = L.near_surface_storage_mm + add_near;
rem = infil_mm - add_near;

root_space = max(L.root_zone_capacity_mm - L.root_zone_storage_mm, 0);
add_root = min(rem, root_space);
L.root_zone_storage_mm = L.root_zone_storage_mm + add_root;
rem = rem - add_root;

trans_space = max(L.transmission_capacity_mm - L.transmission_storage_mm, 0);
add_trans = min(rem, trans_space);
L.transmission_storage_mm = L.transmission_storage_mm + add_trans;
recharge_mm = rem - add_trans;

near_excess = max(L.near_surface_storage_mm - 0.70 .* L.near_surface_capacity_mm, 0);
near_to_root = min(near_excess, ksat_rootzone_mm_h * dt_h);
root_space = max(L.root_zone_capacity_mm - L.root_zone_storage_mm, 0);
near_to_root = min(near_to_root, root_space);
L.near_surface_storage_mm = L.near_surface_storage_mm - near_to_root;
L.root_zone_storage_mm = L.root_zone_storage_mm + near_to_root;

root_excess = max(L.root_zone_storage_mm - 0.80 .* L.root_zone_capacity_mm, 0);
root_to_trans = min(root_excess, 0.5 .* ksat_rootzone_mm_h * dt_h);
trans_space = max(L.transmission_capacity_mm - L.transmission_storage_mm, 0);
root_to_trans = min(root_to_trans, trans_space);
L.root_zone_storage_mm = L.root_zone_storage_mm - root_to_trans;
L.transmission_storage_mm = L.transmission_storage_mm + root_to_trans;

trans_to_gw = min(L.transmission_storage_mm, 0.25 .* ksat_rootzone_mm_h * dt_h);
L.transmission_storage_mm = L.transmission_storage_mm - trans_to_gw;
recharge_mm = recharge_mm + trans_to_gw;
end

function [state, exfiltration_mm] = update_groundwater(PF, Domain, state, recharge_mm, dt_h)
Sy = max(state.GW_States.Sy, 0.01);
wtd = state.BC_States.gw_depth_m;
wtd = wtd - (recharge_mm / 1000) ./ Sy;

sat_thickness = max(state.Soil_Properties.Soil_Depth - wtd, 0);
available_mm = sat_thickness .* Sy * 1000;
baseflow_factor = sat_thickness ./ max(state.Soil_Properties.Soil_Depth, 0.1);
baseflow_mm = min(available_mm, state.GW_States.ksat_mm_h .* baseflow_factor * dt_h);
wtd = wtd + (baseflow_mm / 1000) ./ Sy;

exfiltration_mm = zeros(size(wtd));
idx_exfil = wtd < 0;
exfiltration_mm(idx_exfil) = -wtd(idx_exfil) .* Sy(idx_exfil) * 1000;
wtd(idx_exfil) = 0;

% Treat groundwater drainage as shallow return flow in this synthetic case.
exfiltration_mm = exfiltration_mm + baseflow_mm;
exfiltration_mm(~Domain.valid) = 0;

state.BC_States.gw_depth_m = wtd;
state.BC_States.h_t = Domain.DEM - wtd;
state.GW_States.baseflow_mm = baseflow_mm;
end

function [state, outlet_volume_m3, mean_velocity] = route_surface(PF, Domain, state, dt_s)
depth = state.depths.d_t;
roughness = state.Soil_Properties.roughness;
area = Domain.cell_area_m2;
outlet_volume_m3 = 0;
velocity_sum = 0;
velocity_count = 0;

for ii = 1:numel(Domain.route_order)
    idx = Domain.route_order(ii);
    h = depth(idx);
    if ~isfinite(h) || h <= PF.routing.surface_depth_epsilon_m
        continue;
    end

    recv = Domain.receiver(idx);
    dist = Domain.receiver_distance_m(idx);
    if ~isfinite(dist) || dist <= 0
        dist = Domain.dx;
    end

    eta = Domain.DEM(idx) + h;
    if recv > 0
        slope = (eta - (Domain.DEM(recv) + depth(recv))) / dist;
        if slope <= 0
            continue;
        end
    else
        slope = PF.routing.min_slope;
    end
    slope = max(slope, PF.routing.min_slope);

    n = max(roughness(idx), 0.01);
    velocity = (1 / n) * h^(2/3) * sqrt(slope);
    fraction = min(PF.routing.max_drain_fraction, max(0, velocity * dt_s / dist));
    volume = h * area * fraction;
    if volume <= 0
        continue;
    end

    depth(idx) = depth(idx) - volume / area;
    if recv > 0
        depth(recv) = depth(recv) + volume / area;
    else
        outlet_volume_m3 = outlet_volume_m3 + volume;
    end
    velocity_sum = velocity_sum + velocity;
    velocity_count = velocity_count + 1;
end

state.depths.d_t = max(depth, 0);
if velocity_count > 0
    mean_velocity = velocity_sum / velocity_count;
else
    mean_velocity = 0;
end
end

function y = extract_observation_vector(PF, Domain, state, obs_table)
y = nan(height(obs_table), 1);
for j = 1:height(obs_table)
    obs_id = obs_table.obs_id(j);
    loc_idx = find([Domain.obs_locations.obs_id] == obs_id, 1);
    obs_type = obs_table.type(j);

    if obs_type == "discharge"
        y(j) = state.diagnostics.last_outlet_q_m3s;
    elseif obs_type == "soil_moisture"
        idx = Domain.obs_locations(loc_idx).idx;
        y(j) = soil_moisture_at_index(state, idx);
    elseif obs_type == "groundwater_depth"
        idx = Domain.obs_locations(loc_idx).idx;
        y(j) = state.BC_States.gw_depth_m(idx);
    else
        error('Unknown observation type: %s', obs_type);
    end
end
end

function sm = soil_moisture_at_index(state, idx)
L = state.Soil_Properties.Layers;
storage = L.near_surface_storage_mm(idx) + L.root_zone_storage_mm(idx);
capacity = L.near_surface_capacity_mm(idx) + L.root_zone_capacity_mm(idx);
saturation = min(max(storage / max(capacity, eps), 0), 1);
sm = state.Soil_Properties.theta_r(idx) + saturation * ...
    (state.Soil_Properties.theta_s(idx) - state.Soil_Properties.theta_r(idx));
end

function [ParamRows, PosteriorRows] = append_parameter_rows(ParamRows, PosteriorRows, PF, Particles, truth_theta, window, time_min, weights)
[names, truth_values] = theta_values(truth_theta, PF.parameters);
values = nan(numel(Particles), numel(names));
for i = 1:numel(Particles)
    [~, values(i, :)] = theta_values(Particles(i).theta, PF.parameters);
    for p = 1:numel(names)
        ParamRows(end+1).window = window; %#ok<AGROW>
        ParamRows(end).time_min = time_min;
        ParamRows(end).particle_id = Particles(i).id;
        ParamRows(end).parameter = names(p);
        ParamRows(end).value = values(i, p);
        ParamRows(end).weight = weights(i);
        ParamRows(end).truth_value = truth_values(p);
    end
end

for p = 1:numel(names)
    PosteriorRows(end+1).window = window; %#ok<AGROW>
    PosteriorRows(end).time_min = time_min;
    PosteriorRows(end).parameter = names(p);
    PosteriorRows(end).truth_value = truth_values(p);
    PosteriorRows(end).weighted_mean = sum(weights(:) .* values(:, p), 'omitnan');
    PosteriorRows(end).p05 = weighted_quantile(values(:, p), weights, 0.05);
    PosteriorRows(end).p50 = weighted_quantile(values(:, p), weights, 0.50);
    PosteriorRows(end).p95 = weighted_quantile(values(:, p), weights, 0.95);
end
end

function ObjectiveRows = append_objective_rows(ObjectiveRows, PF, Particles, window, time_min, weights, J, logL, WeightInfo)
for i = 1:PF.n_particles
    ObjectiveRows(end+1).window = window; %#ok<AGROW>
    ObjectiveRows(end).time_min = time_min;
    ObjectiveRows(end).particle_id = Particles(i).id;
    ObjectiveRows(end).objective = J(i);
    ObjectiveRows(end).log_likelihood = logL(i);
    ObjectiveRows(end).paper_likelihood = WeightInfo.total_likelihood(i);
    ObjectiveRows(end).weight = weights(i);
    ObjectiveRows(end).likelihood_discharge = group_likelihood_for_type(WeightInfo, i, "discharge");
    ObjectiveRows(end).likelihood_groundwater_depth = group_likelihood_for_type(WeightInfo, i, "groundwater_depth");
    ObjectiveRows(end).likelihood_soil_moisture = group_likelihood_for_type(WeightInfo, i, "soil_moisture");
end
end

function value = group_likelihood_for_type(WeightInfo, particle_idx, type)
idx = find(WeightInfo.types == type, 1);
if isempty(idx)
    value = 0;
else
    value = WeightInfo.group_likelihood(particle_idx, idx);
end
end

function ObsRows = append_observation_rows(ObsRows, PF, obs_w, sim_matrix, weights, window, time_min)
for j = 1:height(obs_w)
    weighted_sim = sum(weights(:) .* sim_matrix(:, j), 'omitnan');
    for i = 1:PF.n_particles
        ObsRows(end+1).window = window; %#ok<AGROW>
        ObsRows(end).time_min = time_min;
        ObsRows(end).particle_id = i;
        ObsRows(end).obs_id = obs_w.obs_id(j);
        ObsRows(end).type = obs_w.type(j);
        ObsRows(end).observed = obs_w.value(j);
        ObsRows(end).sample_interval_min = obs_w.sample_interval_min(j);
        ObsRows(end).sigma = obs_w.sigma(j);
        ObsRows(end).simulated = sim_matrix(i, j);
        ObsRows(end).weighted_simulated = weighted_sim;
        ObsRows(end).weight = weights(i);
        ObsRows(end).residual = obs_w.value(j) - sim_matrix(i, j);
        ObsRows(end).units = obs_w.units(j);
    end
end
end

function EssRows = append_ess_row(EssRows, window, time_min, Neff, did_resample, weighted_objective)
EssRows(end+1).window = window; %#ok<AGROW>
EssRows(end).time_min = time_min;
EssRows(end).effective_sample_size = Neff;
EssRows(end).resampled = logical(did_resample);
EssRows(end).weighted_objective = weighted_objective;
end

function StateRows = append_state_rows(StateRows, PF, Particles, window, time_min, weights)
surface = zeros(numel(Particles), 1);
soil = zeros(numel(Particles), 1);
gw = zeros(numel(Particles), 1);
for i = 1:numel(Particles)
    S = Particles(i).state;
    surface(i) = sum(S.depths.d_t(:), 'omitnan');
    L = S.Soil_Properties.Layers;
    soil(i) = mean(L.near_surface_storage_mm(:) + L.root_zone_storage_mm(:), 'omitnan');
    gw(i) = mean(S.BC_States.gw_depth_m(:), 'omitnan');
end
StateRows(end+1).window = window; %#ok<AGROW>
StateRows(end).time_min = time_min;
StateRows(end).weighted_surface_storage_index = sum(weights(:) .* surface, 'omitnan');
StateRows(end).weighted_mean_soil_storage_mm = sum(weights(:) .* soil, 'omitnan');
StateRows(end).weighted_mean_groundwater_depth_m = sum(weights(:) .* gw, 'omitnan');
end

function [names, values] = theta_values(theta, specs)
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

function sigma_by_name = inflated_parameter_sigma(PF, Particles, weights)
weights = weights(:);
weights = weights ./ max(sum(weights), eps);
N = numel(Particles);
[names, first_values] = theta_values(Particles(1).theta, PF.parameters);
values = nan(N, numel(first_values));
values(1, :) = first_values;
for i = 2:N
    [~, values(i, :)] = theta_values(Particles(i).theta, PF.parameters);
end

denom = ((N - 1) / N) * sum(weights);
lambda = PF.resampling.inflation_lambda;
sigma_by_name = struct();
for p = 1:numel(names)
    mu = sum(weights .* values(:, p), 'omitnan');
    sigma_hat = sqrt(sum(weights .* (values(:, p) - mu) .^ 2, 'omitnan') ./ max(denom, eps));
    sigma_by_name.(matlab.lang.makeValidName(char(names(p)))) = lambda * sigma_hat;
end
end

function q = weighted_quantile(values, weights, prob)
values = double(values(:));
weights = double(weights(:));
idx = isfinite(values) & isfinite(weights) & weights >= 0;
values = values(idx);
weights = weights(idx);
if isempty(values) || sum(weights) <= 0
    q = nan;
    return;
end
[values, order] = sort(values);
weights = weights(order) ./ sum(weights);
cdf = cumsum(weights);
q = values(find(cdf >= prob, 1, 'first'));
end

function Summary = observation_ensemble_summary(O)
if isempty(O)
    Summary = table();
    return;
end

[G, window, time_min, obs_id, type] = findgroups(O.window, O.time_min, O.obs_id, O.type);
n = max(G);
Rows = repmat(struct( ...
    'window', [], 'time_min', [], 'obs_id', "", 'type', "", ...
    'observed', [], 'weighted_mean', [], 'p05', [], 'p50', [], 'p95', [], ...
    'sigma', [], 'units', ""), n, 1);

for g = 1:n
    idx = G == g;
    w = O.weight(idx);
    sim = O.simulated(idx);
    Rows(g).window = window(g);
    Rows(g).time_min = time_min(g);
    Rows(g).obs_id = obs_id(g);
    Rows(g).type = type(g);
    Rows(g).observed = O.observed(find(idx, 1, 'first'));
    Rows(g).weighted_mean = sum(w .* sim, 'omitnan') ./ max(sum(w, 'omitnan'), eps);
    Rows(g).p05 = weighted_quantile(sim, w, 0.05);
    Rows(g).p50 = weighted_quantile(sim, w, 0.50);
    Rows(g).p95 = weighted_quantile(sim, w, 0.95);
    Rows(g).sigma = O.sigma(find(idx, 1, 'first'));
    Rows(g).units = O.units(find(idx, 1, 'first'));
end

Summary = struct2table(Rows);
end

function write_results(PF, Results)
if ~exist(PF.outputs_dir, 'dir')
    mkdir(PF.outputs_dir);
end
save(fullfile(PF.outputs_dir, 'PF_Results.mat'), 'Results', '-v7.3');
writetable(Results.Parameter_Evolution, fullfile(PF.outputs_dir, 'Parameter_Evolution.csv'));
writetable(Results.Objective_Function, fullfile(PF.outputs_dir, 'Objective_Function.csv'));
writetable(Results.Effective_Sample_Size, fullfile(PF.outputs_dir, 'Effective_Sample_Size.csv'));
writetable(Results.Observation_Fit, fullfile(PF.outputs_dir, 'Observation_Fit.csv'));
writetable(Results.Observation_Ensemble_Summary, fullfile(PF.outputs_dir, 'Observation_Ensemble_Summary.csv'));
writetable(Results.Posterior_Summary, fullfile(PF.outputs_dir, 'Posterior_Summary.csv'));
writetable(Results.Posterior_Summary_Physical, fullfile(PF.outputs_dir, 'Posterior_Summary_Physical.csv'));
writetable(Results.State_Summary, fullfile(PF.outputs_dir, 'State_Summary.csv'));
writetable(Results.Truth, fullfile(PF.outputs_dir, 'Synthetic_Truth.csv'));
writetable(Results.Observations, fullfile(PF.outputs_dir, 'Synthetic_Observations.csv'));
end

function make_pf_plots(PF, Results)
if ~exist(PF.figures_dir, 'dir')
    mkdir(PF.figures_dir);
end

P = Results.Posterior_Summary_Physical;
params = unique(P.parameter, 'stable');
n = numel(params);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1800 1200]);
tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:min(n, 16)
    nexttile;
    idx = P.parameter == params(i);
    T = sortrows(P(idx, :), 'time_min');
    fill([T.time_min; flipud(T.time_min)], [T.p05; flipud(T.p95)], ...
        [0.82 0.89 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.75); hold on;
    plot(T.time_min, T.weighted_mean, 'b-', 'LineWidth', 1.5);
    plot(T.time_min, T.truth_value, 'k--', 'LineWidth', 1);
    title(T.display_name(1), 'Interpreter', 'none', 'FontSize', 8);
    ylabel(T.units(1), 'Interpreter', 'none', 'FontSize', 7);
    xlabel('Time [min]', 'FontSize', 7);
    grid on;
end
save_fig(fig, fullfile(PF.figures_dir, 'parameter_evolution.png'));
make_initial_parameter_spread(PF, Results);

E = Results.Effective_Sample_Size;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
yyaxis left;
plot(E.time_min, E.weighted_objective, 'k-o', 'LineWidth', 1.5);
ylabel('Weighted objective J');
yyaxis right;
plot(E.time_min, E.effective_sample_size, 'b-s', 'LineWidth', 1.5);
ylabel('Effective sample size');
xlabel('Time [min]');
grid on;
title('Objective Function and Ensemble Degeneracy');
save_fig(fig, fullfile(PF.figures_dir, 'objective_function.png'));

O = Results.Observation_Ensemble_Summary;
obs_ids = unique(O.obs_id, 'stable');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 900]);
tiledlayout(numel(obs_ids), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(obs_ids)
    nexttile;
    idx = O.obs_id == obs_ids(i);
    T = sortrows(O(idx, :), 'time_min');
    fill([T.time_min; flipud(T.time_min)], [T.p05; flipud(T.p95)], ...
        [0.85 0.90 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.75); hold on;
    plot(T.time_min, T.weighted_mean, 'b-', 'LineWidth', 1.5);
    plot(T.time_min, T.observed, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
    ylabel(sprintf('%s %s', char(obs_ids(i)), char(T.units(1))), 'Interpreter', 'none');
    grid on;
end
xlabel('Time [min]');
save_fig(fig, { ...
    fullfile(PF.figures_dir, 'observation_fit.png'), ...
    fullfile(PF.figures_dir, 'observation_ensemble_fit.png')});

W = Results.Objective_Function;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
scatter(W.time_min, W.weight, 18, W.objective, 'filled');
xlabel('Time [min]');
ylabel('Particle weight');
cb = colorbar;
cb.Label.String = 'Objective J';
grid on;
title('Particle Weights');
save_fig(fig, fullfile(PF.figures_dir, 'particle_weights.png'));

S = Results.State_Summary;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 700]);
plot(S.time_min, S.weighted_mean_soil_storage_mm, 'g-', 'LineWidth', 1.5); hold on;
plot(S.time_min, S.weighted_mean_groundwater_depth_m, 'b-', 'LineWidth', 1.5);
legend({'Mean soil storage [mm]', 'Mean groundwater depth [m]'}, 'Location', 'best');
xlabel('Time [min]');
grid on;
title('Weighted State Summary');
save_fig(fig, fullfile(PF.figures_dir, 'state_summary.png'));
end

function make_initial_parameter_spread(PF, Results)
I = Results.Parameter_Evolution(Results.Parameter_Evolution.window == 0, :);
if isempty(I)
    return;
end
params = unique(I.parameter, 'stable');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1800 1200]);
tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:min(numel(params), 16)
    nexttile;
    idx = I.parameter == params(i);
    spec = PF.parameters([PF.parameters.name] == params(i));
    [factor, display_name, units] = physical_parameter_metadata(PF, spec(1));
    histogram(I.value(idx) * factor, 14, 'FaceColor', [0.55 0.72 0.95], ...
        'EdgeColor', 'none'); hold on;
    xline(I.truth_value(find(idx, 1, 'first')) * factor, 'k--', 'LineWidth', 1);
    title(display_name, 'Interpreter', 'none', 'FontSize', 8);
    xlabel(units, 'Interpreter', 'none', 'FontSize', 7);
    ylabel('Particles', 'FontSize', 7);
    grid on;
end
save_fig(fig, fullfile(PF.figures_dir, 'initial_parameter_spread.png'));
end

function Tphys = physical_posterior_table(PF, T)
if isempty(T)
    Tphys = T;
    return;
end

Tphys = T;
Tphys.display_name = strings(height(Tphys), 1);
Tphys.units = strings(height(Tphys), 1);

for i = 1:height(Tphys)
    spec = PF.parameters([PF.parameters.name] == Tphys.parameter(i));
    if isempty(spec)
        factor = 1;
        display_name = char(Tphys.parameter(i));
        units = "[-]";
    else
        [factor, display_name, units] = physical_parameter_metadata(PF, spec(1));
    end

    Tphys.truth_value(i) = Tphys.truth_value(i) * factor;
    Tphys.weighted_mean(i) = Tphys.weighted_mean(i) * factor;
    Tphys.p05(i) = Tphys.p05(i) * factor;
    Tphys.p50(i) = Tphys.p50(i) * factor;
    Tphys.p95(i) = Tphys.p95(i) * factor;
    Tphys.display_name(i) = display_name;
    Tphys.units(i) = units;
end
end

function [factor, display_name, units] = physical_parameter_metadata(PF, spec)
class_id = spec.class_id;
property = char(spec.property);

switch char(lower(spec.target_type))
    case 'lulc'
        if strcmp(property, 'manning_mult')
            factor = PF.base.lulc(class_id).manning_n;
            display_name = sprintf('LULC %d Manning n [-]', class_id);
            units = "[-]";
        elseif strcmp(property, 'root_depth_mult')
            factor = PF.base.lulc(class_id).root_depth_m;
            display_name = sprintf('LULC %d root depth [m]', class_id);
            units = "[m]";
        else
            factor = 1;
            display_name = sprintf('LULC %d %s [-]', class_id, property);
            units = "[-]";
        end
    case 'soil'
        if strcmp(property, 'ksat_surface_mult')
            factor = PF.base.soil(class_id).ksat_surface_mm_h;
            display_name = sprintf('SOIL %d near-surface Ksat [mm/h]', class_id);
            units = "[mm/h]";
        elseif strcmp(property, 'ksat_rootzone_mult')
            factor = PF.base.soil(class_id).ksat_rootzone_mm_h;
            display_name = sprintf('SOIL %d root-zone Ksat [mm/h]', class_id);
            units = "[mm/h]";
        elseif strcmp(property, 'storage_mult')
            factor = PF.base.soil(class_id).storage_capacity_mm;
            display_name = sprintf('SOIL %d storage capacity [mm]', class_id);
            units = "[mm]";
        elseif strcmp(property, 'initial_sm_mult')
            factor = PF.base.soil(class_id).initial_sm_fraction;
            display_name = sprintf('SOIL %d initial saturation [-]', class_id);
            units = "[-]";
        else
            factor = 1;
            display_name = sprintf('SOIL %d %s [-]', class_id, property);
            units = "[-]";
        end
    case 'gw'
        if strcmp(property, 'ksat_mult')
            factor = PF.base.gw(class_id).ksat_mm_h;
            display_name = sprintf('GW %d Ksat [mm/h]', class_id);
            units = "[mm/h]";
        elseif strcmp(property, 'specific_yield_mult')
            factor = PF.base.gw(class_id).specific_yield;
            display_name = sprintf('GW %d specific yield [-]', class_id);
            units = "[-]";
        elseif strcmp(property, 'initial_wtd_mult')
            factor = PF.base.gw(class_id).initial_wtd_m;
            display_name = sprintf('GW %d initial water-table depth [m]', class_id);
            units = "[m]";
        else
            factor = 1;
            display_name = sprintf('GW %d %s [-]', class_id, property);
            units = "[-]";
        end
    case 'ic'
        if strcmp(property, 'initial_wtd_mult')
            factor = PF.base.gw(class_id).initial_wtd_m;
            display_name = sprintf('IC %d initial water-table depth [m]', class_id);
            units = "[m]";
        else
            factor = 1;
            display_name = sprintf('IC %d %s [-]', class_id, property);
            units = "[-]";
        end
    otherwise
        factor = 1;
        display_name = char(spec.name);
        units = "[-]";
end
end

function save_fig(fig, paths)
if ischar(paths) || isstring(paths)
    paths = {char(paths)};
end

for i = 1:numel(paths)
    path = paths{i};
    try
        exportgraphics(fig, path, 'Resolution', 200);
    catch
        saveas(fig, path);
    end
end
close(fig);
end

function T = struct2table_safe(S)
if isempty(S)
    T = table();
else
    T = struct2table(S);
end
end

function D = compact_domain_for_save(Domain)
D = rmfield(Domain, {'R'});
end
