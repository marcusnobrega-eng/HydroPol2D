function Results = test_generic_ga_known_optimum()
%TEST_GENERIC_GA_KNOWN_OPTIMUM Verify GA selection, mutation, and objective logic.
%
% This deterministic integration test uses two calibration multipliers and
% a synthetic forward model with a known optimum. The target is not an
% initial ensemble member. The runner returns metrics only, so the generic
% normalized hydrograph objective is also exercised.

this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);

output_root = tempname;
cleanup = onCleanup(@() remove_if_present(output_root)); %#ok<NASGU>

roughness_spec = hydropol2d_make_parameter_spec( ...
    'roughness_multiplier', 'LULC', 'roughness', 1, ...
    'multiplier', 0.4, 1.6, 0.55, 'linear', true);
moisture_spec = hydropol2d_make_parameter_spec( ...
    'initial_moisture_multiplier', 'SOIL', 'theta_i', 1, ...
    'multiplier', 0.4, 1.6, 0.55, 'linear', true);

Cal = struct();
Cal.name = 'generic_ga_known_optimum';
Cal.random_seed = 20260720;
Cal.output_root = output_root;
Cal.ga = struct( ...
    'n_generations', 12, ...
    'n_individuals', 20, ...
    'elite_fraction', 0.15, ...
    'tournament_size', 3, ...
    'crossover_fraction', 0.85, ...
    'mutation_probability', 0.35, ...
    'mutation_scale_fraction', 0.12);
Cal.base_overrides = struct( ...
    'LULC', struct('Index', 1, 'roughness', 0.035), ...
    'SOIL', struct('Index', 1, 'theta_i', 0.25));
Cal.input_data_overrides = struct();
Cal.parameter_specs = hydropol2d_finalize_parameter_specs( ...
    [roughness_spec moisture_spec]);
Cal.cases = struct('name', 'synthetic_known_optimum');
Cal.model_runner = @known_optimum_runner;

Results = run_hydropol2d_ga_calibration(Cal);

first_generation_best = min(Results.History.Objective( ...
    Results.History.Generation == 1));
target = [0.83 1.27];

assert(Results.Best.objective < first_generation_best, ...
    'The GA did not improve on the initial population.');
assert(Results.Best.objective < 0.01, ...
    'The GA did not reach the known objective tolerance.');
assert(max(abs(Results.Best.candidate.values - target)) < 0.06, ...
    'The GA did not recover the known parameter values.');

fprintf(['Generic GA known-optimum test passed: initial best %.6g, ' ...
    'final best %.6g, parameters [%.6f %.6f].\n'], ...
    first_generation_best, Results.Best.objective, ...
    Results.Best.candidate.values(1), Results.Best.candidate.values(2));
end

function Eval = known_optimum_runner(~, ~, candidate, ~)
% Produce hydrograph-like metrics whose unique optimum is known exactly.
target = [0.83 1.27];
error_vector = candidate.values - target;
total_error = sum(error_vector .^ 2);

Metrics = struct();
Metrics.RMSE_m3_s = 100 * total_error;
Metrics.MAE_m3_s = 0;
Metrics.Bias_m3_s = 100 * total_error;
Metrics.NSE = 1;
Metrics.Observed_Volume_m3 = 10000;
Metrics.Modeled_Volume_m3 = 10000 * (1 + total_error);
Metrics.Volume_Error_pct = 100 * total_error;
Metrics.Observed_Peak_m3_s = 100;
Metrics.Modeled_Peak_m3_s = 100 * (1 + error_vector(1) ^ 2);
Metrics.Observed_Peak_Time_min = 60;
Metrics.Modeled_Peak_Time_min = 60 + 120 * error_vector(2) ^ 2;
Metrics.Peak_Timing_Error_min = 120 * error_vector(2) ^ 2;

% Do not set Eval.objective: the engine must combine these metrics itself.
Eval = struct('metrics', Metrics);
end

function remove_if_present(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
