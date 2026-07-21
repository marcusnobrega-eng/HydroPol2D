function Results = run_hydropol2d_ga_calibration(Cal)
%RUN_HYDROPOL2D_GA_CALIBRATION Generic single-objective GA calibration.
%
% The model-specific pieces are supplied by Cal:
%   Cal.parameter_specs  - struct array defining bounds/targets
%   Cal.base_overrides   - baseline ParameterOverrides struct
%   Cal.cases            - rainfall/event case definitions
%   Cal.model_runner     - function handle: Eval = f(Cal, caseDef, candidate, info)
%
% Each candidate is converted to HydroPol2D-compatible overrides, then
% evaluated by Cal.model_runner. The GA minimizes a single scalar objective.

arguments
    Cal struct
end

Cal = fill_default_calibration_options(Cal);
if ~exist(Cal.output_root, 'dir')
    mkdir(Cal.output_root);
end
if ~exist(fullfile(Cal.output_root, 'Runs'), 'dir')
    mkdir(fullfile(Cal.output_root, 'Runs'));
end

rng(Cal.random_seed, 'twister');
specs = Cal.parameter_specs([Cal.parameter_specs.enabled]);
assert(~isempty(specs), 'At least one enabled calibration parameter is required.');

n_params = numel(specs);
n_ind = Cal.ga.n_individuals;
n_gen = Cal.ga.n_generations;

Population = initialize_population(specs, n_ind);
History = table();
Best = struct('objective', inf, 'x', [], 'candidate', [], 'eval', []);

fprintf('\n============================================================\n');
fprintf('HydroPol2D generic GA calibration\n');
fprintf('Output root  : %s\n', Cal.output_root);
fprintf('Parameters   : %d\n', n_params);
fprintf('Individuals  : %d\n', n_ind);
fprintf('Generations  : %d\n', n_gen);
fprintf('Cases        : %d\n', numel(Cal.cases));
fprintf('============================================================\n\n');

for gen = 1:n_gen
    fprintf('GA generation %d/%d\n', gen, n_gen);
    Evaluated = repmat(empty_eval_record(n_params), n_ind, 1);

    for ind = 1:n_ind
        info = struct();
        info.generation = gen;
        info.individual = ind;
        info.run_id = sprintf('gen_%03d_ind_%03d', gen, ind);
        info.run_dir = fullfile(Cal.output_root, 'Runs', info.run_id);
        if ~exist(info.run_dir, 'dir')
            mkdir(info.run_dir);
        end

        candidate = build_candidate_from_vector(Population(ind, :), specs, Cal);
        tic_eval = tic;
        try
            Eval = evaluate_candidate(Cal, candidate, info);
        catch ME
            Eval = failed_evaluation(ME);
        end
        Eval.runtime_seconds = toc(tic_eval);
        Evaluated(ind) = pack_eval_record(gen, ind, Population(ind, :), candidate, Eval);

        fprintf('  ind %02d | J = %.6g | runtime %.1fs | %s\n', ...
            ind, Evaluated(ind).objective, Eval.runtime_seconds, Eval.status);
        if isfield(Eval, 'error_message') && strlength(Eval.error_message) > 0
            fprintf('          error: %s\n', Eval.error_message);
        end

        if Evaluated(ind).objective < Best.objective
            Best.objective = Evaluated(ind).objective;
            Best.x = Population(ind, :);
            Best.candidate = candidate;
            Best.eval = Eval;
            save(fullfile(Cal.output_root, 'BestSoFar.mat'), 'Best', 'Cal');
        end
    end

    GenTable = eval_records_to_table(Evaluated, specs);
    writetable(GenTable, fullfile(Cal.output_root, sprintf('Generation_%03d.csv', gen)));
    History = [History; GenTable]; %#ok<AGROW>
    writetable(History, fullfile(Cal.output_root, 'GA_History.csv'));
    save(fullfile(Cal.output_root, 'GA_Checkpoint.mat'), ...
        'Cal', 'Population', 'History', 'Best', 'gen');
    PartialResults = struct('Cal', Cal, 'History', History, ...
        'Best', Best, 'specs', specs);
    plot_ga_diagnostics(Cal, PartialResults);

    if gen < n_gen
        Population = make_next_generation(Population, Evaluated, specs, Cal.ga);
    end
end

Results = struct();
Results.Cal = Cal;
Results.History = History;
Results.Best = Best;
Results.specs = specs;
save(fullfile(Cal.output_root, 'GA_Results.mat'), 'Results');
write_best_summary(Cal, Results);
plot_ga_diagnostics(Cal, Results);
end

function Cal = fill_default_calibration_options(Cal)
    if ~isfield(Cal, 'random_seed'), Cal.random_seed = 42; end
    if ~isfield(Cal, 'output_root') || isempty(Cal.output_root)
        Cal.output_root = fullfile(pwd, 'Calibration_Output');
    end
    if ~isfield(Cal, 'base_overrides'), Cal.base_overrides = struct(); end
    if ~isfield(Cal, 'input_data_overrides'), Cal.input_data_overrides = struct(); end
    assert(isfield(Cal, 'parameter_specs'), 'Cal.parameter_specs is required.');
    assert(isfield(Cal, 'cases') && ~isempty(Cal.cases), 'Cal.cases is required.');
    assert(isfield(Cal, 'model_runner') && isa(Cal.model_runner, 'function_handle'), ...
        'Cal.model_runner function handle is required.');

    if ~isfield(Cal, 'ga'), Cal.ga = struct(); end
    if ~isfield(Cal.ga, 'n_generations'), Cal.ga.n_generations = 5; end
    if ~isfield(Cal.ga, 'n_individuals'), Cal.ga.n_individuals = 10; end
    if ~isfield(Cal.ga, 'elite_fraction'), Cal.ga.elite_fraction = 0.20; end
    if ~isfield(Cal.ga, 'tournament_size'), Cal.ga.tournament_size = 3; end
    if ~isfield(Cal.ga, 'crossover_fraction'), Cal.ga.crossover_fraction = 0.80; end
    if ~isfield(Cal.ga, 'mutation_probability'), Cal.ga.mutation_probability = 0.20; end
    if ~isfield(Cal.ga, 'mutation_scale_fraction'), Cal.ga.mutation_scale_fraction = 0.12; end

    if ~isfield(Cal, 'objective'), Cal.objective = struct(); end
    if ~isfield(Cal.objective, 'weights')
        Cal.objective.weights = struct('rmse', 1.0, 'volume', 1.0, ...
            'peak', 0.75, 'timing', 0.25, 'bias', 0.25);
    end
    if ~isfield(Cal.objective, 'timing_scale_min')
        Cal.objective.timing_scale_min = 120;
    end
end

function P = initialize_population(specs, n_ind)
    n_params = numel(specs);
    P = zeros(n_ind, n_params);
    for j = 1:n_params
        lo = specs(j).lower_internal;
        hi = specs(j).upper_internal;
        P(:, j) = lo + (hi - lo) .* rand(n_ind, 1);
        if isfield(specs(j), 'initial_internal') && isfinite(specs(j).initial_internal)
            P(1, j) = min(max(specs(j).initial_internal, lo), hi);
        end
    end
end

function candidate = build_candidate_from_vector(x, specs, Cal)
    candidate = struct();
    candidate.x_internal = x(:)';
    candidate.values = zeros(1, numel(specs));
    candidate.names = strings(1, numel(specs));
    candidate.ParameterOverrides = Cal.base_overrides;
    candidate.InputDataOverrides = Cal.input_data_overrides;

    for j = 1:numel(specs)
        value = internal_to_physical(x(j), specs(j));
        candidate.values(j) = value;
        candidate.names(j) = string(specs(j).name);
        candidate = apply_spec_value(candidate, specs(j), value, Cal);
    end
end

function candidate = apply_spec_value(candidate, spec, value, Cal)
    target = lower(char(spec.target));
    switch target
        case {'lulc', 'soil'}
            group = upper(target);
            if ~isfield(candidate.ParameterOverrides, group) || ...
                    ~isstruct(candidate.ParameterOverrides.(group))
                candidate.ParameterOverrides.(group) = struct();
            end
            if ~isfield(candidate.ParameterOverrides.(group), 'Index')
                candidate.ParameterOverrides.(group).Index = spec.index(:);
            end

            base_values = get_base_table_values(Cal.base_overrides, group, spec.field, spec.index);
            if strcmpi(spec.operation, 'multiplier')
                new_values = base_values .* value;
            elseif strcmpi(spec.operation, 'offset')
                new_values = base_values + value;
            else
                new_values = repmat(value, numel(spec.index), 1);
            end

            candidate.ParameterOverrides.(group).Index = union_preserve( ...
                candidate.ParameterOverrides.(group).Index(:), spec.index(:));
            candidate.ParameterOverrides.(group) = assign_override_values( ...
                candidate.ParameterOverrides.(group), spec.field, spec.index, new_values);

        case {'inputdata', 'input_data'}
            candidate.InputDataOverrides = set_nested_field( ...
                candidate.InputDataOverrides, spec.path, value);

        otherwise
            error('Unsupported calibration target "%s" for parameter "%s".', ...
                spec.target, spec.name);
    end
end

function values = get_base_table_values(baseOverrides, group, field, idx)
    if ~isfield(baseOverrides, group) || ~isfield(baseOverrides.(group), 'Index') || ...
            ~isfield(baseOverrides.(group), field)
        values = ones(numel(idx), 1);
        return;
    end
    baseIdx = baseOverrides.(group).Index(:);
    baseVal = baseOverrides.(group).(field);
    values = nan(numel(idx), 1);
    for i = 1:numel(idx)
        row = find(baseIdx == idx(i), 1);
        if isempty(row)
            values(i) = 1;
        elseif numel(baseVal) == 1
            values(i) = baseVal;
        else
            values(i) = baseVal(row);
        end
    end
end

function S = assign_override_values(S, field, idx, values)
    if ~isfield(S, 'Index')
        S.Index = idx(:);
    end
    if ~isfield(S, field) || numel(S.(field)) ~= numel(S.Index)
        S.(field) = nan(numel(S.Index), 1);
    end
    for i = 1:numel(idx)
        row = find(S.Index == idx(i), 1);
        if isempty(row)
            S.Index(end + 1, 1) = idx(i);
            S.(field)(end + 1, 1) = values(i);
        else
            S.(field)(row, 1) = values(i);
        end
    end
end

function out = union_preserve(a, b)
    out = a(:);
    for i = 1:numel(b)
        if ~ismember(b(i), out)
            out(end + 1, 1) = b(i); %#ok<AGROW>
        end
    end
end

function S = set_nested_field(S, path, value)
    if ischar(path) || isstring(path)
        parts = split(string(path), ".");
    else
        parts = string(path);
    end
    if numel(parts) == 1
        S.(parts(1)) = value;
        return;
    end
    name = char(parts(1));
    if ~isfield(S, name) || ~isstruct(S.(name))
        S.(name) = struct();
    end
    S.(name) = set_nested_field(S.(name), parts(2:end), value);
end

function Eval = evaluate_candidate(Cal, candidate, info)
    caseRecords = cell(numel(Cal.cases), 1);
    objectives = zeros(numel(Cal.cases), 1);
    for c = 1:numel(Cal.cases)
        caseEval = Cal.model_runner(Cal, Cal.cases(c), candidate, info);
        if ~isfield(caseEval, 'objective') || ~isfinite(caseEval.objective)
            caseEval.objective = objective_from_metrics(caseEval.metrics, Cal.objective);
        end
        caseRecords{c} = caseEval;
        objectives(c) = caseEval.objective;
    end
    Eval = struct();
    Eval.status = "ok";
    Eval.case = caseRecords;
    Eval.objective = mean(objectives, 'omitnan');
    Eval.metrics = aggregate_case_metrics(caseRecords);
end

function J = objective_from_metrics(M, Obj)
    w = Obj.weights;
    obsPeak = max(abs(M.Observed_Peak_m3_s), eps);
    obsVol = max(abs(M.Observed_Volume_m3), eps);
    rmseTerm = M.RMSE_m3_s / obsPeak;
    volumeTerm = abs(M.Modeled_Volume_m3 - M.Observed_Volume_m3) / obsVol;
    peakTerm = abs(M.Modeled_Peak_m3_s - M.Observed_Peak_m3_s) / obsPeak;
    timingTerm = abs(M.Peak_Timing_Error_min) / max(Obj.timing_scale_min, eps);
    biasTerm = abs(M.Bias_m3_s) / obsPeak;
    J = w.rmse * rmseTerm + w.volume * volumeTerm + ...
        w.peak * peakTerm + w.timing * timingTerm + w.bias * biasTerm;
    if ~isfinite(J)
        J = inf;
    end
end

function metrics = aggregate_case_metrics(caseRecords)
    metrics = caseRecords{1}.metrics;
    if numel(caseRecords) == 1
        return;
    end
    names = fieldnames(metrics);
    for i = 1:numel(names)
        vals = cellfun(@(c) c.metrics.(names{i}), caseRecords);
        metrics.(names{i}) = mean(vals, 'omitnan');
    end
end

function Eval = failed_evaluation(ME)
    Eval = struct();
    Eval.status = "failed";
    Eval.objective = inf;
    Eval.metrics = empty_metrics();
    Eval.error_identifier = string(ME.identifier);
    Eval.error_message = string(ME.message);
end

function M = empty_metrics()
    M = struct();
    M.RMSE_m3_s = nan;
    M.MAE_m3_s = nan;
    M.Bias_m3_s = nan;
    M.NSE = nan;
    M.Observed_Volume_m3 = nan;
    M.Modeled_Volume_m3 = nan;
    M.Volume_Error_pct = nan;
    M.Observed_Peak_m3_s = nan;
    M.Modeled_Peak_m3_s = nan;
    M.Observed_Peak_Time_min = nan;
    M.Modeled_Peak_Time_min = nan;
    M.Peak_Timing_Error_min = nan;
end

function rec = empty_eval_record(n_params)
    rec = struct();
    rec.generation = nan;
    rec.individual = nan;
    rec.x = nan(1, n_params);
    rec.values = nan(1, n_params);
    rec.objective = inf;
    rec.status = "empty";
    rec.runtime_seconds = nan;
    rec.candidate = [];
    rec.eval = [];
end

function rec = pack_eval_record(gen, ind, x, candidate, Eval)
    rec = struct();
    rec.generation = gen;
    rec.individual = ind;
    rec.x = x;
    rec.values = candidate.values;
    rec.objective = Eval.objective;
    rec.status = Eval.status;
    rec.runtime_seconds = Eval.runtime_seconds;
    rec.candidate = candidate;
    rec.eval = Eval;
end

function T = eval_records_to_table(Evaluated, specs)
    n = numel(Evaluated);
    T = table();
    T.Generation = arrayfun(@(r) r.generation, Evaluated(:));
    T.Individual = arrayfun(@(r) r.individual, Evaluated(:));
    T.Objective = arrayfun(@(r) r.objective, Evaluated(:));
    T.Status = arrayfun(@(r) string(r.status), Evaluated(:));
    T.Runtime_s = arrayfun(@(r) r.runtime_seconds, Evaluated(:));
    for j = 1:numel(specs)
        cleanName = matlab.lang.makeValidName(char(specs(j).name));
        T.(cleanName) = arrayfun(@(r) r.values(j), Evaluated(:));
    end
    metricNames = {};
    for i = 1:n
        if isfield(Evaluated(i).eval, 'metrics') && isstruct(Evaluated(i).eval.metrics)
            metricNames = union(metricNames, fieldnames(Evaluated(i).eval.metrics));
        end
    end
    for j = 1:numel(metricNames)
        name = metricNames{j};
        col = nan(n, 1);
        for i = 1:n
            if isfield(Evaluated(i).eval, 'metrics') && isfield(Evaluated(i).eval.metrics, name)
                col(i) = Evaluated(i).eval.metrics.(name);
            end
        end
        T.(matlab.lang.makeValidName(name)) = col;
    end
end

function Next = make_next_generation(Population, Evaluated, specs, ga)
    n_ind = size(Population, 1);
    n_params = size(Population, 2);
    objectives = [Evaluated.objective]';
    [~, order] = sort(objectives, 'ascend', 'MissingPlacement', 'last');
    n_elite = max(1, round(ga.elite_fraction * n_ind));
    Next = zeros(size(Population));
    Next(1:n_elite, :) = Population(order(1:n_elite), :);
    for i = (n_elite + 1):n_ind
        p1 = tournament_select(Population, objectives, ga.tournament_size);
        p2 = tournament_select(Population, objectives, ga.tournament_size);
        if rand < ga.crossover_fraction
            alpha = rand(1, n_params);
            child = alpha .* p1 + (1 - alpha) .* p2;
        else
            child = p1;
        end
        child = mutate_child(child, specs, ga);
        Next(i, :) = child;
    end
end

function parent = tournament_select(Population, objectives, k)
    n = size(Population, 1);
    idx = randi(n, [max(1, k), 1]);
    [~, bestLocal] = min(objectives(idx));
    parent = Population(idx(bestLocal), :);
end

function child = mutate_child(child, specs, ga)
    for j = 1:numel(specs)
        if rand < ga.mutation_probability
            span = specs(j).upper_internal - specs(j).lower_internal;
            child(j) = child(j) + randn * ga.mutation_scale_fraction * span;
        end
        child(j) = min(max(child(j), specs(j).lower_internal), specs(j).upper_internal);
    end
end

function y = internal_to_physical(x, spec)
    if strcmpi(spec.transform, 'log')
        y = exp(x);
    else
        y = x;
    end
end

function write_best_summary(Cal, Results)
    Best = Results.Best;
    specs = Results.specs;
    if isempty(Best.candidate) || ~isfield(Best.candidate, 'values')
        warning('No successful GA candidate was available for Best_Parameters.csv.');
        return;
    end
    T = table();
    T.Parameter = string({specs.name})';
    T.Value = Best.candidate.values(:);
    T.Lower = arrayfun(@(s) s.lower, specs(:));
    T.Upper = arrayfun(@(s) s.upper, specs(:));
    writetable(T, fullfile(Cal.output_root, 'Best_Parameters.csv'));
end

function plot_ga_diagnostics(Cal, Results)
    try
        H = Results.History;
        if isempty(H)
            return;
        end
        fig = figure('Color', 'w', 'Visible', 'off');
        tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        nexttile;
        scatter(H.Generation, H.Objective, 35, 'filled', 'MarkerFaceAlpha', 0.45);
        hold on;
        G = groupsummary(H, 'Generation', 'min', 'Objective');
        plot(G.Generation, G.min_Objective, 'k-', 'LineWidth', 1.8);
        xlabel('Generation');
        ylabel('Objective');
        grid on;
        title('GA objective evolution');
        nexttile;
        runtimeByGen = groupsummary(H, 'Generation', 'sum', 'Runtime_s');
        bar(runtimeByGen.Generation, runtimeByGen.sum_Runtime_s / 60);
        xlabel('Generation');
        ylabel('Runtime (min)');
        grid on;
        title('Evaluation runtime');
        exportgraphics(fig, fullfile(Cal.output_root, 'GA_Diagnostics.png'), 'Resolution', 180);
        close(fig);
    catch ME
        warning('GA diagnostic plot failed: %s', ME.message);
    end
end
