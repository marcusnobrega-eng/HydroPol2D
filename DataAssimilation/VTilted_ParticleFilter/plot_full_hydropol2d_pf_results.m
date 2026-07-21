function plot_full_hydropol2d_pf_results()
%PLOT_FULL_HYDROPOL2D_PF_RESULTS Figures for the full HydroPol2D PF run.

case_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(case_dir, 'Outputs');
fig_dir = fullfile(case_dir, 'Figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

P = readtable(fullfile(out_dir, 'FullHydroPol2D_Parameter_Evolution.csv'));
F = readtable(fullfile(out_dir, 'FullHydroPol2D_Observation_Fit.csv'));
O = readtable(fullfile(out_dir, 'FullHydroPol2D_Synthetic_Observations.csv'));
E = readtable(fullfile(out_dir, 'FullHydroPol2D_Effective_Sample_Size.csv'));

plot_parameters(P, fig_dir);
plot_observations(F, O, fig_dir);
plot_ess(E, fig_dir);
end

function plot_parameters(P, fig_dir)
params = unique(string(P.parameter), 'stable');
n = numel(params);
rows = ceil(n / 3);
fig = figure('Color', 'w', 'Position', [100 100 1500 240 * rows]);
tiledlayout(rows, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:n
    nexttile;
    idx = string(P.parameter) == params(i);
    T = P(idx, :);
    times = unique(T.time_min);
    q05 = nan(size(times)); q50 = q05; q95 = q05; truth = q05; wmean = q05;
    for k = 1:numel(times)
        Tk = T(T.time_min == times(k), :);
        [q05(k), q50(k), q95(k)] = weighted_quantiles(Tk.value, Tk.weight, [0.05 0.50 0.95]);
        wmean(k) = sum(Tk.value .* Tk.weight) / sum(Tk.weight);
        truth(k) = Tk.truth_value(1);
    end
    fill([times; flipud(times)], [q05; flipud(q95)], [0.70 0.82 0.95], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
    plot(times, q50, 'b-', 'LineWidth', 1.2);
    plot(times, wmean, 'Color', [0 0.35 0.65], 'LineWidth', 1.0);
    plot(times, truth, 'k--', 'LineWidth', 1.1);
    title(strrep(params(i), '_', '\_'), 'FontSize', 8);
    xlabel('Time [min]'); ylabel('Multiplier [-]');
    grid on; box on;
end
legend({'5-95%', 'median', 'weighted mean', 'truth'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
exportgraphics(fig, fullfile(fig_dir, 'FullHydroPol2D_parameter_evolution.png'), 'Resolution', 200);
close(fig);
end

function plot_observations(F, O, fig_dir)
obs_ids = unique(string(F.obs_id), 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1400 320 * numel(obs_ids)]);
tiledlayout(numel(obs_ids), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numel(obs_ids)
    nexttile;
    idx = string(F.obs_id) == obs_ids(i);
    T = F(idx, :);
    times = unique(T.time_min);
    q05 = nan(size(times)); q50 = q05; q95 = q05; wmean = q05; obs = q05; truth = q05;
    for k = 1:numel(times)
        Tk = T(T.time_min == times(k), :);
        [q05(k), q50(k), q95(k)] = weighted_quantiles(Tk.simulated, Tk.weight, [0.05 0.50 0.95]);
        wmean(k) = sum(Tk.simulated .* Tk.weight) / sum(Tk.weight);
        Ok = O(string(O.obs_id) == obs_ids(i) & O.time_min == times(k), :);
        obs(k) = Ok.observed(1);
        truth(k) = Ok.truth_value(1);
    end
    fill([times; flipud(times)], [q05; flipud(q95)], [0.72 0.85 0.72], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on;
    plot(times, q50, 'g-', 'LineWidth', 1.2);
    plot(times, wmean, 'Color', [0 0.45 0], 'LineWidth', 1.0);
    plot(times, truth, 'k--', 'LineWidth', 1.1);
    plot(times, obs, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    ylabel(sprintf('%s [%s]', char(obs_ids(i)), char(T.units(1))), 'Interpreter', 'none');
    xlabel('Time [min]');
    grid on; box on;
end
legend({'5-95%', 'median', 'weighted mean', 'truth', 'observation'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
exportgraphics(fig, fullfile(fig_dir, 'FullHydroPol2D_observation_fit.png'), 'Resolution', 200);
close(fig);
end

function plot_ess(E, fig_dir)
fig = figure('Color', 'w', 'Position', [100 100 900 420]);
yyaxis left;
plot(E.time_min, E.effective_sample_size, 'o-', 'LineWidth', 1.5);
ylabel('Effective sample size [-]');
yyaxis right;
plot(E.time_min, E.weighted_objective, 's-', 'LineWidth', 1.5);
ylabel('Weighted objective [-]');
xlabel('Time [min]');
grid on; box on;
title('Full HydroPol2D PF Diagnostics');
exportgraphics(fig, fullfile(fig_dir, 'FullHydroPol2D_ess_objective.png'), 'Resolution', 200);
close(fig);
end

function varargout = weighted_quantiles(x, w, probs)
x = x(:); w = w(:);
valid = isfinite(x) & isfinite(w) & w > 0;
x = x(valid); w = w(valid);
if isempty(x) || sum(w) <= 0
    varargout = num2cell(nan(size(probs)));
    return;
end
[x, order] = sort(x);
w = w(order) ./ sum(w);
cdf = cumsum(w);
if numel(x) == 1
    varargout = num2cell(repmat(x, size(probs)));
    return;
end
[cdf, unique_idx] = unique(cdf, 'stable');
x = x(unique_idx);
q = interp1(cdf, x, probs, 'linear', 'extrap');
q = min(max(q, min(x)), max(x));
varargout = num2cell(q);
end
