% VAL-UD-001: urban drainage runtime against the HydroPol2D-Python kernel.
%
% 1. VAL-UD-001  network_kernel.gpkg (written by HydroPol2D-Python) is read
%    with each available GeoPackage reader and run for 400 surface steps of
%    5 s through hp2d_urban_drainage_advance, with the sources of
%    python_reference_kernel.csv. With free_outfall_depth = "invert" (the
%    Python rule) node heads and volumes, pipe and inlet discharges, inlet
%    cell volumes, outfall volumes and sub-step counts must match the Python
%    engine (_compiled_urban_drainage_advance), and every step must conserve
%    water exactly. The run has capture, surcharge back to the surface at
%    five inlets, pipe-flow reversals and backflow through a fixed-stage
%    outfall. Two altered runs must miss the reference: reverse flow off,
%    and the default free-outfall rule.
% 2. VAL-UD-001R the default free-outfall rule (min of critical and normal
%    depth): a chain of four 50 m pipes to a free outfall carrying a steady
%    0.12 m3/s must reach normal depth upstream of the outfall; the Python
%    rule (outfall at its invert) leaves that node far deeper (review F2).
%    The same 400-step network run with the default rule conserves water.
% VAL-UD-001F, the full grid run, is run_urban_drainage_grid_run.m (its own
% driver, because preprocessing clears the workspace). Both write their rows
% to Outputs/Validation/Pass_Fail.csv and keep the other driver's rows.

clear; clc;
case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
addpath(fullfile(repo_root, 'HydroPol2D_Functions'), '-begin');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
pass_path = fullfile(out_dir, 'Pass_Fail.csv');

Ref = readtable(fullfile(case_dir, 'python_reference_kernel.csv'), 'CommentStyle', '#', ...
    'VariableNamingRule', 'preserve');
terrain = readmatrix(fullfile(case_dir, 'python_reference_inputs.csv'), 'CommentStyle', '#');
fid = fopen(fullfile(case_dir, 'python_reference_kernel.csv'));
python_header = strtrim(erase(fgetl(fid), '#'));
fclose(fid);
gpkg = fullfile(case_dir, 'network_kernel.gpkg');

% ---------------------------------------------------------------------
% Readers
% ---------------------------------------------------------------------
backends = {'mapping', 'sqlite'};
Networks = struct();
read_ok = false(1, numel(backends));
for i = 1:numel(backends)
    try
        Networks.(backends{i}) = hp2d_urban_drainage_read(gpkg, backends{i});
        read_ok(i) = true;
    catch read_error
        fprintf('GeoPackage reader %s unavailable: %s\n', backends{i}, read_error.message);
    end
end
if ~any(read_ok)
    error('No GeoPackage reader is available.');
end
readers_agree = true;
if all(read_ok)
    for layer = {'nodes', 'pipes', 'inlets', 'outfalls'}
        readers_agree = readers_agree && tables_equal(Networks.mapping.(layer{1}), Networks.sqlite.(layer{1}));
    end
end
Network = Networks.(backends{find(read_ok, 1)});

% ---------------------------------------------------------------------
% 1. Python rule against the Python kernel
% ---------------------------------------------------------------------
python_rule = struct('free_outfall_depth', "invert");
[Run, N] = replay(Network, terrain, python_rule, Ref);
Err = compare(Run, Ref, N);
coverage = min(Run.inlet_q(:)) < -1e-3 && max(Run.inlet_q(:)) > 1e-3 && ...
    sum(any(Run.inlet_q < -1e-3, 1)) >= 4 && Run.final.pipe_reversal_count' * ones(numel(N.pipe_ids), 1) > 0 && ...
    sum(Run.inflow) > 1 && sum(Run.outflow) > 1 && max(Run.substeps) >= 2 && ...
    all(ismember([0 1 2], N.node_boundary_type));
% Altered runs that must miss the reference.
NoReverse = Network;
NoReverse.inlets.reverse_flow_allowed(:) = 0;
RunNoReverse = replay(NoReverse, terrain, python_rule, Ref);
ErrNoReverse = compare(RunNoReverse, Ref, N);
RunDefault = replay(Network, terrain, struct(), Ref);
ErrDefault = compare(RunDefault, Ref, N);

tol = struct('head_m', 1e-8, 'volume_m3', 1e-8, 'pipe_q', 1e-9, 'inlet_q', 1e-9, ...
    'surface_m3', 1e-8, 'boundary_m3', 1e-9, 'mass_m3', 1e-9);
parity_ok = Err.head_m < tol.head_m && Err.volume_m3 < tol.volume_m3 && ...
    Err.pipe_q < tol.pipe_q && Err.inlet_q < tol.inlet_q && Err.surface_m3 < tol.surface_m3 && ...
    Err.boundary_m3 < tol.boundary_m3 && Err.substeps == 0 && Run.max_mass_residual_m3 < tol.mass_m3;
controls_ok = ErrNoReverse.head_m > 1e-2 && ErrDefault.pipe_q > 1e-4;
passed_parity = read_ok(1) && readers_agree && parity_ok && coverage && controls_ok;

Metrics = table("VAL-UD-001", string(python_header), height(Ref), Err.head_m, Err.volume_m3, ...
    Err.pipe_q, Err.inlet_q, Err.surface_m3, Err.boundary_m3, Err.substeps, ...
    Run.max_mass_residual_m3, min(Run.inlet_q(:)), sum(Run.outflow), sum(Run.inflow), ...
    sum(Run.final.pipe_reversal_count), max(Run.substeps), ErrNoReverse.head_m, ...
    ErrDefault.pipe_q, readers_agree, string(strjoin(backends(read_ok), '+')), passed_parity, ...
    'VariableNames', {'case_id', 'python_reference', 'steps', 'max_head_error_m', ...
    'max_volume_error_m3', 'max_pipe_q_error_m3_s', 'max_inlet_q_error_m3_s', ...
    'max_surface_volume_error_m3', 'max_outfall_volume_error_m3', 'substep_count_mismatches', ...
    'max_step_mass_residual_m3', 'min_inlet_q_m3_s', 'outfall_outflow_m3', 'outfall_inflow_m3', ...
    'pipe_reversals', 'max_substeps', 'no_reverse_flow_head_error_m', ...
    'default_free_outfall_pipe_q_difference_m3_s', 'readers_agree', 'readers', 'passed'});

% ---------------------------------------------------------------------
% 2. Default free-outfall rule
% ---------------------------------------------------------------------
[chain_default, normal_depth] = free_outfall_chain(struct());
chain_python = free_outfall_chain(python_rule);
default_error = chain_default.depth(end) / normal_depth - 1;
python_error = chain_python.depth(end) / normal_depth - 1;
rule_ok = abs(default_error) < 0.05 && abs(chain_default.depth(end - 1) / normal_depth - 1) < 0.05 && ...
    python_error > 0.2 && abs(chain_default.outflow_m3_s / 0.12 - 1) < 1e-6 && ...
    RunDefault.max_mass_residual_m3 < tol.mass_m3;
RuleMetrics = table("VAL-UD-001R", normal_depth, chain_default.depth(end), chain_python.depth(end), ...
    default_error, python_error, chain_default.outflow_m3_s, RunDefault.max_mass_residual_m3, ...
    ErrDefault.head_m, rule_ok, 'VariableNames', {'case_id', 'normal_depth_m', ...
    'default_rule_depth_above_outfall_m', 'python_rule_depth_above_outfall_m', ...
    'default_rule_relative_error', 'python_rule_relative_error', 'steady_outflow_m3_s', ...
    'default_rule_max_step_mass_residual_m3', 'default_rule_max_head_difference_from_python_m', 'passed'});

writetable(Metrics, fullfile(out_dir, 'UrbanDrainage_Kernel_Metrics.csv'));
writetable(RuleMetrics, fullfile(out_dir, 'UrbanDrainage_FreeOutfall_Metrics.csv'));
writetable(series_table(Run, N), fullfile(out_dir, 'UrbanDrainage_Kernel_MATLAB_Series.csv'));
make_figure(Run, Ref, N, fullfile(out_dir, 'urban_drainage_kernel_vs_python.png'));
disp(Metrics); disp(RuleMetrics);

Pass = [pass_row("VAL-UD-001", passed_parity, Err.head_m, ...
    "Kernel with the Python free-outfall rule against HydroPol2D-Python, 400 steps; exact step mass balance."); ...
    pass_row("VAL-UD-001R", rule_ok, abs(default_error), ...
    "Free outfall at min(critical, normal) depth: normal depth above the outfall; Python rule leaves it deeper.")];
if isfile(pass_path)
    Old = readtable(pass_path, 'TextType', 'string', 'Delimiter', ',');
    Old(ismember(Old.case_id, Pass.case_id), :) = [];
    Pass = [Pass; Old];
end
writetable(Pass, pass_path);
disp(Pass);

% =====================================================================
function [Run, N] = replay(Network, terrain, hydraulics, Ref)
N = hp2d_urban_drainage_prepare(Network, terrain, true(size(terrain)), 100, hydraulics);
S = hp2d_urban_drainage_init(N);
cells = N.inlet_unique_cell;
surface = zeros(size(terrain));
surface(4, 5) = 30;  % the initial pond on cell (3, 4), 0-based
n = height(Ref);
nn = numel(N.node_ids); np = numel(N.pipe_ids); ni = numel(N.inlet_ids);
Run = struct('head', zeros(n, nn), 'volume', zeros(n, nn), 'pipe_q', zeros(n, np), ...
    'inlet_q', zeros(n, ni), 'surface', zeros(n, numel(cells)), 'outflow', zeros(n, 1), ...
    'inflow', zeros(n, 1), 'substeps', zeros(n, 1), 'mass_residual', zeros(n, 1));
for k = 1:n
    before = sum(surface(:)) + sum(S.node_volume_m3);
    surface(1, 1) = surface(1, 1) + Ref.source_00_m3(k);
    surface(2, 3) = surface(2, 3) + Ref.source_12_m3(k);
    [S, values, step] = hp2d_urban_drainage_advance(N, S, surface(cells), 5);
    surface(cells) = values;
    after = sum(surface(:)) + sum(S.node_volume_m3);
    Run.mass_residual(k) = after - before - (Ref.source_00_m3(k) + Ref.source_12_m3(k) + ...
        step.inflow_m3 - step.outflow_m3);
    Run.head(k, :) = S.node_head_m';
    Run.volume(k, :) = S.node_volume_m3';
    Run.pipe_q(k, :) = S.pipe_q_m3_s';
    Run.inlet_q(k, :) = S.inlet_q_m3_s';
    Run.surface(k, :) = surface(cells)';
    Run.outflow(k) = step.outflow_m3;
    Run.inflow(k) = step.inflow_m3;
    Run.substeps(k) = step.substeps;
end
Run.max_mass_residual_m3 = max(abs(Run.mass_residual));
Run.final = S;
Run.cells = cells;
Run.grid_size = size(terrain);
end

function Err = compare(Run, Ref, N)
column = @(prefix, ids) cell2mat(cellfun(@(id) Ref.([prefix char(id)]), cellstr(ids(:))', 'UniformOutput', false));
[r, c] = ind2sub(Run.grid_size, Run.cells);
surface_names = arrayfun(@(i) sprintf('r%dc%d', r(i) - 1, c(i) - 1), 1:numel(r), 'UniformOutput', false);
Err.head_m = max(abs(Run.head - column('h_', N.node_ids)), [], 'all');
Err.volume_m3 = max(abs(Run.volume - column('V_', N.node_ids)), [], 'all');
Err.pipe_q = max(abs(Run.pipe_q - column('q_', N.pipe_ids)), [], 'all');
Err.inlet_q = max(abs(Run.inlet_q - column('qi_', N.inlet_ids)), [], 'all');
Err.surface_m3 = max(abs(Run.surface - column('S_', string(surface_names))), [], 'all');
Err.boundary_m3 = max([abs(Run.outflow - Ref.outflow_m3); abs(Run.inflow - Ref.inflow_m3)]);
Err.substeps = sum(Run.substeps ~= Ref.substeps);
end

function [Result, normal_depth] = free_outfall_chain(hydraulics)
% Four 50 m pipes (0.5 m, n = 0.013, slope 0.004) to a free outfall, with
% 0.12 m3/s added to the top node; 2 h of 5 s steps reach a steady state.
L = 50; D = 0.5; slope = 0.004; n = 0.013; Q0 = 0.12; np = 4;
ids = ["J1"; "J2"; "J3"; "J4"; "OF"];
invert = 10 - slope * L * (0:np)';
nodes = table(ids, invert + 2.5, invert, ones(np + 1, 1), NaN(np + 1, 1), 'VariableNames', ...
    {'node_id', 'rim_elevation_m', 'invert_elevation_m', 'shaft_area_m2', 'initial_head_m'});
pipes = table("P" + (1:np)', ids(1:np), ids(2:end), L * ones(np, 1), D * ones(np, 1), ...
    n * ones(np, 1), invert(1:np), invert(2:end), 'VariableNames', {'pipe_id', 'upstream_node', ...
    'downstream_node', 'length_m', 'diameter_m', 'manning_n', 'upstream_invert_m', 'downstream_invert_m'});
inlets = table(strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'identifier', 'node_id', 'row', 'column', 'crest_elevation_m'});
outfalls = table("OF", "free", NaN, 'VariableNames', {'node_id', 'boundary_type', 'fixed_stage_m'});
Net = struct('nodes', nodes, 'pipes', pipes, 'inlets', inlets, 'outfalls', outfalls);
N = hp2d_urban_drainage_prepare(Net, zeros(2), true(2), 100, hydraulics);
S = hp2d_urban_drainage_init(N);
for k = 1:1440
    S.node_volume_m3(1) = S.node_volume_m3(1) + Q0 * 5;
    S.node_head_m = hp2d_urban_drainage_head(N, S.node_volume_m3, S.node_head_m);
    [S, ~, step] = hp2d_urban_drainage_advance(N, S, zeros(0, 1), 5);
end
Result.depth = S.node_head_m(1:np) - invert(1:np);
Result.outflow_m3_s = step.outflow_m3 / 5;
% Manning normal depth of Q0 by bisection below the peak-conveyance depth.
target = Q0 * n / sqrt(slope);
low = 0; high = 0.9382 * D;
for iteration = 1:60
    middle = 0.5 * (low + high);
    [A, R] = hp2d_urban_drainage_section(middle, D);
    if A * R^(2/3) < target, low = middle; else, high = middle; end
end
normal_depth = 0.5 * (low + high);
end

function tf = tables_equal(A, B)
tf = isequal(sort(A.Properties.VariableNames), sort(B.Properties.VariableNames)) && height(A) == height(B);
if ~tf; return; end
names = A.Properties.VariableNames;
for i = 1:numel(names)
    a = A.(names{i}); b = B.(names{i});
    if isnumeric(a)
        tf = tf && isequaln(double(a), double(b));
    else
        tf = tf && isequal(string(a), string(b));
    end
end
end

function T = series_table(Run, N)
names = [("h_" + N.node_ids(:))', ("q_" + N.pipe_ids(:))', ("qi_" + N.inlet_ids(:))'];
T = array2table([(1:numel(Run.outflow))' * 5, Run.head, Run.pipe_q, Run.inlet_q, Run.outflow, Run.inflow], ...
    'VariableNames', [{'t_end_s'}, cellstr(names), {'outflow_m3', 'inflow_m3'}]);
end

function make_figure(Run, Ref, N, path)
t = Ref.t_end_s / 60;
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 600]);
tiledlayout(fig, 2, 1);
ax = nexttile; hold(ax, 'on');
for i = find(N.dynamic)'
    plot(ax, t, Ref.(['h_' char(N.node_ids(i))]), '-', 'LineWidth', 1.6, 'DisplayName', char(N.node_ids(i)));
    plot(ax, t, Run.head(:, i), 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
ylabel(ax, 'Node head (m)'); legend(ax, 'Location', 'eastoutside'); grid(ax, 'on');
title(ax, 'HydroPol2D-Python (lines) and MATLAB (dotted), Python free-outfall rule');
ax = nexttile; hold(ax, 'on');
for i = 1:numel(N.inlet_ids)
    plot(ax, t, Ref.(['qi_' char(N.inlet_ids(i))]), '-', 'LineWidth', 1.6, 'DisplayName', char(N.inlet_ids(i)));
    plot(ax, t, Run.inlet_q(:, i), 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
xlabel(ax, 'Time (min)'); ylabel(ax, 'Inlet exchange (m^3/s, + into pipes)');
legend(ax, 'Location', 'eastoutside'); grid(ax, 'on');
exportgraphics(fig, path, 'Resolution', 150);
close(fig);
end

function row = pass_row(case_id, passed, max_error, notes)
status = "fail";
if passed; status = "pass"; end
row = table(string(case_id), status, double(passed), max_error, NaN, NaN, string(notes), ...
    'VariableNames', {'case_id', 'status', 'report_ready', 'max_error', 'outfall_outflow_m3', ...
    'outfall_inflow_m3', 'notes'});
end
