% VAL-D8LI-001: D8 local-inertial kernels, MATLAB against HydroPol2D-Python.
%
% The 9 x 11, 20 m V-tilted catchment of the Python D8 tests, in both modes
% (a D8 river graph inside a D4 floodplain, and the full-domain D8 graph).
% make_python_reference.py (part 1) wrote the Python numbers in this folder.
%
% The D8 graph, the river geometry and the hydraulic tables built by
% hp2d_d8_setup must equal Python's, and the tables read from Python's
% netCDF file must equal the tables built in MATLAB. Then
% D8_Local_Inertial_Model replays Python's steps (47 mm/h for 10 min,
% 30 min in all): link discharges, cell volumes, surface and channel
% depths, outlet flow and the link CFL step must match Python every 2 min,
% and rain - outflow - storage must close. A third run adds two reservoir
% controls on the river (one release leaves the model, one goes to a
% downstream river cell) and must still close the water balance.
% run_d8_full_model_check.m (VAL-D8LI-002) runs the whole model.
%
% Results: Outputs/Validation/Pass_Fail.csv (case_id, status, report_ready
% and the metrics).

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(case_dir, '..', '..', '..', 'HydroPol2D_Functions'), '-begin');
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end

RunSummary = d8_kernel_rows(case_dir);
writetable(RunSummary, fullfile(summary_dir, 'Pass_Fail.csv'));
disp(RunSummary(:, {'case_id', 'status', 'report_ready', 'metric', 'value', 'threshold'}));
failed = unique(RunSummary.case_id(RunSummary.status ~= "PASS"));
if ~isempty(failed)
    error('VAL-D8LI-001 failed: %s.', strjoin(failed, ', '));
end
fprintf('VAL-D8LI-001 passed.\n');

%% ------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function Rows = d8_kernel_rows(case_dir)
[z, n_land] = vtilted_dem(0.05);
Rows = table();
for mode = ["hybrid", "full"]
    Options = d8_kernel_options(mode == "full");
    D8 = hp2d_d8_setup(z, n_land, 20, Options);
    Options.table_path = fullfile(case_dir, 'python_river_tables.nc');
    D8_file = hp2d_d8_setup(z, n_land, 20, Options);

    % Graph and river geometry (Python's terrain products).
    G = readtable(fullfile(case_dir, 'reference_graph.csv'));
    idx = sub2ind(size(z), G.row, G.col);
    target = zeros(size(idx));
    has = G.downstream_row > 0;
    target(has) = sub2ind(size(z), G.downstream_row(has), G.downstream_col(has));
    graph_ok = isequal(D8.downstream(idx), target) && ...
        isequal(D8.stream(idx), G.river == 1) && ...
        isequal(D8.contributing_cells(idx), G.contributing_cells);
    geometry_error = max([abs(D8.river_width(idx) - G.river_width_m); ...
        abs(D8.river_depth(idx) - G.river_depth_m)]);

    % Tables: MATLAB-built against Python's (file and CSV).
    T = readtable(fullfile(case_dir, "reference_tables_" + mode + ".csv"));
    links_ok = isequal(D8.Links.source, sub2ind(size(z), T.source_row, T.source_col)) && ...
        isequal(D8.Links.target, sub2ind(size(z), T.target_row, T.target_col)) && ...
        isequal(double(D8.Links.diagonal), T.diagonal) && max(abs(D8.Links.length - T.length_m)) < 1e-12;
    top = sub2ind(size(D8.Faces.area), (1:numel(D8.Links.source)).', D8.Faces.count);
    table_error = max([relative(D8.Faces.datum, T.face_datum_m); ...
        relative(D8.Faces.area(:, 61), T.area_at_level_60_m2); ...
        relative(D8.Faces.conveyance(:, 61), T.conveyance_at_level_60); ...
        relative(D8.Faces.conveyance(top), T.conveyance_at_top); ...
        relative(D8.Cells.volume(:), D8_file.Cells.volume(:)); ...
        relative(D8.Cells.zeta(:), D8_file.Cells.zeta(:)); ...
        relative(D8.Cells.area(:), D8_file.Cells.area(:)); ...
        relative(D8.Cells.datum, D8_file.Cells.datum); ...
        relative(D8.Faces.conveyance(:), D8_file.Faces.conveyance(:))]);
    tables_ok = links_ok && isequal(D8.Cells.is_table, D8_file.Cells.is_table);

    % Step-by-step replay of Python's run.
    [M, P] = replay(case_dir, mode, D8, []);
    q_error = max(abs(M.discharge(:) - P.discharge(:))) / max(abs(P.discharge(:)));
    v_error = max(abs(M.volume(:) - P.volume(:))) / max(P.volume(:));
    depth_error = max([abs(M.surface_depth(:) - P.surface_depth(:)); ...
        abs(M.channel_depth(:) - P.channel_depth(:))]);
    outlet_error = max(abs(M.outlet - P.outlet)) / max(P.outlet);
    later = P.time > 0;
    cfl_error = max(abs(M.cfl(later) - P.cfl(later)) ./ P.cfl(later));
    mass_residual = abs(M.rain(end) - M.outflow(end) - M.storage(end)) / M.rain(end);
    peak_q = max(abs(P.discharge(:)));
    diagonal_share = D8_States_share(M);

    name = "VAL-D8LI-001-" + upper(mode) + "-KERNEL";
    Rows = [Rows; result_rows(name, "D8 " + mode + " kernel vs Python (" + numel(D8.Links.source) + " links)", ...
        ["graph_and_geometry_identical", "river_geometry_abs_error_m", "tables_identical", ...
         "table_rel_error", "link_discharge_rel_error", "cell_volume_rel_error", ...
         "depth_abs_error_m", "outlet_rel_error", "cfl_rel_error", "mass_residual_rel", ...
         "peak_link_discharge_m3_s", "diagonal_transfer_share"], ...
        [double(graph_ok), geometry_error, double(tables_ok), table_error, q_error, v_error, ...
         depth_error, outlet_error, cfl_error, mass_residual, peak_q, diagonal_share], ...
        [1, 1e-12, 1, 1e-12, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-12, NaN, NaN], ...
        ["eq", "le", "eq", "le", "le", "le", "le", "le", "le", "le", "info", "info"])]; %#ok<AGROW>
end

% Reservoir controls on the D8 river (hybrid). MATLAB-only: Python refuses
% reservoirs with D8 routing.
Options = d8_kernel_options(false);
D8 = hp2d_d8_setup(z, n_land, 20, Options);
R = struct('x_index', 6, 'y_index', 5, 'k1', 0.5, 'h1', 0.05, 'k2', 1.5, ...
    'y_ds1_index', 0, 'x_ds1_index', 0, 'k3', 0.3, 'h2', 0.1, 'k4', 1.5, ...
    'y_ds2_index', 7, 'x_ds2_index', 6);
M = replay(case_dir, "hybrid", D8, R);
mass_residual = abs(M.rain(end) - M.outflow(end) - M.storage(end)) / M.rain(end);
Rows = [Rows; result_rows("VAL-D8LI-001-RESERVOIR", "D8 hybrid with two reservoir controls", ...
    ["mass_residual_rel", "external_release_m3", "outlet_volume_m3"], ...
    [mass_residual, M.external_release, M.outflow(end) - M.external_release], ...
    [1e-12, 0, 0], ["le", "gt", "gt"])];
end

function [M, P] = replay(case_dir, mode, D8, Reservoir)
steps = readtable(fullfile(case_dir, "reference_" + mode + "_steps.csv"));
L = readtable(fullfile(case_dir, "reference_" + mode + "_links.csv"));
C = readtable(fullfile(case_dir, "reference_" + mode + "_cells.csv"));
S = readtable(fullfile(case_dir, "reference_" + mode + "_series.csv"));
times = S.time_s;
n_links = numel(D8.Links.source);
P = struct('time', times, 'outlet', S.outlet_m3_s, 'cfl', S.cfl_dt_unclipped_s);
P.discharge = reshape(L.discharge_m3_s, n_links, []);
cells = sub2ind([D8.ny D8.nx], C.row, C.col);
P.volume = zeros(D8.N, numel(times));
P.surface_depth = zeros(D8.N, numel(times));
P.channel_depth = zeros(D8.N, numel(times));
for r = 1:numel(times)
    rows = C.time_s == times(r);
    P.volume(cells(rows), r) = C.volume_m3(rows);
    P.surface_depth(cells(rows), r) = C.surface_depth_m(rows);
    P.channel_depth(cells(rows), r) = C.channel_depth_m(rows);
end

A = D8.cell_area;
idx_nan = ~D8.valid_grid;
d = zeros(D8.ny, D8.nx);
outflow = zeros(D8.ny, D8.nx, 3);
States = struct('discharge_m3_s', zeros(n_links, 1), 'cumulative_cardinal_transfer_m3', 0, ...
    'cumulative_diagonal_transfer_m3', 0, 'last_outlet_m3_s', 0);
use_reservoir = ~isempty(Reservoir);
M = struct('discharge', zeros(n_links, numel(times)), 'volume', zeros(D8.N, numel(times)), ...
    'surface_depth', zeros(D8.N, numel(times)), 'channel_depth', zeros(D8.N, numel(times)), ...
    'outlet', zeros(numel(times), 1), 'cfl', zeros(numel(times), 1), ...
    'rain', zeros(numel(times), 1), 'outflow', zeros(numel(times), 1), ...
    'storage', zeros(numel(times), 1), 'external_release', 0);
rain_m3 = 0;
out_m3 = 0;
record = 1;
store(1, 0);
for s = 1:height(steps)
    dt = steps.dt_s(s);
    rain_mm = steps.rain_mm_h(s) * dt / 3600;
    d = d + rain_mm;
    rain_m3 = rain_m3 + rain_mm / 1000 * A * D8.N;
    [d, outlet_flow, outflow, ~, ~, ~, ~, ~, States] = D8_Local_Inertial_Model( ...
        d, outflow, D8, States, dt / 60, idx_nan, 0, double(use_reservoir), Reservoir);
    out_m3 = out_m3 + sum(outlet_flow(:)) / 1000 / 3600 * A * dt;
    if use_reservoir
        % The D8 outlet is at the bottom-centre cell; outflow booked at the
        % control's cell is the release that left the model.
        M.external_release = M.external_release + ...
            outlet_flow(Reservoir.y_index, Reservoir.x_index) / 1000 / 3600 * A * dt;
    end
    t_end = steps.t_start_s(s) + dt;
    if record < numel(times) && abs(t_end - times(record + 1)) < 1e-6
        record = record + 1;
        store(record, sum(outlet_flow(:)) / 1000 / 3600 * A);
    end
end
M.diagonal = States.cumulative_diagonal_transfer_m3;
M.cardinal = States.cumulative_cardinal_transfer_m3;

    function store(r, outlet_m3_s)
        volume = d(:) / 1000 * A;
        M.discharge(:, r) = States.discharge_m3_s;
        M.volume(:, r) = volume;
        [surface_mm, channel_mm] = hp2d_d8_surface_depth(D8, d);
        M.surface_depth(:, r) = surface_mm(:) / 1000;
        M.channel_depth(:, r) = channel_mm(:) / 1000;
        M.outlet(r) = outlet_m3_s;
        M.cfl(r) = hp2d_d8_stable_timestep(volume, D8, 0.5, 1e9);
        M.rain(r) = rain_m3;
        M.outflow(r) = out_m3;
        M.storage(r) = sum(volume);
    end
end

function share = D8_States_share(M)
share = M.diagonal / max(M.diagonal + M.cardinal, eps);
end

function [z, n_land] = vtilted_dem(manning_n)
[col, row] = meshgrid(0:10, 0:8);
z = 0.004 * abs(col - 5) * 20 + 0.004 * (8 - row) * 20;
n_land = manning_n * ones(size(z));
end

function Options = d8_kernel_options(full_domain)
outlet = false(9, 11);
outlet(9, 6) = true;
Options = struct('full_domain', double(full_domain), 'min_depth_m', 1e-3, ...
    'min_area_km2', 0.0012, 'width_coefficient', 5, 'width_exponent', 0.2, ...
    'depth_coefficient', 0.15, 'depth_exponent', 0.2, 'channel_manning_n', 0.03, ...
    'river_width', [], 'river_depth', [], 'outlet_mask', outlet, 'outlet_slope', 0.004, ...
    'outlet_type', 1, 'table_path', '', 'max_stage_above_bed_m', 100);
end

function e = relative(a, b)
e = max(abs(a(:) - b(:))) / max(max(abs(b(:))), eps);
end

function Rows = result_rows(case_id, case_name, metrics, values, thresholds, rules)
passed = true(numel(metrics), 1);
for k = 1:numel(metrics)
    switch rules(k)
        case "eq"
            passed(k) = values(k) == thresholds(k);
        case "le"
            passed(k) = values(k) <= thresholds(k);
        case "gt"
            passed(k) = values(k) > thresholds(k);
        case "ge"
            passed(k) = values(k) >= thresholds(k);
        otherwise
            passed(k) = isfinite(values(k));
    end
end
status = "PASS";
if ~all(passed); status = "FAIL"; end
Rows = table(repmat(case_id, numel(metrics), 1), repmat(case_name, numel(metrics), 1), ...
    repmat(status, numel(metrics), 1), true(numel(metrics), 1), metrics(:), values(:), ...
    thresholds(:), rules(:), passed, ...
    'VariableNames', {'case_id', 'case_name', 'status', 'report_ready', 'metric', ...
    'value', 'threshold', 'rule', 'metric_passed'});
end
