% VAL-GI-001: green infrastructure in a full MATLAB grid run, against the
% HydroPol2D-Python kernel.
%
% A 10 x 7 tilted V (10 m cells) gets 45 min of 60 mm/h rain and 55 min of
% reference ET, in bypass mode, twice: legacy soil and the layered column
% (flag_layered_subsurface). A rain garden in the valley (3 compartments,
% lateral flow, underdrain to the outlet cell) and permeable pavement on the
% left hillslope (2 compartments, underdrain out of the model) capture
% run-on, overflow, infiltrate, percolate, evaporate and exfiltrate into the
% soil. For both runs the check asks:
%   1. every GI step books its transfers (host + soil + GI stores + water
%      out + ET = GI rain, in double), and the model states agree to their
%      precision (single here: the states take the class of the inputs);
%   2. the GI stores balance their own fluxes;
%   3. the event water balance closes to 1e-5 of the rain (system ledger
%      with GI storage, GI ET and the external underdrain as outlet flow;
%      the same case with GI off closes to 2e-4 to 2e-3 m3, 5e-6 of the rain);
%   4. capture, overflow, native exfiltration, lateral flow and both
%      underdrains are active, and the gross LAI is back after the run;
%   5. the compartment geometry matches Python's prepare step, and the
%      MATLAB kernel replayed on the GI inputs recorded from the run
%      (Reference/replay_inputs_<variant>.csv) gives the states and
%      cumulative fluxes of the Python kernel on the same inputs
%      (Reference/python_reference_<variant>_*.csv).
% A 4 x 5 limiter case (Reference/limiter_case.json) checks the one rule
% MATLAB does not copy (review C14): Python caps each lateral edge by its
% donor alone and creates water; MATLAB must conserve and match Python run
% with MATLAB's cap. Controls that must fail: the ledger without the outlet
% booking or without the GI storage, and a replay with the overflow
% coefficient 10% low.
%
% Regenerate the reference: run with HYDROPOL_GI_WRITE_REPLAY=1 (writes the
% replay inputs and gi_case.json), then Reference/make_python_reference.py.
% HYDROPOL_GI_VARIANTS=legacy (or layered) runs one variant only.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
addpath(fullfile(repo_root, 'HydroPol2D_Functions'), '-begin');
hydropol2d_add_runtime_paths(hydropol2d_find_root(case_dir));
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
reference_dir = fullfile(case_dir, 'Reference');
static_dir = fullfile(case_dir, 'Outputs', 'Static');
forcing_dir = fullfile(case_dir, 'Outputs', 'Forcing');
for d = {summary_dir, static_dir, forcing_dir}
    if ~exist(d{1}, 'dir'); mkdir(d{1}); end
end
write_reference = strcmp(getenv('HYDROPOL_GI_WRITE_REPLAY'), '1');

write_static_rasters(static_dir);
etp_path = fullfile(forcing_dir, 'ETP_input_data.xlsx');
write_station_forcing(etp_path);

variants = {'legacy', 'layered'};
if ~isempty(getenv('HYDROPOL_GI_VARIANTS'))   % e.g. 'layered' while debugging
    variants = strsplit(getenv('HYDROPOL_GI_VARIANTS'), ',');
end
Checks = table();
Summary = table();
for iv = 1:numel(variants)
    variant = variants{iv};
    R = run_variant(case_dir, variant, static_dir, etp_path);
    log_table = replay_inputs_from_log(R.G, R.S.step_log);
    writetable(log_table, fullfile(summary_dir, ['live_replay_inputs_' variant '.csv']));
    if write_reference
        writetable(log_table, fullfile(reference_dir, ['replay_inputs_' variant '.csv']));
        if iv == 1
            write_case_json(fullfile(reference_dir, 'gi_case.json'), R);
        end
    end
    writetable(R.history, fullfile(summary_dir, ['GI_Facilities_' variant '.csv']));
    save(fullfile(summary_dir, ['GI_Run_' variant '.mat']), 'R');

    tag = upper(variant);
    G = R.G; S = R.S;
    gi_storage = sum(S.surface_m3) + sum(S.media_m3) + sum(S.gravel_m3);
    gi_storage0 = sum(G.initial_media_water_volume_m3);
    gi_residual = (gi_storage - gi_storage0) - (sum(S.cumulative_rain_m3) + ...
        sum(S.cumulative_capture_m3) - sum(S.cumulative_overflow_m3) - sum(S.cumulative_et_m3) - ...
        sum(S.cumulative_native_exfiltration_m3) - sum(S.cumulative_underdrain_m3));
    L = R.ledger;
    net = L.cumulative_precipitation_m3 + L.cumulative_boundary_inflow_m3 + ...
        L.cumulative_prescribed_recharge_m3 - L.cumulative_canopy_evaporation_m3 - ...
        L.cumulative_etr_m3 - L.cumulative_open_water_evaporation_m3 - ...
        L.cumulative_snow_sublimation_m3 - L.cumulative_outlet_m3 - ...
        L.cumulative_green_infrastructure_et_m3;
    system_residual = (L.final_storage_m3 - L.initial_storage_m3) - net;
    rain = L.cumulative_precipitation_m3;
    rg = G.facility_index == 1;
    pv = G.facility_index == 2;
    underdrain_surface = sum(S.cumulative_underdrain_m3(rg));
    underdrain_external = S.cumulative_external_m3;

    [replay_error, replay_note, geometry_error, mutated_error] = compare_with_python( ...
        reference_dir, variant, G, write_reference);

    coupling_error = cellfun(@(x) x.coupling_error_m3, S.step_log);
    bookkeeping_error = cellfun(@(x) x.bookkeeping_error_m3, S.step_log);
    Checks = [Checks; ...
        check_row(tag, "GI step books every transfer: host + soil + GI + out - rain [m3]", ...
            max(abs(bookkeeping_error)), 1e-9); ...
        check_row(tag, "GI step on the model states (single precision here) [m3]", ...
            max(abs(coupling_error)), 1e-4); ...
        check_row(tag, "GI stores balance their fluxes [m3]", abs(gi_residual), 1e-9); ...
        check_row(tag, "event water balance (system ledger) [m3]", abs(system_residual), 1e-5 * rain); ...
        check_row(tag, "control: ledger without the outlet booking [m3]", ...
            abs(system_residual - underdrain_external), 1e-5 * rain, true); ...
        check_row(tag, "control: ledger without the GI storage [m3]", ...
            abs(system_residual - (gi_storage - gi_storage0)), 1e-5 * rain, true); ...
        check_row(tag, "capture into the GI [m3]", sum(S.cumulative_capture_m3), 1, true); ...
        check_row(tag, "overflow from the rain garden [m3]", sum(S.cumulative_overflow_m3(rg)), 0.1, true); ...
        check_row(tag, "underdrain to the outlet cell [m3]", underdrain_surface, 0.1, true); ...
        check_row(tag, "underdrain out of the model (booked as outlet) [m3]", underdrain_external, 0.01, true); ...
        check_row(tag, "native exfiltration into the soil [m3]", ...
            sum(S.cumulative_native_exfiltration_m3), 0.01, true); ...
        check_row(tag, "lateral flow between GI compartments [m3/s]", R.max_edge_q, 1e-6, true); ...
        check_row(tag, "gross LAI restored after the run [-]", R.lai_error, 0); ...
        check_row(tag, "geometry vs Python prepare [m2 or m]", geometry_error, 1e-9); ...
        check_row(tag, "kernel replay vs Python kernel [m3]", replay_error, 1e-9); ...
        check_row(tag, "control: replay with overflow coefficient -10% [m3]", mutated_error, 1e-6, true)];

    Summary = [Summary; table(string(variant), rain, system_residual, sum(coupling_error), gi_residual, ...
        sum(S.cumulative_rain_m3), sum(S.cumulative_capture_m3), sum(S.cumulative_infiltration_m3), ...
        sum(S.cumulative_overflow_m3), sum(S.cumulative_et_m3), sum(S.cumulative_native_exfiltration_m3), ...
        underdrain_surface, underdrain_external, gi_storage, L.cumulative_outlet_m3, ...
        replay_error, string(replay_note), numel(S.step_log), ...
        'VariableNames', {'variant', 'rain_m3', 'system_residual_m3', 'gi_coupling_error_m3', 'gi_residual_m3', ...
        'gi_rain_m3', 'capture_m3', 'media_infiltration_m3', 'overflow_m3', 'gi_et_m3', ...
        'native_exfiltration_m3', 'underdrain_to_surface_m3', 'underdrain_external_m3', ...
        'gi_final_storage_m3', 'outlet_m3', 'replay_max_error_m3', 'replay_note', 'steps'})]; %#ok<AGROW>
end

% Limiter case (review C14): Python's per-edge cap overdraws a compartment;
% MATLAB must match Python run with MATLAB's cap, and conserve.
[lim_error, lim_residual, lim_min_store, lim_python_as_is] = limiter_case_check(reference_dir);
Checks = [Checks; ...
    check_row("LIMITER", "MATLAB vs Python with the MATLAB cap [m3]", lim_error, 1e-9); ...
    check_row("LIMITER", "MATLAB GI mass residual [m3]", abs(lim_residual), 1e-9); ...
    check_row("LIMITER", "smallest MATLAB store [m3] (>= 0)", -lim_min_store, 1e-12)];
fprintf('Limiter case: unpatched Python kernel GI mass error %.4g m3 (MATLAB %.2g m3).\n', ...
    lim_python_as_is, lim_residual);

as_is_file = fullfile(reference_dir, 'python_limiter_as_is.csv');
if isfile(as_is_file)
    AsIs = readtable(as_is_file, 'TextType', 'string');
    writetable(AsIs, fullfile(summary_dir, 'Python_Limiter_As_Is.csv'));
    disp(AsIs);
end

passed = all(Checks.passed);
status = "FAIL";
if passed; status = "PASS"; end
PassFail = table("VAL-GI-001", status, true, height(Checks), sum(Checks.passed), ...
    max(abs(Summary.system_residual_m3)), max(Summary.replay_max_error_m3), ...
    'VariableNames', {'case_id', 'status', 'report_ready', 'checks', 'checks_passed', ...
    'max_system_residual_m3', 'max_replay_error_m3'});
writetable(Checks, fullfile(summary_dir, 'GI_Checks.csv'));
writetable(Summary, fullfile(summary_dir, 'GI_Summary.csv'));
writetable(PassFail, fullfile(summary_dir, 'Pass_Fail.csv'));
disp(Checks);
disp(Summary(:, {'variant', 'rain_m3', 'system_residual_m3', 'capture_m3', 'overflow_m3', ...
    'native_exfiltration_m3', 'underdrain_to_surface_m3', 'underdrain_external_m3', ...
    'replay_max_error_m3', 'steps'}));
disp(PassFail);
if ~passed
    error('VAL-GI-001 failed.');
end
fprintf('VAL-GI-001 passed.\n');

%% ------------------------------------------------------------------------
function R = run_variant(case_dir, variant, static_dir, etp_path)
% Runs the normal preprocessing and main loop in this function's workspace
% (they clear every variable outside their keep list; case_dir, summary_dir,
% static_dir and ValidationCase are on it).
ValidationCase = struct('variant', variant, 'record_steps', true, 'etp_path', etp_path); %#ok<NASGU>
summary_dir = fullfile(case_dir, 'Outputs', 'Validation'); %#ok<NASGU>
output_root = fullfile(case_dir, 'Outputs', ['Run_' variant]);
if exist(output_root, 'dir'); rmdir(output_root, 's'); end
Paths = make_output_paths(output_root); %#ok<NASGU>
InputPaths = make_input_paths(static_dir, etp_path); %#ok<NASGU>
input_data_bypass_script_path = fullfile(case_dir, 'Config', 'input_data_bypass_script.m'); %#ok<NASGU>
use_inputpaths_bypass = 1; %#ok<NASGU>
use_inputdata_bypass = 1; %#ok<NASGU>
clean_output_folder = true; %#ok<NASGU>
run_postprocessing = false; %#ok<NASGU>
enable_logging = false; %#ok<NASGU>
model_folder = ''; %#ok<NASGU>
resultsDir = Paths.Results; %#ok<NASGU>
export_root_dir = output_root; %#ok<NASGU>

HydroPol2D_preprocessing;
running_control.flag_system_mass_ledger = true;
HydroPol2D_Main_While;

R = struct();
R.variant = ValidationCase.variant;
R.G = Green_Infrastructure;
R.S = GI_States;
R.ledger = system_mass_ledger;
R.history = Green_Infrastructure_History;
R.max_edge_q = 0;
for i = 1:numel(GI_States.step_log)
    R.max_edge_q = max([R.max_edge_q; abs(GI_States.step_log{i}.edge_q_m3_s(:))]);
end
lai_input = double(GRIDobj(fullfile(case_dir, 'Outputs', 'Static', 'LAI.tif')).Z);
lai_now = double(gather(LAI_raster.Z));
valid = ~idx_nan;
R.lai_error = max(abs(lai_now(valid) - lai_input(valid)));
[gi_xc, gi_yc] = getcoordinates(DEM_raster);
R.dx = double(DEM_raster.cellsize);
R.dem_x_left = double(gi_xc(1)) - R.dx / 2;
R.dem_y_top = double(gi_yc(1)) + R.dx / 2;
R.elevation = double(gather(Elevation_Properties.elevation_cell));
R.valid = ~gather(idx_nan);
R.native_ksat_mm_h = R.G.native_ksat_m_s * 3.6e6;
R.gi_config = read_case_gi(ValidationCase, case_dir);
R.state_class = class(depths.d_t);
end

function Cfg = read_case_gi(ValidationCase, case_dir)
% The case's GI inputs, for gi_case.json (preprocessing clears them).
InputPaths = make_input_paths(fullfile(case_dir, 'Outputs', 'Static'), ValidationCase.etp_path); %#ok<NASGU>
run(fullfile(case_dir, 'Config', 'input_data_bypass_script.m'));
Cfg = InputData_Bypass.Green_Infrastructure;
end

function [max_error, residual, min_store, python_as_is] = limiter_case_check(reference_dir)
% The 4 x 5 case of make_python_reference.py: default designs, 8 m3 in every
% host cell at each step, 60 mm/h rain, 4 mm/day ET, 20 steps of 30 s.
max_error = inf; residual = inf; min_store = -inf; python_as_is = NaN;
json_file = fullfile(reference_dir, 'limiter_case.json');
csv_file = fullfile(reference_dir, 'python_limiter_case.csv');
if ~isfile(json_file) || ~isfile(csv_file)
    return
end
J = jsondecode(fileread(json_file));
dx = J.dx_m;
x = J.x_left_m + ((1:J.ncols) - 0.5) * dx;
y = J.y_top_m - ((1:J.nrows)' - 0.5) * dx;
Z = J.elevation_m;
DEM = GRIDobj(x, y, Z);
P = J.placements;
for f = 1:numel(P)
    d = P(f).underdrain_destination;
    if strcmp(d.kind, 'surface_cell')
        P(f).underdrain_destination = struct('kind', 'surface_cell', 'row', d.row + 1, 'column', d.column + 1);
    else
        P(f).underdrain_destination = struct('kind', d.kind);
    end
end
Cfg = struct('designs', J.designs, 'placements', P, 'courant', J.courant, ...
    'minimum_substep_s', J.minimum_substep_s, 'maximum_substeps', J.maximum_substeps);
G = hp2d_prepare_green_infrastructure(Cfg, DEM, Z, false(size(Z)), J.native_ksat_mm_h);
S = hp2d_green_infrastructure_initial_state(G);
n = G.compartment_count;
Ref = readtable(csv_file);
residual = 0; max_error = 0; min_store = inf;
for step = 1:J.steps
    before = sum(S.surface_m3) + sum(S.media_m3) + sum(S.gravel_m3);
    host = J.host_m3 * ones(numel(G.host_cells), 1);
    [S, F] = hp2d_green_infrastructure_step(G, S, host, J.rain_mm_h * ones(n, 1), ...
        J.et_mm_day * ones(n, 1), true, J.dt_s);
    after = sum(S.surface_m3) + sum(S.media_m3) + sum(S.gravel_m3);
    residual = residual + after - before - (sum(F.rain) + sum(F.capture) - sum(F.overflow) - ...
        sum(F.et) - sum(F.native_request) - sum(F.underdrain));
    min_store = min([min_store; S.surface_m3; S.media_m3; S.gravel_m3]);
    r = Ref(Ref.step == step, :);
    max_error = max([max_error; abs(r.surface_m3 - S.surface_m3); abs(r.media_m3 - S.media_m3); ...
        abs(r.gravel_m3 - S.gravel_m3)]);
end
as_is_file = fullfile(reference_dir, 'python_limiter_as_is.csv');
if isfile(as_is_file)
    A = readtable(as_is_file, 'TextType', 'string');
    python_as_is = A.python_as_is_gi_mass_residual_m3(A.variant == "limiter_case");
end
end

function T = replay_inputs_from_log(G, step_log)
% One row per accepted step: dt, host volumes before the GI step, gross rain
% and ET rates per compartment, and the fraction of the native exfiltration
% request the soil took.
n = numel(step_log);
nh = numel(G.host_cells);
nc = G.compartment_count;
M = zeros(n, 2 + nh + 3 * nc);
for i = 1:n
    s = step_log{i};
    M(i, :) = [s.t_min, s.dt_s, s.host_m3(:)', s.rain_rate_mm_h(:)', s.et_rate_mm_day(:)', ...
        s.native_accepted_fraction(:)'];
end
names = [{'t_min', 'dt_s'}, compose('host_%d_m3', 1:nh), compose('rain_%d_mm_h', 1:nc), ...
    compose('et_%d_mm_day', 1:nc), compose('native_fraction_%d', 1:nc)];
T = array2table(M, 'VariableNames', names);
end

function [max_error, note, geometry_error, mutated_error] = compare_with_python( ...
    reference_dir, variant, G, write_reference)
max_error = inf; geometry_error = inf; mutated_error = 0;
inputs_file = fullfile(reference_dir, ['replay_inputs_' variant '.csv']);
comp_file = fullfile(reference_dir, ['python_reference_' variant '_compartments.csv']);
host_file = fullfile(reference_dir, ['python_reference_' variant '_hosts.csv']);
geometry_file = fullfile(reference_dir, 'python_geometry.csv');
if ~isfile(inputs_file) || ~isfile(comp_file) || ~isfile(host_file) || ~isfile(geometry_file)
    note = 'Python reference missing: run make_python_reference.py';
    if write_reference
        note = 'replay inputs written; now run make_python_reference.py';
    end
    return
end
Geo = readtable(geometry_file, 'TextType', 'string');
cells = Geo(Geo.kind == "compartment", :);
edges = Geo(Geo.kind == "edge", :);
geometry_error = max([abs(cells.row - G.rows); abs(cells.column - G.columns); ...
    abs(cells.facility - G.facility_index); abs(cells.value - G.area_m2); ...
    abs(edges.row - G.edge_start); abs(edges.column - G.edge_end); abs(edges.value - G.edge_width_m)]);
if height(cells) ~= G.compartment_count || height(edges) ~= G.connection_count
    geometry_error = inf;
end

In = readtable(inputs_file);
C = readtable(comp_file);
H = readtable(host_file);
[Mine, MineHost] = replay(G, In, 1);
if height(Mine) ~= height(C) || height(MineHost) ~= height(H)
    note = 'replay and Python reference have different steps: regenerate the reference';
    return
end
max_error = max([max(abs(Mine{:, 3:end} - C{:, 3:end}), [], 'all'), ...
    max(abs(MineHost.host_after_m3 - H.host_after_m3)), ...
    max(abs(MineHost.cumulative_external_m3 - H.cumulative_external_m3))]);
[Mutated, ~] = replay(G, In, 0.9);
mutated_error = max(abs(Mutated{:, 3:end} - C{:, 3:end}), [], 'all');
note = sprintf('%d steps replayed, %d compared', height(In), numel(unique(C.step)));
end

function [T, Hosts] = replay(G, In, overflow_scale)
% hp2d_green_infrastructure_step on the recorded inputs, with the engine's
% native coupling (the soil takes the recorded fraction, the rest returns to
% the gravel). Rows every 10th step and the last, as the Python reference.
G.overflow_coefficient = G.overflow_coefficient * overflow_scale;
S = hp2d_green_infrastructure_initial_state(G);
nh = numel(G.host_cells);
nc = G.compartment_count;
n = height(In);
M = In{:, :};
keep = unique([10:10:n, n]);
names = {'surface_m3', 'media_m3', 'gravel_m3', 'cum_rain_m3', 'cum_capture_m3', ...
    'cum_infiltration_m3', 'cum_overflow_m3', 'cum_et_m3', 'cum_native_request_m3', ...
    'cum_native_m3', 'cum_underdrain_m3'};
cum = zeros(nc, 8);
external = 0;
rows = zeros(numel(keep) * nc, 2 + numel(names));
host_rows = zeros(numel(keep) * nh, 4);
r = 0; rh = 0;
for i = 1:n
    dt = M(i, 2);
    host = M(i, 3:2 + nh)';
    rain = M(i, 3 + nh:2 + nh + nc)';
    et = M(i, 3 + nh + nc:2 + nh + 2 * nc)';
    fraction = M(i, 3 + nh + 2 * nc:2 + nh + 3 * nc)';
    [S, F, host_after] = hp2d_green_infrastructure_step(G, S, host, rain, et, true, dt);
    native = F.native_request .* fraction;
    S.gravel_m3 = S.gravel_m3 + max(F.native_request - native, 0);
    cum = cum + [F.rain, F.capture, F.infiltration, F.overflow, F.et, F.native_request, native, F.underdrain];
    external = external + F.external_m3;
    if any(keep == i)
        for c = 1:nc
            r = r + 1;
            rows(r, :) = [i, c, S.surface_m3(c), S.media_m3(c), S.gravel_m3(c), cum(c, :)];
        end
        for h = 1:nh
            rh = rh + 1;
            host_rows(rh, :) = [i, h, host_after(h), external];
        end
    end
end
T = array2table(rows, 'VariableNames', [{'step', 'compartment'}, names]);
Hosts = array2table(host_rows, 'VariableNames', {'step', 'host', 'host_after_m3', 'cumulative_external_m3'});
end

function write_case_json(path, R)
% Everything make_python_reference.py needs to build the same case in Python.
G = R.G;
[ny, nx] = size(R.elevation);
[host_row, host_col] = ind2sub([ny nx], G.host_cells);
placements = cell(G.facility_count, 1);
cfg = hp2d_green_infrastructure_designs(R.gi_config.designs);
P = R.gi_config.placements;
for f = 1:G.facility_count
    p = P(f);
    destination = struct('kind', 'none', 'row', -1, 'column', -1);
    if G.target_kind(f) == 1
        [tr, tc] = ind2sub([ny nx], G.target_cell(f));
        destination = struct('kind', 'surface_cell', 'row', tr - 1, 'column', tc - 1);
    elseif G.target_kind(f) == 3
        destination = struct('kind', 'outside_domain', 'row', -1, 'column', -1);
    end
    placements{f} = struct('identifier', p.identifier, 'design_identifier', p.design_identifier, ...
        'vertices_m', p.vertices_m, 'underdrain_destination', destination);
end
J = struct();
J.note = 'Written by run_green_infrastructure_check.m (HYDROPOL_GI_WRITE_REPLAY=1). Rows and columns are 0-based.';
J.nrows = ny; J.ncols = nx; J.dx_m = R.dx; J.x_left_m = R.dem_x_left; J.y_top_m = R.dem_y_top;
J.elevation_m = R.elevation;
J.valid = R.valid;
native = zeros(ny, nx);
native(G.cell_index) = R.native_ksat_mm_h;
J.native_ksat_mm_h = native;
J.designs = cfg;
J.placements = placements;
J.courant = G.courant; J.minimum_substep_s = G.minimum_substep_s; J.maximum_substeps = G.maximum_substeps;
J.host_rows = host_row - 1; J.host_columns = host_col - 1;
fid = fopen(path, 'w');
fprintf(fid, '%s', jsonencode(J, 'PrettyPrint', true));
fclose(fid);
end

function Row = check_row(variant, name, value, tolerance, must_exceed)
if nargin < 5; must_exceed = false; end
if must_exceed
    ok = value > tolerance;
    rule = ">";
else
    ok = value <= tolerance;
    rule = "<=";
end
Row = table(string(variant), string(name), value, rule, tolerance, ok, ...
    'VariableNames', {'variant', 'check', 'value', 'rule', 'tolerance', 'passed'});
end

function write_static_rasters(static_dir)
% 10 x 7 cells of 10 m; x from 0 to 70 m, y from 0 to 100 m. The valley is
% column 4 (x = 35 m); slopes 3% across and 1% down to the outlet (row 10).
dx = 10; ny = 10; nx = 7;
xc = ((1:nx) - 0.5) * dx;
yc = 100 - ((1:ny)' - 0.5) * dx;
DEM = 1 + 0.01 * repmat(yc, 1, nx) + 0.03 * abs(repmat(xc, ny, 1) - 35);
R = maprefcells([0 nx * dx], [0 ny * dx], [ny nx], 'ColumnsStartFrom', 'north');
layers = struct('DEM', DEM, 'LULC', ones(ny, nx), 'SOIL', ones(ny, nx), ...
    'DTB', 3 * ones(ny, nx), 'GW_table', DEM - 1.5, 'Initial_SM', 60 * ones(ny, nx), ...
    'LAI', ones(ny, nx), 'Albedo', 0.23 * ones(ny, nx));
names = fieldnames(layers);
for i = 1:numel(names)
    % Double rasters keep the model states in double, so the GI step can be
    % checked to round-off (single rasters give single states).
    geotiffwrite(fullfile(static_dir, [names{i} '.tif']), double(layers.(names{i})), R, ...
        'CoordRefSysCode', 3857);
end
end

function InputPaths = make_input_paths(static_dir, etp_path)
InputPaths = struct();
InputPaths.DEM_path = fullfile(static_dir, 'DEM.tif');
InputPaths.LULC_path = fullfile(static_dir, 'LULC.tif');
InputPaths.SOIL_path = fullfile(static_dir, 'SOIL.tif');
InputPaths.DTB_path = fullfile(static_dir, 'DTB.tif');
InputPaths.GW_table_path = fullfile(static_dir, 'GW_table.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_dir, 'Initial_SM.tif');
InputPaths.LAI_path = fullfile(static_dir, 'LAI.tif');
InputPaths.Albedo_path = fullfile(static_dir, 'Albedo.tif');
InputPaths.Initial_SWE_path = '';
InputPaths.Initial_Snow_Depth_path = '';
InputPaths.ETP_input_spreadsheet = etp_path;
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
InputPaths.Inflow_Hydrograph_CSV = '';
InputPaths.Stage_Hydrograph_CSV = '';
InputPaths.Observed_Gauges_CSV = '';
InputPaths.Rainfall_Timeseries_File = '';
InputPaths.Outlet_Cells_CSV = '';
end

function Paths = make_output_paths(root_dir)
Paths = struct();
Paths.Root = root_dir;
Paths.Results = fullfile(root_dir, 'Modeling_Results');
Paths.Temp = fullfile(root_dir, 'Temporary_Files');
Paths.Logs = fullfile(root_dir, 'Logs');
Paths.FigPDF = fullfile(Paths.Results, 'Figures_PDF');
Paths.FigFIG = fullfile(Paths.Results, 'Figures_FIG');
Paths.Tables = fullfile(Paths.Results, 'Tables_CSV');
Paths.RastersWD = fullfile(Paths.Results, 'Rasters_Water_Depths');
Paths.RastersWSE = fullfile(Paths.Results, 'Rasters_WSE');
Paths.RastersVelocity = fullfile(Paths.Results, 'Rasters_Velocity');
Paths.RastersHazard = fullfile(Paths.Results, 'Rasters_Hazard');
Paths.RastersStatic = fullfile(Paths.Results, 'Rasters_Static');
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

function write_station_forcing(path)
% Hourly station records over the event: warm and dry (reference ET ~5 mm/day).
times = datetime(2025, 6, 1, 0, 0, 0) + minutes([0; 60; 120; 180]);
data = repmat([30 18 24 2 50 0], numel(times), 1);
cells = cell(3 + numel(times), 8);
cells(1,:) = {'Index','Date','Tmax_C','Tmin_C','Tavg_C','U2_m_s','RH_pct','G_MJ_m2_day'};
cells(2,:) = {NaN,NaT,NaN,NaN,NaN,0,NaN,0};
cells(3,:) = {NaN,NaT,NaN,NaN,NaN,NaN,NaN,NaN};
for i = 1:numel(times)
    cells{3+i,1} = i;
    cells{3+i,2} = times(i);
    cells(3+i,3:8) = num2cell(data(i,:));
end
if isfile(path); delete(path); end
writecell(cells, path);
end
