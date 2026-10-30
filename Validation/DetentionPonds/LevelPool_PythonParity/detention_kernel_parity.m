function K = detention_kernel_parity(case_dir)
%DETENTION_KERNEL_PARITY MATLAB pond kernels against HydroPol2D-Python on the same case.
%   Reads kernel_case.json (written by make_python_reference.py), prepares the two
%   ponds with hp2d_prepare_detention_ponds, runs the loop of the Python script with
%   hp2d_detention_pond_forcing, the pond-cell surface take-up and
%   hp2d_detention_pond_exchange, and compares every step with
%   python_reference_steps.csv. Also closes the water balance of the loop:
%   surface + ponds = initial + rain + exfiltration + poured - evaporation - seepage
%   - drained (the routing stand-in of the Python script).
%   K: struct of metrics; K.steps is the MATLAB step table.

C = jsondecode(fileread(fullfile(case_dir, 'kernel_case.json')));
ny = C.nrows; nx = C.ncols; dx = C.dx_m; dt = C.dt_s; A = dx * dx;
z = C.elevation_m;
x_centres = C.x_origin_m + ((1:nx) - 0.5) * dx;
y_centres = C.y_top_m - ((1:ny)' - 0.5) * dx;

% Python rows and columns are 0-based.
Ponds = C.ponds;
if ~iscell(Ponds), Ponds = num2cell(Ponds); end
for p = 1:numel(Ponds)
    devices = Ponds{p}.devices;
    if ~iscell(devices), devices = num2cell(devices); end
    for d = 1:numel(devices)
        for name = {'destination_row', 'destination_column', 'inlet_row', 'inlet_column'}
            if isfield(devices{d}, name{1}) && ~isempty(devices{d}.(name{1}))
                devices{d}.(name{1}) = devices{d}.(name{1}) + 1;
            end
        end
    end
    Ponds{p}.devices = devices;
end
P = hp2d_prepare_detention_ponds(Ponds, z, false(ny, nx), x_centres, y_centres, dx, dx);

% Geometry against Python.
G = readtable(fullfile(case_dir, 'python_reference_geometry.csv'), 'TextType', 'string');
geometry_error = 0;
stage0 = hp2d_detention_pond_stage(P.initial_storage_m3, P);
for p = 1:P.n_ponds
    faces = P.face_pond == p;
    geometry_error = max([geometry_error, ...
        abs(P.footprint_count(p) - G.cells(p)), abs(nnz(faces) - G.faces(p)), ...
        abs(sum(P.face_width_m(faces)) - G.face_width_sum_m(p)), ...
        abs(sum(P.face_sill_m(faces)) - G.face_sill_sum_m(p)), ...
        abs(P.initial_storage_m3(p) - G.initial_storage_m3(p)), ...
        abs(stage0(p) - G.initial_stage_m(p))]);
end

% The loop of make_python_reference.py.
F = C.forcing;
rain_pattern = C.rain_pattern;
evaporation_pattern = C.evaporation_pattern;
ordinary = ~P.mask;
surface = zeros(ny, nx);
surface(ordinary) = C.initial_surface_depth_m * A;
storage = P.initial_storage_m3;
pump = P.pump_active;
cumulative = P.cumulative_device_m3;
pour_cell = sub2ind([ny nx], F.pour_row + 1, F.pour_column + 1);
initial_water = sum(surface(:)) + sum(storage);
added = 0; removed = 0;
n_steps = C.n_steps;
nd = P.n_devices; np_ = P.n_ponds;
S = struct('storage', zeros(n_steps, np_), 'stage', zeros(n_steps, np_), ...
    'rainfall', zeros(n_steps, np_), 'evaporation', zeros(n_steps, np_), ...
    'seepage', zeros(n_steps, np_), 'recharge', zeros(n_steps, np_), ...
    'exfiltration', zeros(n_steps, np_), 'lateral', zeros(n_steps, np_), ...
    'q', zeros(n_steps, nd), 'pump', zeros(n_steps, nd), 'cumulative', zeros(n_steps, nd), ...
    'surface', zeros(n_steps, 1), 'reverse', zeros(n_steps, nd));
for step = 1:n_steps
    rain = F.rain_mm_h(step) * rain_pattern;
    evaporation = F.evaporation_mm_day(step) * evaporation_pattern;
    R = hp2d_detention_pond_forcing(storage, rain(P.cells), evaporation(P.cells), P, dt);
    storage = R.storage_m3;
    ordinary_rain = zeros(ny, nx);
    ordinary_rain(ordinary) = rain(ordinary) / 1000 / 3600 * A * dt;
    surface = surface + ordinary_rain;

    % Groundwater exfiltration: the pond cells' share fills the ponds first.
    exfiltration_mm = F.exfiltration_mm(step) * ones(ny, nx);
    exfiltration_into_ponds = accumarray(P.cell_pond, max(exfiltration_mm(P.cells), 0) / 1000 * A, [np_ 1]);
    storage = storage + exfiltration_into_ponds;
    surface(ordinary) = surface(ordinary) + exfiltration_mm(ordinary) / 1000 * A;
    surface(pour_cell) = surface(pour_cell) + F.pour_m3(step);

    E = hp2d_detention_pond_exchange(surface(P.exchange_cells), storage, pump, cumulative, P, dt);
    surface(P.exchange_cells) = E.surface_m3;
    % Stand-in for routing (make_python_reference.py): drain the ordinary cells.
    drained = zeros(ny, nx);
    drained(ordinary) = surface(ordinary) * C.drain_fraction;
    surface = surface - drained;
    storage = E.storage_m3;
    pump = E.pump_active;
    cumulative = E.cumulative_m3;

    added = added + sum(R.rainfall_m3) + sum(ordinary_rain(:)) + sum(exfiltration_into_ponds) + ...
        sum(exfiltration_mm(ordinary)) / 1000 * A + F.pour_m3(step);
    removed = removed + sum(R.evaporation_m3) + sum(R.seepage_m3) + sum(E.external_m3) + sum(drained(:));
    recharge = accumarray(P.cell_pond, R.recharge_rate_m_s, [np_ 1]) * A * dt;
    S.storage(step, :) = storage';
    S.stage(step, :) = hp2d_detention_pond_stage(storage, P)';
    S.rainfall(step, :) = R.rainfall_m3';
    S.evaporation(step, :) = R.evaporation_m3';
    S.seepage(step, :) = R.seepage_m3';
    S.recharge(step, :) = recharge';
    S.exfiltration(step, :) = exfiltration_into_ponds';
    S.lateral(step, :) = E.lateral_inflow_m3';
    S.q(step, :) = E.device_discharge_m3_s';
    S.pump(step, :) = pump';
    S.cumulative(step, :) = cumulative';
    S.surface(step) = sum(surface(:));
end
final_water = sum(surface(:)) + sum(storage);
K.loop_water_balance_residual_m3 = final_water - initial_water - added + removed;
K.loop_volume_scale_m3 = added;
K.seepage_recharge_mismatch_m3 = max(abs(S.seepage(:) - S.recharge(:)));

% Step-by-step comparison with Python.
Ref = readtable(fullfile(case_dir, 'python_reference_steps.csv'), 'VariableNamingRule', 'preserve');
pond_fields = {'storage', 'storage_m3'; 'stage', 'stage_m'; 'rainfall', 'rainfall_m3'; ...
    'evaporation', 'evaporation_m3'; 'seepage', 'seepage_m3'; 'recharge', 'recharge_m3'; ...
    'exfiltration', 'exfiltration_m3'; 'lateral', 'lateral_m3'};
max_volume_error = 0; max_stage_error = 0; max_q_error = 0; pump_mismatch = 0;
for p = 1:np_
    for f = 1:size(pond_fields, 1)
        ref = Ref.(sprintf('%s_%s', P.identifier{p}, pond_fields{f, 2}));
        err = max(abs(S.(pond_fields{f, 1})(:, p) - ref));
        if strcmp(pond_fields{f, 1}, 'stage')
            max_stage_error = max(max_stage_error, err);
        else
            max_volume_error = max(max_volume_error, err);
        end
    end
end
for d = 1:nd
    name = strrep(P.device_identifier{d}, '/', '__');
    max_q_error = max(max_q_error, max(abs(S.q(:, d) - Ref.([name '_q_m3_s']))));
    max_volume_error = max(max_volume_error, max(abs(S.cumulative(:, d) - Ref.([name '_cumulative_m3']))));
    if P.dev_kind(d) == 2
        pump_mismatch = pump_mismatch + nnz(S.pump(:, d) ~= Ref.([name '_pump_on']));
    end
end
max_volume_error = max(max_volume_error, max(abs(S.surface - Ref.surface_m3)));

K.python_commit = string(C.python_commit);
K.geometry_max_error = geometry_error;
K.max_volume_error_m3 = max_volume_error;
K.max_stage_error_m = max_stage_error;
K.max_discharge_error_m3_s = max_q_error;
K.pump_switch_mismatches = pump_mismatch;
K.peak_discharge_m3_s = max(abs(S.q(:)));
K.total_lateral_inflow_m3 = sum(S.lateral(:));
K.total_seepage_m3 = sum(S.seepage(:));
K.pump_on_steps = sum(S.pump(:, P.dev_kind == 2), 'all');
K.pump_switches = sum(abs(diff(S.pump(:, P.dev_kind == 2))), 'all');
K.reverse_flow_steps = nnz(S.q(:, P.dev_kind < 2) < 0);
K.spillway_flow_steps = nnz(S.q(:, P.dev_kind == 1) > 0);
K.steps = S;
K.prepared = P;
end
