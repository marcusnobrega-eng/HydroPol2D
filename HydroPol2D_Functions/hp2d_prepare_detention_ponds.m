function P = hp2d_prepare_detention_ponds(Ponds, elevation, idx_nan, x_centres, y_centres, dx, dy)
%HP2D_PREPARE_DETENTION_PONDS Footprints, level-pool geometry and devices of the ponds.
%   Ponds: struct array (InputData_Bypass.Detention_Ponds), one element per pond:
%     identifier          name (unique)
%     footprint, one of:
%       vertices_m        [x y] polygon vertices in the DEM's projected CRS; a cell
%                         belongs to the pond when its centre is inside (rasterio
%                         geometry_mask with all_touched=False, as in Python)
%       mask              logical matrix of the DEM's size
%       cells             [row col] list (1-based)
%       raster_path, raster_value   raster of pond ids; cells equal to the value
%     storage_method      'dtm' (default: stage-storage from the DEM cells of the
%                         footprint) or 'power_law' (V = C (stage - invert)^b)
%     initial_stage_m     water level at t = 0 [m]
%     invert_elevation_m, storage_coefficient, storage_exponent   power law only
%     lined               true (default): no bottom seepage
%     bottom_ksat_mm_h    seepage conductivity of an unlined pond [mm/h] (> 0)
%     devices             struct array, each with:
%       kind              'orifice', 'spillway' or 'pump'
%       identifier        name (unique within the pond)
%       destination_row, destination_column   1-based destination cell, or
%       destination_x, destination_y          its map coordinates; (0, 0) sends the
%                                             discharge out of the model (booked as
%                                             outlet flow), as for Reservoir_Data
%       inlet_row, inlet_column (optional)    a cell inside the pond
%       orifice:  invert_elevation_m, opening_area_m2, discharge_coefficient (0.62)
%       spillway: crest_elevation_m, crest_length_m, discharge_coefficient (1.7)
%       pump:     maximum_flow_m3_s, shutoff_head_m, additional_head_loss_m (0),
%                 start_stage_m, stop_stage_m
%       enabled           true (default)
%   elevation: routing DEM [m] (NaN outside the domain); idx_nan: outside cells.
%   x_centres (1 x nx), y_centres (ny x 1): cell-centre coordinates; dx, dy [m].
%
%   P holds the fixed geometry (cells, sorted elevations, critical faces, devices)
%   and the pond state (storage, pump switches, cumulative volumes). This is
%   prepare_detention_ponds in HydroPol2D-Python (routing/detention.py), with the
%   same checks; the (0, 0) destination and the mask/cells/raster footprints are
%   MATLAB additions.

if isempty(Ponds)
    error('HydroPol2D:DetentionPonds', ...
        'flag_detention_ponds = 1 requires at least one pond in InputData_Bypass.Detention_Ponds.');
end
[ny, nx] = size(elevation);
domain = ~idx_nan & isfinite(elevation);
cell_area = dx * dy;
n_ponds = numel(Ponds);

pond_index = zeros(ny, nx);
cells = cell(n_ponds, 1);
identifiers = cell(n_ponds, 1);
dtm = false(n_ponds, 1);
invert = nan(n_ponds, 1);
coefficient = nan(n_ponds, 1);
exponent = nan(n_ponds, 1);
ksat_m_s = zeros(n_ponds, 1);
initial_stage = zeros(n_ponds, 1);

for p = 1:n_ponds
    pond = element(Ponds, p);
    identifiers{p} = char(string(field_or(pond, 'identifier', sprintf('pond_%d', p))));
    mask = pond_footprint(pond, identifiers{p}, ny, nx, x_centres, y_centres);
    mask = mask & domain;
    if ~any(mask(:))
        error('HydroPol2D:DetentionPonds', ...
            'Detention pond "%s" contains no active cell centres.', identifiers{p});
    end
    if any(mask(:) & pond_index(:) > 0)
        error('HydroPol2D:DetentionPonds', 'Detention-pond footprints cannot overlap.');
    end
    pond_index(mask) = p;
    cells{p} = find(mask);

    method = lower(char(string(field_or(pond, 'storage_method', 'dtm'))));
    initial_stage(p) = required_number(pond, 'initial_stage_m', identifiers{p});
    switch method
        case 'dtm'
            dtm(p) = true;
            if initial_stage(p) < min(elevation(cells{p}))
                error('HydroPol2D:DetentionPonds', ...
                    'Detention pond "%s" has an initial stage below its DTM invert (%g m).', ...
                    identifiers{p}, min(elevation(cells{p})));
            end
        case {'power_law', 'power'}
            invert(p) = required_number(pond, 'invert_elevation_m', identifiers{p});
            coefficient(p) = required_number(pond, 'storage_coefficient', identifiers{p});
            exponent(p) = required_number(pond, 'storage_exponent', identifiers{p});
            if ~(coefficient(p) > 0) || ~(exponent(p) >= 1)
                error('HydroPol2D:DetentionPonds', ...
                    'Pond "%s": storage_coefficient must be > 0 and storage_exponent >= 1.', identifiers{p});
            end
            if initial_stage(p) < invert(p)
                error('HydroPol2D:DetentionPonds', ...
                    'Pond "%s": initial_stage_m cannot be below its invert elevation.', identifiers{p});
            end
        otherwise
            error('HydroPol2D:DetentionPonds', ...
                'Pond "%s": storage_method must be ''dtm'' or ''power_law''.', identifiers{p});
    end

    lined = logical(field_or(pond, 'lined', true));
    ksat_mm_h = double(field_or(pond, 'bottom_ksat_mm_h', 0));
    if lined && ksat_mm_h ~= 0
        error('HydroPol2D:DetentionPonds', 'A lined pond ("%s") must have bottom_ksat_mm_h = 0.', identifiers{p});
    end
    if ~lined && ~(ksat_mm_h > 0)
        error('HydroPol2D:DetentionPonds', 'An unlined pond ("%s") requires a positive bottom_ksat_mm_h.', identifiers{p});
    end
    ksat_m_s(p) = ksat_mm_h / 1000 / 3600;
end
pond_mask = pond_index > 0;

% Pond cells, pond by pond.
P = struct();
P.n_ponds = n_ponds;
P.identifier = identifiers;
P.mask = pond_mask;
P.index = pond_index;
P.cells = vertcat(cells{:});
P.cell_pond = pond_index(P.cells);
P.cell_elevation_m = double(elevation(P.cells));
P.footprint_count = accumarray(P.cell_pond, 1, [n_ponds 1]);
P.dtm = dtm;
P.invert_elevation_m = invert;
P.storage_coefficient = coefficient;
P.storage_exponent = exponent;
P.bottom_ksat_m_s = ksat_m_s;
P.cell_area_m2 = cell_area;
P.sorted_elevation_m = cell(n_ponds, 1);
P.cumulative_elevation_m = cell(n_ponds, 1);
for p = 1:n_ponds
    z = sort(double(elevation(cells{p})));
    P.sorted_elevation_m{p} = z(:);
    P.cumulative_elevation_m{p} = cumsum(z(:));
end
P.initial_storage_m3 = zeros(n_ponds, 1);
for p = 1:n_ponds
    if dtm(p)
        P.initial_storage_m3(p) = sum(max(initial_stage(p) - double(elevation(cells{p})), 0)) * cell_area;
    else
        P.initial_storage_m3(p) = coefficient(p) * max(initial_stage(p) - invert(p), 0)^exponent(p);
    end
end

% Critical faces: every face between a pond cell and an active cell outside all
% ponds drains that cell into the pond (critical flow over the higher bed).
face_pond = []; face_cell = []; face_width = []; face_sill = [];
for p = 1:n_ponds
    [rows, cols] = ind2sub([ny nx], cells{p});
    for c = 1:numel(rows)
        for delta = [-1 0; 1 0; 0 -1; 0 1]'
            r = rows(c) + delta(1);
            q = cols(c) + delta(2);
            if r < 1 || q < 1 || r > ny || q > nx || ~domain(r, q) || pond_mask(r, q)
                continue
            end
            face_pond(end + 1, 1) = p; %#ok<AGROW>
            face_cell(end + 1, 1) = sub2ind([ny nx], r, q); %#ok<AGROW>
            if delta(1) ~= 0
                face_width(end + 1, 1) = dx; %#ok<AGROW>
            else
                face_width(end + 1, 1) = dy; %#ok<AGROW>
            end
            face_sill(end + 1, 1) = max(double(elevation(r, q)), double(elevation(rows(c), cols(c)))); %#ok<AGROW>
        end
    end
end

% Devices, in pond order.
dev = struct('pond', [], 'kind', [], 'enabled', [], 'level', [], 'coefficient', [], ...
    'geometry', [], 'start', [], 'stop', [], 'qmax', [], 'loss', [], 'shutoff', [], ...
    'destination', [], 'book_cell', []);
dev_identifier = {};
for p = 1:n_ponds
    devices = field_or(element(Ponds, p), 'devices', []);
    names = {};
    for d = 1:numel(devices)
        D = element(devices, d);
        name = char(string(field_or(D, 'identifier', sprintf('device_%d', d))));
        if any(strcmpi(names, name))
            error('HydroPol2D:DetentionPonds', ...
                'Hydraulic-device identifiers must be unique within each pond ("%s").', identifiers{p});
        end
        names{end + 1} = name; %#ok<AGROW>
        label = sprintf('%s/%s', identifiers{p}, name);
        kind = lower(char(string(required_field(D, 'kind', label))));
        switch kind
            case 'orifice'
                k = 0;
                level = required_number(D, 'invert_elevation_m', label);
                geometry = positive_number(D, 'opening_area_m2', label);
                cd = double(field_or(D, 'discharge_coefficient', 0.62));
            case 'spillway'
                k = 1;
                level = required_number(D, 'crest_elevation_m', label);
                geometry = positive_number(D, 'crest_length_m', label);
                cd = double(field_or(D, 'discharge_coefficient', 1.7));
            case 'pump'
                k = 2;
                level = 0;
                geometry = 0;
                cd = 0;
            otherwise
                error('HydroPol2D:DetentionPonds', ...
                    'Device "%s": kind must be orifice, spillway or pump.', label);
        end
        if k < 2 && ~(cd > 0)
            error('HydroPol2D:DetentionPonds', 'Device "%s": discharge_coefficient must be > 0.', label);
        end
        start_stage = 0; stop_stage = 0; qmax = 0; loss = 0; shutoff = 1;
        if k == 2
            qmax = positive_number(D, 'maximum_flow_m3_s', label);
            shutoff = positive_number(D, 'shutoff_head_m', label);
            loss = double(field_or(D, 'additional_head_loss_m', 0));
            start_stage = required_number(D, 'start_stage_m', label);
            stop_stage = required_number(D, 'stop_stage_m', label);
            if loss < 0
                error('HydroPol2D:DetentionPonds', 'Pump "%s": additional_head_loss_m must be >= 0.', label);
            end
            if stop_stage > start_stage
                error('HydroPol2D:DetentionPonds', 'Pump "%s": stop_stage_m cannot exceed start_stage_m.', label);
            end
        end

        [row, col] = device_destination(D, label, x_centres, y_centres);
        if row == 0 && col == 0
            destination = 0;   % leaves the model
        else
            if row < 1 || col < 1 || row > ny || col > nx || ~domain(row, col)
                error('HydroPol2D:DetentionPonds', ...
                    'Device "%s" must connect to an active surface cell.', label);
            end
            if pond_mask(row, col)
                error('HydroPol2D:DetentionPonds', ...
                    'Device "%s" cannot discharge to a detention-pond cell.', label);
            end
            destination = sub2ind([ny nx], row, col);
        end
        inlet_row = field_or(D, 'inlet_row', []);
        inlet_col = field_or(D, 'inlet_column', []);
        if xor(isempty(inlet_row), isempty(inlet_col))
            error('HydroPol2D:DetentionPonds', ...
                'Device "%s": an inlet requires both inlet_row and inlet_column.', label);
        end
        book_cell = cells{p}(1);
        if ~isempty(inlet_row)
            if inlet_row < 1 || inlet_col < 1 || inlet_row > ny || inlet_col > nx || ...
                    pond_index(inlet_row, inlet_col) ~= p
                error('HydroPol2D:DetentionPonds', ...
                    'Inlet for device "%s" must lie inside pond "%s".', label, identifiers{p});
            end
            book_cell = sub2ind([ny nx], inlet_row, inlet_col);
        end

        dev.pond(end + 1, 1) = p;
        dev.kind(end + 1, 1) = k;
        dev.enabled(end + 1, 1) = logical(field_or(D, 'enabled', true));
        dev.level(end + 1, 1) = level;
        dev.coefficient(end + 1, 1) = cd;
        dev.geometry(end + 1, 1) = geometry;
        dev.start(end + 1, 1) = start_stage;
        dev.stop(end + 1, 1) = stop_stage;
        dev.qmax(end + 1, 1) = qmax;
        dev.loss(end + 1, 1) = loss;
        dev.shutoff(end + 1, 1) = shutoff;
        dev.destination(end + 1, 1) = destination;
        dev.book_cell(end + 1, 1) = book_cell;
        dev_identifier{end + 1, 1} = label; %#ok<AGROW>
    end
end

% The surface cells the exchange reads and writes: face donors and destinations.
P.exchange_cells = unique([face_cell(:); dev.destination(dev.destination > 0)]);
P.exchange_terrain_m = double(elevation(P.exchange_cells));
[~, P.face_local] = ismember(face_cell(:), P.exchange_cells);
P.face_pond = face_pond(:);
P.face_cell = face_cell(:);
P.face_width_m = face_width(:);
P.face_sill_m = face_sill(:);
P.n_devices = numel(dev.pond);
P.device_identifier = dev_identifier;
P.dev_pond = dev.pond(:);
P.dev_kind = dev.kind(:);
P.dev_enabled = logical(dev.enabled(:));
P.dev_level_m = dev.level(:);
P.dev_coefficient = dev.coefficient(:);
P.dev_geometry = dev.geometry(:);
P.dev_start_stage_m = dev.start(:);
P.dev_stop_stage_m = dev.stop(:);
P.dev_maximum_flow_m3_s = dev.qmax(:);
P.dev_head_loss_m = dev.loss(:);
P.dev_shutoff_head_m = dev.shutoff(:);
P.dev_destination = dev.destination(:);
P.dev_book_cell = dev.book_cell(:);
P.dev_local = zeros(P.n_devices, 1);
inside = P.dev_destination > 0;
[~, P.dev_local(inside)] = ismember(P.dev_destination(inside), P.exchange_cells);

% Routing: a pond cell is not surface. A wall above any water level closes its
% faces in the local-inertial solver (the routing mask of HydroPol2D-Python).
P.wall_elevation_m = max(double(elevation(domain))) + 1000;

% State.
P.storage_m3 = P.initial_storage_m3;
P.pump_active = false(P.n_devices, 1);
P.cumulative_device_m3 = zeros(P.n_devices, 1);
P.device_discharge_m3_s = zeros(P.n_devices, 1);
P.cumulative_lateral_m3 = zeros(n_ponds, 1);
P.cumulative_rainfall_m3 = zeros(n_ponds, 1);
P.cumulative_evaporation_m3 = zeros(n_ponds, 1);
P.cumulative_seepage_m3 = zeros(n_ponds, 1);
P.cumulative_seepage_out_m3 = zeros(n_ponds, 1);
P.cumulative_surface_in_m3 = zeros(n_ponds, 1);
P.cumulative_external_m3 = zeros(n_ponds, 1);
P.recharge_rate_m_s = zeros(ny, nx);
P.last_evaporation_m3 = 0;
P.last_seepage_out_m3 = 0;
P.history = struct('time_min', zeros(0, 1), 'stage_m', zeros(0, n_ponds), ...
    'storage_m3', zeros(0, n_ponds), 'device_discharge_m3_s', zeros(0, P.n_devices), ...
    'device_cumulative_m3', zeros(0, P.n_devices), 'lateral_inflow_m3', zeros(0, n_ponds));
end

function mask = pond_footprint(pond, name, ny, nx, x_centres, y_centres)
if has_value(pond, 'vertices_m')
    v = double(pond.vertices_m);
    if size(v, 2) ~= 2 || size(v, 1) < 3
        error('HydroPol2D:DetentionPonds', 'Pond "%s": vertices_m must be an N x 2 [x y] list (N >= 3).', name);
    end
    if polyarea(v(:, 1), v(:, 2)) <= 0
        error('HydroPol2D:DetentionPonds', 'Detention pond "%s" has an invalid polygon.', name);
    end
    [X, Y] = meshgrid(double(x_centres(:)'), double(y_centres(:)));
    mask = inpolygon(X, Y, v(:, 1), v(:, 2));
elseif has_value(pond, 'mask')
    mask = logical(pond.mask);
    if ~isequal(size(mask), [ny nx])
        error('HydroPol2D:DetentionPonds', 'Pond "%s": mask must have the size of the DEM.', name);
    end
elseif has_value(pond, 'cells')
    rc = double(pond.cells);
    mask = false(ny, nx);
    mask(sub2ind([ny nx], rc(:, 1), rc(:, 2))) = true;
elseif has_value(pond, 'raster_path')
    R = GRIDobj(char(pond.raster_path));
    if ~isequal(size(R.Z), [ny nx])
        error('HydroPol2D:DetentionPonds', 'Pond "%s": raster_path must be aligned with the DEM.', name);
    end
    mask = double(R.Z) == double(field_or(pond, 'raster_value', 1));
else
    error('HydroPol2D:DetentionPonds', ...
        'Pond "%s" needs a footprint: vertices_m, mask, cells or raster_path.', name);
end
end

function [row, col] = device_destination(D, label, x_centres, y_centres)
if has_value(D, 'destination_row') && has_value(D, 'destination_column')
    row = double(D.destination_row);
    col = double(D.destination_column);
elseif has_value(D, 'destination_x') && has_value(D, 'destination_y')
    if D.destination_x == 0 && D.destination_y == 0
        row = 0; col = 0;
    else
        col = round((double(D.destination_x) - x_centres(1)) / (x_centres(2) - x_centres(1))) + 1;
        row = round((double(D.destination_y) - y_centres(1)) / (y_centres(2) - y_centres(1))) + 1;
    end
else
    error('HydroPol2D:DetentionPonds', ...
        'Device "%s" needs destination_row/destination_column or destination_x/destination_y.', label);
end
end

function S = element(list, i)
% Ponds and devices may come as a struct array or as a cell array of structs.
if iscell(list)
    S = list{i};
else
    S = list(i);
end
end

function tf = has_value(S, name)
tf = isfield(S, name) && ~isempty(S.(name));
end

function value = field_or(S, name, default)
if has_value(S, name)
    value = S.(name);
else
    value = default;
end
end

function value = required_field(S, name, label)
if ~has_value(S, name)
    error('HydroPol2D:DetentionPonds', '"%s" requires %s.', label, name);
end
value = S.(name);
end

function value = required_number(S, name, label)
value = double(required_field(S, name, label));
if ~isscalar(value) || ~isfinite(value)
    error('HydroPol2D:DetentionPonds', '"%s": %s must be a finite number.', label, name);
end
end

function value = positive_number(S, name, label)
value = required_number(S, name, label);
if ~(value > 0)
    error('HydroPol2D:DetentionPonds', '"%s": %s must be positive.', label, name);
end
end
