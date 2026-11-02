function G = hp2d_prepare_green_infrastructure(Config, DEM_raster, elevation, idx_nan, ...
    native_ksat_mm_h)
%HP2D_PREPARE_GREEN_INFRASTRUCTURE Intersect GI polygons with the grid and build the compartments.
%   Port of prepare_green_infrastructure in HydroPol2D-Python
%   (modules/green_infrastructure.py). Each enabled placement is a polygon in
%   the map coordinates of the model DEM; its part inside each active cell
%   is one compartment, ordered by facility, then row, then column. Two
%   compartments of the same facility in D4-adjacent cells are connected by
%   an edge whose width is the length of the shared cell side inside the
%   polygon (or min(sqrt(A_i), sqrt(A_j)) when that length is 0).
%
%   Config (InputData_Bypass.Green_Infrastructure):
%     designs     design library (hp2d_green_infrastructure_designs); empty
%                 or absent: the three Python templates.
%     placements  struct array (or cell array of structs), one per facility:
%       identifier, design_identifier, vertices_m ([x y] rows, >= 3, in the
%       DEM's map coordinates), enabled (default true), and
%       underdrain_destination: empty (no underdrain), or a struct with kind
%       'surface_cell' (x_m and y_m of a point in the cell, or row and column
%       of the model grid, 1-based) or 'outside_domain' (the water leaves the
%       model and is booked as outlet flow). 'urban_drainage_node' is refused:
%       MATLAB HydroPol2D has no urban drainage network.
%     courant (0.5), minimum_substep_s (0.01), maximum_substeps (10000):
%       substeps of the lateral flow between compartments.
%   DEM_raster: model GRIDobj (cell centres from getcoordinates).
%   elevation: model terrain [m]; idx_nan: inactive cells.
%   native_ksat_mm_h: saturated conductivity of the soil below the GI [mm/h]
%     (Python uses the transmission-zone Ks).

if ~isstruct(Config) || ~isfield(Config, 'placements')
    error('flag_green_infrastructure = 1 requires InputData_Bypass.Green_Infrastructure.placements.');
end
designs = hp2d_green_infrastructure_designs(get_field(Config, 'designs', []));
design_ids = lower({designs.identifier});
placements = as_cell(Config.placements);
enabled = cellfun(@(p) logical(get_field(p, 'enabled', true)), placements);
placements = placements(enabled);
if isempty(placements)
    error('Green infrastructure is enabled but has no enabled placements.');
end
ids = cellfun(@(p) lower(char(p.identifier)), placements, 'UniformOutput', false);
if numel(unique(ids)) ~= numel(ids)
    error('GI placement identifiers must be unique.');
end

[ny, nx] = size(elevation);
valid = ~idx_nan & isfinite(elevation);
dx = double(DEM_raster.cellsize);
cell_area = dx^2;
[xc, yc] = getcoordinates(DEM_raster);
xc = double(xc(:));
yc = double(yc(:));
x_left = xc - dx / 2;
x_right = xc + dx / 2;
y_top = yc + dx / 2;
y_bottom = yc - dx / 2;

warning_state = warning('off', 'MATLAB:polyshape:repairedBySimplify');
restore_warning = onCleanup(@() warning(warning_state));

n_fac = numel(placements);
polygons = cell(n_fac, 1);
for f = 1:n_fac
    p = placements{f};
    vertices = double(p.vertices_m);
    if size(vertices, 2) ~= 2 || size(vertices, 1) < 3
        error('GI placement %s: vertices_m must hold at least three [x y] rows.', char(p.identifier));
    end
    polygon = polyshape(vertices(:, 1), vertices(:, 2));
    if polygon.NumRegions == 0 || area(polygon) <= 1.0e-12
        error('GI polygon %s has no valid area.', char(p.identifier));
    end
    for other = 1:f - 1
        if area(intersect(polygon, polygons{other})) > 1.0e-8
            error('GI polygons %s and %s overlap.', char(placements{other}.identifier), ...
                char(p.identifier));
        end
    end
    polygons{f} = polygon;
end

rows = []; cols = []; facility = []; piece_area = []; design_index = [];
fraction_map = zeros(ny, nx);
for f = 1:n_fac
    p = placements{f};
    d = find(strcmp(design_ids, lower(char(p.design_identifier))), 1);
    if isempty(d)
        error('GI placement %s references the unavailable design %s.', char(p.identifier), ...
            char(p.design_identifier));
    end
    [bx, by] = boundingbox(polygons{f});
    col_range = find(x_right > bx(1) & x_left < bx(2));
    row_range = find(y_top > by(1) & y_bottom < by(2));
    covered = 0;
    for r = row_range(:)'
        for c = col_range(:)'
            cell_polygon = polyshape([x_left(c) x_right(c) x_right(c) x_left(c)], ...
                [y_top(r) y_top(r) y_bottom(r) y_bottom(r)]);
            a = area(intersect(polygons{f}, cell_polygon));
            if a <= 1.0e-8 || ~valid(r, c)
                continue
            end
            rows(end + 1, 1) = r; %#ok<AGROW>
            cols(end + 1, 1) = c; %#ok<AGROW>
            facility(end + 1, 1) = f; %#ok<AGROW>
            piece_area(end + 1, 1) = a; %#ok<AGROW>
            design_index(end + 1, 1) = d; %#ok<AGROW>
            fraction_map(r, c) = fraction_map(r, c) + a / cell_area;
            covered = covered + a;
        end
    end
    if covered < area(polygons{f}) * (1 - 1.0e-6)
        error('GI polygon %s extends outside the active domain.', char(p.identifier));
    end
end
if any(fraction_map(:) > 1 + 1.0e-8)
    error('Combined GI area exceeds one full cell in at least one cell.');
end

% D4 connections inside each facility: the east, then the south neighbour.
n = numel(rows);
key = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:n
    key(sprintf('%d_%d_%d', facility(i), rows(i), cols(i))) = i;
end
edge_start = []; edge_end = []; edge_width = []; edge_length = [];
for i = 1:n
    neighbours = [0 1; 1 0];
    for k = 1:2
        name = sprintf('%d_%d_%d', facility(i), rows(i) + neighbours(k, 1), cols(i) + neighbours(k, 2));
        if ~isKey(key, name)
            continue
        end
        j = key(name);
        r = rows(i); c = cols(i);
        if k == 1   % shared side x = x_right(c)
            segment = [x_right(c) y_bottom(r); x_right(c) y_top(r)];
        else        % shared side y = y_bottom(r)
            segment = [x_left(c) y_bottom(r); x_right(c) y_bottom(r)];
        end
        width = inside_length(polygons{facility(i)}, segment);
        if width <= 1.0e-8
            width = min(sqrt(piece_area(i)), sqrt(piece_area(j)));
        end
        edge_start(end + 1, 1) = i; %#ok<AGROW>
        edge_end(end + 1, 1) = j; %#ok<AGROW>
        edge_width(end + 1, 1) = width; %#ok<AGROW>
        edge_length(end + 1, 1) = dx; %#ok<AGROW>
    end
end

% Underdrain destinations, one per facility.
target_kind = zeros(n_fac, 1);
target_cell = zeros(n_fac, 1);
for f = 1:n_fac
    destination = get_field(placements{f}, 'underdrain_destination', []);
    if isempty(destination) || (isstruct(destination) && strcmpi(char(get_field(destination, 'kind', 'none')), 'none'))
        continue
    end
    kind = lower(char(destination.kind));
    switch kind
        case 'surface_cell'
            if isfield(destination, 'x_m') && ~isempty(destination.x_m)
                [~, c] = min(abs(xc - double(destination.x_m)));
                [~, r] = min(abs(yc - double(destination.y_m)));
                inside = double(destination.x_m) >= x_left(c) && double(destination.x_m) <= x_right(c) && ...
                    double(destination.y_m) >= y_bottom(r) && double(destination.y_m) <= y_top(r);
            else
                r = double(destination.row);
                c = double(destination.column);
                inside = r >= 1 && r <= ny && c >= 1 && c <= nx && r == round(r) && c == round(c);
            end
            if ~inside || ~valid(r, c)
                error('GI underdrain for %s targets an inactive cell.', char(placements{f}.identifier));
            end
            target_kind(f) = 1;
            target_cell(f) = sub2ind([ny nx], r, c);
        case 'outside_domain'
            target_kind(f) = 3;
        case 'urban_drainage_node'
            error(['GI underdrain for %s targets an urban drainage node: MATLAB HydroPol2D has ', ...
                'no urban drainage network. Use surface_cell or outside_domain.'], ...
                char(placements{f}.identifier));
        otherwise
            error('GI underdrain for %s: unknown destination kind %s.', ...
                char(placements{f}.identifier), kind);
    end
end

cell_index = sub2ind([ny nx], rows, cols);
[host_cells, ~, compartment_host] = unique([cell_index; target_cell(target_cell > 0)], 'stable');
compartment_host = compartment_host(1:n);
target_host = zeros(n_fac, 1);
for f = 1:n_fac
    if target_cell(f) > 0
        target_host(f) = find(host_cells == target_cell(f), 1);
    end
end

D = designs(design_index);
param = @(name) reshape([D.(name)], [], 1);
G = struct();
G.shape = [ny nx];
G.cell_area_m2 = cell_area;
G.dx_m = dx;
G.facility_identifiers = cellfun(@(p) char(p.identifier), placements, 'UniformOutput', false);
G.design_names = cellfun(@(p) designs(strcmp(design_ids, lower(char(p.design_identifier)))).name, ...
    placements, 'UniformOutput', false);
G.rows = rows;
G.columns = cols;
G.cell_index = cell_index;
G.facility_index = facility;
G.area_m2 = piece_area;
G.area_fraction = piece_area / cell_area;
G.area_fraction_map = min(max(fraction_map, 0), 1);
G.host_fraction_map = 1 - G.area_fraction_map;
G.terrain_elevation_m = double(elevation(cell_index));
G.surface_storage_depth_m = param('surface_storage_depth_m');
G.bottom_elevation_m = G.terrain_elevation_m - G.surface_storage_depth_m;
G.surface_manning_n = param('surface_manning_n');
G.media_depth_m = param('media_depth_m');
G.media_porosity = param('media_porosity');
G.media_residual_water_content = param('media_residual_water_content');
G.media_initial_water_content = param('media_initial_water_content');
G.media_ksat_m_s = param('media_ksat_mm_h') / 3.6e6;
G.media_alpha_vg_m_inverse = param('media_alpha_vg_m_inverse');
G.media_n_vg = param('media_n_vg');
G.gravel_depth_m = param('gravel_depth_m');
G.gravel_void_fraction = param('gravel_void_fraction');
G.gravel_ksat_m_s = param('gravel_ksat_mm_h') / 3.6e6;
G.et_crop_coefficient = param('et_crop_coefficient');
G.surface_exchange_coefficient = param('surface_exchange_coefficient');
G.inlet_width_m = param('inlet_width_m');
G.overflow_crest_depth_m = param('overflow_crest_depth_m');
G.overflow_coefficient = param('overflow_coefficient');
G.underdrain_coefficient = param('underdrain_coefficient_m3_s_m_exponent');
G.underdrain_exponent = param('underdrain_exponent');
G.underdrain_sill_depth_m = param('underdrain_sill_depth_m');
native = double(native_ksat_mm_h);
if isscalar(native), native = native * ones(ny, nx); end
native(~isfinite(native)) = 0;
G.native_ksat_m_s = max(native(cell_index), 0) / 3.6e6;
G.edge_start = edge_start;
G.edge_end = edge_end;
G.edge_width_m = edge_width;
G.edge_length_m = edge_length;
G.connection_count = numel(edge_start);
G.target_kind = target_kind;
G.target_cell = target_cell;
G.target_host = target_host;
G.host_cells = host_cells;
G.compartment_host = compartment_host;
G.initial_media_water_volume_m3 = piece_area .* G.media_depth_m .* G.media_initial_water_content;
G.courant = double(get_field(Config, 'courant', 0.5));
G.minimum_substep_s = double(get_field(Config, 'minimum_substep_s', 0.01));
G.maximum_substeps = double(get_field(Config, 'maximum_substeps', 10000));
if ~(G.courant > 0 && G.courant <= 1) || ~(G.minimum_substep_s > 0) || ...
        ~(G.maximum_substeps >= 1 && G.maximum_substeps == round(G.maximum_substeps))
    error('GI courant must be in (0, 1], minimum_substep_s > 0, maximum_substeps a positive integer.');
end
G.compartment_count = n;
G.facility_count = n_fac;
if isempty(edge_start)
    G.edge_start = zeros(0, 1); G.edge_end = zeros(0, 1);
    G.edge_width_m = zeros(0, 1); G.edge_length_m = zeros(0, 1);
end
end

function len = inside_length(polygon, segment)
%INSIDE_LENGTH Length of the part of a segment inside (or on) a polygon.
inside = intersect(polygon, segment);
len = 0;
for i = 1:size(inside, 1) - 1
    if all(isfinite(inside(i, :))) && all(isfinite(inside(i + 1, :)))
        len = len + hypot(inside(i + 1, 1) - inside(i, 1), inside(i + 1, 2) - inside(i, 2));
    end
end
end

function value = get_field(s, name, default)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end

function items = as_cell(value)
if iscell(value)
    items = value(:);
elseif isstruct(value)
    items = num2cell(value(:));
else
    error('GI placements must be a struct array or a cell array of structs.');
end
end
