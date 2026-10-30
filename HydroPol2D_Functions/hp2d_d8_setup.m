function D8 = hp2d_d8_setup(elevation, roughness, Resolution, Options)
%HP2D_D8_SETUP Build the D8 local-inertial graph, tables and outlets.
%   D8 = hp2d_d8_setup(elevation, roughness, Resolution, Options) prepares
%   the structures used by D8_Local_Inertial_Model, the Wflow-style D8
%   scheme of HydroPol2D-Python (routing/d8_river_local_inertial.py):
%
%   1. One D8 receiver per cell on the routing DEM (hp2d_d8_receivers) and
%      the contributing area along that graph.
%   2. River cells: contributing area >= Options.min_area_km2. Their width
%      and bank depth follow the area power laws
%          W = width_coefficient * A^width_exponent,
%          H = depth_coefficient * A^depth_exponent   (A in km2),
%      unless positive river rasters are given (Options.river_width/depth).
%   3. The links: river cell -> downstream river cell (hybrid mode), or every
%      valid cell -> its receiver (Options.full_domain = 1). Link length is
%      dx (cardinal) or dx*sqrt(2) (diagonal). Links are ordered as in Python
%      (source cells in row-major order).
%   4. Hydraulic tables. Cell tables give the stored volume against the
%      stage; link tables give flow area and conveyance against the depth.
%      With no Options.table_path, prismatic river tables are built here
%      (write_prismatic_d8_river_tables): a rectangular channel of width W
%      and depth H below the bank on river cells, the whole cell above the
%      bank, conveyance K = A R^(2/3) / n with the channel Manning n. A table
%      file in HydroPol2D-Python's netCDF layout can be read instead. In
%      full-domain mode the other links get flat rectangular rows of width
%      cell_area / L and the mean Manning n of their two cells
%      (expand_d8_tables_to_full_domain).
%   5. The D4 faces that cardinal river links replace (hybrid mode) and the
%      outlets: the outlet cells (Options.outlet_mask) that are terminal D8
%      nodes, each with a normal-depth (slope Options.outlet_slope) or
%      critical-depth (Options.outlet_type = 2) control.
%
%   Options fields: full_domain, min_depth_m, min_area_km2,
%   width_coefficient, width_exponent, depth_coefficient, depth_exponent,
%   channel_manning_n, river_width, river_depth (rasters or []),
%   outlet_mask, outlet_slope, outlet_type, table_path,
%   max_stage_above_bed_m.

z = double(elevation);
[ny, nx] = size(z);
N = ny * nx;
cell_area = Resolution^2;
valid = isfinite(z);
n_land = double(roughness);
n_land(~valid) = NaN;

if ~(Options.min_area_km2 > 0)
    error('HydroPol2D:D8MinArea', ...
        'D8 local-inertial routing needs a positive stream threshold (min_area).');
end
if ~(Options.min_depth_m >= 0 && Options.min_depth_m <= 0.01)
    error('HydroPol2D:D8MinDepth', 'd8_min_depth_m must lie in [0, 0.01] m.');
end

%% Receivers, river cells and river geometry
[downstream, contributing] = hp2d_d8_receivers(z);
area_km2 = contributing * cell_area / 1e6;
stream = isfinite(area_km2) & area_km2 >= Options.min_area_km2 & valid;
if ~any(stream(:))
    error('HydroPol2D:D8NoRivers', ...
        'D8 local-inertial routing needs river cells; lower min_area.');
end
width = zeros(ny, nx);
depth = zeros(ny, nx);
width(stream) = Options.width_coefficient * area_km2(stream).^Options.width_exponent;
depth(stream) = Options.depth_coefficient * area_km2(stream).^Options.depth_exponent;
width = override_positive(width, Options, 'river_width');
depth = override_positive(depth, Options, 'river_depth');
if any(~isfinite(width(stream)) | width(stream) <= 0)
    error('HydroPol2D:D8RiverWidth', ...
        'D8 river width must be finite and positive on river cells (check beta_1, beta_2).');
end
if any(~isfinite(depth(stream)) | depth(stream) <= 0)
    error('HydroPol2D:D8RiverDepth', ...
        'D8 river depth must be finite and positive on river cells (check alfa_1, alfa_2).');
end

%% Topology
RiverLinks = build_links(stream, downstream, Resolution);
if isempty(RiverLinks.source)
    error('HydroPol2D:D8NoRiverLinks', 'The river network has no D8 links.');
end
full_domain = isfield(Options, 'full_domain') && Options.full_domain == 1;

%% Tables
channel_n = Options.channel_manning_n;
if isfield(Options, 'table_path') && strlength(strtrim(string(Options.table_path))) > 0
    [Cells, RiverFaces] = hp2d_d8_read_tables(char(Options.table_path), ny, nx, RiverLinks);
else
    [Cells, RiverFaces] = prismatic_tables(z, stream, width, depth, n_land, channel_n, ...
        RiverLinks, cell_area, Options.max_stage_above_bed_m);
end
Cells.terrain = z(:);
Cells.cell_area = cell_area;
river_cells = unique([RiverLinks.source; RiverLinks.target]);
if ~all(Cells.is_table(river_cells))
    error('HydroPol2D:D8MissingCellTable', 'Every D8 river cell must have a storage table.');
end

if full_domain
    Links = build_links(valid, downstream, Resolution);
    Faces = expand_full_domain(Links, RiverLinks, RiverFaces, z, n_land, cell_area);
else
    Links = RiverLinks;
    Faces = RiverFaces;
end
Faces = share_zeta_if_uniform(Faces);

%% D4 faces replaced by cardinal river links (hybrid floodplain)
inactive_x = false(ny, nx);
inactive_y = false(ny, nx);
if ~full_domain
    [sr, sc] = ind2sub([ny nx], Links.source);
    [tr, tc] = ind2sub([ny nx], Links.target);
    horizontal = sr == tr;
    vertical = sc == tc;
    inactive_x(sub2ind([ny nx], sr(horizontal), min(sc(horizontal), tc(horizontal)))) = true;
    inactive_y(sub2ind([ny nx], max(sr(vertical), tr(vertical)), sc(vertical))) = true;
end

%% Outlets
if full_domain
    nodes = valid;
else
    nodes = stream;
end
terminal = nodes & downstream == 0;
outlet_mask = logical(Options.outlet_mask);
outlet_cells = find(outlet_mask(:) & terminal(:));
if any(outlet_mask(:)) && isempty(outlet_cells)
    error('HydroPol2D:D8NoTerminalOutlet', ...
        'No outlet cell is a terminal D8 node (a river cell with no downstream cell).');
end
Outlets = struct();
Outlets.cell = outlet_cells;
if isscalar(Options.outlet_slope)
    Outlets.slope = Options.outlet_slope * ones(size(outlet_cells));
else
    Outlets.slope = double(Options.outlet_slope(outlet_cells));
end
Outlets.slope = max(Outlets.slope, 0);
Outlets.critical = repmat(Options.outlet_type == 2, size(outlet_cells));
outlet_width = width;
if full_domain
    outlet_width(~stream) = Resolution;
end
Outlets.width = outlet_width(outlet_cells);
outlet_n = n_land;
outlet_n(stream) = channel_n;
Outlets.manning_n = outlet_n(outlet_cells);
Outlets.bed = Cells.datum(outlet_cells);

%% Bankfull channel volume (kept out of infiltration and evaporation)
bankfull = zeros(N, 1);
table_cells = find(Cells.is_table & valid(:));
bankfull(table_cells) = hp2d_d8_cell_volume(Cells, Cells.terrain(table_cells), table_cells);

D8 = struct();
D8.full_domain = full_domain;
D8.ny = ny;
D8.nx = nx;
D8.N = N;
D8.Resolution = Resolution;
D8.cell_area = cell_area;
D8.min_depth_m = Options.min_depth_m;
D8.valid_grid = valid;
D8.terrain_grid = z;
D8.roughness_grid = n_land;
D8.downstream = downstream;
D8.contributing_cells = contributing;
D8.stream = stream;
D8.river_width = width;
D8.river_depth = depth;
D8.channel_manning_n = channel_n;
D8.Links = Links;
D8.RiverLinks = RiverLinks;
D8.Cells = Cells;
D8.Faces = Faces;
D8.inactive_x = inactive_x;
D8.inactive_y = inactive_y;
D8.Outlets = Outlets;
D8.bankfull_volume_m3 = reshape(bankfull, ny, nx);

end

function values = override_positive(values, Options, name)
% A supplied raster replaces the power law where it is finite and positive.
if isfield(Options, name) && ~isempty(Options.(name))
    supplied = double(Options.(name));
    use = isfinite(supplied) & supplied > 0;
    values(use) = supplied(use);
end
end

function Links = build_links(node_mask, downstream, Resolution)
% One link per node with a downstream cell, sources in row-major order.
[ny, nx] = size(node_mask);
source = find(node_mask(:) & downstream(:) > 0);
[sr, sc] = ind2sub([ny nx], source);
[~, order] = sort((sr - 1) * nx + sc);
source = source(order);
target = downstream(source);
if any(~node_mask(target))
    error('HydroPol2D:D8LinkTarget', 'Every D8 link must end in a D8 node.');
end
[sr, sc] = ind2sub([ny nx], source);
[tr, tc] = ind2sub([ny nx], target);
row_step = abs(tr - sr);
col_step = abs(tc - sc);
if any(row_step > 1 | col_step > 1 | (row_step == 0 & col_step == 0))
    error('HydroPol2D:D8LinkNeighbour', 'Every D8 link must join two neighbouring cells.');
end
Links = struct();
Links.source = source;
Links.target = target;
Links.length = hypot(Resolution * (col_step > 0), Resolution * (row_step > 0));
Links.diagonal = row_step == 1 & col_step == 1;
% Node length: the link length at a source; the mean incoming length at a
% terminal node, or dx when nothing flows in (build_d8_river_topology).
node_length = nan(ny * nx, 1);
node_length(source) = Links.length;
terminal = find(node_mask(:) & downstream(:) == 0);
for ii = 1:numel(terminal)
    incoming = Links.length(target == terminal(ii));
    if isempty(incoming)
        node_length(terminal(ii)) = Resolution;
    else
        node_length(terminal(ii)) = mean(incoming);
    end
end
Links.node_length = node_length;
end

function [Cells, Faces] = prismatic_tables(z, stream, width, depth, n_land, channel_n, ...
    Links, cell_area, max_stage)
% write_prismatic_d8_river_tables of HydroPol2D-Python, held in memory.
N = numel(z);
river = stream(:);
land = z(:);
datum = land;
datum(river) = land(river) - depth(river);
plan_channel = min(width(:) .* Links.node_length, cell_area);

zeta = repmat([0 1 10 max_stage], N, 1);
volume = cell_area * zeta;
wet_area = cell_area * ones(N, 4);
bank = depth(river);
river_area = plan_channel(river);
zeta(river, 2) = bank;
zeta(river, 3) = bank + 1;
zeta(river, 4) = bank + max_stage;
volume(river, 2) = river_area .* bank;
volume(river, 3) = volume(river, 2) + cell_area;
volume(river, 4) = volume(river, 3) + cell_area * (max_stage - 1);
wet_area(river, 1) = river_area;

Cells = struct();
Cells.datum = datum;
Cells.zeta = zeta;
Cells.volume = volume;
Cells.area = wet_area;
Cells.count = 4 * ones(N, 1);
Cells.plan_area = cell_area * ones(N, 1);
Cells.is_table = river;

source = Links.source;
target = Links.target;
edge_width = min(width(source), width(target));
edge_n = channel_n * ones(size(source));
invalid = ~isfinite(edge_n) | edge_n <= 0;
edge_n(invalid) = 0.5 * (n_land(source(invalid)) + n_land(target(invalid)));
levels = [0, geometric_levels(1e-4, max_stage, 121)];
flow_area = edge_width .* levels;
perimeter = edge_width + 2 * levels;
perimeter(:, 1) = edge_width;
radius = flow_area ./ max(perimeter, 1e-30);
Faces = struct();
Faces.datum = max(datum(source), datum(target));
Faces.zeta = repmat(levels, numel(source), 1);
Faces.area = flow_area;
Faces.conveyance = flow_area .* radius.^(2/3) ./ edge_n;
Faces.count = numel(levels) * ones(numel(source), 1);
Faces.length = Links.length;
Faces.shared_zeta = [];
end

function Faces = expand_full_domain(Links, RiverLinks, RiverFaces, z, n_land, cell_area)
% expand_d8_tables_to_full_domain of HydroPol2D-Python: river links keep
% their rows; every other link gets a flat rectangular section.
point_count = size(RiverFaces.area, 2);
if isempty(RiverFaces.shared_zeta)
    max_stage = max(RiverFaces.zeta(:));
else
    max_stage = max(RiverFaces.shared_zeta);
end
levels = [0, geometric_levels(1e-4, max_stage, point_count - 1)];
source = Links.source;
target = Links.target;
w = cell_area ./ Links.length;
edge_n = 0.5 * (n_land(source) + n_land(target));
area = w .* levels;
perimeter = w + 2 * levels;
radius = area ./ max(perimeter, 1e-30);
Faces = struct();
Faces.datum = max(z(source), z(target));
Faces.zeta = repmat(levels, numel(source), 1);
Faces.area = area;
Faces.conveyance = area .* radius.^(2/3) ./ max(edge_n, 1e-30);
Faces.count = point_count * ones(numel(source), 1);
Faces.length = Links.length;
Faces.shared_zeta = [];

[is_river, row] = ismember([source target], [RiverLinks.source RiverLinks.target], 'rows');
full_rows = find(is_river);
river_rows = row(is_river);
Faces.datum(full_rows) = RiverFaces.datum(river_rows);
if isempty(RiverFaces.shared_zeta)
    Faces.zeta(full_rows, :) = RiverFaces.zeta(river_rows, :);
else
    Faces.zeta(full_rows, :) = repmat(RiverFaces.shared_zeta, numel(river_rows), 1);
end
Faces.area(full_rows, :) = RiverFaces.area(river_rows, :);
Faces.conveyance(full_rows, :) = RiverFaces.conveyance(river_rows, :);
Faces.count(full_rows) = RiverFaces.count(river_rows);
end

function Faces = share_zeta_if_uniform(Faces)
% One breakpoint row for every link (the prismatic and flat tables) is kept
% once, so the lookup does not compare every column of every row.
if isempty(Faces.shared_zeta) && ~isempty(Faces.zeta)
    first = Faces.zeta(1, :);
    if all(diff(first) > 0) && all(Faces.count == numel(first)) && ...
            all(Faces.zeta == first, 'all')
        Faces.shared_zeta = first;
        Faces.zeta = [];
    end
end
end

function values = geometric_levels(first, last, count)
% numpy.geomspace(first, last, count): evenly spaced log10 levels with the
% end points exact.
exponents = log10(first) + (0:count-1) * ((log10(last) - log10(first)) / (count - 1));
values = 10 .^ exponents;
values(1) = first;
values(end) = last;
end
