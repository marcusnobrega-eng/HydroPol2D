function Network = hp2d_urban_drainage_read(path, backend)
%HP2D_URBAN_DRAINAGE_READ Read a storm-drain network GeoPackage.
%   Network = hp2d_urban_drainage_read(path) reads the layers nodes, pipes,
%   inlets and outfalls of the GeoPackage that HydroPol2D-Python writes
%   (write_drainage_geopackage), so one network file serves both codes.
%   Each field of Network is a table with the layer's attribute columns
%   (numbers as double, text as string) plus x and y, the point geometry in
%   the network CRS (NaN for lines or empty geometry).
%
%   backend: 'auto' (default) tries readgeotable (Mapping Toolbox, which
%   HydroPol2D needs anyway, R2021a or newer) and then sqlite (Database
%   Toolbox). 'mapping' or 'sqlite' forces one of them.
%
%   Defaults of HydroPol2D-Python's read_drainage_geopackage are applied
%   here: crest_elevation_m = surface_elevation_m, opening_area_m2 = 0.10,
%   discharge_coefficient = 0.67 and reverse_flow_allowed = true when an
%   inlet column is missing, and initial_head_m = NaN when the node column
%   is missing.

if nargin < 2 || isempty(backend)
    backend = 'auto';
end
path = char(path);
if exist(path, 'file') ~= 2
    error('HydroPol2D:UrbanDrainage', 'Urban-drainage network not found: %s', path);
end

layers = {'nodes', 'pipes', 'inlets', 'outfalls'};
Network = struct();
switch lower(backend)
    case 'auto'
        try
            for i = 1:numel(layers)
                Network.(layers{i}) = read_mapping(path, layers{i});
            end
        catch mapping_error
            try
                for i = 1:numel(layers)
                    Network.(layers{i}) = read_sqlite(path, layers{i});
                end
            catch sqlite_error
                error('HydroPol2D:UrbanDrainage', ...
                    ['Could not read the drainage GeoPackage %s.\n  readgeotable: %s\n  sqlite: %s\n' ...
                    'Reading it needs readgeotable (Mapping Toolbox, R2021a or newer) or the Database Toolbox.'], ...
                    path, mapping_error.message, sqlite_error.message);
            end
        end
    case 'mapping'
        for i = 1:numel(layers)
            Network.(layers{i}) = read_mapping(path, layers{i});
        end
    case 'sqlite'
        for i = 1:numel(layers)
            Network.(layers{i}) = read_sqlite(path, layers{i});
        end
    otherwise
        error('HydroPol2D:UrbanDrainage', 'Unknown GeoPackage backend %s.', backend);
end

% Column defaults of read_drainage_geopackage (HydroPol2D-Python).
n_nodes = height(Network.nodes);
if ~ismember('initial_head_m', Network.nodes.Properties.VariableNames)
    Network.nodes.initial_head_m = NaN(n_nodes, 1);
end
inlets = Network.inlets;
n_inlets = height(inlets);
names = inlets.Properties.VariableNames;
if ~ismember('crest_elevation_m', names) && ismember('surface_elevation_m', names)
    inlets.crest_elevation_m = inlets.surface_elevation_m;
end
if ~ismember('opening_area_m2', names)
    inlets.opening_area_m2 = 0.10 * ones(n_inlets, 1);
end
if ~ismember('discharge_coefficient', names)
    inlets.discharge_coefficient = 0.67 * ones(n_inlets, 1);
end
if ~ismember('reverse_flow_allowed', names)
    inlets.reverse_flow_allowed = ones(n_inlets, 1);
end
Network.inlets = inlets;
Network.path = path;
end

function T = read_mapping(path, layer)
G = readgeotable(path, Layer = layer);
if isempty(G) || width(G) == 0
    T = table('Size', [0 2], 'VariableTypes', {'double', 'double'}, ...
        'VariableNames', {'x', 'y'});
    return
end
x = NaN(height(G), 1);
y = NaN(height(G), 1);
if ismember('Shape', G.Properties.VariableNames)
    shape = G.Shape;
    if isa(shape, 'mappointshape')
        for i = 1:numel(shape)
            if shape(i).NumPoints >= 1
                x(i) = shape(i).X(1);
                y(i) = shape(i).Y(1);
            end
        end
    elseif isa(shape, 'geopointshape')
        for i = 1:numel(shape)
            if shape(i).NumPoints >= 1
                x(i) = shape(i).Longitude(1);
                y(i) = shape(i).Latitude(1);
            end
        end
    end
    G.Shape = [];
end
T = geotable2table_plain(G);
T.x = x;
T.y = y;
end

function T = geotable2table_plain(G)
% Strip the geotable class and give every column a plain type.
names = G.Properties.VariableNames;
T = table();
for i = 1:numel(names)
    v = G.(names{i});
    if isstring(v) || ischar(v) || iscellstr(v)
        T.(names{i}) = string(v);
    elseif islogical(v)
        T.(names{i}) = double(v);
    elseif isnumeric(v)
        T.(names{i}) = double(v);
    else
        T.(names{i}) = string(v);
    end
end
end

function T = read_sqlite(path, layer)
conn = sqlite(path, 'readonly');
cleanup = onCleanup(@() close(conn));
columns = fetch(conn, sprintf( ...
    'SELECT name, type FROM pragma_table_info(''%s'')', layer));
if isempty(columns)
    error('HydroPol2D:UrbanDrainage', 'The GeoPackage has no layer %s.', layer);
end
geometry_column = '';
geometry = fetch(conn, sprintf( ...
    'SELECT column_name FROM gpkg_geometry_columns WHERE table_name = ''%s''', layer));
if ~isempty(geometry)
    geometry_column = char(geometry.column_name(1));
end
names = cellstr(string(columns.name));
types = upper(cellstr(string(columns.type)));
parts = {};
kinds = {};
keep = {};
for i = 1:numel(names)
    name = names{i};
    if strcmpi(name, 'fid') || strcmp(name, geometry_column)
        continue
    end
    quoted = ['"' strrep(name, '"', '""') '"'];
    if is_real_type(types{i})
        parts{end + 1} = sprintf('(%s IS NULL) AS "%s__null", IFNULL(%s, 0.0) AS "%s"', ...
            quoted, name, quoted, name); %#ok<AGROW>
        kinds{end + 1} = 'number'; %#ok<AGROW>
    elseif is_integer_type(types{i})
        parts{end + 1} = sprintf('(%s IS NULL) AS "%s__null", IFNULL(%s, 0) AS "%s"', ...
            quoted, name, quoted, name); %#ok<AGROW>
        kinds{end + 1} = 'number'; %#ok<AGROW>
    else
        parts{end + 1} = sprintf('(%s IS NULL) AS "%s__null", IFNULL(%s, '''') AS "%s"', ...
            quoted, name, quoted, name); %#ok<AGROW>
        kinds{end + 1} = 'text'; %#ok<AGROW>
    end
    keep{end + 1} = name; %#ok<AGROW>
end
if ~isempty(geometry_column)
    parts{end + 1} = sprintf('IFNULL(hex("%s"), '''') AS "geom__hex"', geometry_column);
end
sql = sprintf('SELECT %s FROM "%s" ORDER BY fid', strjoin(parts, ', '), layer);
R = fetch(conn, sql);
n = height(R);
T = table();
for i = 1:numel(keep)
    name = keep{i};
    is_null = double(R.([name '__null'])) ~= 0;
    if strcmp(kinds{i}, 'number')
        v = double(R.(name));
        v(is_null) = NaN;
    else
        v = string(R.(name));
        v(ismissing(v)) = "";
        v(is_null) = "";
    end
    T.(name) = v;
end
x = NaN(n, 1);
y = NaN(n, 1);
if ~isempty(geometry_column)
    hex = string(R.geom__hex);
    for i = 1:n
        [x(i), y(i)] = gpkg_point(hex(i));
    end
end
T.x = x;
T.y = y;
end

function tf = is_real_type(type)
tf = any(strcmp(type, {'REAL', 'DOUBLE', 'FLOAT', 'NUMERIC', 'DOUBLE PRECISION'}));
end

function tf = is_integer_type(type)
tf = any(strcmp(type, {'INTEGER', 'INT', 'BOOLEAN', 'TINYINT', 'SMALLINT', 'MEDIUMINT', 'BIGINT'}));
end

function [x, y] = gpkg_point(hex)
% First point of a GeoPackage binary geometry (header, envelope, WKB).
x = NaN; y = NaN;
if strlength(hex) < 16
    return
end
bytes = uint8(hex2dec(reshape(char(hex), 2, []).'));
if bytes(1) ~= uint8('G') || bytes(2) ~= uint8('P')
    return
end
flags = bytes(4);
if bitand(flags, uint8(16))
    return  % empty geometry
end
envelope_bytes = [0 32 48 48 64];
envelope = bitand(bitshift(flags, -1), uint8(7));
offset = 8 + envelope_bytes(double(envelope) + 1);
wkb = bytes(offset + 1:end);
little = wkb(1) == 1;
type = typecast(wkb(2:5), 'uint32');
if ~little
    type = swapbytes(type);
end
type = mod(double(type), 1000);
switch type
    case 1        % point
        first = 6;
    case 4        % multipoint: its first point
        first = 15;
    otherwise     % lines and polygons carry no point
        return
end
if numel(wkb) < first + 15
    return
end
xy = typecast(wkb(first:first + 15), 'double');
if ~little
    xy = swapbytes(xy);
end
x = xy(1);
y = xy(2);
end
