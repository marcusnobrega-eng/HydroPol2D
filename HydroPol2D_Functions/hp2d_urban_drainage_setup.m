function Drainage_States = hp2d_urban_drainage_setup(Input, elevation_m, valid, cell_area_m2, x_ul_m, y_ul_m, resolution_m, outlet_mask)
%HP2D_URBAN_DRAINAGE_SETUP Read, check and initialise the storm-drain
%   network of a grid run (flags.flag_urban_drainage = 1).
%
%   Input: InputData_Bypass.Urban_Drainage, with
%       network_path  GeoPackage written by HydroPol2D-Python (layers nodes,
%                     pipes, inlets, outfalls); a relative path is taken
%                     from the current folder.
%       hydraulics    optional struct (see hp2d_urban_drainage_prepare).
%   elevation_m, valid: model terrain and active cells. x_ul_m, y_ul_m:
%   GIS_data.xulcorner / yulcorner; cell (r, c) holds the points with
%   round((x - x_ul)/dx) = c and round((y_ul - y)/dx) = r, the rule of the
%   reservoir controls.
%
%   Inlet cells come from the inlet points, so a cropped or resampled DEM
%   still finds them. An inlet without a point uses its row and column
%   (0-based, as Python writes them). Water that leaves through an outfall
%   is booked as outlet flow at the active cell nearest the outfall point
%   (the model outlet when the outfall has no point).

if ~isstruct(Input) || ~isfield(Input, 'network_path') || isempty(Input.network_path)
    error('HydroPol2D:UrbanDrainage', ...
        'flag_urban_drainage = 1 requires InputData_Bypass.Urban_Drainage.network_path.');
end
hydraulics = struct();
if isfield(Input, 'hydraulics') && isstruct(Input.hydraulics)
    hydraulics = Input.hydraulics;
end
path = char(Input.network_path);
Network = hp2d_urban_drainage_read(path);
elevation_m = double(gather(elevation_m));
valid = logical(gather(valid)) & isfinite(elevation_m);
[ny, nx] = size(elevation_m);

% Inlet cells.
inlets = Network.inlets;
n_inlets = height(inlets);
inlet_rc = zeros(n_inlets, 2);
for i = 1:n_inlets
    [r, c] = point_cell(inlets.x(i), inlets.y(i), x_ul_m, y_ul_m, resolution_m);
    if ~(r >= 1 && r <= ny && c >= 1 && c <= nx)
        if ~all(ismember({'row', 'column'}, inlets.Properties.VariableNames))
            error('HydroPol2D:UrbanDrainage', ...
                'Drainage inlet %s has no point inside the grid and no row/column.', ...
                string(inlets.identifier(i)));
        end
        r = inlets.row(i) + 1;
        c = inlets.column(i) + 1;
    end
    inlet_rc(i, :) = [r c];
end

N = hp2d_urban_drainage_prepare(Network, elevation_m, valid, cell_area_m2, hydraulics, inlet_rc);

% Where outfall water is booked as outlet flow.
[rows, cols] = find(valid);
active = sub2ind([ny nx], rows, cols);
booking = zeros(numel(N.node_ids), 1);
boundary = find(~N.dynamic);
outlet_cells = find(logical(gather(outlet_mask)) & valid);
for k = 1:numel(boundary)
    index = boundary(k);
    x = Network.nodes.x(index);
    y = Network.nodes.y(index);
    if isfinite(x) && isfinite(y)
        xc = x_ul_m + cols * resolution_m;
        yc = y_ul_m - rows * resolution_m;
        [~, nearest] = min((xc - x).^2 + (yc - y).^2);
        booking(index) = active(nearest);
    elseif ~isempty(outlet_cells)
        booking(index) = outlet_cells(1);
    else
        booking(index) = active(1);
    end
end

Drainage_States = struct();
Drainage_States.network_path = path;
Drainage_States.Network = N;
Drainage_States.State = hp2d_urban_drainage_init(N);
Drainage_States.booking_cell = booking;
Drainage_States.step_outflow_m3 = 0;
Drainage_States.step_inflow_m3 = 0;
Drainage_States.History = struct('index', zeros(0, 1), 'time_min', zeros(0, 1));
fprintf(['Urban drainage: %d nodes, %d pipes, %d inlets, %d outfalls from %s ' ...
    '(initial storage %.3f m3).\n'], numel(N.node_ids), numel(N.pipe_ids), ...
    numel(N.inlet_ids), numel(boundary), path, Drainage_States.State.initial_storage_m3);
end

function [r, c] = point_cell(x, y, x_ul, y_ul, dx)
r = NaN; c = NaN;
if isfinite(x) && isfinite(y)
    c = round((x - x_ul) / dx);
    r = round((y_ul - y) / dx);
end
end
