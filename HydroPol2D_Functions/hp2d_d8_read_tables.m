function [Cells, Faces] = hp2d_d8_read_tables(path, ny, nx, Links)
%HP2D_D8_READ_TABLES Read D8 hydraulic tables written by HydroPol2D-Python.
%   The netCDF file holds one storage table per grid cell (cell_datum_m,
%   cell_zeta_m, cell_volume_m3, cell_wet_area_m2, cell_point_count,
%   cell_plan_area_m2, optional cell_is_subgrid) and one conveyance table
%   per D8 river link (face_datum_m, face_zeta_m, face_flow_area_m2,
%   face_conveyance, face_point_count, face_length_m, optional
%   d8_link_source_index / d8_link_target_index / d8_link_length_m), as
%   write_prismatic_d8_river_tables and load_d8_river_subgrid_tables use.
%   Python numbers cells row-major from 0; they are reordered to MATLAB's
%   column-major cells here. The link rows must follow Links (source cells
%   in row-major order), which hp2d_d8_setup builds.

if ~isfile(path)
    error('HydroPol2D:D8TableFile', 'D8 hydraulic table file not found: %s', path);
end
N = ny * nx;
E = numel(Links.source);
% MATLAB cell i (column-major) is Python cell (row - 1) * nx + (col - 1).
[row, col] = ind2sub([ny nx], (1:N).');
python_cell = (row - 1) * nx + col;

Cells = struct();
Cells.datum = column(read_var(path, 'cell_datum_m'), N, python_cell);
Cells.zeta = table_rows(read_var(path, 'cell_zeta_m'), N, python_cell);
Cells.volume = table_rows(read_var(path, 'cell_volume_m3'), N, python_cell);
Cells.area = table_rows(read_var(path, 'cell_wet_area_m2'), N, python_cell);
Cells.count = column(read_var(path, 'cell_point_count'), N, python_cell);
Cells.plan_area = column(read_var(path, 'cell_plan_area_m2'), N, python_cell);
if has_var(path, 'cell_is_subgrid')
    Cells.is_table = column(read_var(path, 'cell_is_subgrid'), N, python_cell) ~= 0;
else
    Cells.is_table = true(N, 1);
end

identity = (1:E).';
Faces = struct();
Faces.datum = column(read_var(path, 'face_datum_m'), E, identity);
Faces.zeta = table_rows(read_var(path, 'face_zeta_m'), E, identity);
Faces.area = table_rows(read_var(path, 'face_flow_area_m2'), E, identity);
Faces.conveyance = table_rows(read_var(path, 'face_conveyance'), E, identity);
Faces.count = column(read_var(path, 'face_point_count'), E, identity);
Faces.length = column(read_var(path, 'face_length_m'), E, identity);
Faces.shared_zeta = [];

if has_var(path, 'd8_link_source_index') && has_var(path, 'd8_link_target_index') && ...
        has_var(path, 'd8_link_length_m')
    source = to_matlab_cell(column(read_var(path, 'd8_link_source_index'), E, identity), ny, nx);
    target = to_matlab_cell(column(read_var(path, 'd8_link_target_index'), E, identity), ny, nx);
    link_length = column(read_var(path, 'd8_link_length_m'), E, identity);
    if ~isequal(source, Links.source) || ~isequal(target, Links.target) || ...
            any(abs(link_length - Links.length) > 1e-9)
        error('HydroPol2D:D8TableLinks', ...
            'The table D8 links do not match the D8 river graph of this DEM.');
    end
elseif any(abs(Faces.length - Links.length) > 1e-9)
    error('HydroPol2D:D8TableLengths', ...
        'The table link lengths do not match the D8 river graph of this DEM.');
end

end

function value = read_var(path, name)
if ~has_var(path, name)
    error('HydroPol2D:D8TableVariable', 'D8 table file %s has no variable %s.', path, name);
end
value = double(ncread(path, name));
end

function found = has_var(path, name)
info = ncinfo(path);
found = any(strcmp({info.Variables.Name}, name));
end

function values = column(values, n, order)
values = values(:);
if numel(values) ~= n
    error('HydroPol2D:D8TableSize', 'D8 table has %d rows; the grid needs %d.', numel(values), n);
end
values = values(order);
end

function values = table_rows(values, n, order)
% ncread returns netCDF (row, point) variables as point x row.
if size(values, 2) == n && size(values, 1) ~= n
    values = values.';
end
if size(values, 1) ~= n
    error('HydroPol2D:D8TableSize', 'D8 table has %d rows; the grid needs %d.', size(values, 1), n);
end
values = values(order, :);
end

function cells = to_matlab_cell(python_index, ny, nx)
% Python row-major index from 0 to MATLAB linear index.
row = floor(python_index / nx) + 1;
col = mod(python_index, nx) + 1;
cells = sub2ind([ny nx], row, col);
end
