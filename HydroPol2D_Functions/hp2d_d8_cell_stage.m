function [stage_m, channel_depth_m] = hp2d_d8_cell_stage(Cells, volume_m3, idx)
%HP2D_D8_CELL_STAGE Water-surface elevation of each cell from its volume.
%   Cells with a storage table (Cells.is_table, the D8 river cells) invert
%   the piecewise-linear volume curve; every other cell is flat:
%       table cell: stage = datum + zeta_j + (V - V_j) / A_j
%       flat cell:  stage = terrain + V / cell_area
%   where j is the bracketing breakpoint and A_j its wet plan area (the cell
%   plan area above the top breakpoint or where A_j <= 0). channel_depth_m is
%   stage - datum on table cells and 0 elsewhere. This is
%   d8_river_water_surface_jax / d8_river_channel_depth_jax in
%   HydroPol2D-Python.
%   idx (optional) restricts the evaluation to those linear cell indices;
%   volume_m3 then holds one value per entry of idx.

if nargin < 3 || isempty(idx)
    idx = (1:numel(Cells.datum)).';
end
idx = idx(:);
query = max(volume_m3(:), 0);

table_volume = Cells.volume(idx, :);
slot = hp2d_d8_row_slot(table_volume, query, Cells.count(idx));
rows = (1:numel(idx)).';
lin = sub2ind(size(table_volume), rows, slot);
base_volume = table_volume(lin);
zeta = Cells.zeta(idx, :);
base_zeta = zeta(lin);
area = Cells.area(idx, :);
wet_area = area(lin);
plan = Cells.plan_area(idx);
use_plan = slot >= Cells.count(idx) | wet_area <= 0;
wet_area(use_plan) = plan(use_plan);

table_stage = Cells.datum(idx) + base_zeta + (query - base_volume) ./ max(wet_area, 1e-30);
stage_m = Cells.terrain(idx) + query ./ Cells.cell_area;
is_table = Cells.is_table(idx);
stage_m(is_table) = table_stage(is_table);

if nargout > 1
    channel_depth_m = zeros(size(stage_m), 'like', stage_m);
    channel_depth_m(is_table) = max(stage_m(is_table) - Cells.datum(idx(is_table)), 0);
end

end
