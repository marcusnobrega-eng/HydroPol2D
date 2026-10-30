function volume_m3 = hp2d_d8_cell_volume(Cells, stage_m, idx)
%HP2D_D8_CELL_VOLUME Stored volume of table cells at a water-surface elevation.
%   Interpolates the piecewise-linear volume curve of each cell table
%   between its breakpoints; above the top breakpoint the cell fills at its
%   plan area. Returns 0 where the stage is at or below the datum. This is
%   CellVolumeTable.volume in HydroPol2D-Python (subgrid.py). It is the
%   inverse of hp2d_d8_cell_stage on table cells.
%   idx (optional) restricts the evaluation to those linear cell indices.

if nargin < 3 || isempty(idx)
    idx = (1:numel(Cells.datum)).';
end
idx = idx(:);
zeta = max(stage_m(:) - Cells.datum(idx), 0);
keys = Cells.zeta(idx, :);
table_volume = Cells.volume(idx, :);
count = Cells.count(idx);
slot = hp2d_d8_row_slot(keys, zeta, count);
next_slot = min(slot + 1, count);
rows = (1:numel(idx)).';
z0 = keys(sub2ind(size(keys), rows, slot));
z1 = keys(sub2ind(size(keys), rows, next_slot));
v0 = table_volume(sub2ind(size(keys), rows, slot));
v1 = table_volume(sub2ind(size(keys), rows, next_slot));
top = slot >= count;
weight = min(max((zeta - z0) ./ max(z1 - z0, 1e-30), 0), 1);
weight(top) = 0;
volume_m3 = v0 + weight .* (v1 - v0);
volume_m3(top) = v0(top) + (zeta(top) - z0(top)) .* Cells.plan_area(idx(top));
volume_m3 = max(volume_m3, 0);
volume_m3(~(zeta > 0)) = 0;

end
