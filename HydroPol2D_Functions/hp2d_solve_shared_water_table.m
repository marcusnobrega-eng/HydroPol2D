function R = hp2d_solve_shared_water_table(column, cell_head_m, water_mm, net_volume_m3, host, column_area_m2, ground_m, base_depth_m, iterations)
%HP2D_SOLVE_SHARED_WATER_TABLE Apply a net volume per cell to columns that share one table.
%   column, water_mm and column_area_m2 are per column (mesh HRUs); cell_head_m,
%   ground_m, base_depth_m and net_volume_m3 are per cell; host (1-based) maps each
%   column to its cell. The common depth is found by bisection on the cell storage
%   curve, which is monotone in depth (default 60 iterations).
if nargin < 9 || isempty(iterations), iterations = 60; end
host = host(:);
area = column_area_m2(:);
n_cells = numel(cell_head_m);
head_columns = cell_head_m(host);
volume = net_volume_m3(:);
    function change = cell_change(depth_cell)
        delta = hp2d_storage_change_to_depth(column, head_columns, water_mm, depth_cell(host));
        change = accumarray(host, delta/1000.*area, [n_cells, 1]);
    end
base = max(base_depth_m(:), 0);
top_change = cell_change(zeros(n_cells, 1));
bottom_change = cell_change(base);
R.exfiltration_m3 = max(volume - top_change, 0);
R.unmet_removal_m3 = max(bottom_change - volume, 0);
target = min(max(volume, bottom_change), top_change);
low = zeros(n_cells, 1);
high = base;
for i = 1:iterations
    middle = 0.5*(low + high);
    above = cell_change(middle) > target;
    low(above) = middle(above);
    high(~above) = middle(~above);
end
depth = 0.5*(low + high);
depth(volume <= bottom_change) = base(volume <= bottom_change);
depth(volume >= top_change) = 0;
[~, R.water_mm] = hp2d_storage_change_to_depth(column, head_columns, water_mm, depth(host));
R.head_m = ground_m(:) - depth;
end
