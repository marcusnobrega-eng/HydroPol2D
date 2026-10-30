function [downstream, contributing] = hp2d_d8_receivers(elevation)
%HP2D_D8_RECEIVERS One D8 receiver per cell and the contributing-cell count.
%   downstream(i) is the linear index of the cell that cell i drains to, and
%   0 at sinks, boundary outlets and NaN cells. The receiver is the neighbour
%   with the steepest drop divided by the cardinal or diagonal distance; a tie
%   goes to the neighbour with the larger linear index. A flat cell (best
%   drop 0) drains to an equal neighbour only when that neighbour has the
%   larger linear index, so the graph has no cycles. This is
%   _hydropol_d8_receivers in HydroPol2D-Python (terrain/pyflwdir_backend.py);
%   MATLAB linear indices are Python's column-major ("Fortran") indices + 1.
%
%   contributing(i) counts the active cells that drain through cell i,
%   itself included (NaN outside the domain). The count follows the graph in
%   topological order.

z = double(elevation);
[ny, nx] = size(z);
active = isfinite(z);

% (row, column) offsets in the order HydroPol2D-Python tests them.
offsets = [-1 0; 0 1; 1 0; 0 -1; -1 1; 1 1; 1 -1; -1 -1];
[col_grid, row_grid] = meshgrid(1:nx, 1:ny);
best_gradient = -inf(ny, nx);
best_index = zeros(ny, nx);
for d = 1:size(offsets, 1)
    nr = row_grid + offsets(d, 1);
    nc = col_grid + offsets(d, 2);
    inside = nr >= 1 & nr <= ny & nc >= 1 & nc <= nx;
    neighbour = zeros(ny, nx);
    neighbour(inside) = sub2ind([ny nx], nr(inside), nc(inside));
    z_neighbour = nan(ny, nx);
    z_neighbour(inside) = z(neighbour(inside));
    gradient = (z - z_neighbour) / hypot(offsets(d, 1), offsets(d, 2));
    take = active & isfinite(z_neighbour) & (gradient > best_gradient | ...
        (gradient == best_gradient & neighbour > best_index));
    best_gradient(take) = gradient(take);
    best_index(take) = neighbour(take);
end

own_index = reshape(1:ny*nx, ny, nx);
drains = active & best_index > 0 & (best_gradient > 0 | ...
    (best_gradient == 0 & best_index > own_index));
downstream = zeros(ny, nx);
downstream(drains) = best_index(drains);

% Contributing cells, one topological wave at a time (Kahn's algorithm).
contributing = double(active);
source_cells = find(drains);
pending = accumarray(downstream(source_cells), 1, [ny*nx 1]);
front = find(active(:) & pending == 0);
while ~isempty(front)
    front = front(downstream(front) > 0);
    targets = downstream(front);
    if isempty(targets)
        break
    end
    contributing(:) = contributing(:) + accumarray(targets, contributing(front), [ny*nx 1]);
    pending = pending - accumarray(targets, 1, [ny*nx 1]);
    released = unique(targets);
    front = released(pending(released) == 0);
end
contributing(~active) = NaN;

end
