function [volume_m3, outlet_m3_s, outlet_depth_m] = hp2d_d8_outlet_flow(volume_m3, D8, dt_s)
%HP2D_D8_OUTLET_FLOW Normal- or critical-depth outflow at the D8 terminal cells.
%   Each outlet is a rectangular section of width W at a terminal D8 node
%   (a node with no downstream cell). With h = max(stage - bed - min_depth, 0):
%       normal:   Q = W h R^(2/3) sqrt(S0) / n,  R = W h / (W + 2 h)
%       critical: Q = W sqrt(g) h^(3/2)
%   capped at the cell volume over the step. The water leaves the model.
%   This is apply_d8_river_outlets_jax in HydroPol2D-Python. The bed is the
%   table datum (channel bed of a river cell) and n is the channel Manning n
%   on river cells (the n of its link tables) and the land-cover n elsewhere.
%   outlet_m3_s is per cell (column, D8.N cells); outlet_depth_m the depth h.

g = 9.81;
volume_m3 = max(volume_m3(:), 0);
n_cells = D8.N;
outlet_m3_s = zeros(n_cells, 1, 'like', volume_m3);
outlet_depth_m = zeros(n_cells, 1, 'like', volume_m3);
O = D8.Outlets;
if isempty(O.cell)
    return
end
stage = hp2d_d8_cell_stage(D8.Cells, volume_m3(O.cell), O.cell);
depth = max(stage - O.bed - D8.min_depth_m, 0);
area = O.width .* depth;
radius = area ./ max(O.width + 2 * depth, 1e-30);
normal = area .* radius.^(2/3) ./ max(O.manning_n, 1e-30) .* sqrt(O.slope);
critical = O.width * sqrt(g) .* depth.^1.5;
requested = normal;
requested(O.critical) = critical(O.critical);
requested_by_cell = accumarray(O.cell, requested, [n_cells 1]);
scale = min(1, volume_m3 ./ max(requested_by_cell * dt_s, 1e-30));
actual = requested .* scale(O.cell);
removed = accumarray(O.cell, actual * dt_s, [n_cells 1]);
volume_m3 = max(volume_m3 - removed, 0);
outlet_m3_s = accumarray(O.cell, actual, [n_cells 1]);
outlet_depth_m(O.cell) = depth;

end
