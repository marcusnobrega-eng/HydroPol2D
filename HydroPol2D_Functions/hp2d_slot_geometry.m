function [z_top, z_bot] = hp2d_slot_geometry(root_m, soil_base_m, column_base_m, near_surface_m, deep_m)
%HP2D_SLOT_GEOMETRY Slot tops and bottoms (m below ground), each [N x (3 + numel(deep_m))].
%   Soil slots: near-surface [0, ltop], root zone [ltop, root], transmission [root, soil
%   base]. Deep slots follow deep_m from the soil base, cut at the column base. A slot
%   thinner than 5 mm is folded into the slot above. Inputs are [N x 1] (near_surface_m
%   may be scalar). Same rule as HydroPol2D-Python build_slot_geometry.
sliver = 0.005;
base = max(column_base_m(:), 0);
soil = min(max(soil_base_m(:), 0), base);
ltop = min(near_surface_m(:) + zeros(size(soil), 'like', soil), soil);
root = min(max(root_m(:), ltop), soil);
n_deep = numel(deep_m);
L = 3 + n_deep;
z_top = zeros(numel(soil), L, 'like', soil);
z_bot = z_top;
z_top(:, 1) = 0;    z_bot(:, 1) = ltop;
z_top(:, 2) = ltop; z_bot(:, 2) = root;
z_top(:, 3) = root; z_bot(:, 3) = soil;
cumulative = 0;
for k = 1:n_deep
    z_top(:, 3+k) = min(soil + cumulative, base);
    cumulative = cumulative + deep_m(k);
    z_bot(:, 3+k) = min(soil + cumulative, base);
end
for k = L:-1:2
    dz = z_bot(:, k) - z_top(:, k);
    tiny = dz > 0 & dz < sliver;
    z_bot(tiny, k-1) = z_bot(tiny, k);
    z_top(tiny, k) = z_bot(tiny, k);
end
end
