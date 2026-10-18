function V = hp2d_soil_view(column, head_m, water_mm, root_depth_m)
%HP2D_SOIL_VIEW Three-layer view of the unsaturated soil slots (full-cell storage).
%   Slots 1-3 above the water table, with capacities (theta_s - theta_r) * unsat and
%   a pervious fraction of 1 (layered storage is per full cell). [N x 1] fields.
unsat = hp2d_unsaturated_thickness(column, head_m);
capacity = (column.theta_s(:, 1:3) - column.theta_r(:, 1:3)).*unsat(:, 1:3)*1000;
V.root_depth_m = root_depth_m(:) + zeros(size(unsat, 1), 1, 'like', unsat);
V.water_table_depth_m = hp2d_water_table_depth(column, head_m);
V.vadose_depth_m = sum(unsat(:, 1:3), 2);
V.pervious_fraction = ones(size(unsat, 1), 1, 'like', unsat);
V.near_surface_thickness_m = unsat(:, 1);
V.root_zone_thickness_m = unsat(:, 2);
V.transmission_thickness_m = unsat(:, 3);
V.capacity_mm = capacity;
V.total_capacity_mm = sum(capacity, 2);
V.water_mm = water_mm(:, 1:3);
end
