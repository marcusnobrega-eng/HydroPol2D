function water = hp2d_initial_slot_water(column, head_m, soil_theta_initial)
%HP2D_INITIAL_SLOT_WATER Initial unsaturated water (mm): given soil theta in slots 1-3
%   ([N x 3], optional) and hydrostatic equilibrium with the table in the deep slots.
unsat = hp2d_unsaturated_thickness(column, head_m);
depth = hp2d_water_table_depth(column, head_m);
middle = column.z_top_m + 0.5*unsat;
head_above_table = -max(depth - middle, 0);
theta = hp2d_theta_from_matric_head(head_above_table, column.theta_r, column.theta_s, ...
    column.alpha_m_inverse, column.n_vg);
if nargin > 2 && ~isempty(soil_theta_initial)
    theta(:, 1:3) = soil_theta_initial + zeros(size(theta(:, 1:3)), 'like', theta);
end
theta = min(max(theta, column.theta_r), column.theta_s);
water = (theta - column.theta_r).*unsat*1000;
water(~column.valid, :) = 0;
end
