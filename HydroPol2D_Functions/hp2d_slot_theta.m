function theta = hp2d_slot_theta(column, head_m, water_mm)
%HP2D_SLOT_THETA Water content of the unsaturated part of each slot (theta_r where empty).
unsat = hp2d_unsaturated_thickness(column, head_m);
theta = column.theta_r + max(water_mm, 0)./max(unsat*1000, 1e-12);
theta = min(max(theta, column.theta_r), column.theta_s);
theta(unsat <= 0) = column.theta_r(unsat <= 0);
end
