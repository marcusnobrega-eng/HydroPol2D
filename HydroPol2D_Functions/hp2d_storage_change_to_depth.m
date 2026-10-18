function [delta_mm, water_mm] = hp2d_storage_change_to_depth(column, head_m, water_mm, target_depth_m)
%HP2D_STORAGE_CHANGE_TO_DEPTH Move each column's table to target_depth_m (MMF rules).
%   Returns the storage change (mm over the column, [N x 1]) and the new water [N x L].
water_in = water_mm;
depth = hp2d_water_table_depth(column, head_m);
target = min(max(target_depth_m(:), 0), column.base_depth_m);
unsat = max(min(depth, column.z_bot_m) - column.z_top_m, 0);
theta = min(max(column.theta_r + water_mm./max(unsat*1000, 1e-12), column.theta_r), column.theta_s);
rise = max(min(depth, column.z_bot_m) - max(target, column.z_top_m), 0).*(target < depth);
fall = max(min(target, column.z_bot_m) - max(depth, column.z_top_m), 0).*(target > depth);
fraction = zeros(size(unsat), 'like', unsat);
has = unsat > 0;
fraction(has) = rise(has)./max(unsat(has), 1e-300);
water_mm = water_mm.*(1 - fraction) + (column.theta_retained - column.theta_r).*fall*1000;
delta_mm = sum(max(column.theta_s - theta, 0).*rise*1000 ...
    - max(column.theta_s - column.theta_retained, 0).*fall*1000, 2);
delta_mm(~column.valid) = 0;
water_mm(~column.valid, :) = water_in(~column.valid, :);
end
