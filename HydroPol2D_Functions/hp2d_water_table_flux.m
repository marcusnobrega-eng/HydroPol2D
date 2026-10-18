function [water_mm, recharge_mm, capillary_mm] = hp2d_water_table_flux(column, head_m, water_mm, dt_s, l_vg)
%HP2D_WATER_TABLE_FLUX Darcy flux between the lowest unsaturated part and the water table.
%   q = K(theta)*(psi/dz + 1), positive downward, with dz = max(unsat/2, 1 mm) of that
%   slot and psi = 0 at the table. Downward flow is recharge, limited by the slot's
%   water; upward flow is capillary rise, limited by its free space and the drainable
%   water below. The caller applies recharge - capillary with hp2d_water_table_update.
if nargin < 5 || isempty(l_vg), l_vg = 0.5; end
water_in = water_mm;
unsat = hp2d_unsaturated_thickness(column, head_m);
lowest = hp2d_lowest_unsaturated(unsat);
theta = hp2d_slot_theta(column, head_m, water_mm);
psi = hp2d_matric_head_from_theta(theta, column.theta_r, column.theta_s, ...
    column.alpha_m_inverse, column.n_vg);
conductivity = column.kv_m_s.*hp2d_relative_conductivity(column, theta, l_vg);
distance = max(0.5*unsat, 1e-3);
depth_mm = conductivity.*(psi./distance + 1)*dt_s*1000;
space = max(column.theta_s - theta, 0).*unsat*1000;
down = min(max(depth_mm, 0), water_mm).*lowest;
up = min(max(-depth_mm, 0), space).*lowest;
available = hp2d_drainable_storage_mm(column, head_m);
up_total = sum(up, 2);
scale = ones(size(up_total), 'like', up_total);
over = up_total > available;
scale(over) = available(over)./max(up_total(over), 1e-300);
up = up.*scale;
water_mm = water_mm - down + up;
water_mm(~column.valid, :) = water_in(~column.valid, :);
recharge_mm = sum(down, 2);
capillary_mm = sum(up, 2);
recharge_mm(~column.valid) = 0;
capillary_mm(~column.valid) = 0;
end
