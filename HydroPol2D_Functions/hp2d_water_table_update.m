function [head_m, water_mm, exfiltration_mm, unmet_mm] = hp2d_water_table_update(column, head_m, water_mm, net_mm)
%HP2D_WATER_TABLE_UPDATE Add (net_mm > 0) or remove (< 0) water at the water table (MMF).
%   Rise: fills theta_s - theta of the unsaturated parts, slot by slot upward; the water
%   there stays. What is left at the ground is exfiltration. Fall: lowers the table and
%   leaves theta_retained in the exposed thickness (releases theta_s - theta_retained per
%   metre); what is left at the column base is unmet removal.
%   A move continues to the next slot only once the current slot is full (rise) or
%   drained (fall); a partial move ends exactly, so round-off cannot skip a slot.
%   Same rules as HydroPol2D-Python modules.subsurface.water_table_update.
head_in = head_m(:);
water_in = water_mm;
valid = column.valid;
net = net_mm(:) + zeros(size(head_in), 'like', head_in);
net(~valid) = 0;
depth = hp2d_water_table_depth(column, head_in);
remaining = max(net, 0);
removal = max(-net, 0);
L = size(water_mm, 2);
for k = L:-1:1
    top = column.z_top_m(:, k); bot = column.z_bot_m(:, k);
    unsat = max(min(depth, bot) - top, 0);
    theta = min(max(column.theta_r(:, k) + water_mm(:, k)./max(unsat*1000, 1e-12), ...
        column.theta_r(:, k)), column.theta_s(:, k));
    cost = max(column.theta_s(:, k) - theta, 0)*1000;
    full = remaining >= cost.*unsat;
    delta = remaining./max(cost, 1e-300);
    delta(full) = unsat(full);
    delta = min(delta, unsat);
    delta(~(remaining > 0)) = 0;
    fraction = zeros(size(unsat), 'like', unsat);
    has = unsat > 0;
    fraction(has) = delta(has)./max(unsat(has), 1e-300);
    water_mm(:, k) = water_mm(:, k).*(1 - fraction);
    leftover = max(remaining - cost.*unsat, 0);
    remaining = leftover.*full;
    depth = depth - delta;
end
exfiltration_mm = remaining;
for k = 1:L
    top = column.z_top_m(:, k); bot = column.z_bot_m(:, k);
    saturated = max(bot - max(depth, top), 0);
    release = max(column.theta_s(:, k) - column.theta_retained(:, k), 1e-9)*1000;
    full = removal >= release.*saturated;
    delta = removal./release;
    delta(full) = saturated(full);
    delta = min(delta, saturated);
    delta(~(removal > 0)) = 0;
    water_mm(:, k) = water_mm(:, k) + (column.theta_retained(:, k) - column.theta_r(:, k)).*delta*1000;
    leftover = max(removal - release.*saturated, 0);
    removal = leftover.*full;
    depth = depth + delta;
end
head_m = column.ground_m - depth;
head_m(~valid) = head_in(~valid);
water_mm(~valid, :) = water_in(~valid, :);
exfiltration_mm(~valid) = 0;
unmet_mm = removal;
unmet_mm(~valid) = 0;
end
