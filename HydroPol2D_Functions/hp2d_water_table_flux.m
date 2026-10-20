function [water_mm, recharge_mm, capillary_mm] = hp2d_water_table_flux(column, head_m, water_mm, dt_s, l_vg, implicit)
%HP2D_WATER_TABLE_FLUX Darcy flux between the lowest unsaturated part and the water table.
%   q = K(theta)*(psi/dz + 1), positive downward, with dz = max(unsat/2, 1 mm) of that
%   slot and psi = 0 at the table. Downward flow is recharge, limited by the slot's
%   water; upward flow is capillary rise, limited by the drainable water below and by
%   the slot's deficit below theta_retained (field capacity, or the rock Sy). A falling
%   table leaves theta_retained behind, so capillary rise refills only what roots took
%   below it (HBV-style, Lindstrom et al. 1997); filling above it circulated water
%   between the flux and the table moves. The caller applies recharge - capillary with
%   hp2d_water_table_update.
%
%   implicit (default true, as HydroPol2D-Python) takes the flux at the end of the step,
%   linearised (one backward-Euler step per cell, no iteration):
%   q = q0 / (1 + dt*(dq/dW - dq/dh/Sy)), keeping only the damping part. dq/dh follows
%   the move hp2d_water_table_update makes (a rise carries the risen slice's water, so
%   the gap keeps its theta; a fall exposes theta_retained), and Sy is that move's cost
%   (theta_s - theta of the gap for recharge, the fall release for capillary rise).
%   Python takes the derivatives with jax.jvp; here they are fourth-order central
%   differences (agreeing to about 1e-9 relative away from the clips). A thin
%   unsaturated gap above a shallow table makes the explicit flux stiff.
if nargin < 5 || isempty(l_vg), l_vg = 0.5; end
if nargin < 6 || isempty(implicit), implicit = true; end
water_in = water_mm;
unsat = hp2d_unsaturated_thickness(column, head_m);
lowest = hp2d_lowest_unsaturated(unsat);
q = darcy_flux_to_table(column, head_m, water_mm, l_vg);
if implicit && any(lowest(:))
    % Only the lowest unsaturated slot's flux reaches the table, so the derivatives are
    % taken on a one-slot column of it (cells x 1; the same arithmetic, L times less work).
    head_m = head_m(:);
    [rows, slot] = find(lowest);
    idx = sub2ind(size(lowest), rows, slot);
    one = struct('z_top_m', column.z_top_m(idx), 'z_bot_m', column.z_bot_m(idx), ...
        'theta_r', column.theta_r(idx), 'theta_s', column.theta_s(idx), ...
        'alpha_m_inverse', column.alpha_m_inverse(idx), 'n_vg', column.n_vg(idx), ...
        'kv_m_s', column.kv_m_s(idx), 'ground_m', column.ground_m(rows), ...
        'base_depth_m', column.base_depth_m(rows));
    h = head_m(rows);
    w = water_mm(idx);
    q0 = q(idx);
    % dq/dW of the slot's water (per mm).
    step_w = 1e-5*max(abs(w), 1);
    dq_dwater = central_difference(@(k) darcy_flux_to_table(one, h, w + k*step_w, l_vg), step_w);
    % dq/dh along the table move (per m of head).
    theta_now = hp2d_slot_theta(one, h, w);
    downward = q0 > 0;
    moved_theta = column.theta_retained(idx);
    moved_theta(downward) = theta_now(downward);
    water_per_metre = -(moved_theta - one.theta_r)*1000;
    step_h = 1e-6;
    dq_dhead = central_difference(@(k) darcy_flux_to_table(one, h + k*step_h, ...
        w + k*water_per_metre*step_h, l_vg), step_h);
    specific_yield = hp2d_effective_specific_yield(column, head_m);
    specific_yield = specific_yield(rows);
    specific_yield(downward) = max(one.theta_s(downward) - theta_now(downward), 0);
    specific_yield = max(specific_yield, 1e-4);
    damping = dq_dwater*1000 - dq_dhead./specific_yield;
    damping(isnan(damping)) = 0;
    damping(damping == Inf) = 1e30;
    damping(damping == -Inf) = 0;
    % A saturated gap drains at Ks undamped, as in Python: the Mualem slope is infinite
    % there, and the table then rises through the gap at no cost.
    damping(theta_now >= one.theta_s) = 0;
    q(idx) = q0./(1 + dt_s*max(damping, 0));
end
depth_mm = q*dt_s*1000;
theta = hp2d_slot_theta(column, head_m, water_mm);
deficit = max(column.theta_retained - theta, 0).*unsat*1000;
down = min(max(depth_mm, 0), water_mm).*lowest;
up = min(max(-depth_mm, 0), deficit).*lowest;
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

function d = central_difference(f, h)
% Fourth-order central difference of f(k) = q(x + k h) at k = 0.
d = (8*(f(1) - f(-1)) - (f(2) - f(-2)))./(12*h);
end

function q = darcy_flux_to_table(column, head_m, water_mm, l_vg)
% Per-slot Darcy flux q = K(theta) (psi / dz + 1) toward the table (m/s, down > 0).
unsat = hp2d_unsaturated_thickness(column, head_m);
theta = hp2d_slot_theta(column, head_m, water_mm);
psi = hp2d_matric_head_from_theta(theta, column.theta_r, column.theta_s, ...
    column.alpha_m_inverse, column.n_vg);
conductivity = column.kv_m_s.*hp2d_relative_conductivity(column, theta, l_vg);
distance = max(0.5*unsat, 1e-3);
q = conductivity.*(psi./distance + 1);
end
