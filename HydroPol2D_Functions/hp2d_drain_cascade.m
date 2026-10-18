function water_mm = hp2d_drain_cascade(column, head_m, water_mm, dt_s, substeps, l_vg)
%HP2D_DRAIN_CASCADE Unit-gradient drainage Kv*Kr(Se) down through the unsaturated slots.
%   Fluxes come from start-of-substep states (Jacobi). Each flux is limited by the
%   donor's water and by the free space of the next deeper slot with an unsaturated
%   part, where it arrives. The lowest unsaturated slot drains through
%   hp2d_water_table_flux instead.
if nargin < 5 || isempty(substeps), substeps = 1; end
if nargin < 6 || isempty(l_vg), l_vg = 0.5; end
water_in = water_mm;
unsat = hp2d_unsaturated_thickness(column, head_m);
has = unsat > 0;
lowest = hp2d_lowest_unsaturated(unsat);
dt_sub = dt_s/substeps;
[N, L] = size(water_mm);
for step = 1:substeps
    theta = min(max(column.theta_r + water_mm./max(unsat*1000, 1e-12), column.theta_r), column.theta_s);
    theta(~has) = column.theta_r(~has);
    space = max(column.theta_s - theta, 0).*unsat*1000;
    conductivity = column.kv_m_s.*hp2d_relative_conductivity(column, theta, l_vg);
    flux = min(conductivity*dt_sub*1000, water_mm);
    receiver_space = zeros(N, L, 'like', water_mm);
    carry = zeros(N, 1, 'like', water_mm);
    for k = L:-1:1
        receiver_space(:, k) = carry;
        carry(has(:, k)) = space(has(:, k), k);
    end
    flux = min(flux, receiver_space);
    flux(~(has & ~lowest)) = 0;
    inflow = zeros(N, L, 'like', water_mm);
    carry = zeros(N, 1, 'like', water_mm);
    for k = 1:L
        inflow(has(:, k), k) = carry(has(:, k));
        carry(has(:, k)) = flux(has(:, k), k);
    end
    water_mm = water_mm - flux + inflow;
end
water_mm(~column.valid, :) = water_in(~column.valid, :);
end
