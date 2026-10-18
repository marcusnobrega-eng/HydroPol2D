function [water3, accepted_mm] = hp2d_add_slot_infiltration(V, depth_mm, valid)
%HP2D_ADD_SLOT_INFILTRATION Place infiltrated water top-down in the three soil slots.
requested = max(depth_mm(:), 0);
remaining = requested;
remaining(~valid) = 0;
water3 = V.water_mm;
for k = 1:3
    add = min(remaining, max(V.capacity_mm(:, k) - water3(:, k), 0));
    water3(:, k) = water3(:, k) + add;
    remaining = remaining - add;
end
water3 = min(max(water3, 0), V.capacity_mm);
accepted_mm = requested - max(remaining, 0);
accepted_mm(~valid) = 0;
end
