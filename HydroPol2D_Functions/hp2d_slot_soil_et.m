function [water3, actual_mm] = hp2d_slot_soil_et(V, demand_mm, extract_mask, valid, feddes, theta_r, theta_s, alpha, n_vg, rate_mm_day)
%HP2D_SLOT_SOIL_ET ET from the rooted part of the near slot, then the root slot.
%   Without Feddes: the near slot (rooted fraction) first, then the root zone, up to the
%   demand. With Feddes (struct, or [] for off): uniform root density; each rooted layer
%   is asked for alpha(h)*demand*share, no compensation. Same rule as HydroPol2D-Python
%   soil.extract_layered_et.
allowed = logical(extract_mask(:)) & valid;
requested = max(demand_mm(:), 0);
remaining = requested;
remaining(~allowed) = 0;
root = max(V.root_depth_m, 0);
near_t = max(V.near_surface_thickness_m, 0);
access = ones(size(near_t), 'like', near_t);
partial = near_t > 0 & root > 0 & root < near_t;
access(partial) = root(partial)./near_t(partial);
access(root <= 0) = 0;
access = min(max(access, 0), 1);
water3 = V.water_mm;
near_available = max(water3(:, 1).*access, 0);
root_available = max(water3(:, 2), 0);
if isempty(feddes)
    take_near = min(remaining, near_available);
    remaining = remaining - take_near;
    take_root = min(remaining, root_available);
    remaining = remaining - take_root;
else
    near_rooted = near_t.*access;
    root_rooted = max(V.root_zone_thickness_m, 0);
    rooted = near_rooted + root_rooted;
    share_near = zeros(size(rooted), 'like', rooted);
    share_root = share_near;
    has = rooted > 0;
    share_near(has) = near_rooted(has)./max(rooted(has), 1e-12);
    share_root(has) = root_rooted(has)./max(rooted(has), 1e-12);
    pervious = max(V.pervious_fraction, 1e-12);
    theta_n = theta_r + max(water3(:, 1), 0)./max(pervious.*near_t*1000, 1e-12);
    theta_z = theta_r + max(water3(:, 2), 0)./max(pervious.*V.root_zone_thickness_m*1000, 1e-12);
    alpha_near = hp2d_feddes_alpha(hp2d_matric_head_from_theta(theta_n, theta_r, theta_s, alpha, n_vg), feddes, rate_mm_day);
    alpha_root = hp2d_feddes_alpha(hp2d_matric_head_from_theta(theta_z, theta_r, theta_s, alpha, n_vg), feddes, rate_mm_day);
    take_near = min(remaining.*share_near.*alpha_near, near_available);
    take_root = min(remaining.*share_root.*alpha_root, root_available);
    remaining = remaining - take_near - take_root;
end
water3(:, 1) = water3(:, 1) - take_near;
water3(:, 2) = water3(:, 2) - take_root;
water3 = min(max(water3, 0), V.capacity_mm);
actual_mm = requested - max(remaining, 0);
actual_mm(~allowed) = 0;
end
