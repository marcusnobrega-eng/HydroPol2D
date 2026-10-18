function [Soil_Properties, actual_et_depth_mm] = extract_layered_et(Soil_Properties, et_demand_depth_mm, idx_extract, idx_nan, potential_rate_mm_day)
%EXTRACT_LAYERED_ET Remove ET from root-accessible soil storage.
%   Without Feddes: the rooted fraction of the near-surface layer first, then the root
%   zone. With Soil_Properties.Feddes (struct from hp2d_layered_subsurface_config):
%   uniform root density, each rooted layer asked for alpha(h)*demand*share with the
%   Feddes factor of its matric head; no compensation (HydroPol2D-Python
%   soil.extract_layered_et). potential_rate_mm_day sets the demand-dependent h3.
if nargin < 5 || isempty(potential_rate_mm_day), potential_rate_mm_day = 0; end

layers = Soil_Properties.Layers;
remaining = max(et_demand_depth_mm, 0);
remaining(~idx_extract) = 0;
remaining(idx_nan) = 0;

root_depth = max(layers.root_depth_m, 0);
near_thick = max(layers.near_surface_thickness_m, 0);

near_access_fraction = ones(size(near_thick), 'like', near_thick);
idx_partial_near = near_thick > 0 & root_depth > 0 & root_depth < near_thick;
near_access_fraction(idx_partial_near) = root_depth(idx_partial_near) ./ max(near_thick(idx_partial_near), eps);
near_access_fraction(root_depth <= 0) = 0;
near_access_fraction = min(max(near_access_fraction, 0), 1);

near_available = max(layers.near_surface_storage_mm .* near_access_fraction, 0);
root_available = max(layers.root_zone_storage_mm, 0);
use_feddes = isfield(Soil_Properties, 'Feddes') && ~isempty(Soil_Properties.Feddes);
if ~use_feddes
    take_near = min(remaining, near_available);
    remaining = remaining - take_near;
    take_root = min(remaining, root_available);
    remaining = remaining - take_root;
else
    near_rooted = near_thick .* near_access_fraction;
    root_rooted = max(layers.root_zone_thickness_m, 0);
    rooted = near_rooted + root_rooted;
    share_near = zeros(size(rooted), 'like', rooted);
    share_root = share_near;
    has_roots = rooted > 0;
    share_near(has_roots) = near_rooted(has_roots) ./ max(rooted(has_roots), 1e-12);
    share_root(has_roots) = root_rooted(has_roots) ./ max(rooted(has_roots), 1e-12);
    pervious = max(layers.pervious_fraction, 1e-12);
    theta_r = Soil_Properties.theta_r; theta_s = Soil_Properties.theta_sat;
    alpha_vg = Soil_Properties.alpha_vg; n_vg = Soil_Properties.n_vg;
    theta_near = theta_r + max(layers.near_surface_storage_mm, 0) ./ max(pervious .* near_thick * 1000, 1e-12);
    theta_root = theta_r + max(layers.root_zone_storage_mm, 0) ./ ...
        max(pervious .* layers.root_zone_thickness_m * 1000, 1e-12);
    alpha_near = hp2d_feddes_alpha(hp2d_matric_head_from_theta(theta_near, theta_r, theta_s, ...
        alpha_vg, n_vg), Soil_Properties.Feddes, potential_rate_mm_day);
    alpha_root = hp2d_feddes_alpha(hp2d_matric_head_from_theta(theta_root, theta_r, theta_s, ...
        alpha_vg, n_vg), Soil_Properties.Feddes, potential_rate_mm_day);
    take_near = min(remaining .* share_near .* alpha_near, near_available);
    take_root = min(remaining .* share_root .* alpha_root, root_available);
    take_near(~isfinite(take_near)) = 0;
    take_root(~isfinite(take_root)) = 0;
    remaining = remaining - take_near - take_root;
end
layers.near_surface_storage_mm = layers.near_surface_storage_mm - take_near;
layers.root_zone_storage_mm = layers.root_zone_storage_mm - take_root;

actual_et_depth_mm = max(et_demand_depth_mm, 0) - max(remaining, 0);
actual_et_depth_mm(~idx_extract) = 0;
actual_et_depth_mm(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end
