function [Soil_Properties, actual_et_depth_mm] = extract_layered_et(Soil_Properties, et_demand_depth_mm, idx_extract, idx_nan)
%EXTRACT_LAYERED_ET Remove ET from root-accessible soil storage.

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

near_available = layers.near_surface_storage_mm .* near_access_fraction;
take_near = min(remaining, max(near_available, 0));
layers.near_surface_storage_mm = layers.near_surface_storage_mm - take_near;
remaining = remaining - take_near;

root_available = layers.root_zone_storage_mm;
take_root = min(remaining, max(root_available, 0));
layers.root_zone_storage_mm = layers.root_zone_storage_mm - take_root;
remaining = remaining - take_root;

actual_et_depth_mm = max(et_demand_depth_mm, 0) - max(remaining, 0);
actual_et_depth_mm(~idx_extract) = 0;
actual_et_depth_mm(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end
