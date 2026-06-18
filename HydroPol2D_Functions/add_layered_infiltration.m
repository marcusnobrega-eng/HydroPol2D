function [Soil_Properties, accepted_depth_mm] = add_layered_infiltration(Soil_Properties, infiltration_depth_mm, idx_nan)
%ADD_LAYERED_INFILTRATION Add infiltrated water top-down into soil layers.

layers = Soil_Properties.Layers;
remaining = max(infiltration_depth_mm, 0);
remaining(idx_nan) = 0;

add_near = min(remaining, max(layers.near_surface_capacity_mm - layers.near_surface_storage_mm, 0));
layers.near_surface_storage_mm = layers.near_surface_storage_mm + add_near;
remaining = remaining - add_near;

add_root = min(remaining, max(layers.root_zone_capacity_mm - layers.root_zone_storage_mm, 0));
layers.root_zone_storage_mm = layers.root_zone_storage_mm + add_root;
remaining = remaining - add_root;

add_trans = min(remaining, max(layers.transmission_capacity_mm - layers.transmission_storage_mm, 0));
layers.transmission_storage_mm = layers.transmission_storage_mm + add_trans;
remaining = remaining - add_trans;

accepted_depth_mm = max(infiltration_depth_mm, 0) - max(remaining, 0);
accepted_depth_mm(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end
