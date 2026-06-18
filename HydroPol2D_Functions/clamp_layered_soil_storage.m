function [Soil_Properties, excess_mm] = clamp_layered_soil_storage(Soil_Properties, idx_nan)
%CLAMP_LAYERED_SOIL_STORAGE Enforce current layer capacities after WT motion.

layers = Soil_Properties.Layers;

excess_near = max(layers.near_surface_storage_mm - layers.near_surface_capacity_mm, 0);
excess_root = max(layers.root_zone_storage_mm - layers.root_zone_capacity_mm, 0);
excess_trans = max(layers.transmission_storage_mm - layers.transmission_capacity_mm, 0);

layers.near_surface_storage_mm = min(max(layers.near_surface_storage_mm, 0), layers.near_surface_capacity_mm);
layers.root_zone_storage_mm = min(max(layers.root_zone_storage_mm, 0), layers.root_zone_capacity_mm);
layers.transmission_storage_mm = min(max(layers.transmission_storage_mm, 0), layers.transmission_capacity_mm);

excess_mm = excess_near + excess_root + excess_trans;
excess_mm(idx_nan) = nan;

layers.saturation_excess_mm = excess_mm;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end
