function remaining_mm = layered_soil_remaining_capacity_mm(Soil_Properties, idx_nan)
%LAYERED_SOIL_REMAINING_CAPACITY_MM Total available vadose storage in layers.

layers = Soil_Properties.Layers;
remaining_mm = ...
    max(layers.near_surface_capacity_mm - layers.near_surface_storage_mm, 0) + ...
    max(layers.root_zone_capacity_mm - layers.root_zone_storage_mm, 0) + ...
    max(layers.transmission_capacity_mm - layers.transmission_storage_mm, 0);
remaining_mm(idx_nan) = nan;
end
