function Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan)
%SYNC_LAYERED_SOIL_STORAGE Keep aggregate I_t equal to layer storages.

if ~isfield(Soil_Properties, 'Layers') || ~isstruct(Soil_Properties.Layers)
    return
end

layers = Soil_Properties.Layers;
required = {'near_surface_storage_mm','root_zone_storage_mm','transmission_storage_mm'};
for i = 1:numel(required)
    if ~isfield(layers, required{i}) || isempty(layers.(required{i}))
        return
    end
end

layers.near_surface_storage_mm = min(max(layers.near_surface_storage_mm, 0), layers.near_surface_capacity_mm);
layers.root_zone_storage_mm = min(max(layers.root_zone_storage_mm, 0), layers.root_zone_capacity_mm);
layers.transmission_storage_mm = min(max(layers.transmission_storage_mm, 0), layers.transmission_capacity_mm);

layers.near_surface_storage_mm(idx_nan) = nan;
layers.root_zone_storage_mm(idx_nan) = nan;
layers.transmission_storage_mm(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties.I_t = layers.near_surface_storage_mm + ...
    layers.root_zone_storage_mm + ...
    layers.transmission_storage_mm;
Soil_Properties.I_t(idx_nan) = nan;
end
