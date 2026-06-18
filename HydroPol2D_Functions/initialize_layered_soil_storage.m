function Soil_Properties = initialize_layered_soil_storage(Soil_Properties, idx_nan)
%INITIALIZE_LAYERED_SOIL_STORAGE Create dynamic layer storages from theta_i.

if ~isfield(Soil_Properties, 'Layers') || ~isstruct(Soil_Properties.Layers)
    return
end

layers = Soil_Properties.Layers;

if isfield(layers, 'near_surface_storage_mm') && ...
        isfield(layers, 'root_zone_storage_mm') && ...
        isfield(layers, 'transmission_storage_mm')
    Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
    return
end

theta_avail = max(Soil_Properties.theta_i - Soil_Properties.theta_r, 0);

layers.near_surface_storage_mm = theta_avail .* layers.near_surface_thickness_m .* 1000;
layers.root_zone_storage_mm = theta_avail .* layers.root_zone_thickness_m .* 1000;
layers.transmission_storage_mm = theta_avail .* layers.transmission_thickness_m .* 1000;

layers.near_surface_storage_mm = min(max(layers.near_surface_storage_mm, 0), layers.near_surface_capacity_mm);
layers.root_zone_storage_mm = min(max(layers.root_zone_storage_mm, 0), layers.root_zone_capacity_mm);
layers.transmission_storage_mm = min(max(layers.transmission_storage_mm, 0), layers.transmission_capacity_mm);

layers.near_surface_storage_mm(idx_nan) = nan;
layers.root_zone_storage_mm(idx_nan) = nan;
layers.transmission_storage_mm(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
Soil_Properties.I_0 = Soil_Properties.I_t;
Soil_Properties.I_p = Soil_Properties.I_t;
end
