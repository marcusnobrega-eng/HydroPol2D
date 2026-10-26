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

% theta_i is the water content of the pervious part of the cell, so each store is
% f_p (theta_i - theta_r) times its thickness, as the capacities are (2026-10-26).
theta_avail = max(Soil_Properties.theta_i - Soil_Properties.theta_r, 0);
pore_mm = max(layers.pervious_fraction, 0) .* 1000;

layers.near_surface_storage_mm = theta_avail .* layers.near_surface_thickness_m .* pore_mm;
layers.root_zone_storage_mm = theta_avail .* layers.root_zone_thickness_m .* pore_mm;
layers.transmission_storage_mm = theta_avail .* layers.transmission_thickness_m .* pore_mm;

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
