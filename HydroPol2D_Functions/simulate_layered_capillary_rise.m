function [capillary_rate, Soil_Properties] = simulate_layered_capillary_rise(Soil_Properties, BC_States, elevation, dt_s, idx_nan, max_capillary_depth_m)
%SIMULATE_LAYERED_CAPILLARY_RISE Move shallow groundwater upward to vadose layers.
%
% capillary_rate is positive upward from groundwater to vadose storage [m/s].

dt_s = max(dt_s, eps);
layers = Soil_Properties.Layers;
template = Soil_Properties.I_t;
capillary_depth_mm = zeros(size(template), 'like', template);

if nargin < 6 || isempty(max_capillary_depth_m)
    max_capillary_depth_m = inf(size(template), 'like', template);
end

remaining_mm = max(max_capillary_depth_m, 0) .* 1000;
remaining_mm(idx_nan) = 0;

water_table_depth_m = max(elevation - BC_States.h_t, 0);
water_table_depth_m(idx_nan) = nan;

capillary_extinction_depth_m = 2.0;

[layers.transmission_storage_mm, moved_trans] = local_fill_layer( ...
    layers.transmission_storage_mm, layers.transmission_capacity_mm, ...
    layers.ksat_transmission, water_table_depth_m, capillary_extinction_depth_m, ...
    dt_s, remaining_mm);
remaining_mm = remaining_mm - moved_trans;
capillary_depth_mm = capillary_depth_mm + moved_trans;

[layers.root_zone_storage_mm, moved_root] = local_fill_layer( ...
    layers.root_zone_storage_mm, layers.root_zone_capacity_mm, ...
    layers.ksat_root_zone, water_table_depth_m, capillary_extinction_depth_m, ...
    dt_s, remaining_mm);
remaining_mm = remaining_mm - moved_root;
capillary_depth_mm = capillary_depth_mm + moved_root;

idx_no_lower_vadose = layers.transmission_capacity_mm <= 0 & layers.root_zone_capacity_mm <= 0;
[near_candidate, moved_near] = local_fill_layer( ...
    layers.near_surface_storage_mm, layers.near_surface_capacity_mm, ...
    layers.ksat_near_surface, water_table_depth_m, capillary_extinction_depth_m, ...
    dt_s, remaining_mm);
moved_near(~idx_no_lower_vadose) = 0;
layers.near_surface_storage_mm = layers.near_surface_storage_mm + moved_near;
layers.near_surface_storage_mm(idx_no_lower_vadose) = near_candidate(idx_no_lower_vadose);
capillary_depth_mm = capillary_depth_mm + moved_near;

capillary_depth_mm(idx_nan) = nan;
capillary_rate = capillary_depth_mm ./ 1000 ./ dt_s;
capillary_rate(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end

function [storage_mm, moved_mm] = local_fill_layer(storage_mm, capacity_mm, ksat_mm_h, water_table_depth_m, extinction_depth_m, dt_s, remaining_mm)
deficit_mm = max(capacity_mm - storage_mm, 0);
storage_fraction = storage_mm ./ max(capacity_mm, eps);
dryness = max(1 - storage_fraction, 0);
proximity = max(1 - water_table_depth_m ./ max(extinction_depth_m, eps), 0);
potential_mm = max(ksat_mm_h, 0) .* proximity .* dryness .* dt_s ./ 3600;
moved_mm = min(deficit_mm, potential_mm);
moved_mm = min(moved_mm, max(remaining_mm, 0));
moved_mm(~isfinite(moved_mm)) = 0;
storage_mm = storage_mm + moved_mm;
end
