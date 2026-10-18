function [capillary_rate, Soil_Properties] = simulate_layered_capillary_rise(Soil_Properties, BC_States, elevation, dt_s, idx_nan, max_capillary_depth_m, field_capacity_head_m)
%SIMULATE_LAYERED_CAPILLARY_RISE Refill root-accessible storage with VG-Mualem flux.
%
% capillary_rate is positive upward from groundwater to vadose storage [m/s].

dt_s = max(dt_s, eps);
layers = Soil_Properties.Layers;
template = Soil_Properties.I_t;
capillary_depth_mm = zeros(size(template), 'like', template);

if nargin < 6 || isempty(max_capillary_depth_m)
    max_capillary_depth_m = inf(size(template), 'like', template);
end
if nargin < 7 || isempty(field_capacity_head_m), field_capacity_head_m = 3.3; end

remaining_mm = max(max_capillary_depth_m, 0) .* 1000;
remaining_mm(idx_nan) = 0;

water_table_depth_m = max(elevation - BC_States.h_t, 0);
water_table_depth_m(idx_nan) = nan;

alpha=max(abs(Soil_Properties.alpha_vg),eps);
n_value=max(Soil_Properties.n_vg,1+eps);
root_depth_m=max(layers.root_depth_m,0);
gap=max(water_table_depth_m-root_depth_m,0);
use_root=layers.root_zone_capacity_mm>0;
rooted_capacity=layers.root_zone_capacity_mm;
rooted_capacity(~use_root)=layers.near_surface_capacity_mm(~use_root).*(root_depth_m(~use_root)>0);
rooted_storage=layers.root_zone_storage_mm;
rooted_storage(~use_root)=layers.near_surface_storage_mm(~use_root);
capillary_depth_mm=1000.*hp2d_vg_mualem_capillary_transfer( ...
    layers.ksat_root_zone./1000./3600,alpha,n_value,0.5,gap, ...
    field_capacity_head_m,rooted_capacity./1000,rooted_storage./1000, ...
    remaining_mm./1000,dt_s);
layers.root_zone_storage_mm=layers.root_zone_storage_mm+capillary_depth_mm.*use_root;
layers.near_surface_storage_mm=layers.near_surface_storage_mm+capillary_depth_mm.*~use_root;

capillary_depth_mm(idx_nan) = nan;
capillary_rate = capillary_depth_mm ./ 1000 ./ dt_s;
capillary_rate(idx_nan) = nan;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end
