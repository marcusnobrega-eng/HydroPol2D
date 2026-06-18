function [Soil_Properties, LULC_Properties] = derive_layered_soil_profile(Soil_Properties, LULC_Properties, BC_States, elevation, idx_nan, options)
%DERIVE_LAYERED_SOIL_PROFILE Build robust soil-layer geometry and parameters.
%
% This helper prepares the layer geometry and layer-specific parameter maps
% for:
%   1) near-surface layer
%   2) root-zone layer
%   3) transmission-zone layer
%   4) groundwater zone
%
% Fallback rules:
%   - If root depth is shallower than the near-surface target, the explicit
%     root-zone layer collapses and roots are represented in the near-surface
%     layer.
%   - If bedrock is shallower than root depth, root depth is truncated to
%     available soil depth.
%   - If the water table is shallower than root depth, vadose root storage is
%     truncated at the water table and the remaining root zone is saturated.
%   - If water table reaches the surface, all vadose layers have zero
%     thickness and the groundwater zone starts at the surface.

if nargin < 6 || isempty(options)
    options = struct();
end

old_layers = struct();
if isfield(Soil_Properties, 'Layers') && isstruct(Soil_Properties.Layers)
    old_layers = Soil_Properties.Layers;
end

near_target = local_getopt(options, 'near_surface_depth_m', 0.10);
min_layer   = local_getopt(options, 'min_layer_thickness_m', 0.005);

soil_depth = max(Soil_Properties.Soil_Depth, 0);
soil_depth(idx_nan) = nan;

if ~isfield(LULC_Properties, 'root_depth_m') || isempty(LULC_Properties.root_depth_m)
    root_depth = ones(size(soil_depth), 'like', soil_depth);
    if isfield(LULC_Properties, 'idx_imp')
        root_depth(LULC_Properties.idx_imp) = 0;
    end
else
    root_depth = max(LULC_Properties.root_depth_m, 0);
end
root_depth(idx_nan) = nan;

zwt = elevation - BC_States.h_t;
zwt = max(zwt, 0);
zwt = min(zwt, soil_depth);
zwt(idx_nan) = nan;

vadose_depth = min(soil_depth, zwt);
vadose_depth = max(vadose_depth, 0);

near_thick = min(near_target, vadose_depth);

root_bottom = min(root_depth, vadose_depth);
root_bottom = max(root_bottom, 0);

root_thick = max(root_bottom - near_thick, 0);

trans_thick = max(vadose_depth - near_thick - root_thick, 0);

% Collapse tiny layers into the layer above to avoid numerically meaningless
% slivers. The near-surface layer may remain tiny only when total vadose
% depth itself is tiny.
idx_tiny_root = root_thick > 0 & root_thick < min_layer;
near_thick(idx_tiny_root) = near_thick(idx_tiny_root) + root_thick(idx_tiny_root);
root_thick(idx_tiny_root) = 0;

idx_tiny_trans = trans_thick > 0 & trans_thick < min_layer;
root_thick(idx_tiny_trans) = root_thick(idx_tiny_trans) + trans_thick(idx_tiny_trans);
trans_thick(idx_tiny_trans) = 0;

root_within_near_surface = root_depth <= near_target & root_depth > 0 & vadose_depth > 0;
root_truncated_by_bedrock = root_depth > soil_depth & soil_depth > 0;
root_truncated_by_water_table = root_depth > zwt & zwt >= 0;
water_table_at_surface = zwt <= 0;

layers = struct();
layers.near_surface_target_m = near_target;
layers.min_layer_thickness_m = min_layer;
layers.root_depth_m = root_depth;
layers.vadose_depth_m = vadose_depth;
layers.water_table_depth_m = zwt;
layers.near_surface_thickness_m = near_thick;
layers.root_zone_thickness_m = root_thick;
layers.transmission_thickness_m = trans_thick;
layers.groundwater_zone_top_depth_m = zwt;
layers.root_within_near_surface = root_within_near_surface;
layers.root_truncated_by_bedrock = root_truncated_by_bedrock;
layers.root_truncated_by_water_table = root_truncated_by_water_table;
layers.water_table_at_surface = water_table_at_surface;

layers.theta_sat_near_surface = Soil_Properties.theta_sat;
layers.theta_sat_root_zone = Soil_Properties.theta_sat;
layers.theta_sat_transmission = Soil_Properties.theta_sat;
layers.theta_r_near_surface = Soil_Properties.theta_r;
layers.theta_r_root_zone = Soil_Properties.theta_r;
layers.theta_r_transmission = Soil_Properties.theta_r;
layers.alpha_vg_near_surface = Soil_Properties.alpha_vg;
layers.alpha_vg_root_zone = Soil_Properties.alpha_vg;
layers.alpha_vg_transmission = Soil_Properties.alpha_vg;
layers.n_vg_near_surface = Soil_Properties.n_vg;
layers.n_vg_root_zone = Soil_Properties.n_vg;
layers.n_vg_transmission = Soil_Properties.n_vg;

layers.ksat_near_surface = Soil_Properties.ksat .* local_multiplier(Soil_Properties, 'Ks_multiplier_near_surface', soil_depth);
layers.ksat_root_zone = Soil_Properties.ksat .* local_multiplier(Soil_Properties, 'Ks_multiplier_root_zone', soil_depth);
layers.ksat_transmission = Soil_Properties.ksat .* local_multiplier(Soil_Properties, 'Ks_multiplier_transmission', soil_depth);

layers.near_surface_capacity_mm = near_thick .* max(Soil_Properties.theta_sat - Soil_Properties.theta_r, 0) .* 1000;
layers.root_zone_capacity_mm = root_thick .* max(Soil_Properties.theta_sat - Soil_Properties.theta_r, 0) .* 1000;
layers.transmission_capacity_mm = trans_thick .* max(Soil_Properties.theta_sat - Soil_Properties.theta_r, 0) .* 1000;
layers.total_vadose_capacity_mm = layers.near_surface_capacity_mm + layers.root_zone_capacity_mm + layers.transmission_capacity_mm;

storage_fields = { ...
    'near_surface_storage_mm', ...
    'root_zone_storage_mm', ...
    'transmission_storage_mm'};

for i = 1:numel(storage_fields)
    field_name = storage_fields{i};
    if isfield(old_layers, field_name) && ~isempty(old_layers.(field_name))
        layers.(field_name) = old_layers.(field_name);
    end
end

Soil_Properties.Layers = layers;
LULC_Properties.root_depth_m = root_depth;
end

function value = local_getopt(options, name, default_value)
if isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = default_value;
end
end

function multiplier = local_multiplier(Soil_Properties, field_name, template)
if isfield(Soil_Properties, field_name) && ~isempty(Soil_Properties.(field_name))
    multiplier = Soil_Properties.(field_name);
else
    multiplier = ones(size(template), 'like', template);
end
multiplier(~isfinite(multiplier)) = 1;
multiplier = max(multiplier, 0);
end
