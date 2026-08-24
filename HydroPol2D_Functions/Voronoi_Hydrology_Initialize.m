function state = Voronoi_Hydrology_Initialize(mesh, hydrology, groundwater_head)
%VORONOI_HYDROLOGY_INITIALIZE Build fractional polygon hydrology columns.

arguments
    mesh struct
    hydrology struct
    groundwater_head double
end

n = mesh.n_cells;
if ~isfield(hydrology, 'hru_fraction') || isempty(hydrology.hru_fraction)
    hydrology.hru_fraction = ones(n, 1);
end
fraction = double(hydrology.hru_fraction);
if isvector(fraction) && numel(fraction) == n, fraction = fraction(:); end
assert(size(fraction,1) == n && all(isfinite(fraction(:)) & fraction(:) >= 0), ...
    'HydroPol2D:InvalidVoronoiHydrology', 'hru_fraction must be a nonnegative n-cell-by-n-HRU matrix.');
assert(all(abs(sum(fraction,2) - 1) <= 1e-10), ...
    'HydroPol2D:InvalidVoronoiHydrology', 'HRU fractions must sum to one in every polygon.');
h = size(fraction,2);

state.enabled = true;
state.fraction = fraction;
state.area_m2 = mesh.cell_area(:) .* fraction;
state.pervious_fraction = field(hydrology, 'pervious_fraction', 1, n, h, 0, 1);
state.lai = field(hydrology, 'lai', 0, n, h, 0, inf);
state.crop_coefficient = field(hydrology, 'crop_coefficient', 1, n, h, 0, inf);
state.canopy_capacity_per_lai_m = field(hydrology, 'canopy_capacity_per_lai_m', 2e-4, n, h, 0, inf);
state.soil_depth_m = field(hydrology, 'soil_depth_m', 2, n, h, 0, inf);
state.root_depth_m = min(field(hydrology, 'root_depth_m', 1, n, h, 0, inf), state.soil_depth_m);
state.theta_r = field(hydrology, 'theta_r', 0.05, n, h, 0, 1);
state.theta_sat = field(hydrology, 'theta_sat', 0.45, n, h, 0, 1);
assert(all(state.theta_sat(:) > state.theta_r(:)), ...
    'HydroPol2D:InvalidVoronoiHydrology', 'theta_sat must exceed theta_r.');
state.alpha_vg_1_m = field(hydrology, 'alpha_vg_1_m', 3.6, n, h, eps, inf);
state.n_vg = field(hydrology, 'n_vg', 1.56, n, h, 1 + eps, inf);
state.l_vg = field(hydrology, 'l_vg', 0.5, n, h, -inf, inf);
state.ksat_m_s = field(hydrology, 'ksat_m_s', 1e-6, n, h, 0, inf);
state.surface_layer_depth_m = field(hydrology, 'surface_layer_depth_m', 0.10, n, h, eps, inf);
state.infiltration_head_cap_m = field(hydrology, 'infiltration_head_cap_m', 1, n, h, 0, inf);
state.surface_conductivity_weight = field(hydrology, 'surface_conductivity_weight', 0.5, n, h, 0, 1);
state.wilting_saturation = field(hydrology, 'wilting_saturation', 0.10, n, h, 0, 1);
state.capillary_extinction_depth_m = field(hydrology, 'capillary_extinction_depth_m', 2, n, h, eps, inf);
state.capillary_rise_enabled = logical(getfield_default(hydrology, 'capillary_rise_enabled', true));

state.canopy_storage_m = field(hydrology, 'initial_canopy_storage_m', 0, n, h, 0, inf);
initial_saturation = field(hydrology, 'initial_soil_saturation', 0.5, n, h, 0, 1);
state = update_geometry(state, mesh.surface_bed(:), groundwater_head(:));
state.near_storage_m = initial_saturation .* state.near_capacity_m;
state.root_storage_m = initial_saturation .* state.root_capacity_m;
state.transmission_storage_m = initial_saturation .* state.transmission_capacity_m;
state.cumulative_infiltration_m = zeros(n,h);
state.cumulative_recharge_m = zeros(n,h);
state.cumulative_actual_et_m = zeros(n,h);
state.cumulative_surface_evaporation_m = zeros(n,1);
state.last_infiltration_rate_m_s = zeros(n,h);
state.last_recharge_rate_m_s = zeros(n,h);
state.last_capillary_rate_m_s = zeros(n,h);
state.last_actual_et_rate_m_s = zeros(n,h);
state.last_canopy_evaporation_rate_m_s = zeros(n,h);
state.last_surface_evaporation_rate_m_s = zeros(n,1);
end

function state = update_geometry(state, surface_bed, groundwater_head)
zwt = max(surface_bed - groundwater_head, 0);
zwt = min(zwt, max(state.soil_depth_m, 0));
near_thickness = min(state.surface_layer_depth_m, zwt);
root_bottom = min(state.root_depth_m, zwt);
root_thickness = max(root_bottom - near_thickness, 0);
transmission_thickness = max(zwt - near_thickness - root_thickness, 0);
water_capacity = state.pervious_fraction .* max(state.theta_sat - state.theta_r, 0);
state.water_table_depth_m = zwt;
state.near_thickness_m = near_thickness;
state.root_thickness_m = root_thickness;
state.transmission_thickness_m = transmission_thickness;
state.near_capacity_m = water_capacity .* near_thickness;
state.root_capacity_m = water_capacity .* root_thickness;
state.transmission_capacity_m = water_capacity .* transmission_thickness;
end

function value = field(source, name, default_value, n, h, lower, upper)
value = getfield_default(source, name, default_value);
value = double(value);
if isscalar(value)
    value = repmat(value,n,h);
elseif isequal(size(value),[n 1]) && h > 1
    value = repmat(value,1,h);
elseif isvector(value) && numel(value) == h
    value = repmat(reshape(value,1,h),n,1);
end
assert(isequal(size(value),[n h]) && all(isfinite(value(:))) && ...
    all(value(:) >= lower & value(:) <= upper), ...
    'HydroPol2D:InvalidVoronoiHydrology', '%s has invalid dimensions or values.', name);
end

function value = getfield_default(source, name, default_value)
if isfield(source,name) && ~isempty(source.(name)), value=source.(name); else, value=default_value; end
end
