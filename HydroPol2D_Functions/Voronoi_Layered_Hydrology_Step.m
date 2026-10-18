function [surface_volume, state, groundwater, diagnostics] = Voronoi_Layered_Hydrology_Step( ...
    mesh, surface_volume, state, groundwater, precipitation_m_s, potential_et_m_s, open_water_evaporation_m_s, dt_s)
%VORONOI_LAYERED_HYDROLOGY_STEP Canopy, infiltration, ET (Feddes), drainage and table flux per HRU.
%   Layered-subsurface counterpart of Voronoi_Hydrology_Step (state.layered from
%   hp2d_voronoi_layered_initialize). Same order as HydroPol2D-Python
%   voronoi_layered.layered_hydrology_step; the net vadose-table exchange of every HRU
%   goes to the cell's pending exchange, applied by Voronoi_Layered_Groundwater_Advance.
n = mesh.n_cells;
L = state.layered;
host = L.host;
area = L.area_m2;
column = L.column;
precipitation_m_s = cell_vector(precipitation_m_s, n);
potential_et_m_s = cell_vector(potential_et_m_s, n);
open_water_evaporation_m_s = cell_vector(open_water_evaporation_m_s, n);
dt_s = max(double(dt_s), eps);
head = groundwater_head(groundwater, L);
head_hru = head(host);
initial_storage = layered_storage(surface_volume, state, L, groundwater, head);

rain = precipitation_m_s(host)*dt_s;
pet = potential_et_m_s(host)*dt_s;
capacity = state.canopy_capacity_per_lai_m(:).*state.lai(:);
canopy = state.canopy_storage_m(:);
beta = zeros(size(capacity));
active = capacity > 0;
beta(active) = canopy(active)./capacity(active);
canopy_evaporation = min(max(beta, 0).*pet, canopy + rain);
canopy_evaporation(~active) = 0;
provisional = canopy + rain - canopy_evaporation;
throughfall = max(provisional - capacity, 0);
state.canopy_storage_m(:) = min(provisional, capacity);
surface_volume = surface_volume + accumarray(host, throughfall.*area, [n 1]);

water = L.water_mm;
pervious = state.pervious_fraction(:);
root = state.root_depth_m(:);
valid = true(numel(host), 1);
V = hp2d_soil_view(column, head_hru, water, root);
available_mm = surface_volume(host)./mesh.cell_area(host)*1000;
rate = hp2d_layered_infiltration_potential(V, available_mm, column.theta_r(:, 1), ...
    column.theta_s(:, 1), column.alpha_m_inverse(:, 1), column.n_vg(:, 1), ...
    state.ksat_m_s(:)*3.6e6, dt_s, state.l_vg(:), state.surface_layer_depth_m(:), ...
    state.infiltration_head_cap_m(:), state.surface_conductivity_weight(:), valid);
rate(~isfinite(rate)) = 0;
request = min(rate*dt_s/3600.*pervious, available_mm);
requested = accumarray(host, request/1000.*area, [n 1]);
scale = min(1, surface_volume./max(requested, 1e-300));
request = request.*scale(host);
[water3, accepted] = hp2d_add_slot_infiltration(V, request, valid);
water(:, 1:3) = water3;
infiltration = accepted/1000;
surface_volume = surface_volume - accumarray(host, infiltration.*area, [n 1]);

ponded = surface_volume > 0;
surface_evaporation = min(surface_volume, open_water_evaporation_m_s*dt_s.*mesh.cell_area(:));
surface_evaporation(~ponded) = 0;
surface_volume = surface_volume - surface_evaporation;
demand_mm = pet.*state.crop_coefficient(:)*1000.*pervious;
V = hp2d_soil_view(column, head_hru, water, root);
[water3, actual] = hp2d_slot_soil_et(V, demand_mm, ~ponded(host), valid, L.feddes, ...
    column.theta_r(:, 2), column.theta_s(:, 2), column.alpha_m_inverse(:, 2), column.n_vg(:, 2), ...
    potential_et_m_s(host)*86.4e6);
water(:, 1:3) = water3;
soil_et = actual/1000;

water = hp2d_drain_cascade(column, head_hru, water, dt_s, L.drainage_substeps);
[water, recharge_mm, capillary_mm] = hp2d_water_table_flux(column, head_hru, water, dt_s);
recharge = recharge_mm/1000;
capillary = capillary_mm/1000;
deep_drainage = 0;
if ~groundwater.enabled
    % Without an aquifer there is no capillary source: undo the upward part.
    lowest = hp2d_lowest_unsaturated(hp2d_unsaturated_thickness(column, head_hru));
    water = water - lowest.*capillary_mm;
    capillary(:) = 0;
end
net_exchange = accumarray(host, (recharge - capillary).*area, [n 1]);
if groundwater.enabled
    groundwater.pending_exchange_m3 = groundwater.pending_exchange_m3 + net_exchange;
else
    deep_drainage = sum(max(net_exchange, 0));
end
L.water_mm = water;
state.layered = L;
state = sync_soil(state, L);

shape = size(state.fraction);
state.cumulative_infiltration_m = state.cumulative_infiltration_m + reshape(infiltration, shape);
state.cumulative_recharge_m = state.cumulative_recharge_m + reshape(recharge, shape);
state.cumulative_capillary_m = state.cumulative_capillary_m + reshape(capillary, shape);
state.cumulative_potential_et_m = state.cumulative_potential_et_m + reshape(pet, shape);
state.cumulative_actual_et_m = state.cumulative_actual_et_m + reshape(soil_et + canopy_evaporation, shape);
state.cumulative_surface_evaporation_m = state.cumulative_surface_evaporation_m + ...
    surface_evaporation./mesh.cell_area(:);
state.last_infiltration_rate_m_s = reshape(infiltration/dt_s, shape);
state.last_recharge_rate_m_s = reshape(recharge/dt_s, shape);
state.last_capillary_rate_m_s = reshape(capillary/dt_s, shape);
state.last_actual_et_rate_m_s = reshape(soil_et/dt_s, shape);
state.last_canopy_evaporation_rate_m_s = reshape(canopy_evaporation/dt_s, shape);
state.last_surface_evaporation_rate_m_s = surface_evaporation./mesh.cell_area(:)/dt_s;

rain_volume = sum(rain.*area);
et_volume = sum((soil_et + canopy_evaporation).*area) + sum(surface_evaporation);
final_storage = layered_storage(surface_volume, state, L, groundwater, head);
diagnostics = struct('precipitation_volume_m3', rain_volume, ...
    'infiltration_volume_m3', sum(infiltration.*area), 'actual_et_volume_m3', et_volume, ...
    'recharge_volume_m3', sum(recharge.*area), 'capillary_volume_m3', sum(capillary.*area), ...
    'deep_drainage_volume_m3', deep_drainage, 'saturation_excess_volume_m3', 0, ...
    'mass_residual_m3', final_storage - (initial_storage + rain_volume - et_volume - deep_drainage));
end

function head = groundwater_head(groundwater, L)
if groundwater.enabled
    head = groundwater.head(:);
else
    head = L.head_m(:);
end
end

function total = layered_storage(surface_volume, state, L, groundwater, head)
column_mm = hp2d_column_storage_mm(L.column, head(L.host), L.water_mm);
total = sum(surface_volume) + sum(state.canopy_storage_m(:).*L.area_m2) + sum(column_mm/1000.*L.area_m2);
if groundwater.enabled
    total = total + sum(groundwater.pending_exchange_m3);
end
end

function state = sync_soil(state, L)
shape = size(state.fraction);
state.near_storage_m = reshape(L.water_mm(:, 1)/1000, shape);
state.root_storage_m = reshape(L.water_mm(:, 2)/1000, shape);
state.transmission_storage_m = reshape(L.water_mm(:, 3)/1000, shape);
end

function value = cell_vector(value, n)
value = double(value);
if isscalar(value), value = repmat(value, n, 1); else, value = value(:); end
end
