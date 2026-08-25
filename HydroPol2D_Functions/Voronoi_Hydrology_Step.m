function [surface_volume, state, groundwater, diagnostics] = Voronoi_Hydrology_Step( ...
    mesh, surface_volume, state, groundwater, precipitation_m_s, potential_et_m_s, open_water_evaporation_m_s, dt_s)
%VORONOI_HYDROLOGY_STEP Conservative canopy, vadose, ET, and recharge step.

n = mesh.n_cells;
precipitation_m_s = cell_vector(precipitation_m_s,n);
potential_et_m_s = cell_vector(potential_et_m_s,n);
open_water_evaporation_m_s = cell_vector(open_water_evaporation_m_s,n);
dt_s = max(double(dt_s),eps);
initial_storage = sum(surface_volume) + hru_storage(state);
if groundwater.enabled, initial_storage = initial_storage + sum(groundwater.pending_exchange_m3); end
saturation_excess = zeros(size(state.fraction));

% A rising water table reduces vadose capacity. Displaced vadose water is
% saturation excess and returns to the polygon surface exactly once.
if groundwater.enabled
    [state, saturation_excess] = update_geometry_and_clamp(state, mesh.surface_bed(:), groundwater.head(:));
    surface_volume = surface_volume + sum(saturation_excess .* state.area_m2,2);
end

rain_depth = precipitation_m_s .* dt_s;
pet_depth = potential_et_m_s .* dt_s;
rain_hru = rain_depth + zeros(size(state.fraction));
pet_hru = pet_depth + zeros(size(state.fraction));

% Rutter canopy bucket, identical storage convention to interceptionModel.
capacity = state.canopy_capacity_per_lai_m .* state.lai;
beta = zeros(size(capacity));
active_canopy = capacity > 0;
beta(active_canopy) = state.canopy_storage_m(active_canopy) ./ capacity(active_canopy);
canopy_evaporation = min(max(beta,0) .* pet_hru, state.canopy_storage_m + rain_hru);
canopy_evaporation(~active_canopy) = 0;
provisional = state.canopy_storage_m + rain_hru - canopy_evaporation;
throughfall = max(provisional - capacity,0);
state.canopy_storage_m = min(provisional,capacity);
surface_volume = surface_volume + sum(throughfall .* state.area_m2,2);

% Ponded cells evaporate from the shared surface store; otherwise ET is
% extracted from root-accessible HRU storage above its wilting threshold.
ponded = surface_volume > 0;
surface_evaporation_volume = min(surface_volume, open_water_evaporation_m_s .* dt_s .* mesh.cell_area(:));
surface_evaporation_volume(~ponded) = 0;
surface_volume = surface_volume - surface_evaporation_volume;
soil_et_demand = pet_hru .* state.crop_coefficient;
soil_et_demand(ponded,:) = 0;
[state.near_storage_m, take_near] = extract_above_wilting( ...
    state.near_storage_m,state.near_capacity_m,state.wilting_saturation,soil_et_demand);
remaining_et = soil_et_demand - take_near;
[state.root_storage_m, take_root] = extract_above_wilting( ...
    state.root_storage_m,state.root_capacity_m,state.wilting_saturation,remaining_et);
actual_soil_et = take_near + take_root;

% Darcy-vG/Mualem infiltration closure from the raster implementation.
surface_depth = surface_volume ./ mesh.cell_area(:);
available = surface_depth + zeros(size(state.fraction));
remaining_capacity = max(state.near_capacity_m-state.near_storage_m,0) + ...
    max(state.root_capacity_m-state.root_storage_m,0) + ...
    max(state.transmission_capacity_m-state.transmission_storage_m,0);
theta = state.theta_r + state.near_storage_m ./ max(state.pervious_fraction .* state.near_thickness_m,eps);
theta = min(max(theta,state.theta_r),state.theta_sat);
se = min(max((theta-state.theta_r) ./ max(state.theta_sat-state.theta_r,eps),1e-6),1);
m_vg = 1 - 1 ./ state.n_vg;
pressure_head = -(1 ./ state.alpha_vg_1_m) .* max(se.^(-1./m_vg)-1,0).^(1./state.n_vg);
pressure_head(se >= 0.999999) = 0;
wetting = min(min(available,remaining_capacity),max(state.near_capacity_m-state.near_storage_m,0));
theta_top = theta + wetting ./ max(state.pervious_fraction .* state.near_thickness_m,eps);
theta_top = min(max(theta_top,state.theta_r),state.theta_sat);
se_top = min(max((theta_top-state.theta_r) ./ max(state.theta_sat-state.theta_r,eps),1e-6),1);
term = min(max(1-se_top.^(1./m_vg),0),1);
relative_k = se_top.^state.l_vg .* (1-term.^m_vg).^2;
k_top = state.ksat_m_s .* min(max(relative_k,0),1);
k_effective = state.surface_conductivity_weight .* state.ksat_m_s + ...
    (1-state.surface_conductivity_weight) .* k_top;
head_difference = min(max(available-pressure_head,0),state.infiltration_head_cap_m);
capacity_rate = k_effective .* (1 + head_difference ./ max(state.near_thickness_m,eps)) .* state.pervious_fraction;
infiltration = min(min(available,capacity_rate.*dt_s),remaining_capacity);
infiltration(state.pervious_fraction <= 1e-12 | state.water_table_depth_m <= 0) = 0;

% All HRUs see the same ponded depth. Fractions guarantee conservation, but
% this final limiter protects roundoff and partially specified HRU matrices.
infiltration_volume = sum(infiltration .* state.area_m2,2);
scale = min(1,surface_volume ./ max(infiltration_volume,eps));
infiltration = infiltration .* scale;
infiltration_volume = sum(infiltration .* state.area_m2,2);
surface_volume = max(surface_volume-infiltration_volume,0);
[state.near_storage_m,state.root_storage_m,state.transmission_storage_m] = ...
    add_top_down(state.near_storage_m,state.root_storage_m,state.transmission_storage_m, ...
    state.near_capacity_m,state.root_capacity_m,state.transmission_capacity_m,infiltration);

% Reuse the same free-drainage ordering as the layered raster model.
[state.near_storage_m,state.root_storage_m,state.transmission_storage_m,recharge] = drain_layers(state,dt_s);
capillary = zeros(size(recharge));
if groundwater.enabled && state.capillary_rise_enabled
    available_groundwater_depth = max((groundwater.storage + groundwater.pending_exchange_m3) ./ mesh.cell_area(:),0);
    [state,capillary] = capillary_rise(state,available_groundwater_depth,dt_s);
end

net_exchange_volume = sum((recharge-capillary) .* state.area_m2,2);
deep_drainage_volume = 0;
if groundwater.enabled
    groundwater.pending_exchange_m3 = groundwater.pending_exchange_m3 + net_exchange_volume;
else
    deep_drainage_volume = sum(max(net_exchange_volume,0));
end

state.cumulative_infiltration_m = state.cumulative_infiltration_m + infiltration;
state.cumulative_recharge_m = state.cumulative_recharge_m + recharge;
state.cumulative_actual_et_m = state.cumulative_actual_et_m + actual_soil_et + canopy_evaporation;
state.cumulative_surface_evaporation_m = state.cumulative_surface_evaporation_m + ...
    surface_evaporation_volume ./ mesh.cell_area(:);
state.last_infiltration_rate_m_s = infiltration ./ dt_s;
state.last_recharge_rate_m_s = recharge ./ dt_s;
state.last_capillary_rate_m_s = capillary ./ dt_s;
state.last_actual_et_rate_m_s = actual_soil_et ./ dt_s;
state.last_canopy_evaporation_rate_m_s = canopy_evaporation ./ dt_s;
state.last_surface_evaporation_rate_m_s = surface_evaporation_volume ./ mesh.cell_area(:) ./ dt_s;

rain_volume = sum(precipitation_m_s .* mesh.cell_area(:)) .* dt_s;
et_volume = sum((actual_soil_et+canopy_evaporation).*state.area_m2,'all') + sum(surface_evaporation_volume);
final_storage = sum(surface_volume) + hru_storage(state);
if groundwater.enabled, final_storage = final_storage + sum(groundwater.pending_exchange_m3); end
diagnostics = struct('precipitation_volume_m3',rain_volume,'infiltration_volume_m3',sum(infiltration_volume), ...
    'actual_et_volume_m3',et_volume,'recharge_volume_m3',sum(recharge.*state.area_m2,'all'), ...
    'capillary_volume_m3',sum(capillary.*state.area_m2,'all'),'deep_drainage_volume_m3',deep_drainage_volume, ...
    'saturation_excess_volume_m3',sum(saturation_excess.*state.area_m2,'all'), ...
    'mass_residual_m3',final_storage-(initial_storage+rain_volume-et_volume-deep_drainage_volume));
end

function [state,excess] = update_geometry_and_clamp(state,surface_bed,groundwater_head)
old_near=state.near_storage_m; old_root=state.root_storage_m; old_trans=state.transmission_storage_m;
zwt=min(max(surface_bed-groundwater_head,0),state.soil_depth_m);
state.water_table_depth_m=zwt;
state.near_thickness_m=min(state.surface_layer_depth_m,zwt);
root_bottom=min(state.root_depth_m,zwt);
state.root_thickness_m=max(root_bottom-state.near_thickness_m,0);
state.transmission_thickness_m=max(zwt-state.near_thickness_m-state.root_thickness_m,0);
water_capacity=state.pervious_fraction.*max(state.theta_sat-state.theta_r,0);
state.near_capacity_m=water_capacity.*state.near_thickness_m;
state.root_capacity_m=water_capacity.*state.root_thickness_m;
state.transmission_capacity_m=water_capacity.*state.transmission_thickness_m;
state.near_storage_m=min(old_near,state.near_capacity_m);
state.root_storage_m=min(old_root,state.root_capacity_m);
state.transmission_storage_m=min(old_trans,state.transmission_capacity_m);
excess=(old_near-state.near_storage_m)+(old_root-state.root_storage_m)+(old_trans-state.transmission_storage_m);
end

function [storage,taken] = extract_above_wilting(storage,capacity,wilting,demand)
available=max(storage-wilting.*capacity,0);
taken=min(max(demand,0),available);
storage=storage-taken;
end

function [near,root,trans] = add_top_down(near,root,trans,near_cap,root_cap,trans_cap,water)
add=min(water,max(near_cap-near,0)); near=near+add; water=water-add;
add=min(water,max(root_cap-root,0)); root=root+add; water=water-add;
add=min(water,max(trans_cap-trans,0)); trans=trans+add;
end

function [near,root,trans,recharge] = drain_layers(state,dt_s)
near=state.near_storage_m; root=state.root_storage_m; trans=state.transmission_storage_m;
near_drain=drainage(near,state.near_capacity_m,state.near_thickness_m,state.ksat_m_s,state.theta_r,state.theta_sat,state.n_vg,dt_s);
[near,root,moved]=move_down(near,root,state.root_capacity_m,near_drain,state.root_capacity_m>0);
remaining=max(near_drain-moved,0);
[near,trans,moved]=move_down(near,trans,state.transmission_capacity_m,remaining,state.root_capacity_m<=0 & state.transmission_capacity_m>0);
near_to_gw=min(max(remaining-moved,0),near).*(state.root_capacity_m<=0 & state.transmission_capacity_m<=0);
near=near-near_to_gw;
root_drain=drainage(root,state.root_capacity_m,state.root_thickness_m,state.ksat_m_s,state.theta_r,state.theta_sat,state.n_vg,dt_s);
[root,trans,moved]=move_down(root,trans,state.transmission_capacity_m,root_drain,state.transmission_capacity_m>0);
root_to_gw=min(max(root_drain-moved,0),root).*(state.root_capacity_m>0 & state.transmission_capacity_m<=0);
root=root-root_to_gw;
trans_drain=drainage(trans,state.transmission_capacity_m,state.transmission_thickness_m,state.ksat_m_s,state.theta_r,state.theta_sat,state.n_vg,dt_s);
trans_drain=min(trans_drain,trans); trans=trans-trans_drain;
recharge=near_to_gw+root_to_gw+trans_drain;
end

function amount = drainage(storage,capacity,thickness,ksat,theta_r,theta_s,n_vg,dt_s)
theta=theta_r+storage./max(thickness,eps);
theta=min(max(theta,theta_r),theta_s);
se=min(max((theta-theta_r)./max(theta_s-theta_r,eps),1e-6),1);
m=1-1./n_vg; term=min(max(1-se.^(1./m),0),1);
kr=se.^0.5.*(1-term.^m).^2;
amount=min(max(ksat.*kr.*dt_s,0),storage);
amount(capacity<=0 | thickness<=0)=0;
end

function [upper,lower,moved] = move_down(upper,lower,lower_capacity,demand,route)
moved=zeros(size(upper));
moved(route)=min(min(demand(route),max(lower_capacity(route)-lower(route),0)),upper(route));
upper=upper-moved; lower=lower+moved;
end

function [state,capillary] = capillary_rise(state,available_cell_depth,dt_s)
remaining=available_cell_depth+zeros(size(state.fraction)); capillary=zeros(size(remaining));
proximity=max(1-state.water_table_depth_m./state.capillary_extinction_depth_m,0);
[state.transmission_storage_m,moved]=fill(state.transmission_storage_m,state.transmission_capacity_m,state.ksat_m_s,proximity,dt_s,remaining);
remaining=remaining-moved; capillary=capillary+moved;
[state.root_storage_m,moved]=fill(state.root_storage_m,state.root_capacity_m,state.ksat_m_s,proximity,dt_s,remaining);
remaining=remaining-moved; capillary=capillary+moved;
eligible=state.transmission_capacity_m<=0 & state.root_capacity_m<=0;
[candidate,moved]=fill(state.near_storage_m,state.near_capacity_m,state.ksat_m_s,proximity,dt_s,remaining);
moved(~eligible)=0; state.near_storage_m(eligible)=candidate(eligible); capillary=capillary+moved;
end

function [storage,moved] = fill(storage,capacity,ksat,proximity,dt_s,remaining)
dryness=max(1-storage./max(capacity,eps),0);
moved=min(min(max(capacity-storage,0),ksat.*proximity.*dryness.*dt_s),max(remaining,0));
moved(capacity<=0)=0; storage=storage+moved;
end

function total = hru_storage(state)
total=sum((state.canopy_storage_m+state.near_storage_m+state.root_storage_m+state.transmission_storage_m).*state.area_m2,'all');
end

function value = cell_vector(value,n)
if isscalar(value), value=repmat(double(value),n,1); else, value=double(value(:)); end
assert(numel(value)==n && all(isfinite(value)));
end
