function [recharge_rate, Soil_Properties, cumulative_recharge] = simulate_layered_groundwater_recharge(Soil_Properties, dt_s, idx_nan, cumulative_recharge)
%SIMULATE_LAYERED_GROUNDWATER_RECHARGE Percolate water through vadose layers.

dt_s = max(dt_s, 1e-12);
layers = Soil_Properties.Layers;

zero_matrix = zeros(size(Soil_Properties.I_t), 'like', Soil_Properties.I_t);
recharge_depth_mm = zero_matrix;

% Near-surface drainage goes first to the explicit root zone, then to the
% transmission zone, and only directly to groundwater if no lower vadose
% layer exists.
near_drain_mm = local_drainage_depth( ...
    layers.near_surface_storage_mm, layers.near_surface_capacity_mm, ...
    layers.near_surface_thickness_m, layers.ksat_near_surface, ...
    layers.theta_r_near_surface, layers.theta_sat_near_surface, ...
    layers.n_vg_near_surface, dt_s);

idx_to_root = layers.root_zone_capacity_mm > 0;
[layers.near_surface_storage_mm, layers.root_zone_storage_mm, near_to_root] = ...
    local_move_down(layers.near_surface_storage_mm, layers.root_zone_storage_mm, ...
    layers.root_zone_capacity_mm, near_drain_mm, idx_to_root);

near_remaining_drain = max(near_drain_mm - near_to_root, 0);
idx_to_trans_from_near = ~idx_to_root & layers.transmission_capacity_mm > 0;
[layers.near_surface_storage_mm, layers.transmission_storage_mm, near_to_trans] = ...
    local_move_down(layers.near_surface_storage_mm, layers.transmission_storage_mm, ...
    layers.transmission_capacity_mm, near_remaining_drain, idx_to_trans_from_near);

near_to_gw = max(near_remaining_drain - near_to_trans, 0) .* ...
    double(~idx_to_root & layers.transmission_capacity_mm <= 0);
near_to_gw = min(near_to_gw, layers.near_surface_storage_mm);
layers.near_surface_storage_mm = layers.near_surface_storage_mm - near_to_gw;
recharge_depth_mm = recharge_depth_mm + near_to_gw;

% Root-zone drainage goes to the transmission zone when present.
root_drain_mm = local_drainage_depth( ...
    layers.root_zone_storage_mm, layers.root_zone_capacity_mm, ...
    layers.root_zone_thickness_m, layers.ksat_root_zone, ...
    layers.theta_r_root_zone, layers.theta_sat_root_zone, ...
    layers.n_vg_root_zone, dt_s);

idx_root_to_trans = layers.transmission_capacity_mm > 0;
[layers.root_zone_storage_mm, layers.transmission_storage_mm, root_to_trans] = ...
    local_move_down(layers.root_zone_storage_mm, layers.transmission_storage_mm, ...
    layers.transmission_capacity_mm, root_drain_mm, idx_root_to_trans);

root_to_gw = max(root_drain_mm - root_to_trans, 0) .* ...
    double(layers.root_zone_capacity_mm > 0 & layers.transmission_capacity_mm <= 0);
root_to_gw = min(root_to_gw, layers.root_zone_storage_mm);
layers.root_zone_storage_mm = layers.root_zone_storage_mm - root_to_gw;
recharge_depth_mm = recharge_depth_mm + root_to_gw;

% Transmission-zone drainage is recharge to groundwater.
trans_drain_mm = local_drainage_depth( ...
    layers.transmission_storage_mm, layers.transmission_capacity_mm, ...
    layers.transmission_thickness_m, layers.ksat_transmission, ...
    layers.theta_r_transmission, layers.theta_sat_transmission, ...
    layers.n_vg_transmission, dt_s);

trans_drain_mm = min(trans_drain_mm, layers.transmission_storage_mm);
layers.transmission_storage_mm = layers.transmission_storage_mm - trans_drain_mm;
recharge_depth_mm = recharge_depth_mm + trans_drain_mm;

recharge_depth_mm(idx_nan) = nan;
recharge_rate = recharge_depth_mm ./ 1000 ./ dt_s;
recharge_rate(idx_nan) = nan;

cumulative_recharge = cumulative_recharge + max(recharge_rate, 0) .* dt_s .* 1000;

Soil_Properties.Layers = layers;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
end

function drainage_mm = local_drainage_depth(storage_mm, capacity_mm, thickness_m, ksat_mm_h, theta_r, theta_s, n_vg, dt_s)

storage_mm = max(storage_mm, 0);
capacity_mm = max(capacity_mm, 0);
thickness_m = max(thickness_m, 0);

drainage_mm = zeros(size(storage_mm), 'like', storage_mm);
idx = capacity_mm > 0 & thickness_m > 0 & storage_mm > 0;
if ~any(idx(:))
    return
end

m_vg = 1 - 1 ./ n_vg;
theta = theta_r + storage_mm ./ max(thickness_m .* 1000, eps);
theta = max(theta, theta_r);
theta = min(theta, theta_s);

Se = (theta - theta_r) ./ max(theta_s - theta_r, 1e-12);
Se = min(max(Se, 1e-6), 1);

term = 1 - Se .^ (1 ./ m_vg);
term = min(max(term, 0), 1);
Kr = Se .^ 0.5 .* (1 - term .^ m_vg) .^ 2;
Kr = min(max(Kr, 0), 1);

K_m_s = ksat_mm_h ./ 1000 ./ 3600 .* Kr;
drainage_mm(idx) = K_m_s(idx) .* dt_s .* 1000;
drainage_mm = min(max(drainage_mm, 0), storage_mm);
end

function [upper_storage, lower_storage, moved] = local_move_down(upper_storage, lower_storage, lower_capacity, demand, idx_route)

moved = zeros(size(upper_storage), 'like', upper_storage);
idx_route = logical(idx_route);

lower_space = max(lower_capacity - lower_storage, 0);
moved(idx_route) = min(demand(idx_route), lower_space(idx_route));
moved(idx_route) = min(moved(idx_route), upper_storage(idx_route));

upper_storage = upper_storage - moved;
lower_storage = lower_storage + moved;
end
