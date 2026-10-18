function Summary = run_feddes_legacy_consistency()
%RUN_FEDDES_LEGACY_CONSISTENCY Legacy three-layer ET with Feddes equals the layered Feddes ET.
%   extract_layered_et (legacy scheme) and hp2d_slot_soil_et (layered scheme) implement the
%   same HydroPol2D-Python rule; for the same storages, thicknesses and parameters they
%   must take the same water. Feddes off must reproduce the near-then-root rule, and Feddes
%   must never take more than without stress.
root = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'));
rng(3);
n = 50;
near_t = 0.05 + zeros(n, 1);
root_t = rand(n, 1)*0.9 + 0.05;
theta_r = 0.05; theta_s = 0.45; alpha = 3.0; n_vg = 1.6;
near_cap = (theta_s - theta_r)*near_t*1000;
root_cap = (theta_s - theta_r)*root_t*1000;
near_storage = rand(n, 1).*near_cap;
root_storage = rand(n, 1).*root_cap;
root_depth = near_t + root_t;
root_depth(1:5) = 0.03;  % roots inside the near-surface layer
demand = rand(n, 1)*3;
rate = demand*48;  % mm/day for a 30 min step
feddes = struct('h1_m', -0.1, 'h2_m', -0.25, 'h3_high_m', -2, 'h3_low_m', -6, 'h4_m', -80);

L = struct('root_depth_m', root_depth, 'near_surface_thickness_m', near_t, ...
    'root_zone_thickness_m', root_t, 'transmission_thickness_m', zeros(n, 1), ...
    'near_surface_storage_mm', near_storage, 'root_zone_storage_mm', root_storage, ...
    'transmission_storage_mm', zeros(n, 1), 'near_surface_capacity_mm', near_cap, ...
    'root_zone_capacity_mm', root_cap, 'transmission_capacity_mm', zeros(n, 1), ...
    'total_vadose_capacity_mm', near_cap + root_cap, 'pervious_fraction', ones(n, 1));
S = struct('Layers', L, 'theta_r', theta_r + zeros(n, 1), 'theta_sat', theta_s + zeros(n, 1), ...
    'alpha_vg', alpha + zeros(n, 1), 'n_vg', n_vg + zeros(n, 1), 'I_t', near_storage + root_storage);
V = struct('root_depth_m', root_depth, 'near_surface_thickness_m', near_t, ...
    'root_zone_thickness_m', root_t, 'pervious_fraction', ones(n, 1), ...
    'capacity_mm', [near_cap, root_cap, zeros(n, 1)], ...
    'water_mm', [near_storage, root_storage, zeros(n, 1)]);
valid = true(n, 1); idx_nan = false(n, 1);

S.Feddes = feddes;
[~, legacy_feddes] = extract_layered_et(S, demand, valid, idx_nan, rate);
[~, layered_feddes] = hp2d_slot_soil_et(V, demand, valid, valid, feddes, theta_r + zeros(n, 1), ...
    theta_s + zeros(n, 1), alpha + zeros(n, 1), n_vg + zeros(n, 1), rate);
S.Feddes = [];
[~, legacy_plain] = extract_layered_et(S, demand, valid, idx_nan, rate);
[~, layered_plain] = hp2d_slot_soil_et(V, demand, valid, valid, [], theta_r + zeros(n, 1), ...
    theta_s + zeros(n, 1), alpha + zeros(n, 1), n_vg + zeros(n, 1), rate);
feddes_difference = max(abs(legacy_feddes - layered_feddes));
plain_difference = max(abs(legacy_plain - layered_plain));
stressed = all(legacy_feddes <= legacy_plain + 1e-12);
Summary = table(feddes_difference, plain_difference, stressed, ...
    'VariableNames', {'feddes_max_difference_mm', 'plain_max_difference_mm', 'feddes_never_exceeds_plain'});
disp(Summary);
assert(feddes_difference < 1e-12 && plain_difference < 1e-12 && stressed, ...
    'HydroPol2D:FeddesConsistency', 'Legacy and layered Feddes ET disagree.');
end
