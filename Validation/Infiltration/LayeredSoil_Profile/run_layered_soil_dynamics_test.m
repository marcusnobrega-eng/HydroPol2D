clear; clc;

this_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(this_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(this_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

idx_nan = false(1,1);
[Soil_Properties, ~] = base_profile(0.60, 0.30);

Soil_Properties = initialize_layered_soil_storage(Soil_Properties, idx_nan);
[Soil_Properties, accepted_mm] = add_layered_infiltration(Soil_Properties, 50, idx_nan);
topdown_pass = abs(accepted_mm - 50) < 1e-10 && ...
    abs(Soil_Properties.Layers.near_surface_storage_mm - 30) < 1e-10 && ...
    abs(Soil_Properties.Layers.root_zone_storage_mm - 20) < 1e-10 && ...
    abs(Soil_Properties.Layers.transmission_storage_mm) < 1e-10;

Soil_Properties.Layers.near_surface_storage_mm = 20;
Soil_Properties.Layers.root_zone_storage_mm = 30;
Soil_Properties.Layers.transmission_storage_mm = 0;
Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
[Soil_Properties, actual_et_mm] = extract_layered_et(Soil_Properties, 35, true(1,1), idx_nan);
et_pass = abs(actual_et_mm - 35) < 1e-10 && ...
    abs(Soil_Properties.Layers.near_surface_storage_mm) < 1e-10 && ...
    abs(Soil_Properties.Layers.root_zone_storage_mm - 15) < 1e-10;

[Soil_Shallow, ~] = base_profile(0.60, 0.05);
Soil_Shallow = initialize_layered_soil_storage(Soil_Shallow, idx_nan);
Soil_Shallow.Layers.near_surface_storage_mm = 20;
Soil_Shallow = sync_layered_soil_storage(Soil_Shallow, idx_nan);
[Soil_Shallow, actual_shallow_et_mm] = extract_layered_et(Soil_Shallow, 20, true(1,1), idx_nan);
shallow_root_pass = abs(actual_shallow_et_mm - 10) < 1e-10 && ...
    abs(Soil_Shallow.Layers.near_surface_storage_mm - 10) < 1e-10;

[Soil_Recharge, ~] = base_profile(0.60, 0.30);
Soil_Recharge.ksat(:) = 1e12;
Soil_Recharge.Layers.ksat_near_surface(:) = 1e12;
Soil_Recharge.Layers.ksat_root_zone(:) = 1e12;
Soil_Recharge.Layers.ksat_transmission(:) = 1e12;
Soil_Recharge = initialize_layered_soil_storage(Soil_Recharge, idx_nan);
Soil_Recharge.Layers.near_surface_storage_mm = 10;
Soil_Recharge.Layers.root_zone_storage_mm = 0;
Soil_Recharge.Layers.transmission_storage_mm = 0;
Soil_Recharge = sync_layered_soil_storage(Soil_Recharge, idx_nan);
initial_storage_mm = Soil_Recharge.I_t;
[recharge_rate, Soil_Recharge, cumulative_recharge] = simulate_layered_groundwater_recharge( ...
    Soil_Recharge, 60, idx_nan, zeros(1,1));
recharged_mm = recharge_rate * 60 * 1000;
percolation_pass = abs((initial_storage_mm - Soil_Recharge.I_t) - recharged_mm) < 1e-8 && ...
    abs(cumulative_recharge - recharged_mm) < 1e-8 && recharged_mm > 0;

[Soil_Clamp, ~] = base_profile(0.60, 0.30);
Soil_Clamp = initialize_layered_soil_storage(Soil_Clamp, idx_nan);
Soil_Clamp.Layers.near_surface_storage_mm = 50;
Soil_Clamp.Layers.root_zone_storage_mm = 70;
Soil_Clamp.Layers.transmission_storage_mm = 100;
[Soil_Clamp, excess_mm] = clamp_layered_soil_storage(Soil_Clamp, idx_nan);
clamp_pass = abs(excess_mm - 40) < 1e-10 && ...
    abs(Soil_Clamp.I_t - Soil_Clamp.Layers.total_vadose_capacity_mm) < 1e-10;

case_id = ["top_down_infiltration"; "root_access_et"; "shallow_root_fallback"; ...
    "layer_percolation_recharge"; "capacity_clamp"];
metric_value = [accepted_mm; actual_et_mm; actual_shallow_et_mm; recharged_mm; excess_mm];
threshold = [1e-10; 1e-10; 1e-10; 1e-8; 1e-10];
passed = [topdown_pass; et_pass; shallow_root_pass; percolation_pass; clamp_pass];

Diagnostics = table(case_id, metric_value, threshold, passed);
writetable(Diagnostics, fullfile(out_dir, 'Layered_Dynamics_Diagnostics.csv'));

status = "pass";
if ~all(passed)
    status = "fail";
end

PassFail = table( ...
    "P1-INFIL-LAYERS-002", "Layered soil dynamics", status, all(passed), ...
    'VariableNames', {'case_id','case_name','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Layered_Dynamics_Pass_Fail.csv'));

disp(Diagnostics);
disp(PassFail);

function [Soil_Properties, LULC_Properties] = base_profile(soil_depth_m, root_depth_m)

Soil_Properties = struct();
Soil_Properties.Soil_Depth = soil_depth_m;
Soil_Properties.theta_sat = 0.40;
Soil_Properties.theta_r = 0.10;
Soil_Properties.theta_i = 0.10;
Soil_Properties.alpha_vg = 1.5;
Soil_Properties.n_vg = 1.6;
Soil_Properties.ksat = 100;
Soil_Properties.Ks_multiplier_near_surface = 1;
Soil_Properties.Ks_multiplier_root_zone = 1;
Soil_Properties.Ks_multiplier_transmission = 1;
Soil_Properties.I_t = 0;
Soil_Properties.I_p = 0;
Soil_Properties.I_0 = 0;

LULC_Properties = struct();
LULC_Properties.root_depth_m = root_depth_m;
LULC_Properties.idx_imp = false(1,1);

elevation = 100;
BC_States = struct();
BC_States.h_t = elevation - soil_depth_m;

idx_nan = false(1,1);
options = struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005);
[Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
    Soil_Properties, LULC_Properties, BC_States, elevation, idx_nan, options);
end
