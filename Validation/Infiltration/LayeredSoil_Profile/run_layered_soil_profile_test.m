% P1-INFIL-LAYERS-001: test layered soil profile fallback rules.

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
output_dir = fullfile(case_dir, 'Outputs', 'Validation');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
addpath(functions_dir);

scenario = ["deep_profile"; "root_shallower_than_surface"; "bedrock_shallower_than_root"; "water_table_at_surface"; "water_table_in_root_zone"];
soil_depth = [2.00; 1.00; 0.06; 2.00; 2.00];
root_depth = [1.20; 0.05; 1.00; 1.00; 1.50];
water_table_depth = [1.50; 1.00; 0.06; 0.00; 0.35];

n = numel(scenario);
template = zeros(n, 1);
idx_nan = false(n, 1);
elevation = 100 * ones(n, 1);

Soil_Properties = struct();
Soil_Properties.Soil_Depth = soil_depth;
Soil_Properties.theta_sat = 0.40 * ones(n, 1);
Soil_Properties.theta_r = 0.10 * ones(n, 1);
Soil_Properties.theta_i = 0.20 * ones(n, 1);
Soil_Properties.ksat = 10.0 * ones(n, 1);
Soil_Properties.alpha_vg = 2.0 * ones(n, 1);
Soil_Properties.n_vg = 1.5 * ones(n, 1);
Soil_Properties.Ks_multiplier_near_surface = [1.0; 0.8; 1.0; 1.0; 1.0];
Soil_Properties.Ks_multiplier_root_zone = [1.2; 1.1; 1.0; 1.0; 1.0];
Soil_Properties.Ks_multiplier_transmission = [0.7; 0.9; 1.0; 1.0; 1.0];

LULC_Properties = struct();
LULC_Properties.root_depth_m = root_depth;
LULC_Properties.idx_imp = false(n, 1);

BC_States = struct();
BC_States.h_t = elevation - water_table_depth;

options = struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005);
[Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
    Soil_Properties, LULC_Properties, BC_States, elevation, idx_nan, options);

L = Soil_Properties.Layers;
capacity_sum_error_mm = L.total_vadose_capacity_mm - ...
    (L.near_surface_capacity_mm + L.root_zone_capacity_mm + L.transmission_capacity_mm);
thickness_sum_error_m = L.vadose_depth_m - ...
    (L.near_surface_thickness_m + L.root_zone_thickness_m + L.transmission_thickness_m);

expected_near = [0.10; 0.10; 0.06; 0.00; 0.10];
expected_root = [1.10; 0.00; 0.00; 0.00; 0.25];
expected_trans = [0.30; 0.90; 0.00; 0.00; 0.00];

tol = 1e-12;
pass_thickness = abs(L.near_surface_thickness_m - expected_near) <= tol & ...
    abs(L.root_zone_thickness_m - expected_root) <= tol & ...
    abs(L.transmission_thickness_m - expected_trans) <= tol;
pass_nonnegative = L.near_surface_thickness_m >= 0 & ...
    L.root_zone_thickness_m >= 0 & ...
    L.transmission_thickness_m >= 0 & ...
    L.total_vadose_capacity_mm >= 0;
pass_sum = abs(thickness_sum_error_m) <= tol & abs(capacity_sum_error_mm) <= tol;

Diagnostics = table( ...
    scenario, soil_depth, root_depth, water_table_depth, ...
    L.near_surface_thickness_m, L.root_zone_thickness_m, L.transmission_thickness_m, ...
    L.total_vadose_capacity_mm, thickness_sum_error_m, capacity_sum_error_mm, ...
    L.root_within_near_surface, L.root_truncated_by_bedrock, ...
    L.root_truncated_by_water_table, L.water_table_at_surface, ...
    pass_thickness, pass_nonnegative, pass_sum, ...
    'VariableNames', {'scenario','soil_depth_m','root_depth_m','water_table_depth_m', ...
    'near_surface_thickness_m','root_zone_thickness_m','transmission_thickness_m', ...
    'total_vadose_capacity_mm','thickness_sum_error_m','capacity_sum_error_mm', ...
    'root_within_near_surface','root_truncated_by_bedrock', ...
    'root_truncated_by_water_table','water_table_at_surface', ...
    'pass_thickness','pass_nonnegative','pass_sum'});

writetable(Diagnostics, fullfile(output_dir, 'Layered_Profile_Diagnostics.csv'));

all_pass = all(pass_thickness) && all(pass_nonnegative) && all(pass_sum);
PassFail = table( ...
    "P1-INFIL-LAYERS-001", "Infiltration layered profile", ...
    string(local_status(all_pass)), ...
    all(pass_thickness), all(pass_nonnegative), all(pass_sum), all_pass, ...
    "run_layered_soil_profile_test.m", "2026-06-15", ...
    'VariableNames', {'case_id','module','status','thickness_pass', ...
    'nonnegative_pass','closure_pass','report_ready','reviewer','date_utc'});
writetable(PassFail, fullfile(output_dir, 'Pass_Fail.csv'));

disp('Wrote layered soil profile diagnostics to:');
disp(fullfile(output_dir, 'Layered_Profile_Diagnostics.csv'));

function status = local_status(pass)
if pass
    status = 'pass';
else
    status = 'fail';
end
end
