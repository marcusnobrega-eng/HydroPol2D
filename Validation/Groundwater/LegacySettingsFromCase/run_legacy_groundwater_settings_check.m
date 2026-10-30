% VAL-GW-SETTINGS-001: legacy groundwater settings come from the case.
%
% HydroPol2D-Python's legacy groundwater takes from the configuration:
% the Boussinesq Courant number (groundwater_courant, L17), the adaptive
% solver floor (groundwater_min_dt_s) and a stable step of at least 1 s
% (L15), and the capillary-rise field-capacity head and the soil's Mualem
% l_vg (L16) (core/engine.py, modules/groundwater.py). MATLAB fixed the
% Courant number at 0.25, the floor at min(dt, 1 s), had no 1 s floor on
% the stable step, and fixed the head at 3.3 m and l at 0.5.
%   A. estimate_groundwater_stable_timestep scales with the Courant number
%      and is never below 1 s.
%   B. simulate_layered_capillary_rise uses the soil's l_vg and the given
%      field-capacity head: one cell matches hp2d_vg_mualem_capillary_transfer
%      with l = -1 and a 1.0 m head.
%   C. Groundwater_Module reads groundwater_courant,
%      capillary_field_capacity_head_m and groundwater_min_dt_s.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

% A. stable step
soil = struct('Sy', 0.2 * ones(2), 'ksat_gw', 36 * ones(2));   % 1e-5 m/s
bc = struct('h_t', 5 * ones(2));
dt_025 = estimate_groundwater_stable_timestep(soil, bc, zeros(2), 10, 10, false(2), 0.25);
dt_050 = estimate_groundwater_stable_timestep(soil, bc, zeros(2), 10, 10, false(2), 0.5);
dt_default = estimate_groundwater_stable_timestep(soil, bc, zeros(2), 10, 10, false(2));
fast = soil; fast.ksat_gw = 3.6e7 * ones(2);                   % 10 m/s
dt_floor = estimate_groundwater_stable_timestep(fast, bc, zeros(2), 1, 1, false(2), 0.25);
expected_025 = 0.25 * 100 / (1e-5 * 5 / 0.2) / 4;
stable_ok = abs(dt_025 - expected_025) < 1e-9 && abs(dt_050 - 2 * dt_025) < 1e-9 && ...
    dt_default == dt_025 && dt_floor == 1;

% B. capillary rise with l_vg and the field-capacity head
L = struct('root_depth_m', 0.5, 'root_zone_capacity_mm', 160, 'near_surface_capacity_mm', 20, ...
    'root_zone_storage_mm', 20, 'near_surface_storage_mm', 5, 'ksat_root_zone', 20, ...
    'near_surface_thickness_m', 0.05, 'root_zone_thickness_m', 0.45, 'transmission_thickness_m', 0, ...
    'transmission_capacity_mm', 0, 'transmission_storage_mm', 0, 'pervious_fraction', 1);
S = struct('Layers', L, 'alpha_vg', 2.0, 'n_vg', 1.4, 'l_vg', -1.0, 'I_t', 25, 'theta_r', 0.05, ...
    'theta_sat', 0.45);
bc_cap = struct('h_t', 8.5);                                       % table 1.5 m down
dt_s = 3600;
[rate, ~] = simulate_layered_capillary_rise(S, bc_cap, 10, dt_s, false, 1, 1.0);
expected_mm = 1000 * hp2d_vg_mualem_capillary_transfer(20 / 3.6e6, 2.0, 1.4, -1.0, 1.0, 1.0, ...
    0.160, 0.020, 1, dt_s);
capillary_error = abs(rate * dt_s * 1000 - expected_mm);
old_rule_mm = 1000 * hp2d_vg_mualem_capillary_transfer(20 / 3.6e6, 2.0, 1.4, 0.5, 1.0, 3.3, ...
    0.160, 0.020, 1, dt_s);
capillary_ok = capillary_error < 1e-12 && abs(expected_mm - old_rule_mm) > 1e-6;

% C. the module reads the settings
module_text = fileread(fullfile(functions_dir, 'Groundwater_Module.m'));
reads_settings = contains(module_text, 'flags.groundwater_courant') && ...
    contains(module_text, 'flags.capillary_field_capacity_head_m') && ...
    contains(module_text, 'minimum_step_s=groundwater_min_dt_s');

passed = stable_ok && capillary_ok && reads_settings;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-GW-SETTINGS-001", "Legacy groundwater settings from the case", ...
    dt_025, dt_050, dt_floor, expected_mm, capillary_error, old_rule_mm, reads_settings, status, true, ...
    'VariableNames', {'case_id','case_name','stable_dt_courant_025_s','stable_dt_courant_05_s', ...
    'stable_dt_floor_s','capillary_mm','capillary_error_mm','old_rule_capillary_mm', ...
    'module_reads_settings','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-GW-SETTINGS-001 failed.');
end
fprintf('VAL-GW-SETTINGS-001 passed.\n');
