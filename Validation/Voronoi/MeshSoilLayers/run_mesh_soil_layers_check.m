% VAL-MESH-SOIL-001: mesh soil layers drain at their own Ksat with the soil's l.
%
% HydroPol2D-Python's legacy mesh soil drains the near-surface, root-zone
% and transmission layers with their own saturated conductivities, takes
% capillary rise at the root-zone conductivity, and uses the soil's Mualem
% pore-connectivity l_vg (voronoi.py). MATLAB used one ksat_m_s for every
% layer and fixed l at 0.5 in the drainage (audit M19).
% One 100 m2 cell: the root zone holds water at saturation 0.6, the layers
% above and below are empty, the transmission layer is nearly impermeable,
% and there is no rain, ET or groundwater. After one 600 s step the water
% that left the root zone must equal f_p K_root kr(Se; l) dt with
% K_root = 2e-5 m/s and l = -1.2 (the near-surface Ksat is 1e-7 m/s).

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

mesh = struct('n_cells', 1, 'cell_area', 100, 'surface_area', 100, 'surface_bed', 10, 'cell_bed', 10);
fp = 0.8; tr = 0.05; ts = 0.45; n_vg = 1.4; l_vg = -1.2; K_root = 2e-5; dt = 600;
hydrology = struct('pervious_fraction', fp, 'soil_depth_m', 2, 'root_depth_m', 0.6, ...
    'surface_layer_depth_m', 0.1, 'theta_r', tr, 'theta_sat', ts, 'alpha_vg_1_m', 2, ...
    'n_vg', n_vg, 'l_vg', l_vg, 'ksat_m_s', 1e-7, 'ksat_root_m_s', K_root, ...
    'ksat_transmission_m_s', 1e-14, 'initial_soil_saturation', 0, 'capillary_rise_enabled', false);
state = Voronoi_Hydrology_Initialize(mesh, hydrology, mesh.surface_bed - 3);
se = 0.6;
state.root_storage_m = se * state.root_capacity_m;
root_before = state.root_storage_m;
[~, state] = Voronoi_Hydrology_Step(mesh, 0, state, struct('enabled', false), 0, 0, 0, dt);
moved = root_before - state.root_storage_m;
m = 1 - 1 / n_vg;
kr = min(max(se^l_vg * (1 - (1 - se^(1 / m))^m)^2, 0), 1);
expected = min(fp * K_root * kr * dt, root_before);
error_m = abs(moved - expected);
passed = error_m < 1e-12 && expected > 0;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-SOIL-001", "Mesh soil layers drain at their own Ksat with l_vg", ...
    moved, expected, error_m, status, true, 'VariableNames', {'case_id','case_name', ...
    'root_drainage_m','expected_m','error_m','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-SOIL-001 failed.');
end
fprintf('VAL-MESH-SOIL-001 passed.\n');
