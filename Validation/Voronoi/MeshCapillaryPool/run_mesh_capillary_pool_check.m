% VAL-MESH-CAP-001: the HRUs of a mesh cell share its groundwater for
% capillary rise.
%
% HydroPol2D-Python lets each HRU take the capillary rise its deficit and
% conductivity allow, then scales the HRUs of a cell together when they
% ask for more than the cell's groundwater holds (voronoi.py,
% _capillary_rise). MATLAB capped every HRU at the cell's available depth
% instead (audit L27): with a dry HRU next to a wet one, the dry one could
% take only part of a store the other did not need.
% One 100 m2 cell, two HRUs of 50 %: HRU 1 is dry with a fast soil, HRU 2
% is at field capacity. The cell holds 2 mm of groundwater, less than HRU 1
% asks for. The rise must take the whole 2 mm (0.2 m3), all of it into HRU 1.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

mesh = struct('n_cells', 1, 'cell_area', 100, 'surface_area', 100, 'surface_bed', 10, 'cell_bed', 10);
head = mesh.surface_bed - 1.0;                      % water table 1 m down
hydrology = struct('hru_fraction', [0.5 0.5], 'pervious_fraction', 1, 'soil_depth_m', 2, ...
    'root_depth_m', 0.6, 'surface_layer_depth_m', 0.1, 'theta_r', 0.05, 'theta_sat', 0.45, ...
    'alpha_vg_1_m', 1.0, 'n_vg', 1.5, 'ksat_m_s', [1e-4 1e-4], 'initial_soil_saturation', [0 0.9], ...
    'capillary_rise_enabled', true);
state = Voronoi_Hydrology_Initialize(mesh, hydrology, head);
available_m3 = 0.002 * mesh.cell_area;              % 2 mm of groundwater
groundwater = struct('enabled', true, 'head', head, 'storage', available_m3, ...
    'pending_exchange_m3', 0, 'bottom', mesh.surface_bed - 5);
before = state.root_storage_m .* state.area_m2;
[~, state] = Voronoi_Hydrology_Step(mesh, 0, state, groundwater, 0, 0, 0, 3600);
taken = state.cumulative_capillary_m .* state.area_m2;
total_taken = sum(taken);
pool_error = abs(total_taken - available_m3);
passed = pool_error < 1e-12 && taken(2) < 1e-12;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-MESH-CAP-001", "HRUs of a mesh cell share its groundwater", ...
    available_m3, total_taken, taken(1), taken(2), pool_error, status, true, ...
    'VariableNames', {'case_id','case_name','available_m3','capillary_m3', ...
    'hru1_m3','hru2_m3','pool_error_m3','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-MESH-CAP-001 failed.');
end
fprintf('VAL-MESH-CAP-001 passed.\n');
