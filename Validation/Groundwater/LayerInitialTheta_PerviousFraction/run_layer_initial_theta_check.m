% VAL-SOIL-THETA0-001: per-layer initial moisture rasters and the pervious
% fraction.
%
% HydroPol2D-Python starts each soil layer with f_p (theta - theta_r)
% thickness, where theta is that layer's initial-moisture raster or, when
% the raster is absent, the cell's initial theta (modules/soil.py,
% core/preprocess.py). MATLAB used the three layer rasters only when all
% three existed, and then left out f_p: a cell that is half paved started
% with twice its water, often full.
% The check builds the storages of six cells (pervious fractions 1 to 0)
% with only the root-zone raster given, through
% initialize_layered_soil_storage, and compares them with Python's
% initialize_layered_soil_storage on the same numbers
% (python_reference.csv, made by make_python_reference.py).

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

fp = [1; 0.8; 0.5; 0.25; 0.05; 0];
n = numel(fp);
theta_r = 0.05 * ones(n, 1); theta_s = 0.45 * ones(n, 1);
theta_i = [0.30; 0.25; 0.40; 0.20; 0.35; 0.30];
theta_root = [0.40; 0.10; 0.30; 0.44; 0.20; 0.30];          % the only layer raster
near = 0.05 * ones(n, 1); root = [0.45; 0.30; 0.60; 0.40; 0.25; 0.5]; trans = [1.5; 1.2; 0.9; 1.0; 1.6; 1.5];
capacity = @(thickness) fp .* (theta_s - theta_r) .* thickness * 1000;
Layers = struct('pervious_fraction', fp, 'near_surface_thickness_m', near, ...
    'root_zone_thickness_m', root, 'transmission_thickness_m', trans, ...
    'near_surface_capacity_mm', capacity(near), 'root_zone_capacity_mm', capacity(root), ...
    'transmission_capacity_mm', capacity(trans));
Soil_Properties = struct('theta_i', theta_i, 'theta_r', theta_r, 'theta_sat', theta_s, ...
    'Layers', Layers, 'initial_layer_theta', struct('root_zone', theta_root));
idx_nan = false(n, 1);
Soil_Properties = initialize_layered_soil_storage(Soil_Properties, idx_nan);
L = Soil_Properties.Layers;

writetable(table(fp, near, root, trans, Layers.near_surface_capacity_mm, ...
    Layers.root_zone_capacity_mm, Layers.transmission_capacity_mm, theta_r, theta_i, ...
    theta_root, theta_i, 'VariableNames', {'pervious_fraction','near_thickness_m', ...
    'root_thickness_m','transmission_thickness_m','near_capacity_mm','root_capacity_mm', ...
    'transmission_capacity_mm','theta_r','theta_near','theta_root','theta_transmission'}), ...
    fullfile(case_dir, 'layer_inputs.csv'));
reference = readtable(fullfile(case_dir, 'python_reference.csv'));
error_mm = max(abs([L.near_surface_storage_mm, L.root_zone_storage_mm, L.transmission_storage_mm] - ...
    [reference.near_storage_mm, reference.root_storage_mm, reference.transmission_storage_mm]), [], 'all');
passed = error_mm < 1e-9;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-SOIL-THETA0-001", "Per-layer initial theta with the pervious fraction", ...
    error_mm, status, true, 'VariableNames', {'case_id','case_name', ...
    'max_storage_error_vs_python_mm','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(table(fp, L.near_surface_storage_mm, reference.near_storage_mm, L.root_zone_storage_mm, ...
    reference.root_storage_mm, 'VariableNames', {'fp','near_matlab','near_python','root_matlab','root_python'}));
disp(PassFail);
if ~passed
    error('VAL-SOIL-THETA0-001 failed.');
end
fprintf('VAL-SOIL-THETA0-001 passed.\n');
