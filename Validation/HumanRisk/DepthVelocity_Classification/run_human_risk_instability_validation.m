clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Model', 'HydroPol2D_Functions');
addpath(functions_dir);

out_dir = fullfile(case_dir, 'Outputs', 'Validation');
fig_dir = fullfile(out_dir, 'Figures');
table_dir = fullfile(out_dir, 'Tables');
ensure_dir(out_dir);
ensure_dir(fig_dir);
ensure_dir(table_dir);

[Diagnostics, MassBalance, PassFail, SimpleGrid, ForceGrid, BoundaryTable, ConfusionTable, Params] = ...
    run_human_risk_cases();

writetable(Diagnostics, fullfile(out_dir, 'Metric_Summary.csv'));
writetable(MassBalance, fullfile(out_dir, 'Mass_Balance.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
writetable(SimpleGrid, fullfile(table_dir, 'P1-HR-001_simple_instability_grid.csv'));
writetable(ForceGrid, fullfile(table_dir, 'P1-HR-001_force_balance_grid.csv'));
writetable(BoundaryTable, fullfile(table_dir, 'P1-HR-001_boundary_checks.csv'));
writetable(ConfusionTable, fullfile(table_dir, 'P1-HR-001_confusion_table.csv'));
writetable(struct2table(Params, 'AsArray', true), fullfile(table_dir, 'P1-HR-001_parameters.csv'));

make_figures(SimpleGrid, ForceGrid, fig_dir);

disp(Diagnostics);
disp(MassBalance);
disp(PassFail);

function [Diagnostics, MassBalance, PassFail, SimpleGrid, ForceGrid, BoundaryTable, ConfusionTable, Params] = ...
        run_human_risk_cases()
Params = base_parameters();

[diag_simple, mass_simple, pass_simple, SimpleGrid] = simple_instability_case(Params);
[diag_force, mass_force, pass_force, ForceGrid, BoundaryTable, ConfusionTable] = force_balance_case(Params);

Diagnostics = [diag_simple; diag_force];
MassBalance = [mass_simple; mass_force];
PassFail = [pass_simple; pass_force];
end

function Params = base_parameters()
Params = struct();
Params.case_id = "P1-HR-001";
Params.mu = 0.6;
Params.Cd = 1.0;
Params.ro_water = 1000;
Params.gravity = 9.81;
Params.weight_person_kg = 75;
Params.width1_person_m = 0.35;
Params.width2_person_m = 0.20;
Params.slope_deg = 10;
Params.mass_kg = 75;
Params.height_m = 1.75;
Params.body_width_m = 0.35;
Params.body_depth_m = 0.20;
Params.tolerance = 1e-12;
end

function [Diag, Mass, Pass, Grid] = simple_instability_case(P)
case_id = "P1-HR-001A";
case_name = "Simple drag-friction instability index";

depth_values_m = [0, 0.05, 0.15, 0.30, 0.60, 1.10, 1.30];
velocity_values_m_s = [0, 0.5, 1.0, 1.5, 2.5];
[H, V] = ndgrid(depth_values_m, velocity_values_m_s);

model_risk = call_simple_instability_module(H, V, P);
expected_risk = independent_simple_instability(H, V, P);
err = model_risk - expected_risk;
finite_mask = isfinite(expected_risk) & isfinite(model_risk);

rmse = sqrt(mean(err(finite_mask).^2));
max_error = max(abs(err(finite_mask)));
nan_count = sum(~isfinite(model_risk), 'all');
classification_accuracy_pct = mean((model_risk(finite_mask) >= 1 - P.tolerance) == ...
    (expected_risk(finite_mask) >= 1 - P.tolerance)) * 100;

boundary_depth_m = P.weight_person_kg / ...
    (P.width1_person_m * P.width2_person_m * P.ro_water);
boundary_velocity_m_s = sqrt((2 * P.mu * ...
    (P.weight_person_kg * P.gravity - P.width1_person_m * P.width2_person_m * 0.30 * P.ro_water * P.gravity)) / ...
    (P.Cd * P.ro_water * P.width1_person_m * 0.30));
boundary_model = call_simple_instability_module([0, boundary_depth_m, boundary_depth_m + 1e-6, 0.30], ...
    [0, 0, 0, boundary_velocity_m_s], P);
boundary_expected = independent_simple_instability([0, boundary_depth_m, boundary_depth_m + 1e-6, 0.30], ...
    [0, 0, 0, boundary_velocity_m_s], P);
threshold_boundary_error_count = sum(abs(boundary_model - boundary_expected) > 1e-10);

passed = rmse < 1e-12 && max_error < 1e-10 && classification_accuracy_pct == 100 && ...
    threshold_boundary_error_count == 0 && nan_count == 0;

Grid = table(H(:), V(:), expected_risk(:), model_risk(:), err(:), ...
    expected_risk(:) >= 1 - P.tolerance, model_risk(:) >= 1 - P.tolerance, ...
    'VariableNames', {'depth_m','velocity_m_s','expected_risk_index','model_risk_index', ...
    'risk_error','expected_unstable','model_unstable'});

Diag = diagnostic_row(case_id, case_name, "independent_drag_friction_formula", rmse, max_error, ...
    classification_accuracy_pct, threshold_boundary_error_count, nan_count, passed, ...
    "HydroPol2D flag 1 risk index matches independent drag, buoyancy, and friction reference.");
Mass = mass_row(case_id, case_name, 0, 0, passed);
Pass = pass_row(case_id, case_name, passed);
end

function [Diag, Mass, Pass, Grid, BoundaryTable, ConfusionTable] = force_balance_case(P)
case_id = "P1-HR-001B";
case_name = "Detailed slide-topple-drowning classifier";

depth_values_m = unique([0, 0.05, 0.25, 0.50, 1.05, P.height_m * 13 / 16, ...
    P.height_m * 13 / 16 + 1e-6, P.height_m * cosd(P.slope_deg) / 2, 1.60]);
velocity_values_m_s = [0, 0.5, 1.0, 2.0, 3.0, 6.0];
[H, V] = ndgrid(depth_values_m, velocity_values_m_s);

model_class = Human_risk(3, V, H, P.ro_water, P.gravity, P.mu, P.Cd, P.slope_deg, ...
    P.mass_kg, P.height_m, P.body_width_m, P.body_depth_m);
expected_class = independent_force_balance_class(V, H, P);
class_error = model_class - expected_class;
classification_accuracy_pct = mean(model_class(:) == expected_class(:)) * 100;
threshold_boundary_error_count = boundary_error_count(P);
nan_count = sum(~isfinite(model_class), 'all');
rmse = sqrt(mean(class_error(:).^2));
max_error = max(abs(class_error(:)));
passed = classification_accuracy_pct == 100 && threshold_boundary_error_count == 0 && ...
    nan_count == 0 && rmse == 0 && max_error == 0;

Grid = table(H(:), V(:), expected_class(:), model_class(:), class_error(:), ...
    class_label(expected_class(:)), class_label(model_class(:)), ...
    'VariableNames', {'depth_m','velocity_m_s','expected_class_code','model_class_code', ...
    'class_error','expected_class','model_class'});

BoundaryTable = build_boundary_table(P);
ConfusionTable = build_confusion_table(expected_class(:), model_class(:));

Diag = diagnostic_row(case_id, case_name, "independent_force_balance_classifier", rmse, max_error, ...
    classification_accuracy_pct, threshold_boundary_error_count, nan_count, passed, ...
    "HydroPol2D flag 3 slide/topple/drowning map matches independent force-balance classifier.");
Mass = mass_row(case_id, case_name, 0, 0, passed);
Pass = pass_row(case_id, case_name, passed);
end

function model_risk = call_simple_instability_module(H, V, P)
flags.flag_human_instability = 1;
k = 1;
DEM_raster.Z = zeros(size(H));
depths.d_t = H * 1000;
velocities.total_velocity = V;
idx_nan = false(size(H));
Human_Instability.mu = P.mu;
Human_Instability.Cd = P.Cd;
Human_Instability.ro_water = P.ro_water;
Human_Instability.gravity = P.gravity;
Human_Instability.weight_person = P.weight_person_kg;
Human_Instability.width1_person = P.width1_person_m;
Human_Instability.width2_person = P.width2_person_m;
Human_Instability_Module;
model_risk = Human_Instability.risk_t;
end

function risk = independent_simple_instability(H, V, P)
F_person = P.weight_person_kg * P.gravity;
F_buoy = P.width1_person_m * P.width2_person_m .* H * P.ro_water * P.gravity;
available_friction = P.mu * max(F_person - F_buoy, 0);
hydro_force = max(0.5 * P.Cd * P.ro_water * P.width1_person_m .* H .* V.^2, 0);
risk = zeros(size(H));
idx_supported = available_friction > 0;
risk(idx_supported) = hydro_force(idx_supported) ./ available_friction(idx_supported);
risk(~idx_supported & H > 0) = 1;
risk = min(max(risk, 0), 1);
end

function map = independent_force_balance_class(v, h, P)
slope = P.slope_deg;
m = P.mass_kg;
y = P.height_m;
w = P.body_width_m;
d = P.body_depth_m;
ro = P.ro_water;
g = P.gravity;
mu = P.mu;
cd = P.Cd;

sub_temp_1 = h ./ cosd(slope) <= y ./ 2;
V_s1 = sub_temp_1 .* (((2 * h .* d^2) .* pi()) ./ (cosd(slope) .* 4));
A_s1 = sub_temp_1 .* ((2 * h .* d^2) ./ cosd(slope));
X_gs1 = sub_temp_1 .* (d / 2);
Y_gs1 = sub_temp_1 .* (h ./ (2 * cosd(slope)));

sub_temp_2 = h ./ cosd(slope) > y / 2;
V_s2 = sub_temp_2 .* (((y * pi() * d^2) / 4) + ((h ./ cosd(slope) - y / 2) .* pi * ((w^2) / 4)));
A_s2 = sub_temp_2 .* ((y * d) + ((h ./ cosd(slope)) - (y / 2)) * w);
X_gs2 = sub_temp_2 .* (((2 * ((pi() * d^2) / 4) * (y / 2) * (d / 2)) + ...
    ((pi() * w^2) / 4) * ((h ./ cosd(slope)) - (y / 2)) * (w / 2)) ./ ...
    (2 * ((pi() * d^2) / 4) * (y / 2) + ((pi() * w^2) / 4) * ((h ./ cosd(slope)) - (y / 2))));
Y_gs2 = sub_temp_2 .* ((2 * ((pi() * d^2) / 4) * ((y^2) / 8) + ...
    ((pi() * w^2) / 4) * ((h ./ cosd(slope)) - (y / 2)) .* ...
    ((0.5 * (h ./ cosd(slope) - (y / 2))) + (y / 2))) ./ ...
    (2 * ((pi() * d^2) / 4) * (y / 2) + ((pi() * w^2) / 4) * ((h ./ cosd(slope)) - (y / 2))));

L1 = 0.5 * ro * cd * (sind(90 - slope).^2) .* cosd(90 - slope) .* (v.^2) .* A_s1;
L2 = 0.5 * ro * cd * (sind(90 - slope).^2) .* cosd(90 - slope) .* (v.^2) .* A_s2;
L = L1 + L2;
D1 = 0.5 * ro * cd * (sind(90 - slope).^3) .* (v.^2) .* A_s1;
D2 = 0.5 * ro * cd * (sind(90 - slope).^3) .* (v.^2) .* A_s2;
Drag = D1 + D2;

Wp = m * g * sind(slope);
Wn = m * g * cosd(slope);
Bn = ro * g * V_s1 .* cosd(slope) + ro * g * V_s2 .* cosd(slope);
X_gs = X_gs1 + X_gs2;
Y_gs = Y_gs1 + Y_gs2;
T = mu * (Wn - Bn - L);
M = Drag .* (h ./ 2) + (Wp .* (7 / 12) .* y .* cosd(slope)) + ...
    (Bn .* ((X_gs ./ cosd(slope)) + (Y_gs .* sind(slope)))) + ...
    L .* ((X_gs ./ cosd(slope)) + ((h ./ 2) .* tan(slope)));

slide = double(Drag + Wp > T) * 1;
topple = double(M > Wn .* ((((5 / 6) * d) ./ cosd(slope)) + (((7 / 12) * y) ./ sind(slope)))) * 2;
drowning = double(h > ((13 / 16) * y)) * 3;
map = max(max(slide, topple), drowning);
end

function count = boundary_error_count(P)
T = build_boundary_table(P);
count = sum(T.model_class_code ~= T.expected_class_code);
end

function T = build_boundary_table(P)
geometry_depth = P.height_m * cosd(P.slope_deg) / 2;
drowning_depth = P.height_m * 13 / 16;
depth_m = [0; geometry_depth; drowning_depth; drowning_depth + 1e-6];
velocity_m_s = [0; 1.0; 0; 0];
expected_class_code = independent_force_balance_class(velocity_m_s, depth_m, P);
model_class_code = Human_risk(3, velocity_m_s, depth_m, P.ro_water, P.gravity, P.mu, P.Cd, ...
    P.slope_deg, P.mass_kg, P.height_m, P.body_width_m, P.body_depth_m);
boundary_name = ["dry_zero"; "geometry_split_inclusive"; "drowning_exact_open"; "drowning_above"];
expected_class = class_label(expected_class_code);
model_class = class_label(model_class_code);
T = table(boundary_name, depth_m, velocity_m_s, expected_class_code, model_class_code, ...
    expected_class, model_class);
end

function T = build_confusion_table(expected_class, model_class)
classes = (0:3)';
expected_count = zeros(numel(classes), 1);
model_count = zeros(numel(classes), 1);
correct_count = zeros(numel(classes), 1);
for i = 1:numel(classes)
    c = classes(i);
    expected_count(i) = sum(expected_class == c);
    model_count(i) = sum(model_class == c);
    correct_count(i) = sum(expected_class == c & model_class == c);
end
class_code = classes;
class_name = class_label(classes);
T = table(class_code, class_name, expected_count, model_count, correct_count);
end

function labels = class_label(codes)
labels = strings(size(codes));
labels(codes == 0) = "stable";
labels(codes == 1) = "sliding";
labels(codes == 2) = "toppling";
labels(codes == 3) = "drowning";
end

function row = diagnostic_row(case_id, case_name, evidence_type, rmse, max_error, ...
        classification_accuracy_pct, threshold_boundary_error_count, nan_count, passed, metric_note)
row = table(case_id, case_name, evidence_type, rmse, max_error, classification_accuracy_pct, ...
    threshold_boundary_error_count, nan_count, passed, metric_note, ...
    'VariableNames', {'case_id','case_name','evidence_type','rmse','max_error', ...
    'classification_accuracy_pct','threshold_boundary_error_count','nan_count','passed','metric_note'});
end

function row = mass_row(case_id, case_name, diagnostic_input_volume_m3, residual_m3, passed)
row = table(case_id, case_name, diagnostic_input_volume_m3, residual_m3, passed, ...
    'VariableNames', {'case_id','case_name','diagnostic_input_volume_m3','residual_m3','passed'});
end

function row = pass_row(case_id, case_name, passed)
status = "fail";
if passed
    status = "pass";
end
report_ready = passed;
row = table(case_id, case_name, status, report_ready, ...
    'VariableNames', {'case_id','case_name','status','report_ready'});
end

function make_figures(SimpleGrid, ForceGrid, fig_dir)
simple_depths = unique(SimpleGrid.depth_m);
simple_velocities = unique(SimpleGrid.velocity_m_s);
simple_risk = reshape(SimpleGrid.model_risk_index, numel(simple_depths), numel(simple_velocities));
simple_error = reshape(SimpleGrid.risk_error, numel(simple_depths), numel(simple_velocities));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 430]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imagesc(simple_velocities, simple_depths, simple_risk);
set(gca, 'YDir', 'normal');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
title('Model instability index');
colorbar;
caxis([0 1]);
nexttile;
imagesc(simple_velocities, simple_depths, simple_error);
set(gca, 'YDir', 'normal');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
title('Model minus reference');
colorbar;
exportgraphics(fig, fullfile(fig_dir, 'P1_HR_001_SIMPLE_INSTABILITY.png'), 'Resolution', 200);
close(fig);

force_depths = unique(ForceGrid.depth_m);
force_velocities = unique(ForceGrid.velocity_m_s);
force_class = reshape(ForceGrid.model_class_code, numel(force_depths), numel(force_velocities));
force_error = reshape(ForceGrid.class_error, numel(force_depths), numel(force_velocities));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 430]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imagesc(force_velocities, force_depths, force_class);
set(gca, 'YDir', 'normal');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
title('Model class: 0 stable, 1 slide, 2 topple, 3 drowning');
colorbar;
caxis([0 3]);
nexttile;
imagesc(force_velocities, force_depths, force_error);
set(gca, 'YDir', 'normal');
xlabel('Velocity (m/s)');
ylabel('Depth (m)');
title('Model minus reference class');
colorbar;
exportgraphics(fig, fullfile(fig_dir, 'P1_HR_001_FORCE_BALANCE.png'), 'Resolution', 200);
close(fig);
end

function ensure_dir(path_in)
if ~exist(path_in, 'dir')
    mkdir(path_in);
end
end
