% VAL-HR-VEL-001: human risk uses the cell-centre speed of the step.
%
% HydroPol2D-Python evaluates both human-risk modes with the cell-centre
% speed: the mean of the net discharges through each pair of opposite faces,
% over the face width and the cell depth, in cells deeper than 1 mm, every
% step (core/engine.py, _cell_velocity_m_s). MATLAB used the largest face
% speed (mode 3) or the cell velocity of refreshing_timestep (mode 1), both
% refreshed only when the time step was recomputed, and mode 1 without the
% west face (audit M14, L4).
% A 1 x 3 strip of 10 m cells, 0.2 m deep: 0.2 m3/s enters the middle cell
% through its west face and 0.4 m3/s leaves through its east face. Its
% cell-centre speed is 0.5 (0.2 + 0.4) / 10 / 0.2 = 0.15 m/s; the largest
% face speed is 0.2 m/s. Human_Instability_Module (mode 1) must give the
% risk index of 0.15 m/s.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

dx = 10; h = 0.2;
to_mm_h = @(q) q / dx^2 * 1000 * 3600;
flow_rate = struct('qout_right_t', to_mm_h([0.2 0.4 0]), 'qout_left_t', zeros(1, 3), ...
    'qout_up_t', zeros(1, 3), 'qout_down_t', zeros(1, 3));
speed = hp2d_cell_centre_velocity(flow_rate, h * 1000 * ones(1, 3), dx, false(1, 3));
expected = [0.5 * (0 + 0.2), 0.5 * (0.2 + 0.4), 0.5 * (0.4 + 0)] / dx / h;
speed_error = max(abs(speed - expected));

P = struct('mu', 0.6, 'Cd', 1.0, 'ro_water', 1000, 'gravity', 9.81, 'weight', 75, 'w1', 0.35, 'w2', 0.2);
risk = run_module(flow_rate, h, dx, P);
F = P.weight * P.gravity; B = P.w1 * P.w2 * h * P.ro_water * P.gravity;
expected_risk = 0.5 * P.Cd * P.ro_water * P.w1 * h * expected.^2 / (P.mu * (F - B));
risk_error = max(abs(risk - expected_risk));
passed = speed_error < 1e-12 && risk_error < 1e-12;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-HR-VEL-001", "Human risk uses the cell-centre speed", speed(2), ...
    expected(2), speed_error, risk_error, status, true, 'VariableNames', {'case_id','case_name', ...
    'middle_cell_speed_m_s','expected_speed_m_s','speed_error_m_s','risk_error','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-HR-VEL-001 failed.');
end
fprintf('VAL-HR-VEL-001 passed.\n');

function risk = run_module(flow_rate, h, dx, P)
flags.flag_human_instability = 1; %#ok<STRNU>
k = 1; %#ok<NASGU>
DEM_raster.Z = zeros(1, 3); %#ok<STRNU>
depths.d_t = h * 1000 * ones(1, 3); %#ok<STRNU>
velocities.total_velocity = zeros(1, 3); %#ok<STRNU>   % the old mode-1 input, not used now
idx_nan = false(1, 3); %#ok<NASGU>
Wshed_Properties.Resolution = dx; %#ok<STRNU>
Human_Instability = struct('mu', P.mu, 'Cd', P.Cd, 'ro_water', P.ro_water, 'gravity', P.gravity, ...
    'weight_person', P.weight, 'width1_person', P.w1, 'width2_person', P.w2);
Human_Instability_Module;
risk = Human_Instability.risk_t;
end
