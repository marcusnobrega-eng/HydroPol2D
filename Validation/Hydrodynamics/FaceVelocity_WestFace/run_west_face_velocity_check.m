% VAL-HYDRO-VELFACE-001: every cell sees the velocity on all four faces.
%
% refreshing_timestep builds the face velocities of each cell from the
% east-face (qout_right_t) and north-face (qout_up_t) discharges. The west
% face of cell j is the east face of cell j-1. MATLAB took the east face of
% cell j itself (vel_right(:,2:end) instead of vel_right(:,1:end-1)), so the
% west face dropped out of the cell velocity, the velocity raster, the wave
% speed and human risk mode 1. HydroPol2D-Python uses all four faces
% (routing/timestep.py).
% Here only the east faces of column 1 carry water, into column 2: column 2
% must see that flow on its west face.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

ny = 3; nx = 4; dx = 10; depth_m = 0.2; q_m3_s = 0.5;
velocities = run_refresh(functions_dir, ny, nx, dx, depth_m, q_m3_s);
v = q_m3_s / (dx * depth_m);                       % face speed, m/s
expected_raster = zeros(ny, nx); expected_raster(:, 1:2) = v;
expected_east = zeros(ny, nx); expected_east(:, 1:2) = v / 2;   % mean of the two x faces
raster_error = max(abs(velocities.velocity_raster - expected_raster), [], 'all');
component_error = max(abs(velocities.right_component - expected_east), [], 'all');
west_face_error = max(abs(velocities.vel_left(:, 2) + v));
passed = raster_error < 1e-12 && component_error < 1e-12 && west_face_error < 1e-12;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-HYDRO-VELFACE-001", "Cell velocity uses all four faces", ...
    raster_error, component_error, west_face_error, status, true, ...
    'VariableNames', {'case_id','case_name','velocity_raster_error_m_s', ...
    'x_component_error_m_s','west_face_error_m_s','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(velocities.velocity_raster); disp(PassFail);
if ~passed
    error('VAL-HYDRO-VELFACE-001 failed.');
end
fprintf('VAL-HYDRO-VELFACE-001 passed.\n');

function velocities = run_refresh(functions_dir, ny, nx, dx, depth_m, q_m3_s)
k = 1; %#ok<NASGU>
flags = struct('flag_CA', 0, 'flag_D8', 0, 'flag_inertial', 1, 'flag_diffusive', 0, ...
    'flag_kinematic', 0, 'flag_full_momentum', 0, 'flag_subgrid', 0, 'flag_overbanks', 0, ...
    'flag_waterquality', 0, 'flag_human_instability', 0, 'flag_adaptive_timestepping', 0); %#ok<NASGU>
Hf = depth_m * ones(ny, nx, 2); %#ok<NASGU>
q_mm_h = q_m3_s / dx^2 * 1000 * 3600;
flow_rate = struct('qout_right_t', [q_mm_h * ones(ny, 1), zeros(ny, nx - 1)], ...
    'qout_up_t', zeros(ny, nx), 'qout_left_t', zeros(ny, nx), 'qout_down_t', zeros(ny, nx)); %#ok<NASGU>
CA_States = struct('depth_tolerance', 1); %#ok<NASGU>
Wshed_Properties = struct('Resolution', dx); %#ok<NASGU>
running_control = struct('delta_time_save', 1, 'max_time_step', 60, 'min_time_step', 0.01, ...
    'time_step_change', 1, 'time_step_increments', 0.001); %#ok<NASGU>
Courant_Parameters = struct('alfa_max', 0.7, 'alfa_min', 0.3, 'v_threshold', 0.1, 'slope_alfa', 2); %#ok<NASGU>
depths = struct('d_t', depth_m * 1000 * ones(ny, nx)); %#ok<NASGU>
idx_nan = false(ny, nx); %#ok<NASGU>
time_step = 1/60; %#ok<NASGU>
velocities = struct();
try
    run(fullfile(functions_dir, 'refreshing_timestep.m'));
catch caught
    % Later parts of the script need more of the run's workspace; the face
    % velocities are complete by then.
    if ~isfield(velocities, 'velocity_raster') || ~isfield(velocities, 'right_component')
        rethrow(caught);
    end
end
end
