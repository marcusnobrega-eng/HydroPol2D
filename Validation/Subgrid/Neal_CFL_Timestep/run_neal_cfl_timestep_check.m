% VAL-SUBGRID-NEAL-CFL-001: the Neal channel scheme gets a CFL time step.
%
% A row of ten 50 m cells with a 5 m wide, 3 m deep channel. Cells 1-5 flow
% in bank (2 m above the bed, channel discharge 20 m3/s); cells 6-10 flow over
% the banks (3.5 m above the bed, 20 m3/s in the channel and 15 m3/s on the
% floodplain). refreshing_timestep must return
%     dt = alfa_min * dx / max(|v| + sqrt(g h)),
% with the cell discharge the mean of its two faces, divided by the channel
% width x depth in bank and by the cell width x depth over the banks
% (estimate_neal_subgrid_dt_s in HydroPol2D-Python). The reference below is
% computed cell by cell, independently of the model helper.
% Before the fix the Neal solver returned no face depth (Hf = 0), so every
% rule saw still water and the step stayed at max_time_step.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(functions_dir, '-begin');

[dt_model_s, dt_reference_s, max_dt_s] = run_case(functions_dir);

% refreshing_timestep floors the step to time_step_increments (0.001 s).
error_s = abs(dt_model_s - dt_reference_s);
passed = dt_model_s < max_dt_s && error_s <= 0.001 + 1e-12;
status = "fail";
if passed; status = "pass"; end

Results = table("VAL-SUBGRID-NEAL-CFL-001", "Neal channel CFL time step", ...
    dt_model_s, dt_reference_s, max_dt_s, error_s, status, true, ...
    'VariableNames', {'case_id','case_name','dt_model_s','dt_reference_s', ...
    'max_time_step_s','abs_error_s','status','report_ready'});
writetable(Results, fullfile(out_dir, 'Pass_Fail.csv'));
disp(Results);
if ~passed
    error('VAL-SUBGRID-NEAL-CFL-001 failed.');
end
fprintf('VAL-SUBGRID-NEAL-CFL-001 passed.\n');

function [dt_model_s, dt_reference_s, max_dt_s] = run_case(functions_dir)
ny = 1; nx = 10;
dx = 50;
width = 5;
bank = 3;
g = 9.81;
alfa_min = 0.4;
max_dt_s = 60;

h = [2 * ones(1, 5), 3.5 * ones(1, 5)];          % depth above the bed [m]
Qc = zeros(ny, nx, 2);
Qf = zeros(ny, nx, 2);
Qc(1, 1:nx-1, 1) = 20;                             % east faces [m3/s]
Qf(1, 6:nx-1, 1) = 15;

flags = struct('flag_CA', 0, 'flag_D8', 0, 'flag_inertial', 1, ...
    'flag_diffusive', 0, 'flag_kinematic', 0, 'flag_full_momentum', 0, ...
    'flag_timestep', 1, 'flag_waterquality', 0, 'flag_GPU', 0, ...
    'flag_subgrid', 1, 'flag_overbanks', 1); %#ok<NASGU>
running_control = struct('delta_time_save', 1, 'max_time_step', max_dt_s, ...
    'min_time_step', 0.1, 'time_step_increments', 0.001);
Courant_Parameters = struct('alfa_max', 0.7, 'alfa_min', alfa_min); %#ok<NASGU>
CA_States = struct('depth_tolerance', 1e-3); %#ok<NASGU>
Wshed_Properties = struct('Resolution', dx, 'cell_area', dx^2, ...
    'River_Width', width * ones(ny, nx), 'River_Depth', bank * ones(ny, nx)); %#ok<NASGU>
depths = struct('d_t', 1000 * h); %#ok<NASGU>
idx_nan = false(ny, nx); %#ok<NASGU>
velocities = struct(); %#ok<NASGU>
k = 2; %#ok<NASGU>
time_step = 1 / 60; %#ok<NASGU>

% What the Neal solver hands to the main loop: face outflows in mm/h and no
% face depth.
Q = Qc + Qf;
outflow = Q ./ dx^2 * 1000 * 3600;
flow_rate = struct();
flow_rate.qout_right_t = outflow(:,:,1);
flow_rate.qout_left_t = -[zeros(ny,1), outflow(:,1:end-1,1)];
flow_rate.qout_up_t = outflow(:,:,2);
flow_rate.qout_down_t = -[outflow(2:end,:,2); zeros(1,nx)]; %#ok<STRNU>
Hf = zeros(ny, nx, 3); %#ok<NASGU>

run(fullfile(functions_dir, 'refreshing_timestep.m'));
dt_model_s = running_control.time_calculation_routing;

% Independent reference, cell by cell.
max_speed = 0;
for j = 1:nx
    q_east = Q(1, j, 1);
    if j > 1
        q_west = Q(1, j-1, 1);
    else
        q_west = 0;
    end
    q_cell = 0.5 * (q_east + q_west);
    if h(j) <= bank
        section = width * h(j);
    else
        section = dx * h(j);
    end
    speed = abs(q_cell) / section + sqrt(g * h(j));
    max_speed = max(max_speed, speed);
end
dt_reference_s = min(max_dt_s, alfa_min * dx / max_speed);
end
