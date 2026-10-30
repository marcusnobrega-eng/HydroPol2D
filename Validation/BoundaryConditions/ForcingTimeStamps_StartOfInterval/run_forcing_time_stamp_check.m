% VAL-FORCING-STAMP-001: forcing value i is the rate for [t_i, t_i+1).
%
% update_spatial_BC integrates the rain and the inflow over each model step.
% The forcing tables start at 0 min, as date-based tables do:
%     time (min)      0    10    20
%     rain (mm/h)     6    12     3
%     inflow (mm/h)  60     0    30
% Value i is the rate from t_i to t_i+1; the last value lasts one more
% interval (to 30 min), then the rate is 0. This is the convention of
% HydroPol2D-Python (forcing/raster_forcing.py, forcing/gauge_rainfall.py).
% Model steps of 7 min from 0 to 42 min cross every interval edge. The
% reference integrates the same rule exactly.
% Before the fix MATLAB read value i as the rate for (t_i-1, t_i]: row 1
% covered (-10, 0] and was never used, and every rate came one interval
% early.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

addpath(functions_dir, '-begin');

times_min = [0; 10; 20];
rain_mm_h = [6; 12; 3];
inflow_mm_h = [60; 0; 30];
step_edges = (0:7:42)';

[rain_model, inflow_model] = run_steps(functions_dir, times_min, rain_mm_h, inflow_mm_h, step_edges);
rain_reference = reference_depths(times_min, rain_mm_h, step_edges);
inflow_reference = reference_depths(times_min, inflow_mm_h, step_edges);

Steps = table(step_edges(1:end-1), step_edges(2:end), rain_model, rain_reference, ...
    inflow_model, inflow_reference, ...
    'VariableNames', {'t0_min','t1_min','rain_model_mm','rain_reference_mm', ...
    'inflow_model_mm','inflow_reference_mm'});

rain_error = max(abs(rain_model - rain_reference));
inflow_error = max(abs(inflow_model - inflow_reference));
total_rain_error = abs(sum(rain_model) - (6 + 12 + 3) * 10 / 60);
passed = rain_error < 1e-12 && inflow_error < 1e-12 && total_rain_error < 1e-12;
status = "fail";
if passed; status = "pass"; end

PassFail = table("VAL-FORCING-STAMP-001", "Forcing value i covers [t_i, t_i+1)", ...
    rain_error, inflow_error, total_rain_error, status, true, ...
    'VariableNames', {'case_id','case_name','max_rain_error_mm', ...
    'max_inflow_error_mm','total_rain_error_mm','status','report_ready'});

writetable(Steps, fullfile(out_dir, 'Forcing_Time_Stamp_Steps.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(Steps);
disp(PassFail);
if ~passed
    error('VAL-FORCING-STAMP-001 failed.');
end
fprintf('VAL-FORCING-STAMP-001 passed.\n');

function [rain_mm, inflow_mm] = run_steps(functions_dir, times_min, rain_mm_h, inflow_mm_h, step_edges)
ny = 2; nx = 2;
flags = struct('flag_rainfall', 1, 'flag_spatial_rainfall', 0, ...
    'flag_input_rainfall_map', 0, 'flag_satellite_rainfall', 0, ...
    'flag_real_time_satellite_rainfall', 0, 'flag_inflow', 1, ...
    'flag_stage_hydrograph', 0, 'flag_ETP', 0, 'flag_input_ETP_map', 0); %#ok<NASGU>
Elevation_Properties = struct('elevation_cell', zeros(ny, nx)); %#ok<NASGU>
idx_nan = false(ny, nx); %#ok<NASGU>
Hydro_States = struct(); %#ok<NASGU>
BC_States = struct();
Rainfall_Parameters = struct('time_rainfall', times_min, 'intensity_rainfall', rain_mm_h); %#ok<NASGU>
Inflow_Parameters = struct('n_stream_gauges', 1, 'time_inflow', times_min, ...
    'inflow_hydrograph_rate', inflow_mm_h); %#ok<NASGU>
Wshed_Properties = struct('inflow_cells', ones(ny, nx, 1)); %#ok<NASGU>
script = fullfile(functions_dir, 'update_spatial_BC.m');

n_steps = numel(step_edges) - 1;
rain_mm = zeros(n_steps, 1);
inflow_mm = zeros(n_steps, 1);
for k = 1:n_steps
    t_previous = step_edges(k); %#ok<NASGU>
    t = step_edges(k+1); %#ok<NASGU>
    run(script);
    rain_mm(k) = BC_States.delta_p_agg;
    inflow_mm(k) = BC_States.delta_inflow_agg(1);
end
end

function depth_mm = reference_depths(times_min, rate_mm_h, step_edges)
dt_force = times_min(2) - times_min(1);
starts = times_min;
ends = [times_min(2:end); times_min(end) + dt_force];
n_steps = numel(step_edges) - 1;
depth_mm = zeros(n_steps, 1);
for k = 1:n_steps
    t0 = step_edges(k); t1 = step_edges(k+1);
    for i = 1:numel(times_min)
        overlap = max(0, min(t1, ends(i)) - max(t0, starts(i)));
        depth_mm(k) = depth_mm(k) + rate_mm_h(i) * overlap / 60;
    end
end
end
