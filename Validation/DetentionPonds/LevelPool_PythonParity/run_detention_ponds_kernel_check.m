% VAL-POND-001K: the MATLAB detention-pond kernels against HydroPol2D-Python.
%
% The case of make_python_reference.py (kernel_case.json) runs through
% hp2d_prepare_detention_ponds, hp2d_detention_pond_forcing and
% hp2d_detention_pond_exchange for 600 steps of 10 s and is compared, step by step,
% with python_reference_steps.csv (see README.md). A few seconds; the release
% runner's fast profile. run_detention_ponds_check.m adds the full-model runs.
%
% Writes Outputs/Validation/Kernel_Parity.csv and Kernel_Pass_Fail.csv.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
addpath(fullfile(repo_root, 'HydroPol2D_Functions'), '-begin');
addpath(case_dir, '-begin');
summary_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(summary_dir, 'dir'); mkdir(summary_dir); end

K = detention_kernel_parity(case_dir);
Kernel = struct2table(rmfield(K, {'steps', 'prepared'}));
Kernel.passed = K.geometry_max_error <= 1e-10 && K.max_volume_error_m3 <= 1e-9 && ...
    K.max_stage_error_m <= 1e-10 && K.max_discharge_error_m3_s <= 1e-10 && ...
    K.pump_switch_mismatches == 0 && abs(K.loop_water_balance_residual_m3) <= 1e-9 && ...
    K.seepage_recharge_mismatch_m3 <= 1e-12 && K.pump_switches >= 2 && ...
    K.reverse_flow_steps > 0 && K.spillway_flow_steps > 0;
writetable(Kernel, fullfile(summary_dir, 'Kernel_Parity.csv'));

note = sprintf(['MATLAB vs Python kernel (commit %s), 600 steps: max volume error %.2e m3, ', ...
    'stage %.2e m, discharge %.2e m3/s, %d pump-switch mismatches; loop balance %.2e m3.'], ...
    extractBefore(K.python_commit, 9), K.max_volume_error_m3, K.max_stage_error_m, ...
    K.max_discharge_error_m3_s, K.pump_switch_mismatches, K.loop_water_balance_residual_m3);
status = "FAIL";
if Kernel.passed, status = "PASS"; end
Pass = table("VAL-POND-001K", status, double(Kernel.passed), double(Kernel.passed), ...
    K.max_volume_error_m3, string(note), ...
    'VariableNames', {'case_id', 'status', 'passed', 'report_ready', 'metric', 'notes'});
writetable(Pass, fullfile(summary_dir, 'Kernel_Pass_Fail.csv'));
disp(Kernel);
disp(Pass);
