% VAL-GW-LATERAL-DT-001: when the lateral groundwater update runs.
%
% HydroPol2D-Python runs the lateral groundwater update after
% max(min(target, stable, suggested), min_dt) seconds, with target
% groundwater_target_dt_s = 86400 s and min_dt groundwater_min_dt_s = 60 s,
% or as soon as the head has moved by groundwater_max_head_change_m
% (0.25 m), or on the last step (modules/groundwater.py,
% core/layered_vertical.py). MATLAB's legacy scheme used min(stable,
% suggested), with no one-day cap and no 60 s floor (audit M9); the layered
% scheme defaulted to 3600 s (M10) and had no head-change trigger (M11).
% hp2d_groundwater_lateral_due now holds the rule for both schemes.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

% name, elapsed, target, stable, suggested, min, head change, max change, final, due, interval
cases = {
    'long stable step is capped at one day',      86399, 86400, 1e6, inf, 60, [], [], false, false, 86400
    'one day reached',                             86400, 86400, 1e6, inf, 60, [], [], false, true,  86400
    'stability limit shortens the interval',        3600, 86400, 3600, inf, 60, [], [], false, true,  3600
    'a tiny stable step is floored at 60 s',          59, 86400,   10, inf, 60, [], [], false, false, 60
    'the adaptive suggestion shortens it',          7200, 86400,  inf, 7200, 60, [], [], false, true, 7200
    'a 0.3 m head change triggers early',             10, 86400,  inf, inf, 60, [0; -0.3], 0.25, false, true, 86400
    'a 0.2 m head change does not',                   10, 86400,  inf, inf, 60, [0; 0.2], 0.25, false, false, 86400
    'the last step always updates',                   10, 86400,  inf, inf, 60, [], [], true, true, 86400};
n = size(cases, 1);
due = false(n, 1); interval = zeros(n, 1); ok = false(n, 1);
for k = 1:n
    c = cases(k, :);
    [due(k), interval(k)] = hp2d_groundwater_lateral_due(c{2}, c{3}, c{4}, c{5}, c{6}, c{7}, c{8}, c{9});
    ok(k) = due(k) == c{10} && abs(interval(k) - c{11}) < 1e-9;
end
cfg = hp2d_layered_subsurface_config(struct('flag_layered_subsurface', 1));
default_target = cfg.groundwater_target_dt_s;
default_head_change = cfg.groundwater_max_head_change_m;
legacy_uses_rule = contains(fileread(fullfile(functions_dir, 'Groundwater_Module.m')), 'hp2d_groundwater_lateral_due');
layered_uses_rule = contains(fileread(fullfile(functions_dir, 'Layered_Subsurface_Module.m')), 'hp2d_groundwater_lateral_due');

Cases = table(string(cases(:, 1)), due, interval, ok, 'VariableNames', {'case','due','interval_s','ok'});
passed = all(ok) && default_target == 86400 && default_head_change == 0.25 && ...
    legacy_uses_rule && layered_uses_rule;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-GW-LATERAL-DT-001", "Lateral groundwater update interval", nnz(ok), n, ...
    default_target, default_head_change, legacy_uses_rule, layered_uses_rule, status, true, ...
    'VariableNames', {'case_id','case_name','cases_passed','cases','layered_default_target_s', ...
    'layered_default_max_head_change_m','legacy_uses_rule','layered_uses_rule','status','report_ready'});
writetable(Cases, fullfile(out_dir, 'Lateral_Update_Cases.csv'));
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(Cases); disp(PassFail);
if ~passed
    error('VAL-GW-LATERAL-DT-001 failed.');
end
fprintf('VAL-GW-LATERAL-DT-001 passed.\n');
