% VAL-GW-BOUSS-STEP-001: the Boussinesq step ends at the requested time and
% reports all of its exfiltration.
%
% Boussinesq_2D_explicit splits a step where the Darcy-velocity Courant
% limit asks for it. Its counter was reset on every inner step, so the last
% inner step was never shortened: the loop ran past the requested step, and
% the exfiltration it returned covered only the last inner step
% (audit L14). HydroPol2D-Python takes one predictor-corrector per
% requested step (modules/groundwater.py). Now the last inner step ends
% exactly at T, and exfiltration and seepage are averaged over T.
% The case is a 3 x 8 strip of 10 m cells with a fast aquifer (K = 5e-3
% m/s), so the Courant limit splits a two-day step. Recharge pushes the
% table above the soil surface. The water balance over T must close: recharge =
% storage change + exfiltration.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

ny = 3; nx = 8; dx = 10; T = 2 * 86400;
z0 = zeros(ny, nx);                         % aquifer bottom
soil = 2 * ones(ny, nx);                    % surface 2 m above the bottom
head = repmat(linspace(1.9, 0.5, nx), ny, 1);  % steep gradient, near the surface upstream
Sy = 0.2 * ones(ny, nx); K = 5e-3 * ones(ny, nx);
R = 2e-7 * ones(ny, nx); R(:, 1:3) = 1e-6;  % m/s, enough to exfiltrate
mask = true(ny, nx); none = false(ny, nx);
[h_new, ~, ~, q_exf, ~, mass_error] = Boussinesq_2D_explicit(T, dx, dx, head, z0, Sy, R, K, ...
    none, zeros(ny, nx), z0 + soil, z0 + soil, 0.5, soil, mask, none, [], none, [], none);
area = dx^2;
recharge = sum(R, 'all') * T * area;
storage_change = sum(Sy .* (max(h_new - z0, 0) - max(head - z0, 0)), 'all') * area;
exfiltration = sum(q_exf, 'all') * T * area;
residual = abs(recharge - storage_change - exfiltration) / recharge;
passed = residual < 1e-9 && exfiltration > 0 && abs(mass_error) < 1e-6 * recharge;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-GW-BOUSS-STEP-001", "Boussinesq inner steps end at T", recharge, ...
    storage_change, exfiltration, residual, mass_error, status, true, ...
    'VariableNames', {'case_id','case_name','recharge_m3','storage_change_m3', ...
    'exfiltration_m3','balance_residual_fraction','reported_mass_error_m3','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-GW-BOUSS-STEP-001 failed.');
end
fprintf('VAL-GW-BOUSS-STEP-001 passed.\n');
