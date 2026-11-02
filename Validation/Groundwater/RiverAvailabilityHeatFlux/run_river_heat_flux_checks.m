% VAL-GW-RIVAVAIL-001, VAL-ET-GHEAT-001, VAL-GW-SINGLE-001.
%
% Three small parity and robustness checks of the grid groundwater and ET:
%   VAL-GW-RIVAVAIL-001 (audit L18). HydroPol2D-Python lets a gaining river
%     take at most the groundwater stored before this step's recharge
%     (core/engine.py); MATLAB counted the recharge too. One cell holds 2 mm
%     of groundwater, of which 1.5 mm is this step's recharge, under a fast
%     riverbed: the river may take 0.5 mm.
%   VAL-ET-GHEAT-001 (audit L11). Python reads the stations' ground heat
%     flux G; ETP_model forced G = 0. One station with G = 1.5 MJ/m2/day must
%     give the reference ET of Evapotranspiration with that G.
%   VAL-GW-SINGLE-001. hp2d_boussinesq_adaptive carried a single-precision
%     step when the heads were single, so an interval could end a rounding
%     error short and the loop never closed (found by the detention-pond
%     port). With single heads over a 135.127 s interval the call must
%     return, and match the double-precision call to single precision.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

% L18
area = 100; Sy = 0.2;
soil = struct('riverbed_hydraulic_conductivity_m_s', 1e-2, 'cell_riverbed_area_m2', area, ...
    'riverbed_thickness_m', 0.1, 'Sy', Sy);
bedrock = 0; head = bedrock + 0.01;             % 2 mm of drainable water
[~, surface_after] = hp2d_structured_river_exchange(head, 0, bedrock, -1, soil, true, 3600, area, ...
    false, 0.0015);
taken_mm = surface_after;                         % surface gain [mm]
[~, surface_old] = hp2d_structured_river_exchange(head, 0, bedrock, -1, soil, true, 3600, area, ...
    false, 0);
river_ok = abs(taken_mm - 0.5) < 1e-9 && abs(surface_old - 2) < 1e-9;

% L11
grid_x = [5 15]; grid_y = [15 5]; dem = 100 * ones(2); mask = false(2);   % coordinate vectors
[etp, ~] = ETP_model(1, 180, 5, 5, grid_x, grid_y, 30, 18, 24, 2, 60, 1.5, dem, 45 * ones(2), ...
    0.16, 0.23, mask);
[expected_etp, ~] = Evapotranspiration(dem, 24 * ones(2), 30 * ones(2), 18 * ones(2), 180, ...
    45 * ones(2), 2 * ones(2), 60 * ones(2), 0.16, 0.23, 1.5 * ones(2));
[zero_g_etp, ~] = Evapotranspiration(dem, 24 * ones(2), 30 * ones(2), 18 * ones(2), 180, ...
    45 * ones(2), 2 * ones(2), 60 * ones(2), 0.16, 0.23, zeros(2));
heat_ok = max(abs(etp(:) - expected_etp(:))) < 1e-12 && max(abs(expected_etp(:) - zero_g_etp(:))) > 1e-3;

% Single-precision interval
n = 6; head0 = single(repmat(linspace(9.8, 9.2, n), 3, 1)); bed = single(zeros(3, n));
sy = single(0.2 * ones(3, n)); k = single(1e-4 * ones(3, n)); surface = single(10 * ones(3, n));
soil_depth = single(10 * ones(3, n)); domain = true(3, n); none = false(3, n);
args = @(h, b, s, kk, z, d) {135.127, single(135.127), 1, 1e-3, 1e-3, 10, 10, h, b, s, ...
    zeros(3, n, 'like', h), kk, z, d, domain, none, [], none, [], none, 0.25};
single_args = args(head0, bed, sy, k, surface, soil_depth);
started = tic;
h_single = hp2d_boussinesq_adaptive(single_args{:});
seconds_single = toc(started);
double_args = args(double(head0), double(bed), double(sy), double(k), double(surface), double(soil_depth));
h_double = hp2d_boussinesq_adaptive(double_args{:});
single_ok = seconds_single < 60 && max(abs(double(h_single(:)) - h_double(:))) < 1e-4;

checks = ["VAL-GW-RIVAVAIL-001"; "VAL-ET-GHEAT-001"; "VAL-GW-SINGLE-001"];
names = ["River takes groundwater stored before this step's recharge"; ...
    "Reference ET uses the stations' ground heat flux"; ...
    "Adaptive groundwater interval closes with single heads"];
ok = [river_ok; heat_ok; single_ok];
status = repmat("fail", 3, 1); status(ok) = "pass";
PassFail = table(checks, names, status, true(3, 1), 'VariableNames', ...
    {'case_id','case_name','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
fprintf('river taken %.4f mm (old rule %.4f mm); ETP with G %.5f vs %.5f mm/day; single run %.2f s\n', ...
    taken_mm, surface_old, etp(1), expected_etp(1), seconds_single);
disp(PassFail);
if ~all(ok)
    error('River availability / ground heat flux / single-precision checks failed.');
end
fprintf('All three checks passed.\n');
