% VAL-INT-LAIGAP-001: active cells without an LAI value keep their rain.
%
% A NaN LAI in an active cell made interceptionModel's canopy store NaN, so
% the cell's throughfall was 0 and its rain never reached the ground.
% HydroPol2D-Python fills such cells with the configured LAI fallback
% (static.lai, default 0) and warns (core/preprocess.py). MATLAB now does
% the same in preprocessing with hp2d_fill_lai_gaps (lai_fallback in the
% flags, default 0) (audit L10).
% A 2 x 2 LAI map with one gap in an active cell and one cell outside the
% domain, 5 mm of rain on an empty canopy: the gap cell gets all 5 mm as
% throughfall; the outside cell stays NaN.

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(case_dir, '..', '..', '..');
functions_dir = fullfile(repo_root, 'HydroPol2D_Functions');
out_dir = fullfile(case_dir, 'Outputs', 'Validation');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(functions_dir, '-begin');

lai = [2 NaN; 1 NaN];
idx_nan = logical([0 0; 0 1]);
warning('off', 'HydroPol2D:LAIFallback');
[filled_lai, filled] = hp2d_fill_lai_gaps(lai, idx_nan, 0);
warning('on', 'HydroPol2D:LAIFallback');
rain = 5 * ones(2);
[~, throughfall_old] = interceptionModel(rain, zeros(2), lai, zeros(2), 0.2, zeros(2));
[~, throughfall_new] = interceptionModel(rain, zeros(2), filled_lai, zeros(2), 0.2, zeros(2));
gap_throughfall = throughfall_new(1, 2);
passed = filled == 1 && filled_lai(1, 2) == 0 && isnan(filled_lai(2, 2)) && ...
    abs(gap_throughfall - 5) < 1e-12 && throughfall_old(1, 2) == 0;
status = "fail";
if passed; status = "pass"; end
PassFail = table("VAL-INT-LAIGAP-001", "Active cells without an LAI value keep their rain", ...
    filled, throughfall_old(1, 2), gap_throughfall, status, true, 'VariableNames', {'case_id', ...
    'case_name','filled_cells','gap_throughfall_unfilled_mm','gap_throughfall_mm','status','report_ready'});
writetable(PassFail, fullfile(out_dir, 'Pass_Fail.csv'));
disp(PassFail);
if ~passed
    error('VAL-INT-LAIGAP-001 failed.');
end
fprintf('VAL-INT-LAIGAP-001 passed.\n');
