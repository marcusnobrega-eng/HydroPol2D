function Summary = run_single_workbook_checks()
%RUN_SINGLE_WORKBOOK_CHECKS The single case workbook against the older multi-workbook case.
%   Examples/VTilted_Excel holds the same case twice: HydroPol2D_Inputs.xlsx (one
%   workbook, built by tools/workbook_builder) and General_Data.xlsx with its companion
%   workbooks. Checks:
%   1) tables: each class-table and series tab reads (hp2d_read_case_table) exactly as
%      its old workbook does (readtable);
%   2) runs: the case run from the single workbook, with every old workbook deleted,
%      gives the same outlet hydrograph, final depths and volume error as the old layout;
%   3) the tab is really read: doubling the rain on the Rain series tab doubles the rain
%      volume of the run.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
addpath(fullfile(root, 'HydroPol2D_Functions'));
example = fullfile(root, 'Examples', 'VTilted_Excel');
single = fullfile(example, 'Input_Data_Sheets', 'HydroPol2D_Inputs.xlsx');
sheets = fullfile(example, 'Input_Data_Sheets');
names = {}; details = {};

% ---- 1) tables ---------------------------------------------------------------------
pairs = {'Soil classes', 'SOIL_parameters.xlsx', {}; ...
    'Land use classes', 'LULC_parameters.xlsx', {'Range', 'A2', 'VariableNamingRule', 'preserve'}; ...
    'Rain series', 'Rainfall_Intensity_Data.xlsx', {}};
for i = 1:size(pairs, 1)
    old = readtable(fullfile(sheets, pairs{i, 2}), pairs{i, 3}{:});
    new = hp2d_read_case_table(single, pairs{i, 1}, '', pairs{i, 3}{:});
    assert(isequaln(old, new), '%s: the tab does not read as %s', pairs{i, 1}, pairs{i, 2});
end
names{end + 1} = 'tables'; details{end + 1} = 'Soil classes, Land use classes, Rain series read as their old workbooks';

% ---- 2) runs -------------------------------------------------------------------------
A = run_case(example, 'General_Data.xlsx', 'old', {}, []);
B = run_case(example, 'HydroPol2D_Inputs.xlsx', 'single', ...
    {'General_Data.xlsx', 'SOIL_parameters.xlsx', 'LULC_parameters.xlsx', ...
    'Rainfall_Intensity_Data.xlsx', 'Voronoi_Settings.xlsx'}, []);
dq = max(abs(A.outlet_m3_s - B.outlet_m3_s));
dd = max(abs(A.depth_mm(:) - B.depth_mm(:)), [], 'omitnan');
assert(numel(A.outlet_m3_s) == numel(B.outlet_m3_s) && dq == 0 && dd == 0 && ...
    A.volume_error_m3 == B.volume_error_m3, ...
    'single workbook run differs: outlet %.3g m3/s, depth %.3g mm, volume error %.4g vs %.4g m3', ...
    dq, dd, A.volume_error_m3, B.volume_error_m3);
names{end + 1} = 'same run'; details{end + 1} = sprintf(['old workbooks deleted; outlet peak %.4f m3/s, ', ...
    'volume error %.3g m3, identical in both layouts'], max(A.outlet_m3_s), A.volume_error_m3);

% ---- 3) the tab is read ----------------------------------------------------------------
C = run_case(example, 'HydroPol2D_Inputs.xlsx', 'single_x2', ...
    {'General_Data.xlsx', 'SOIL_parameters.xlsx', 'LULC_parameters.xlsx', ...
    'Rainfall_Intensity_Data.xlsx', 'Voronoi_Settings.xlsx'}, 2);
ratio = C.rain_m3/B.rain_m3;
assert(abs(ratio - 2) < 1e-9, 'doubling the Rain series tab gave %.6f x the rain', ratio);
names{end + 1} = 'tab is read'; details{end + 1} = sprintf('Rain series x2 -> rain volume x%.6f', ratio);

for i = 1:numel(names)
    fprintf('[PASS] %s: %s\n', names{i}, details{i});
end
Summary = table(string(names(:)), string(details(:)), 'VariableNames', {'check', 'detail'});
fprintf('Single-workbook checks: %d of 3 passed\n', height(Summary));
end

function R = run_case(example, workbook_name, tag, delete_files, rain_factor)
%RUN_CASE Copy the example, optionally delete old workbooks or scale the rain, run it.
work = fullfile(tempdir, ['hp2d_single_workbook_' tag]);
if isfolder(work), rmdir(work, 's'); end
copyfile(example, work);
if isfolder(fullfile(work, 'Outputs')), rmdir(fullfile(work, 'Outputs'), 's'); end
sheets = fullfile(work, 'Input_Data_Sheets');
for k = 1:numel(delete_files)
    delete(fullfile(sheets, delete_files{k}));
end
workbook = fullfile(sheets, workbook_name);
if ~isempty(rain_factor)
    C = readcell(workbook, 'Sheet', 'Rain series', 'Range', 'A1:C2000');   % column C = intensity
    last = 2;
    for r = 3:size(C, 1)
        if isnumeric(C{r, 3}) && isscalar(C{r, 3})
            C{r, 3} = rain_factor*C{r, 3};
            last = r;
        end
    end
    writecell(C(3:last, 3), workbook, 'Sheet', 'Rain series', 'Range', 'C3');
end
setenv('HYDROPOL_RUN_MODE', 'excel');
setenv('HYDROPOL_INPUT_EXCEL_FILE', workbook);
setenv('HYDROPOL_EXPORT_ROOT_DIR', fullfile(work, 'Outputs'));
setenv('HYDROPOL_SKIP_POSTPROCESS', '1');
setenv('HP2D_SINGLE_RESULT', fullfile(work, 'result.mat'));
run(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'HydroPol2D_V115.m'));
R = struct();
R.outlet_m3_s = double(gather(outlet_states.outlet_hydrograph));
R.depth_mm = double(gather(depths.d_t));
R.volume_error_m3 = double(gather(volume_error));
R.rain_m3 = double(gather(BC_States.inflow_volume));
save(getenv('HP2D_SINGLE_RESULT'), 'R');
end
