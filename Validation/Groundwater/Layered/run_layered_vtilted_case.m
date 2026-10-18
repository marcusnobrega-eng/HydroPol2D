function Result = run_layered_vtilted_case(variant)
%RUN_LAYERED_VTILTED_CASE Full-model run of the V-tilted example with the layered subsurface.
%   variant: 'local' (no lateral flow), 'lateral' (flag_baseflow = 1) or 'feddes'
%   (local + Feddes + ET off). Copies Examples/VTilted_Excel to tempdir, enables
%   infiltration, groundwater, capillary rise and the layered scheme in the copy's
%   Flags sheet (double precision), runs HydroPol2D_V115.m without post-processing and
%   returns the mass-balance diagnostics and the final layered state.
arguments
    variant (1,:) char {mustBeMember(variant, {'local', 'lateral', 'feddes', 'legacy'})} = 'local'
end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..', '..');
work = fullfile(tempdir, ['hp2d_layered_vtilted_' variant]);
if isfolder(work), rmdir(work, 's'); end
copyfile(fullfile(root, 'Examples', 'VTilted_Excel'), work);
workbook = fullfile(work, 'Input_Data_Sheets', 'General_Data.xlsx');
% Groundwater runs need DTB_path in Excel mode: point it at a 6 m DTB raster in the
% temporary copy (keys are added if a workbook lacks them).
G = readcell(workbook, 'Sheet', 'General_Data');
% A 6 m depth to bedrock gives the column deep (rock) slots below a 1.2 m soil base.
addpath(fullfile(root, 'third_party', 'topotoolbox_lite'));
dtb_file = fullfile(work, 'Static', 'DTB_layered.tif');
D = GRIDobj(fullfile(work, 'Static', 'DTB.tif'));
D.Z(isfinite(D.Z)) = 6;
GRIDobj2geotiff(D, dtb_file);
% Numeric layered settings are General_Data keys (the Flags sheet reads only flag_* rows).
keys = {'DTB_path', dtb_file; 'GW_table_path', ''; 'subsurface_soil_base_depth_m', 1.2};
for i = 1:size(keys, 1)
    [r, c] = find(strcmpi(string(G), keys{i, 1}), 1);
    if isempty(r)
        G(end + 1, :) = {missing}; %#ok<AGROW>
        r = size(G, 1); c = 1;
        G{r, 1} = keys{i, 1};
    end
    G{r, c + 1} = keys{i, 2};
end
G(cellfun(@(x) isa(x, 'missing'), G)) = {''};
writecell(G, workbook, 'Sheet', 'General_Data');
C = readcell(workbook, 'Sheet', 'Flags');
settings = {'flag_infiltration', 1; 'flag_groundwater_modeling', 1; 'flag_capillary_rise', 1; ...
    'flag_single', 0; 'flag_GPU', 0; 'flag_layered_subsurface', double(~strcmp(variant, 'legacy')); ...
    'flag_baseflow', double(strcmp(variant, 'lateral')); 'flag_feddes', double(strcmp(variant, 'feddes'))};
for i = 1:size(settings, 1)
    C = set_flag(C, settings{i, 1}, settings{i, 2});
end
C(cellfun(@(x) isa(x, 'missing'), C)) = {''};
writecell(C, workbook, 'Sheet', 'Flags');
setenv('HYDROPOL_RUN_MODE', 'excel');
setenv('HYDROPOL_INPUT_EXCEL_FILE', workbook);
setenv('HYDROPOL_EXPORT_ROOT_DIR', fullfile(work, 'Outputs'));
setenv('HYDROPOL_SKIP_POSTPROCESS', '1');
setenv('HYDROPOL_LAYERED_RESULT_FILE', fullfile(work, 'layered_result.mat'));
run(fullfile(root, 'HydroPol2D_V115.m'));
% The model script clears the workspace; recover the paths from the environment.
result_file = getenv('HYDROPOL_LAYERED_RESULT_FILE');
Result = struct();
Result.volume_error_m3 = double(gather(volume_error));
Result.module_errors_m3 = double(gather(cum_errors_m3));
Result.head_m = double(gather(BC_States.h_t));
Result.z0_m = double(gather(BC_States.z0));
Result.n_slots = 0;
if exist('Subsurface', 'var') && isstruct(Subsurface)
    Result.n_slots = Subsurface.n_slots;
    Result.storage_m3 = Wshed_Properties.Resolution^2*sum(hp2d_column_storage_mm(Subsurface.column, ...
        double(gather(BC_States.h_t(:))), double(gather(Subsurface.water_mm))), 'omitnan')/1000;
end
Result.total_recharge_m3 = GW_States.total_recharge_m3;
Result.total_capillary_m3 = GW_States.total_capillary_m3;
Result.total_exfiltration_m3 = GW_States.total_exfiltration_m3;
Result.mass_balance_history = mass_balance_history;
save(result_file, 'Result');
fprintf('Layered V-tilted (%s): %d slots, volume error %.3e m3, module errors %s m3\n', ...
    getenv('HYDROPOL_EXPORT_ROOT_DIR'), Result.n_slots, Result.volume_error_m3, ...
    mat2str(Result.module_errors_m3, 3));
end

function C = set_flag(C, name, value)
%SET_FLAG Set a name/value pair in a Flags sheet cell array (append a row if absent).
[r, c] = find(strcmpi(string(C), name), 1);
if isempty(r)
    r = size(C, 1) + 1; c = 1;
    C(r, :) = {missing};
    C{r, 1} = name;
end
C{r, c + 1} = value;
end
