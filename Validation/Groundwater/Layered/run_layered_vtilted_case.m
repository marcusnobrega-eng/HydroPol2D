function Result = run_layered_vtilted_case(variant)
%RUN_LAYERED_VTILTED_CASE Full-model run of the V-tilted example with the layered subsurface.
%   variant: 'local' (no lateral flow), 'lateral' (flag_baseflow = 1), 'feddes'
%   (local + Feddes + ET off), 'legacy' (old groundwater scheme), or 'deep' /
%   'deep_legacy' (deep aquifer on, one subbasin, started from a 0.5 mm/day baseflow),
%   'deep_empty' (as 'deep', starting empty), 'lateral_scaled' (as 'lateral' with soil and
%   regolith lateral-conductivity multipliers 2 and 0.5), 'no_dtb' (as 'local' with DTB_path
%   empty and a 4.5 m 'DTB (m)' for every class in SOIL_parameters, the last fallback for
%   the depth to bedrock), 'explicit_flux' (as 'local' with
%   flag_implicit_water_table_flux = 0), or 'soil_rasters' (as 'local' with ThetaS_path =
%   0.41 and KsatNearSurface_path = 25 mm/h maps over the SOIL class values), or 'gauges'
%   (surface only, with four gauges and their discharge cross-sections; see
%   Validation/Gauges/run_gauge_cross_section_checks.m).
%   Copies Examples/VTilted_Excel to tempdir, enables
%   infiltration, groundwater, capillary rise and the layered scheme in the copy's
%   Flags sheet (double precision), runs run_HydroPol2D.m without post-processing and
%   returns the mass-balance diagnostics and the final layered state.
arguments
    variant (1,:) char {mustBeMember(variant, {'local', 'lateral', 'feddes', 'legacy', ...
        'deep', 'deep_legacy', 'deep_empty', 'lateral_scaled', 'no_dtb', 'explicit_flux', 'soil_rasters', 'gauges'})} = 'local'
end
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..', '..');
work = fullfile(tempdir, ['hp2d_layered_vtilted_' variant]);
if isfolder(work), rmdir(work, 's'); end
copyfile(fullfile(root, 'Examples', 'VTilted_Excel'), work);
workbook = fullfile(work, 'Input_Data_Sheets', 'General_Data.xlsx');
% Groundwater runs need DTB_path in Excel mode: point it at a 6 m DTB raster in the
% temporary copy (keys are added if a workbook lacks them).
G = readcell(workbook, 'Sheet', 'General_Data', 'Range', 'A1:AZ400');
% A 6 m depth to bedrock gives the column deep (rock) slots below a 1.2 m soil base.
addpath(fullfile(root, 'third_party', 'topotoolbox_lite'));
dtb_file = fullfile(work, 'Static', 'DTB_layered.tif');
D = GRIDobj(fullfile(work, 'Static', 'DTB.tif'));
D.Z(isfinite(D.Z)) = 6;
GRIDobj2geotiff(D, dtb_file);
% Numeric layered settings are General_Data keys (the Flags sheet reads only flag_* rows).
keys = {'DTB_path', dtb_file; 'GW_table_path', ''; 'subsurface_soil_base_depth_m', 1.2};
deep = startsWith(variant, 'deep');
if strcmp(variant, 'no_dtb')
    % No DTB raster and no subsurface_depth_to_bedrock_m: the SOIL sheet's 'DTB (m)'.
    keys{1, 2} = '';
    soil_book = fullfile(work, 'Input_Data_Sheets', 'SOIL_parameters.xlsx');
    S = readcell(soil_book, 'Range', 'A1:Z200');
    [hr, hc] = find(strcmpi(string(S), 'DTB (m)'), 1);
    for r = hr + 1:size(S, 1)
        if isnumeric(S{r, hc}) && isscalar(S{r, hc}), S{r, hc} = 4.5; end
    end
    S(cellfun(@(x) isa(x, 'missing'), S)) = {''};
    writecell(S, soil_book, 'Range', 'A1');
end
if strcmp(variant, 'soil_rasters')
    theta_file = fullfile(work, 'Static', 'ThetaS.tif');
    ksat_file = fullfile(work, 'Static', 'KsatNearSurface.tif');
    T = D; T.Z(isfinite(T.Z)) = 0.41; GRIDobj2geotiff(T, theta_file);
    K = D; K.Z(isfinite(K.Z)) = 25; GRIDobj2geotiff(K, ksat_file);
    keys = [keys; {'ThetaS_path', theta_file; 'KsatNearSurface_path', ksat_file}];
end
if strcmp(variant, 'lateral_scaled')
    keys = [keys; {'soil_lateral_conductivity_multiplier', 2; ...
        'regolith_lateral_conductivity_multiplier', 0.5}];
end
if deep
    % One deep-aquifer subbasin returning to every active cell, started at q0.
    B = D; B.Z(:) = 1;  % no NaN: bilinear resampling would spread it into the domain
    basin_file = fullfile(work, 'Static', 'DeepAquiferSubbasins.tif');
    weight_file = fullfile(work, 'Static', 'DeepAquiferChannelWeight.tif');
    GRIDobj2geotiff(B, basin_file);
    GRIDobj2geotiff(B, weight_file);
    q0 = 0.5;
    if strcmp(variant, 'deep_empty'), q0 = ''; end
    keys = [keys; {'DeepAquiferSubbasins_path', basin_file; ...
        'DeepAquiferChannelWeight_path', weight_file; ...
        'deep_aquifer_initial_baseflow_mm_day', q0; 'deep_aquifer_tau_days', 100}];
end
for i = 1:size(keys, 1)
    [r, c] = find(strcmpi(string(G), keys{i, 1}), 1);
    if isempty(r)
        G(end + 1, :) = {missing}; %#ok<AGROW>
        r = size(G, 1); c = 1;
        G{r, 1} = keys{i, 1};
    end
    G{r, c + 1} = keys{i, 2};
end
if strcmp(variant, 'gauges')
    % Channel (column 21) flows south to the outlet (row 31). Cell centres in the DEM CRS:
    % x = (col - 0.5)*20, y = 620 - (row - 0.5)*20.
    % The old layout has the Design Storms block right of the gauge table, so the two
    % optional columns go in free columns of the same rows (the reader finds them
    % anywhere in the gauge header row).
    [gr, gc] = find(strcmpi(string(G), 'Gauge'), 1);
    table_rows = {'Gauge', 'Easting (m)', 'Northing (m)', 'Label Name';
        1, 410, 310, 'mid valley, full width';
        2, 410, 30, 'above the outlet, full width';
        3, 410, 310, 'mid valley, channel only';
        4, 410, 310, 'mid valley, even width'};
    section_cols = {'Cross-section width (cells)', 'Downstream direction';
        41, 'auto'; 41, 'auto'; 5, 'auto'; 4, 'S'};
    G(gr:end, gc:gc + 3) = {''};   % the example's 351 old gauges below must go too
    G(gr:gr + size(table_rows, 1) - 1, gc:gc + 3) = table_rows;
    G(gr:gr + size(section_cols, 1) - 1, 44:45) = section_cols;
end
G(cellfun(@(x) isa(x, 'missing'), G)) = {''};
writecell(G, workbook, 'Sheet', 'General_Data', 'Range', 'A1');
% Read and write the sheet from A1: its data start at B2, and readcell's detected range
% would shift every flag up and left on write-back, leaving the old last row behind
% (which then overrode the new value, e.g. flag_feddes).
C = readcell(workbook, 'Sheet', 'Flags', 'Range', 'A1:Z200');
% The gauge run is surface-only: the case's rain would otherwise all infiltrate and the
% outlet would carry nothing to compare the sections against.
subsurface = double(~strcmp(variant, 'gauges'));
settings = {'flag_infiltration', subsurface; 'flag_groundwater_modeling', subsurface; ...
    'flag_capillary_rise', subsurface; 'flag_single', 0; 'flag_GPU', 0; ...
    'flag_layered_subsurface', subsurface*double(~any(strcmp(variant, {'legacy', 'deep_legacy'}))); ...
    'flag_baseflow', double(any(strcmp(variant, {'lateral', 'lateral_scaled'}))); ...
    'flag_feddes', double(strcmp(variant, 'feddes')); ...
    'flag_deep_aquifer', double(deep); ...
    'flag_implicit_water_table_flux', double(~strcmp(variant, 'explicit_flux')); ...
    'flag_obs_gauges', double(strcmp(variant, 'gauges'))};
for i = 1:size(settings, 1)
    C = set_flag(C, settings{i, 1}, settings{i, 2});
end
C(cellfun(@(x) isa(x, 'missing'), C)) = {''};
writecell(C, workbook, 'Sheet', 'Flags', 'Range', 'A1');
setenv('HYDROPOL_RUN_MODE', 'excel');
setenv('HYDROPOL_INPUT_EXCEL_FILE', workbook);
setenv('HYDROPOL_EXPORT_ROOT_DIR', fullfile(work, 'Outputs'));
setenv('HYDROPOL_SKIP_POSTPROCESS', '1');
setenv('HYDROPOL_LAYERED_RESULT_FILE', fullfile(work, 'layered_result.mat'));
run(fullfile(root, 'run_HydroPol2D.m'));
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
    Result.kh_m_s = double(gather(Subsurface.column.kh_m_s));
    Result.column_base_depth_m = double(gather(Subsurface.column_base_depth_m));
    Result.implicit_water_table_flux = Subsurface.implicit_water_table_flux;
    Result.feddes_on = ~isempty(Subsurface.feddes);
    Result.slot_theta_s = double(gather(Subsurface.column.theta_s));
    Result.slot_kv_m_s = double(gather(Subsurface.column.kv_m_s));
    Result.storage_m3 = Wshed_Properties.Resolution^2*sum(hp2d_column_storage_mm(Subsurface.column, ...
        double(gather(BC_States.h_t(:))), double(gather(Subsurface.water_mm))), 'omitnan')/1000;
end
Result.total_recharge_m3 = 0;   % surface-only runs have no ledger
if exist('GW_States', 'var') && isfield(GW_States, 'total_recharge_m3'), Result.total_recharge_m3 = GW_States.total_recharge_m3; end
if exist('gauges', 'var') && isstruct(gauges) && isfield(gauges, 'section_discharge')
    Result.gauge_sections = gauges.sections;
    Result.gauge_section_discharge_m3_s = gauges.section_discharge;
    Result.gauge_cell_outflow_m3_s = double(gather(gauges.hydrograph_cell));
    Result.time_hydrograph_min = double(gather(running_control.time_hydrograph));
    Result.outlet_hydrograph_m3_s = double(gather(outlet_states.outlet_hydrograph));
end
Result.total_capillary_m3 = 0;   % surface-only runs have no ledger
if exist('GW_States', 'var') && isfield(GW_States, 'total_capillary_m3'), Result.total_capillary_m3 = GW_States.total_capillary_m3; end
Result.total_exfiltration_m3 = 0;   % surface-only runs have no ledger
if exist('GW_States', 'var') && isfield(GW_States, 'total_exfiltration_m3'), Result.total_exfiltration_m3 = GW_States.total_exfiltration_m3; end
Result.mass_balance_history = mass_balance_history;
Result.deep_initial_storage_m3 = [];
Result.deep_final_storage_m3 = [];
if exist('BC_States', 'var') && isfield(BC_States, 'deep_aquifer_initial_storage_m3')
    Result.deep_initial_storage_m3 = BC_States.deep_aquifer_initial_storage_m3;
    Result.deep_active_area_m2 = sum(isfinite(BC_States.deep_aquifer_subbasin_id(:)) & ...
        BC_States.deep_aquifer_subbasin_id(:) > 0)*Wshed_Properties.cell_area;
end
if exist('GW_States', 'var') && isfield(GW_States, 'deep_aquifer_storage_m3')
    Result.deep_final_storage_m3 = GW_States.deep_aquifer_storage_m3;
end
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
