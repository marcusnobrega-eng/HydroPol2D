function Results = advanced_case3_saturation_excess_analysis(Opt)
%% ========================================================================
% ADVANCED CASE 3 SATURATION-EXCESS ANALYSIS
% ========================================================================
%
% Updated for the revised Case 3:
%   Case3_shallow_groundwater_saturation_excess
%
% Revised physical setup:
%   - DTB / soil-aquifer thickness = 0.30 m
%   - initial water-table depth    = 0.05 m below land surface
%   - rainfall                     = 20 mm/h for 48 h
%   - simulation duration          = 48 h
%
% Purpose:
%   Diagnose whether near-surface groundwater causes rapid saturation-excess
%   overland flow and whether the outlet hydrograph approaches the
%   no-infiltration rainfall-runoff benchmark.
%
% This script reads HydroPol2D temporary map chunks:
%   Paths.Temp/save_map_hydro_1.mat
%   Paths.Temp/save_map_hydro_2.mat
%   ...
%
% Minimum required map:
%   Maps.Hydro.d
%
% Strongly recommended maps:
%   Maps.Hydro.I_t
%   Maps.Hydro.S_rem
%   Maps.Hydro.zwt
%   Maps.Hydro.h_t
%   Maps.Hydro.GW_depth
%   Maps.Hydro.f
%   Maps.Hydro.C
%   Maps.Hydro.recharge_mm_h
%
% If zwt/h_t/GW_depth is missing, the script uses the initial UZ capacity.
% For the updated Case 3, that initial capacity is based on zwt0 = 0.05 m.
%
% The analysis excludes the computational perimeter from vadose/GW metrics
% to avoid the perimeter artifacts observed in previous diagnostics.
% ========================================================================

%% ------------------------------------------------------------------------
% 0) USER OPTIONS
% -------------------------------------------------------------------------
if nargin < 1 || isempty(Opt)
    Opt = struct();
end

Opt = set_default(Opt, 'workspaceFile', fullfile(pwd, 'modeled_results.mat'));
Opt = set_default(Opt, 'caseRoot', pwd);
Opt = set_default(Opt, 'outputRoot', fullfile(pwd, 'Advanced_Case03_NearSurfaceGW_Analysis'));

Opt = set_default(Opt, 'caseFolderName', 'Case3_shallow_groundwater_saturation_excess');
Opt = set_default(Opt, 'rain_mm_h', []);
Opt = set_default(Opt, 'rain_duration_h', 48);
Opt = set_default(Opt, 'sim_duration_h', 48);

Opt = set_default(Opt, 'DTB_m', 0.30);
Opt = set_default(Opt, 'zwt0_m', []);
Opt = set_default(Opt, 'expected_fill_time_h', []);

Opt = set_default(Opt, 'depth_threshold_mm', 0.5);
Opt = set_default(Opt, 'Srem_threshold_mm', 0.5);
Opt = set_default(Opt, 'zwt_threshold_m', 0.005);
Opt = set_default(Opt, 'fill_threshold_pct', 99.0);
Opt = set_default(Opt, 'f_threshold_mm_h', 0.25);

Opt = set_default(Opt, 'snapshot_hours', [0 0.5 1 1.5 2 3 4 6 8 12 24 36 48]);
Opt = set_default(Opt, 'write_geotiff', true);
Opt = set_default(Opt, 'show_figures', false);
Opt = set_default(Opt, 'exclude_perimeter', true);

fprintf('\n============================================================\n');
fprintf('Advanced Case 3 Near-Surface Groundwater Analysis\n');
fprintf('Case root: %s\n', Opt.caseRoot);
fprintf('Workspace: %s\n', Opt.workspaceFile);
fprintf('Output: %s\n', Opt.outputRoot);
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% 1) LOAD MAIN WORKSPACE METADATA
% -------------------------------------------------------------------------
if ~isfile(Opt.workspaceFile)
    error('Workspace file not found:\n%s', Opt.workspaceFile);
end

varsAvailable = who('-file', Opt.workspaceFile);
varsWanted = { ...
    'Paths', 'flags', 'running_control', 'saver_memory_maps', 'idx_nan', ...
    'DEM_raster', 'Elevation_Properties', 'Soil_Properties', 'Wshed_Properties', ...
    'outlet_states', 'Rainfall_Parameters', 'BC_States', 'depths', ...
    'Hydro_States', 'cumulative_infiltration', 'max_GW_depth'};

varsToLoad = intersect(varsWanted, varsAvailable);
W = load(Opt.workspaceFile, varsToLoad{:});

requiredVars = {'running_control','DEM_raster','Elevation_Properties','Soil_Properties','Wshed_Properties'};
for ii = 1:numel(requiredVars)
    if ~isfield(W, requiredVars{ii})
        error('Required variable "%s" not found in workspace.', requiredVars{ii});
    end
end

if isfield(W, 'Paths') && isfield(W.Paths, 'Temp')
    tempDir = W.Paths.Temp;
else
    tempDir = fullfile(Opt.caseRoot, 'Temporary_Files');
end

if ~isfolder(tempDir)
    error('Temporary map folder not found:\n%s', tempDir);
end

%% ------------------------------------------------------------------------
% 2) OUTPUT FOLDERS
% -------------------------------------------------------------------------
Dirs.Root    = Opt.outputRoot;
Dirs.Figures = fullfile(Dirs.Root, 'Figures');
Dirs.Tables  = fullfile(Dirs.Root, 'Tables');
Dirs.Maps    = fullfile(Dirs.Root, 'Diagnostic_Rasters');
Dirs.Snaps   = fullfile(Dirs.Root, 'Snapshot_Figures');

ensure_folder(Dirs.Root);
ensure_folder(Dirs.Figures);
ensure_folder(Dirs.Tables);
ensure_folder(Dirs.Maps);
ensure_folder(Dirs.Snaps);

%% ------------------------------------------------------------------------
% 3) READ CASE CONFIG IF AVAILABLE
% -------------------------------------------------------------------------
caseConfigPath = fullfile(Opt.caseRoot, 'Case_Config.csv');
CaseConfig = table();

if isfile(caseConfigPath)
    CaseConfig = readtable(caseConfigPath);

    if isempty(Opt.rain_mm_h)
        Opt.rain_mm_h = read_table_var(CaseConfig, 'rain_mm_h', []);
    end
    if isempty(Opt.zwt0_m)
        Opt.zwt0_m = read_table_var(CaseConfig, 'zwt0_m', []);
    end
    Opt.DTB_m = read_table_var(CaseConfig, 'DTB_m', Opt.DTB_m);

    rd_min = read_table_var(CaseConfig, 'rain_duration_min', []);
    if ~isempty(rd_min) && isfinite(rd_min)
        Opt.rain_duration_h = rd_min / 60;
    end

    sd_min = read_table_var(CaseConfig, 'sim_duration_min', []);
    if ~isempty(sd_min) && isfinite(sd_min)
        Opt.sim_duration_h = sd_min / 60;
    end

    eft = read_table_var(CaseConfig, 'expected_fill_time_h', []);
    if isempty(Opt.expected_fill_time_h) && ~isempty(eft) && isfinite(eft)
        Opt.expected_fill_time_h = eft;
    end
end

if isempty(Opt.rain_mm_h)
    Opt.rain_mm_h = infer_rainfall_intensity(W);
end
if isempty(Opt.zwt0_m)
    Opt.zwt0_m = 0.05;
end

%% ------------------------------------------------------------------------
% 4) STATIC GEOMETRY AND SOIL FIELDS
% -------------------------------------------------------------------------
DEM = double(gather_if_needed(W.Elevation_Properties.elevation_cell));
[nRows, nCols] = size(DEM);

if isfield(W, 'idx_nan')
    idx_nan = logical(gather_if_needed(W.idx_nan));
else
    idx_nan = isnan(DEM) | isinf(DEM);
end

validMask = ~idx_nan & isfinite(DEM);

if Opt.exclude_perimeter
    idx_perimeter_diag = false(size(validMask));
    if isfield(W, 'Wshed_Properties') && isfield(W.Wshed_Properties, 'perimeter') && ...
            ~isempty(W.Wshed_Properties.perimeter) && isequal(size(W.Wshed_Properties.perimeter), size(validMask))
        idx_perimeter_diag = logical(gather_if_needed(W.Wshed_Properties.perimeter));
    else
        idx_perimeter_diag(1,:)   = true;
        idx_perimeter_diag(end,:) = true;
        idx_perimeter_diag(:,1)   = true;
        idx_perimeter_diag(:,end) = true;
    end
    validMask = validMask & ~idx_perimeter_diag;
end

SoilDepth = to_grid(W.Soil_Properties.Soil_Depth, nRows, nCols);
theta_s   = to_grid(W.Soil_Properties.theta_sat,  nRows, nCols);
theta_r   = to_grid(W.Soil_Properties.theta_r,    nRows, nCols);

if isfield(W.Soil_Properties, 'theta_i')
    theta_i = to_grid(W.Soil_Properties.theta_i, nRows, nCols);
elseif isfield(W.Soil_Properties, 'teta_i')
    theta_i = to_grid(W.Soil_Properties.teta_i, nRows, nCols);
else
    theta_i = 0.0997 + zeros(nRows, nCols);
end

Sy_local = max(theta_s - theta_r, 1e-6);

if isempty(Opt.expected_fill_time_h)
    mean_deficit = local_nanmean(Opt.zwt0_m .* (theta_s - theta_i) * 1000, validMask);
    Opt.expected_fill_time_h = mean_deficit / 8.5;
end

if isfield(W.Wshed_Properties, 'cell_area')
    cellArea = to_grid(W.Wshed_Properties.cell_area, nRows, nCols);
elseif isfield(W.Wshed_Properties, 'Resolution')
    cellArea = to_grid(W.Wshed_Properties.Resolution.^2, nRows, nCols);
else
    cellArea = ones(nRows, nCols);
end

areaTotal_m2 = nansum_local(cellArea(validMask));
Q_noinf_m3s = Opt.rain_mm_h / 1000 / 3600 * areaTotal_m2;

UZ_capacity_initial_mm = Opt.zwt0_m .* Sy_local * 1000;
I_initial_mm = Opt.zwt0_m .* max(theta_i - theta_r, 0) * 1000;
Srem_initial_mm = max(UZ_capacity_initial_mm - I_initial_mm, 0);

fprintf('Rainfall intensity        = %.4f mm/h\n', Opt.rain_mm_h);
fprintf('Rainfall duration         = %.3f h\n', Opt.rain_duration_h);
fprintf('Simulation duration       = %.3f h\n', Opt.sim_duration_h);
fprintf('DTB / soil thickness      = %.4f m\n', Opt.DTB_m);
fprintf('Initial zwt0              = %.4f m\n', Opt.zwt0_m);
fprintf('No-infiltration Qss       = %.6f m3/s\n', Q_noinf_m3s);
fprintf('Mean initial UZ capacity  = %.3f mm\n', local_nanmean(UZ_capacity_initial_mm, validMask));
fprintf('Mean initial I0           = %.3f mm\n', local_nanmean(I_initial_mm, validMask));
fprintf('Mean initial Srem         = %.3f mm\n', local_nanmean(Srem_initial_mm, validMask));
fprintf('Expected fill marker      = %.3f h\n\n', Opt.expected_fill_time_h);

%% ------------------------------------------------------------------------
% 5) TIME VECTOR AND MAP CHUNK INFO
% -------------------------------------------------------------------------
time_min = get_time_minutes(W.running_control, get_field_or(W, 'flags', struct()));
firstMapFile = find_map_file(tempDir, 1);
Tmp = load(firstMapFile, 'Maps');
Maps0 = Tmp.Maps;

if isfield(W, 'saver_memory_maps') && ~isempty(W.saver_memory_maps)
    nPerChunk = double(gather_if_needed(W.saver_memory_maps));
else
    d0 = get_map_by_candidates(Maps0, {'Hydro.d'}, 1);
    nPerChunk = size(d0, 3);
end

if isempty(time_min)
    nTimesTotal = count_total_saved_frames(tempDir);
    time_min = (0:nTimesTotal-1)' * 10;
else
    nTimesTotal = numel(time_min);
end

time_h_all = time_min(:) / 60;
idxUse = time_h_all <= Opt.sim_duration_h + 1e-9;
time_h = time_h_all(idxUse);
time_min = time_min(idxUse);
nTimes = numel(time_h);

snapHours = Opt.snapshot_hours(:);
snapHours = snapHours(snapHours <= max(time_h) + 1e-9);
snapshotIdx = zeros(numel(snapHours), 1);
for is = 1:numel(snapHours)
    [~, snapshotIdx(is)] = min(abs(time_h - snapHours(is)));
end
snapshotIdx = unique(snapshotIdx(:));

fprintf('Temporary map directory: %s\n', tempDir);
fprintf('Available time records: %d\n', nTimesTotal);
fprintf('Used time records:      %d\n', nTimes);
fprintf('Maps per chunk:         %d\n', nPerChunk);
fprintf('Snapshot hours:         %s\n\n', mat2str(round(time_h(snapshotIdx)',3)));

%% ------------------------------------------------------------------------
% 6) PREALLOCATE DIAGNOSTIC TABLES AND MAPS
% -------------------------------------------------------------------------
nanv = nan(nTimes,1);
TS = table();
TS.time_min = time_min(:);
TS.time_h   = time_h(:);

TS.mean_depth_m = nanv; TS.max_depth_m = nanv; TS.wet_area_pct = nanv;
TS.mean_I_t_mm = nanv; TS.max_I_t_mm = nanv;
TS.mean_fill_pct = nanv; TS.p05_fill_pct = nanv; TS.p50_fill_pct = nanv; TS.p95_fill_pct = nanv;
TS.mean_Srem_mm = nanv; TS.p05_Srem_mm = nanv; TS.p50_Srem_mm = nanv;
TS.sat_excess_area_pct = nanv; TS.storage_full_area_pct = nanv;
TS.mean_zwt_m = nanv; TS.min_zwt_m = nanv; TS.p05_zwt_m = nanv;
TS.mean_f_mm_h = nanv; TS.mean_C_mm_h = nanv; TS.mean_recharge_mm_h = nanv;
TS.surface_storage_m3 = nanv; TS.UZ_storage_m3 = nanv; TS.GW_storage_rel_m3 = nanv;
TS.total_storage_rel_m3 = nanv;
TS.cumulative_rain_m3 = nanv; TS.cumulative_outlet_m3 = nanv;
TS.mass_residual_m3 = nanv; TS.mass_residual_mm = nanv;
TS.outlet_Q_m3s = nanv; TS.noinf_Q_m3s = Q_noinf_m3s + zeros(nTimes,1); TS.Q_ratio_to_noinf = nanv;

MaxDepth_m = -inf(nRows, nCols);
MaxIt_mm = -inf(nRows, nCols);
MaxFill_pct = -inf(nRows, nCols);
MinSrem_mm = inf(nRows, nCols);
MinZwt_m = inf(nRows, nCols);
MaxRecharge_mm_h = -inf(nRows, nCols);
MaxSatMask = false(nRows, nCols);
TimeToSat_h = nan(nRows, nCols);

%% ------------------------------------------------------------------------
% 7) OUTLET AND RAINFALL VOLUMES
% -------------------------------------------------------------------------
if isfield(W, 'outlet_states') && isfield(W.outlet_states, 'outlet_hydrograph')
    q_out_raw = double(gather_if_needed(W.outlet_states.outlet_hydrograph(:)));
    if isfield(W.running_control, 'time_hydrograph')
        t_q_min = get_time_vector_generic(W.running_control.time_hydrograph);
    else
        t_q_min = linspace(min(time_min), max(time_min), numel(q_out_raw))';
    end
    q_out = interp1(t_q_min(:), q_out_raw(:), time_min(:), 'linear', 'extrap');
    TS.outlet_Q_m3s = q_out;
    TS.Q_ratio_to_noinf = q_out ./ max(Q_noinf_m3s, eps);
    TS.cumulative_outlet_m3 = cumtrapz(time_min(:) * 60, q_out);
else
    warning('outlet_states.outlet_hydrograph not found. Outlet volume diagnostics will be NaN.');
end

rain_elapsed_h = min(time_h(:), Opt.rain_duration_h);
TS.cumulative_rain_m3 = Opt.rain_mm_h .* rain_elapsed_h / 1000 * areaTotal_m2;

%% ------------------------------------------------------------------------
% 8) MAIN LOOP OVER SAVED MAPS
% -------------------------------------------------------------------------
currentStore = -1;
Maps = [];

for it = 1:nTimes
    storeID = ceil(it / nPerChunk);
    localID = it - (storeID - 1) * nPerChunk;

    if storeID ~= currentStore
        mapFile = find_map_file(tempDir, storeID);
        fprintf('Loading map chunk %d: %s\n', storeID, mapFile);
        Tmp = load(mapFile, 'Maps');
        Maps = Tmp.Maps;
        currentStore = storeID;
    end

    d_mm = get_map_by_candidates(Maps, {'Hydro.d','Hydro.depth','Hydro.depths'}, localID);
    I_t_mm = get_map_by_candidates(Maps, {'Hydro.I_t','Hydro.It','Hydro.soil_storage','Hydro.Soil_I_t'}, localID);
    f_mm_h = get_map_by_candidates(Maps, {'Hydro.f','Hydro.infiltration_rate','Hydro.Hydro_States_f'}, localID);
    C_mm_h = get_map_by_candidates(Maps, {'Hydro.C','Hydro.infiltration_capacity','Hydro.capacity'}, localID);
    recharge_mm_h = get_map_by_candidates(Maps, {'Hydro.recharge_mm_h','Hydro.recharge','Hydro.recharge_rate_mm_h'}, localID);
    Srem_saved_mm = get_map_by_candidates(Maps, {'Hydro.S_rem','Hydro.Srem','Hydro.remaining_storage_mm'}, localID);
    zwt_m = get_map_by_candidates(Maps, {'Hydro.zwt','Hydro.zwt_m','Hydro.water_table_depth'}, localID);
    GW_depth_m = get_map_by_candidates(Maps, {'Hydro.GW_depth','Hydro.gw_depth','Hydro.saturated_thickness'}, localID);
    h_t_m = get_map_by_candidates(Maps, {'Hydro.h_t','Hydro.GW_head','Hydro.water_table_elevation'}, localID);

    d_mm = clean_map(d_mm, validMask);
    I_t_mm = clean_map(I_t_mm, validMask);

    if all(isnan(d_mm(:)))
        error('Could not find Maps.Hydro.d in stored maps.');
    end
    if all(isnan(I_t_mm(:)))
        I_t_mm = nan(nRows, nCols);
    end

    hasZwt = ~all(isnan(zwt_m(:)));
    hasGWDepth = ~all(isnan(GW_depth_m(:)));
    hasHead = ~all(isnan(h_t_m(:)));

    if hasZwt
        zwt_current_m = clean_map(zwt_m, validMask);
    elseif hasGWDepth
        GW_depth_m = clean_map(GW_depth_m, validMask);
        zwt_current_m = min(max(SoilDepth - GW_depth_m, 0), SoilDepth);
        hasZwt = true;
    elseif hasHead
        h_t_m = clean_map(h_t_m, validMask);
        zwt_current_m = min(max(DEM - h_t_m, 0), SoilDepth);
        hasZwt = true;
    else
        zwt_current_m = nan(nRows, nCols);
    end

    if hasZwt
        UZ_capacity_mm = zwt_current_m .* Sy_local * 1000;
        capacityMode = "dynamic_actual_zwt";
    else
        UZ_capacity_mm = UZ_capacity_initial_mm;
        capacityMode = "initial_capacity_only";
    end

    if ~all(isnan(Srem_saved_mm(:)))
        Srem_mm = clean_map(Srem_saved_mm, validMask);
    else
        Srem_mm = max(UZ_capacity_mm - I_t_mm, 0);
    end

    fill_pct = 100 * I_t_mm ./ max(UZ_capacity_mm, 1e-9);
    fill_pct = min(max(fill_pct, 0), 150);

    wetMask = d_mm > Opt.depth_threshold_mm;

    if hasZwt || ~all(isnan(Srem_saved_mm(:)))
        storageFullMask = Srem_mm <= Opt.Srem_threshold_mm;
    else
        storageFullMask = fill_pct >= Opt.fill_threshold_pct;
    end

    satExcessMask = wetMask & storageFullMask;

    if hasZwt
        satExcessMask = satExcessMask | (wetMask & zwt_current_m <= Opt.zwt_threshold_m);
    end

    if ~all(isnan(f_mm_h(:)))
        f_mm_h = clean_map(f_mm_h, validMask);
        satExcessMask = satExcessMask | (wetMask & f_mm_h <= Opt.f_threshold_mm_h & fill_pct >= 95);
    end

    TS.mean_depth_m(it) = area_weighted_mean(d_mm / 1000, cellArea, validMask);
    TS.max_depth_m(it)  = local_nanmax(d_mm(validMask)) / 1000;
    TS.wet_area_pct(it) = 100 * nansum_local(cellArea(wetMask & validMask)) / areaTotal_m2;

    TS.mean_I_t_mm(it) = area_weighted_mean(I_t_mm, cellArea, validMask);
    TS.max_I_t_mm(it)  = local_nanmax(I_t_mm(validMask));

    TS.mean_fill_pct(it) = area_weighted_mean(fill_pct, cellArea, validMask);
    TS.p05_fill_pct(it) = percentile_valid(fill_pct(validMask), 5);
    TS.p50_fill_pct(it) = percentile_valid(fill_pct(validMask), 50);
    TS.p95_fill_pct(it) = percentile_valid(fill_pct(validMask), 95);

    TS.mean_Srem_mm(it) = area_weighted_mean(Srem_mm, cellArea, validMask);
    TS.p05_Srem_mm(it) = percentile_valid(Srem_mm(validMask), 5);
    TS.p50_Srem_mm(it) = percentile_valid(Srem_mm(validMask), 50);

    TS.storage_full_area_pct(it) = 100 * nansum_local(cellArea(storageFullMask & validMask)) / areaTotal_m2;
    TS.sat_excess_area_pct(it)   = 100 * nansum_local(cellArea(satExcessMask & validMask)) / areaTotal_m2;

    if hasZwt
        TS.mean_zwt_m(it) = area_weighted_mean(zwt_current_m, cellArea, validMask);
        TS.min_zwt_m(it)  = local_nanmin(zwt_current_m(validMask));
        TS.p05_zwt_m(it)  = percentile_valid(zwt_current_m(validMask), 5);
    end

    if ~all(isnan(f_mm_h(:)))
        TS.mean_f_mm_h(it) = area_weighted_mean(f_mm_h, cellArea, validMask);
    end
    if ~all(isnan(C_mm_h(:)))
        C_mm_h = clean_map(C_mm_h, validMask);
        TS.mean_C_mm_h(it) = area_weighted_mean(C_mm_h, cellArea, validMask);
    end
    if ~all(isnan(recharge_mm_h(:)))
        recharge_mm_h = clean_map(recharge_mm_h, validMask);
        TS.mean_recharge_mm_h(it) = area_weighted_mean(recharge_mm_h, cellArea, validMask);
    end

    TS.surface_storage_m3(it) = nansum_local((d_mm(validMask) / 1000) .* cellArea(validMask));
    TS.UZ_storage_m3(it) = nansum_local((I_t_mm(validMask) / 1000) .* cellArea(validMask));

    if hasZwt
        GWrel_m = max(Opt.zwt0_m - zwt_current_m, 0) .* Sy_local;
        TS.GW_storage_rel_m3(it) = nansum_local(GWrel_m(validMask) .* cellArea(validMask));
    end

    if hasZwt
        TS.total_storage_rel_m3(it) = TS.surface_storage_m3(it) + TS.UZ_storage_m3(it) + TS.GW_storage_rel_m3(it);
    else
        TS.total_storage_rel_m3(it) = TS.surface_storage_m3(it) + TS.UZ_storage_m3(it);
    end

    MaxDepth_m(validMask) = max(MaxDepth_m(validMask), d_mm(validMask) / 1000);
    MaxIt_mm(validMask) = max(MaxIt_mm(validMask), I_t_mm(validMask));
    MaxFill_pct(validMask) = max(MaxFill_pct(validMask), fill_pct(validMask));
    MinSrem_mm(validMask) = min(MinSrem_mm(validMask), Srem_mm(validMask));

    if hasZwt
        MinZwt_m(validMask) = min(MinZwt_m(validMask), zwt_current_m(validMask));
    end
    if ~all(isnan(recharge_mm_h(:)))
        MaxRecharge_mm_h(validMask) = max(MaxRecharge_mm_h(validMask), recharge_mm_h(validMask));
    end

    newlySat = satExcessMask & validMask & isnan(TimeToSat_h);
    TimeToSat_h(newlySat) = time_h(it);
    MaxSatMask = MaxSatMask | satExcessMask;

    if ismember(it, snapshotIdx)
        make_snapshot_figure(Dirs.Snaps, time_h(it), d_mm, I_t_mm, fill_pct, Srem_mm, ...
            zwt_current_m, satExcessMask, validMask, hasZwt, capacityMode, Opt);

        if Opt.write_geotiff
            write_snapshot_rasters(Dirs.Maps, W.DEM_raster, time_h(it), d_mm/1000, I_t_mm, ...
                fill_pct, Srem_mm, zwt_current_m, double(satExcessMask), validMask, hasZwt);
        end
    end
end

%% ------------------------------------------------------------------------
% 9) MASS RESIDUAL, SUMMARY, EXPORTS
% -------------------------------------------------------------------------
storage0 = TS.total_storage_rel_m3(1);
dS = TS.total_storage_rel_m3 - storage0;
TS.mass_residual_m3 = TS.cumulative_rain_m3 - TS.cumulative_outlet_m3 - dS;
TS.mass_residual_mm = TS.mass_residual_m3 / areaTotal_m2 * 1000;

Event = table();
Event.metric = strings(0,1); Event.value = zeros(0,1); Event.units = strings(0,1); Event.note = strings(0,1);
Event = add_event_metric(Event, 'Q_noinf', Q_noinf_m3s, 'm3/s', 'Rainfall-only no-infiltration steady discharge');
Event = add_event_metric(Event, 'final_Q', TS.outlet_Q_m3s(end), 'm3/s', 'Final outlet discharge');
Event = add_event_metric(Event, 'final_Q_ratio_to_noinf', TS.Q_ratio_to_noinf(end), '-', 'Final outlet discharge divided by no-infiltration Qss');
Event = add_event_metric(Event, 'max_sat_excess_area', local_nanmax(TS.sat_excess_area_pct), '%', 'Maximum area classified as saturation-excess');
Event = add_event_metric(Event, 'final_sat_excess_area', TS.sat_excess_area_pct(end), '%', 'Final saturation-excess area percentage');
Event = add_event_metric(Event, 'time_sat_area_50pct', first_time_threshold(TS.time_h, TS.sat_excess_area_pct, 50), 'h', 'First time saturation-excess area exceeds 50%');
Event = add_event_metric(Event, 'time_sat_area_90pct', first_time_threshold(TS.time_h, TS.sat_excess_area_pct, 90), 'h', 'First time saturation-excess area exceeds 90%');
Event = add_event_metric(Event, 'time_mean_f_below_threshold', first_time_below(TS.time_h, TS.mean_f_mm_h, Opt.f_threshold_mm_h), 'h', 'First time mean infiltration falls below threshold');
Event = add_event_metric(Event, 'time_Q_above_95pct_noinf', first_time_threshold(TS.time_h, TS.Q_ratio_to_noinf, 0.95), 'h', 'First time outlet Q exceeds 95% of Qnoinf');
Event = add_event_metric(Event, 'final_mean_I_t', TS.mean_I_t_mm(end), 'mm', 'Final mean above-residual vadose storage');
Event = add_event_metric(Event, 'final_mean_fill', TS.mean_fill_pct(end), '%', 'Final mean storage fill percentage');
Event = add_event_metric(Event, 'final_mass_residual', TS.mass_residual_mm(end), 'mm', 'Water-balance residual; depends on saved GW states');

writetable(TS, fullfile(Dirs.Tables, 'Case03_Advanced_TimeSeries.csv'));
writetable(Event, fullfile(Dirs.Tables, 'Case03_Advanced_EventSummary.csv'));
make_event_summary_text_file(Dirs.Tables, Event);

MaxDepth_m(~validMask) = nan; MaxDepth_m(isinf(MaxDepth_m)) = nan;
MaxIt_mm(~validMask) = nan; MaxIt_mm(isinf(MaxIt_mm)) = nan;
MaxFill_pct(~validMask) = nan; MaxFill_pct(isinf(MaxFill_pct)) = nan;
MinSrem_mm(~validMask) = nan; MinSrem_mm(isinf(MinSrem_mm)) = nan;
MinZwt_m(~validMask) = nan; MinZwt_m(isinf(MinZwt_m)) = nan;
MaxRecharge_mm_h(~validMask) = nan; MaxRecharge_mm_h(isinf(MaxRecharge_mm_h)) = nan;
TimeToSat_h(~validMask) = nan;

if Opt.write_geotiff
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_MaxDepth_m.tif'), W.DEM_raster, MaxDepth_m);
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_Max_I_t_mm.tif'), W.DEM_raster, MaxIt_mm);
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_MaxStorageFill_pct.tif'), W.DEM_raster, MaxFill_pct);
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_MinRemainingStorage_mm.tif'), W.DEM_raster, MinSrem_mm);
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_SatExcess_EverMask.tif'), W.DEM_raster, double(MaxSatMask));
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_TimeToSatExcess_h.tif'), W.DEM_raster, TimeToSat_h);
    if any(isfinite(MinZwt_m(:)))
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_MinZWT_m.tif'), W.DEM_raster, MinZwt_m);
    end
    if any(isfinite(MaxRecharge_mm_h(:)))
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case03_MaxRecharge_mm_h.tif'), W.DEM_raster, MaxRecharge_mm_h);
    end
end

make_timeseries_figure(Dirs.Figures, TS, Q_noinf_m3s, Opt);
make_storage_volume_figure(Dirs.Figures, TS, Opt);

Results = struct();
Results.Options = Opt;
Results.TimeSeries = TS;
Results.EventSummary = Event;
Results.StaticMaps.MaxDepth_m = MaxDepth_m;
Results.StaticMaps.MaxIt_mm = MaxIt_mm;
Results.StaticMaps.MaxFill_pct = MaxFill_pct;
Results.StaticMaps.MinSrem_mm = MinSrem_mm;
Results.StaticMaps.MinZwt_m = MinZwt_m;
Results.StaticMaps.TimeToSat_h = TimeToSat_h;
Results.StaticMaps.SatExcessEverMask = MaxSatMask;
Results.OutputDirs = Dirs;

save(fullfile(Dirs.Root, 'Case03_Advanced_Analysis_Results.mat'), 'Results', '-v7.3');

fprintf('\n============================================================\n');
fprintf('Advanced Case 3 analysis complete.\n');
fprintf('Tables:  %s\n', Dirs.Tables);
fprintf('Figures: %s\n', Dirs.Figures);
fprintf('Rasters: %s\n', Dirs.Maps);
fprintf('============================================================\n\n');

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================
function S = set_default(S, fieldName, value)
if ~isfield(S, fieldName) || isempty(S.(fieldName)); S.(fieldName) = value; end
end

function ensure_folder(p)
if ~isfolder(p); mkdir(p); end
end

function x = gather_if_needed(x)
try
    if isa(x, 'gpuArray'); x = gather(x); end
catch
end
end

function val = get_field_or(S, fieldName, defaultVal)
if isfield(S, fieldName); val = S.(fieldName); else; val = defaultVal; end
end

function A = to_grid(x, nRows, nCols)
x = double(gather_if_needed(x));
if isscalar(x)
    A = x + zeros(nRows, nCols);
elseif isequal(size(x), [nRows, nCols])
    A = x;
else
    A = mean(x(:), 'omitnan') + zeros(nRows, nCols);
end
end

function v = read_table_var(T, preferredName, defaultVal)
v = defaultVal;
if isempty(T); return; end
names = T.Properties.VariableNames;
idx = find(strcmpi(names, preferredName), 1, 'first');
if isempty(idx); idx = find(contains(lower(names), lower(preferredName)), 1, 'first'); end
if ~isempty(idx)
    tmp = T.(names{idx});
    if isnumeric(tmp)
        v = tmp(1);
    elseif iscell(tmp)
        v = str2double(tmp{1});
    elseif isstring(tmp)
        v = str2double(tmp(1));
    end
end
end

function rain_mm_h = infer_rainfall_intensity(W)
rain_mm_h = 20;
if isfield(W, 'Rainfall_Parameters') && isfield(W.Rainfall_Parameters, 'intensity_rainfall')
    rr = double(gather_if_needed(W.Rainfall_Parameters.intensity_rainfall(:)));
    rr = rr(isfinite(rr) & rr > 0);
    if ~isempty(rr); rain_mm_h = median(rr); end
end
end

function t_min = get_time_minutes(running_control, flags)
t_min = [];
if ~isfield(running_control, 'time_records'); return; end
t = gather_if_needed(running_control.time_records(:));
if isdatetime(t)
    t_min = minutes(t - t(1));
elseif isduration(t)
    t_min = minutes(t);
else
    t = double(t);
    if isfield(flags, 'flag_elapsed_time') && flags.flag_elapsed_time == 1
        t_min = t(:);
    else
        if max(t) > 10 && max(t) <= 1e6
            t_min = t(:);
        else
            t_min = (t(:) - t(1)) * 24 * 60;
        end
    end
end
end

function t_min = get_time_vector_generic(t)
t = gather_if_needed(t(:));
if isdatetime(t)
    t_min = minutes(t - t(1));
elseif isduration(t)
    t_min = minutes(t);
else
    t = double(t);
    if max(t) <= 30
        t_min = (t - t(1)) * 24 * 60;
    else
        t_min = t;
    end
end
end

function mapFile = find_map_file(tempDir, storeID)
candidates = {fullfile(tempDir, sprintf('save_map_hydro_%d.mat', storeID)), fullfile(tempDir, sprintf('save_map_hydro_%d', storeID))};
mapFile = '';
for i = 1:numel(candidates)
    if isfile(candidates{i}); mapFile = candidates{i}; return; end
end
error('Could not find map chunk %d in:\n%s', storeID, tempDir);
end

function nTotal = count_total_saved_frames(tempDir)
nTotal = 0; store = 1;
while true
    try; f = find_map_file(tempDir, store); catch; break; end
    Tmp = load(f, 'Maps'); Maps = Tmp.Maps;
    if isfield(Maps, 'Hydro') && isfield(Maps.Hydro, 'd')
        nTotal = nTotal + size(Maps.Hydro.d, 3);
    else
        break
    end
    store = store + 1;
end
end

function Z = get_map_by_candidates(Maps, candidates, localID)
Z = nan;
for i = 1:numel(candidates)
    [tf, val] = get_nested_field(Maps, candidates{i});
    if tf
        val = gather_if_needed(val);
        if ndims(val) == 3
            if localID <= size(val,3); Z = double(val(:,:,localID)); return; end
        elseif ismatrix(val)
            Z = double(val); return;
        end
    end
end
if isfield(Maps, 'Hydro') && isfield(Maps.Hydro, 'd')
    d = Maps.Hydro.d; Z = nan(size(d,1), size(d,2));
end
end

function [tf, val] = get_nested_field(S, pathStr)
tf = false; val = S; parts = strsplit(pathStr, '.');
for i = 1:numel(parts)
    if isstruct(val) && isfield(val, parts{i})
        val = val.(parts{i});
    else
        return
    end
end
tf = true;
end

function Z = clean_map(Z, validMask)
if isscalar(Z) && isnan(Z); Z = nan(size(validMask)); end
Z = double(gather_if_needed(Z));
Z(~validMask) = nan;
Z(~isfinite(Z)) = nan;
end

function m = area_weighted_mean(Z, A, validMask)
idx = validMask & isfinite(Z) & isfinite(A);
if ~any(idx(:)); m = nan; else; m = nansum_local(Z(idx) .* A(idx)) / nansum_local(A(idx)); end
end

function p = percentile_valid(x, q)
x = x(isfinite(x));
if isempty(x); p = nan; else; p = prctile(x, q); end
end

function t0 = first_time_threshold(t, y, threshold)
idx = find(isfinite(y) & y >= threshold, 1, 'first');
if isempty(idx); t0 = nan; else; t0 = t(idx); end
end

function t0 = first_time_below(t, y, threshold)
idx = find(isfinite(y) & y <= threshold, 1, 'first');
if isempty(idx); t0 = nan; else; t0 = t(idx); end
end

function T = add_event_metric(T, metric, value, units, note)
newRow = table(string(metric), value, string(units), string(note), 'VariableNames', {'metric','value','units','note'});
T = [T; newRow];
end

function s = nansum_local(x)
x = x(isfinite(x));
if isempty(x); s = 0; else; s = sum(x(:)); end
end

function m = local_nanmean(Z, validMask)
if nargin > 1; x = Z(validMask); else; x = Z(:); end
x = x(isfinite(x));
if isempty(x); m = nan; else; m = mean(x); end
end

function m = local_nanmax(x)
x = x(isfinite(x));
if isempty(x); m = nan; else; m = max(x); end
end

function m = local_nanmin(x)
x = x(isfinite(x));
if isempty(x); m = nan; else; m = min(x); end
end

function make_snapshot_figure(outDir, time_h, d_mm, I_t_mm, fill_pct, Srem_mm, zwt_m, satMask, validMask, hasZwt, capacityMode, Opt)
if Opt.show_figures; vis = 'on'; else; vis = 'off'; end
fig = figure('Color','w','Units','pixels','Position',[50 50 1700 950], 'Visible',vis);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
nexttile; plot_map(d_mm/1000, validMask, sprintf('Depth [m], t = %.2f h', time_h)); colorbar;
nexttile; plot_map(I_t_mm, validMask, 'I_t [mm]'); colorbar;
nexttile; plot_map(fill_pct, validMask, sprintf('Storage fill [%%]\n%s', capacityMode)); colorbar; caxis([0 110]);
nexttile; plot_map(Srem_mm, validMask, 'Remaining storage S_{rem} [mm]'); colorbar;
nexttile;
if hasZwt; plot_map(zwt_m, validMask, 'Water-table depth z_{wt} [m]'); else; plot_map(nan(size(d_mm)), validMask, 'z_{wt} not stored'); end
colorbar;
nexttile; plot_map(double(satMask), validMask, 'Saturation-excess mask'); colorbar; caxis([0 1]);
sgtitle(sprintf('Updated Case 3 diagnostics at t = %.2f h', time_h), 'Interpreter','none');
outPng = fullfile(outDir, sprintf('Case03_snapshot_t_%07.3f_h.png', time_h));
exportgraphics(fig, outPng, 'Resolution', 250);
close(fig);
end

function plot_map(Z, validMask, ttl)
Z = double(Z); Z(~validMask) = nan;
imagesc(Z); axis image tight; set(gca,'YDir','normal');
title(ttl, 'Interpreter','tex'); xlabel('Column'); ylabel('Row');
end

function write_snapshot_rasters(outDir, DEM_raster, time_h, depth_m, I_t_mm, fill_pct, Srem_mm, zwt_m, satMask, validMask, hasZwt)
tag = sprintf('t_%07.3f_h', time_h); tag = strrep(tag, '.', 'p');
depth_m(~validMask) = nan; I_t_mm(~validMask) = nan; fill_pct(~validMask) = nan; Srem_mm(~validMask) = nan; satMask(~validMask) = nan;
write_geotiff_safe(fullfile(outDir, ['Case03_Depth_' tag '.tif']), DEM_raster, depth_m);
write_geotiff_safe(fullfile(outDir, ['Case03_I_t_' tag '.tif']), DEM_raster, I_t_mm);
write_geotiff_safe(fullfile(outDir, ['Case03_FillPct_' tag '.tif']), DEM_raster, fill_pct);
write_geotiff_safe(fullfile(outDir, ['Case03_Srem_' tag '.tif']), DEM_raster, Srem_mm);
write_geotiff_safe(fullfile(outDir, ['Case03_SatMask_' tag '.tif']), DEM_raster, satMask);
if hasZwt
    zwt_m(~validMask) = nan;
    write_geotiff_safe(fullfile(outDir, ['Case03_ZWT_' tag '.tif']), DEM_raster, zwt_m);
end
end

function write_geotiff_safe(outPath, DEM_raster, Z)
try
    if isfield(DEM_raster, 'georef') && isfield(DEM_raster.georef, 'GeoKeyDirectoryTag')
        geotiffwrite(outPath, single(Z), DEM_raster.georef.SpatialRef, 'GeoKeyDirectoryTag', DEM_raster.georef.GeoKeyDirectoryTag);
    else
        geotiffwrite(outPath, single(Z), DEM_raster.georef.SpatialRef);
    end
catch ME
    warning('Could not write GeoTIFF:\n%s\nReason: %s', outPath, ME.message);
end
end

function make_timeseries_figure(outDir, TS, Q_noinf_m3s, Opt)
if Opt.show_figures; vis = 'on'; else; vis = 'off'; end
fig = figure('Color','w','Units','pixels','Position',[50 50 1700 1000], 'Visible',vis);
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
nexttile; plot(TS.time_h, TS.outlet_Q_m3s, 'LineWidth', 2); hold on; yline(Q_noinf_m3s, '--', 'Q_{no-inf}', 'LineWidth', 1.5); xlabel('Time [h]'); ylabel('Outlet Q [m^3/s]'); title('Outlet hydrograph'); grid on;
nexttile; plot(TS.time_h, TS.mean_I_t_mm, 'LineWidth', 2); hold on; plot(TS.time_h, TS.mean_Srem_mm, 'LineWidth', 2); xline(Opt.expected_fill_time_h, ':', 'Expected fill time', 'LineWidth', 1.5); xlabel('Time [h]'); ylabel('Storage [mm]'); legend('Mean I_t','Mean S_{rem}','Location','best'); title('Vadose storage'); grid on;
nexttile; plot(TS.time_h, TS.mean_fill_pct, 'LineWidth', 2); hold on; plot(TS.time_h, TS.storage_full_area_pct, 'LineWidth', 2); plot(TS.time_h, TS.sat_excess_area_pct, 'LineWidth', 2); xlabel('Time [h]'); ylabel('[%]'); legend('Mean fill','Storage-full area','Sat-excess area','Location','best'); title('Saturation-excess expansion'); grid on;
nexttile; plot(TS.time_h, TS.mean_depth_m, 'LineWidth', 2); hold on; plot(TS.time_h, TS.max_depth_m, 'LineWidth', 2); xlabel('Time [h]'); ylabel('Depth [m]'); legend('Mean depth','Max depth','Location','best'); title('Surface water depth'); grid on;
nexttile; plot(TS.time_h, TS.mean_f_mm_h, 'LineWidth', 2); hold on; plot(TS.time_h, TS.mean_C_mm_h, 'LineWidth', 2); plot(TS.time_h, TS.mean_recharge_mm_h, 'LineWidth', 2); xlabel('Time [h]'); ylabel('Rate [mm/h]'); legend('Mean f','Mean C','Mean recharge','Location','best'); title('Infiltration / capacity / recharge'); grid on;
nexttile; plot(TS.time_h, TS.mean_zwt_m, 'LineWidth', 2); hold on; plot(TS.time_h, TS.min_zwt_m, 'LineWidth', 2); yline(0.005,':','z_{wt}=5 mm'); xlabel('Time [h]'); ylabel('z_{wt} [m]'); legend('Mean z_{wt}','Minimum z_{wt}','Location','best'); title('Water-table depth'); grid on;
sgtitle('Updated Case 3: near-surface groundwater saturation-excess diagnostics');
exportgraphics(fig, fullfile(outDir, 'Case03_Advanced_TimeSeries.png'), 'Resolution', 250);
saveas(fig, fullfile(outDir, 'Case03_Advanced_TimeSeries.fig'));
close(fig);
end

function make_storage_volume_figure(outDir, TS, Opt)
if Opt.show_figures; vis = 'on'; else; vis = 'off'; end
fig = figure('Color','w','Units','pixels','Position',[80 80 1600 850], 'Visible',vis);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile; plot(TS.time_h, TS.cumulative_rain_m3, 'LineWidth', 2); hold on; plot(TS.time_h, TS.cumulative_outlet_m3, 'LineWidth', 2); xlabel('Time [h]'); ylabel('Cumulative volume [m^3]'); legend('Rainfall','Outlet','Location','best'); title('Cumulative rainfall and outlet volume'); grid on;
nexttile; plot(TS.time_h, TS.surface_storage_m3, 'LineWidth', 2); hold on; plot(TS.time_h, TS.UZ_storage_m3, 'LineWidth', 2); plot(TS.time_h, TS.GW_storage_rel_m3, 'LineWidth', 2); xlabel('Time [h]'); ylabel('Storage [m^3]'); legend('Surface','UZ','Relative GW','Location','best'); title('Storage components'); grid on;
nexttile; plot(TS.time_h, TS.mass_residual_mm, 'LineWidth', 2); xlabel('Time [h]'); ylabel('Residual [mm]'); title('Water-balance residual'); grid on;
nexttile; plot(TS.time_h, TS.Q_ratio_to_noinf, 'LineWidth', 2); hold on; yline(1, '--', 'Q/Q_{no-inf}=1'); xlabel('Time [h]'); ylabel('Q / Q_{no-inf}'); title('Approach to no-infiltration benchmark'); grid on;
sgtitle('Updated Case 3 volume and benchmark diagnostics');
exportgraphics(fig, fullfile(outDir, 'Case03_Advanced_VolumeDiagnostics.png'), 'Resolution', 250);
saveas(fig, fullfile(outDir, 'Case03_Advanced_VolumeDiagnostics.fig'));
close(fig);
end

function make_event_summary_text_file(outDir, Event)
fid = fopen(fullfile(outDir, 'Case03_Advanced_EventSummary.txt'), 'w');
fprintf(fid, 'Updated Case 3 Advanced Saturation-Excess Event Summary\n');
fprintf(fid, '============================================================\n\n');
for i = 1:height(Event)
    fprintf(fid, '%-35s : %12.5g %-8s | %s\n', Event.metric(i), Event.value(i), Event.units(i), Event.note(i));
end
fclose(fid);
end
