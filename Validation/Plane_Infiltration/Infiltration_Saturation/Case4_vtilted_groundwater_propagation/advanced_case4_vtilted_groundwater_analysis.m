function Results = advanced_case4_vtilted_groundwater_analysis(Opt)
%% ========================================================================
% ADVANCED CASE 4 — V-TILTED 2D GROUNDWATER PROPAGATION ANALYSIS
% ========================================================================
%
% Purpose
% -------------------------------------------------------------------------
% Advanced post-processing for:
%
%   Case4_vtilted_groundwater_propagation
%
% This benchmark is designed to test the 2D Boussinesq groundwater module
% without rainfall:
%
%   - V-shaped tilted catchment
%   - soil/aquifer depth ≈ 1 m
%   - initial water table ≈ 0.25 m below land surface
%   - rainfall off
%   - infiltration off
%   - baseflow / 2D groundwater propagation on
%
% Expected behavior
% -------------------------------------------------------------------------
% With no rainfall, any surface water should come from groundwater
% redistribution and exfiltration. The main question is whether groundwater
% heads redistribute toward the central valley / downstream direction and
% whether exfiltration generates surface water in physically meaningful
% locations.
%
% Minimum stored maps
% -------------------------------------------------------------------------
%   Maps.Hydro.d
%
% Strongly recommended stored maps
% -------------------------------------------------------------------------
%   Maps.Hydro.h_t              groundwater head elevation [m]
%   Maps.Hydro.zwt              depth to water table below surface [m]
%   Maps.Hydro.GW_depth         saturated thickness [m]
%   Maps.Hydro.q_exf            exfiltration flux [m/s]
%   Maps.Hydro.q_river          river exchange flux [m/s], optional
%
% If h_t / zwt / GW_depth are not saved, the script can still analyze
% surface depth and outlet response, but it cannot prove the groundwater
% propagation mechanism.
%
% Usage
% -------------------------------------------------------------------------
% Run from the Case 4 root folder:
%
%   Results = advanced_case4_vtilted_groundwater_analysis;
%
% Or pass options:
%
%   Opt = struct();
%   Opt.workspaceFile = '/path/to/modeled_results.mat';
%   Opt.caseRoot      = '/path/to/Case4_vtilted_groundwater_propagation';
%   Opt.sim_duration_h = 72;
%   Results = advanced_case4_vtilted_groundwater_analysis(Opt);
%
% ========================================================================

%% ------------------------------------------------------------------------
% 0) USER OPTIONS
% -------------------------------------------------------------------------
if nargin < 1 || isempty(Opt)
    Opt = struct();
end

Opt = set_default(Opt, 'workspaceFile', fullfile(pwd, 'modeled_results.mat'));
Opt = set_default(Opt, 'caseRoot', pwd);
Opt = set_default(Opt, 'outputRoot', fullfile(pwd, 'Advanced_Case04_GW_Analysis'));

% Case defaults
Opt = set_default(Opt, 'caseFolderName', 'Case4_vtilted_groundwater_propagation');
Opt = set_default(Opt, 'sim_duration_h', 72);
Opt = set_default(Opt, 'rain_mm_h', 0);

% Geometry / initial condition defaults
Opt = set_default(Opt, 'DTB_m', []);
Opt = set_default(Opt, 'zwt0_m', []);
Opt = set_default(Opt, 'initial_head_mode', 'DEM_minus_zwt0');

% Diagnostic thresholds
Opt = set_default(Opt, 'depth_threshold_mm', 0.5);           % wet surface threshold
Opt = set_default(Opt, 'exfil_threshold_mm_h', 0.05);        % exfiltration threshold
Opt = set_default(Opt, 'surface_sat_zwt_threshold_m', 0.005);% water table within 5 mm
Opt = set_default(Opt, 'exclude_perimeter', true);
Opt = set_default(Opt, 'perimeter_buffer_cells', 1);

% V-catchment diagnostic masks
Opt = set_default(Opt, 'center_band_fraction', 0.15);        % middle 15% of width
Opt = set_default(Opt, 'hillslope_band_fraction', 0.25);     % outer 25% each side

% Snapshot times
Opt = set_default(Opt, 'snapshot_hours', [0 0.5 1 2 4 8 12 24 36 48 72]);

% Exports
Opt = set_default(Opt, 'write_geotiff', true);
Opt = set_default(Opt, 'show_figures', true);

fprintf('\n============================================================\n');
fprintf('Advanced Case 4 V-Tilted Groundwater Analysis\n');
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
    'Paths', ...
    'flags', ...
    'running_control', ...
    'saver_memory_maps', ...
    'idx_nan', ...
    'DEM_raster', ...
    'Elevation_Properties', ...
    'Soil_Properties', ...
    'Wshed_Properties', ...
    'outlet_states', ...
    'BC_States', ...
    'depths', ...
    'max_GW_depth'};

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
Dirs.Root     = Opt.outputRoot;
Dirs.Figures  = fullfile(Dirs.Root, 'Figures');
Dirs.Tables   = fullfile(Dirs.Root, 'Tables');
Dirs.Maps     = fullfile(Dirs.Root, 'Diagnostic_Rasters');
Dirs.Snaps    = fullfile(Dirs.Root, 'Snapshot_Figures');
Dirs.Profiles = fullfile(Dirs.Root, 'Profile_Figures');

ensure_folder(Dirs.Root);
ensure_folder(Dirs.Figures);
ensure_folder(Dirs.Tables);
ensure_folder(Dirs.Maps);
ensure_folder(Dirs.Snaps);
ensure_folder(Dirs.Profiles);

%% ------------------------------------------------------------------------
% 3) READ CASE CONFIG IF AVAILABLE
% -------------------------------------------------------------------------
caseConfigPath = fullfile(Opt.caseRoot, 'Case_Config.csv');

if isfile(caseConfigPath)
    CaseConfig = readtable(caseConfigPath);

    if isempty(Opt.DTB_m)
        Opt.DTB_m = read_table_var(CaseConfig, 'DTB_m', []);
    end

    if isempty(Opt.zwt0_m)
        Opt.zwt0_m = read_table_var(CaseConfig, 'zwt0_m', []);
    end

    if isempty(Opt.sim_duration_h) || isnan(Opt.sim_duration_h)
        sim_min = read_table_var(CaseConfig, 'sim_duration_min', []);
        if ~isempty(sim_min)
            Opt.sim_duration_h = sim_min / 60;
        end
    end
else
    CaseConfig = table();
end

if isempty(Opt.DTB_m)
    Opt.DTB_m = infer_scalar_soil_depth(W);
end

if isempty(Opt.zwt0_m)
    Opt.zwt0_m = 0.25;
end

fprintf('Rainfall intensity used in diagnostics = %.4f mm/h\n', Opt.rain_mm_h);
fprintf('Soil/aquifer depth DTB                 = %.4f m\n', Opt.DTB_m);
fprintf('Initial water-table depth zwt0         = %.4f m\n', Opt.zwt0_m);
fprintf('Simulation duration                    = %.4f h\n\n', Opt.sim_duration_h);

%% ------------------------------------------------------------------------
% 4) STATIC GEOMETRY AND MASKS
% -------------------------------------------------------------------------
DEM = double(gather_if_needed(W.Elevation_Properties.elevation_cell));
[nRows, nCols] = size(DEM);

if isfield(W, 'idx_nan')
    idx_nan = logical(gather_if_needed(W.idx_nan));
else
    idx_nan = isnan(DEM) | isinf(DEM);
end

validMaskFull = ~idx_nan & isfinite(DEM);

% Optional perimeter exclusion for diagnostics only.
diagMask = validMaskFull;

idx_perimeter_diag = false(size(diagMask));
if Opt.exclude_perimeter
    if isfield(W, 'Wshed_Properties') && ...
            isfield(W.Wshed_Properties, 'perimeter') && ...
            ~isempty(W.Wshed_Properties.perimeter)

        idx_perimeter_diag = logical(gather_if_needed(W.Wshed_Properties.perimeter));
    else
        b = max(1, round(Opt.perimeter_buffer_cells));
        idx_perimeter_diag(1:b,:) = true;
        idx_perimeter_diag(end-b+1:end,:) = true;
        idx_perimeter_diag(:,1:b) = true;
        idx_perimeter_diag(:,end-b+1:end) = true;
    end

    diagMask = diagMask & ~idx_perimeter_diag;
end

SoilDepth = to_grid(W.Soil_Properties.Soil_Depth, nRows, nCols);
theta_s   = to_grid(W.Soil_Properties.theta_sat,  nRows, nCols);
theta_r   = to_grid(W.Soil_Properties.theta_r,    nRows, nCols);

if isfield(W.Soil_Properties, 'Sy')
    Sy = to_grid(W.Soil_Properties.Sy, nRows, nCols);
else
    Sy = max(theta_s - theta_r, 1e-6);
end

if isfield(W.Soil_Properties, 'ksat_gw')
    Kgw_mm_h = to_grid(W.Soil_Properties.ksat_gw, nRows, nCols);
elseif isfield(W.Soil_Properties, 'ksat_gw_mm_h')
    Kgw_mm_h = to_grid(W.Soil_Properties.ksat_gw_mm_h, nRows, nCols);
else
    Kgw_mm_h = nan(nRows, nCols);
end

cellArea = infer_cell_area_grid(W, nRows, nCols);

areaTotal_m2 = nansum(cellArea(diagMask));
z_bed = DEM - SoilDepth;

% Initial groundwater head inferred from DEM - zwt0
h_initial = DEM - Opt.zwt0_m;
sat_thick_initial = max(h_initial - z_bed, 0);
sat_thick_initial = min(sat_thick_initial, SoilDepth);
GW_storage_initial_m3 = nansum(sat_thick_initial(diagMask) .* Sy(diagMask) .* cellArea(diagMask));

fprintf('Diagnostic area excluding perimeter = %.3f m2\n', areaTotal_m2);
fprintf('Mean Sy over diagnostic domain       = %.5f\n', area_weighted_mean(Sy, cellArea, diagMask));
fprintf('Mean Kgw over diagnostic domain      = %.5f mm/h\n', area_weighted_mean(Kgw_mm_h, cellArea, diagMask));
fprintf('Initial GW storage estimate          = %.5f m3\n\n', GW_storage_initial_m3);

%% ------------------------------------------------------------------------
% 5) V-CATCHMENT MASKS: CENTER VALLEY VS HILLSLOPES
% -------------------------------------------------------------------------
[colGrid, rowGrid] = meshgrid(1:nCols, 1:nRows);
midCol = (nCols + 1) / 2;

centerHalfWidth = max(1, round(0.5 * Opt.center_band_fraction * nCols));
sideWidth = max(1, round(Opt.hillslope_band_fraction * nCols));

centerMask = abs(colGrid - midCol) <= centerHalfWidth;
leftMask   = colGrid <= sideWidth;
rightMask  = colGrid >= (nCols - sideWidth + 1);
sideMask   = leftMask | rightMask;

centerMask = centerMask & diagMask;
sideMask   = sideMask & diagMask;

% Longitudinal sections
upstreamMask   = rowGrid <= round(0.25*nRows) & diagMask;
midstreamMask  = rowGrid > round(0.45*nRows) & rowGrid < round(0.55*nRows) & diagMask;
downstreamMask = rowGrid >= round(0.75*nRows) & diagMask;

%% ------------------------------------------------------------------------
% 6) TIME VECTOR AND MAP CHUNK INFO
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

nTotalSaved = count_total_saved_frames(tempDir);

if isempty(time_min)
    time_min = (0:nTotalSaved-1)';
else
    time_min = time_min(:);
    nUse = min(numel(time_min), nTotalSaved);
    time_min = time_min(1:nUse);
end

time_h = time_min / 60;

idxUse = time_h <= Opt.sim_duration_h + 1e-9;
time_h = time_h(idxUse);
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
fprintf('Saved map frames available: %d\n', nTotalSaved);
fprintf('Frames used in analysis:   %d\n', nTimes);
fprintf('Maps per chunk:            %d\n', nPerChunk);
fprintf('Snapshot hours:            %s\n\n', mat2str(round(time_h(snapshotIdx)',3)));

%% ------------------------------------------------------------------------
% 7) PREALLOCATE TIME-SERIES DIAGNOSTICS
% -------------------------------------------------------------------------
nanv = nan(nTimes,1);

TS = table();
TS.time_min = time_min(:);
TS.time_h   = time_h(:);

% Surface routing response
TS.outlet_Q_m3s = nanv;
TS.cumulative_outlet_m3 = nanv;

TS.mean_depth_m = nanv;
TS.max_depth_m  = nanv;
TS.wet_area_pct = nanv;
TS.surface_storage_m3 = nanv;

% Groundwater state
TS.mean_h_t_m = nanv;
TS.mean_zwt_m = nanv;
TS.min_zwt_m  = nanv;
TS.p05_zwt_m  = nanv;

TS.mean_saturated_thickness_m = nanv;
TS.max_saturated_thickness_m  = nanv;
TS.surface_saturated_area_pct = nanv;

TS.mean_gw_storage_m3 = nanv;
TS.gw_storage_change_m3 = nanv;

% Exfiltration
TS.mean_q_exf_mm_h = nanv;
TS.max_q_exf_mm_h  = nanv;
TS.exfil_area_pct  = nanv;
TS.exfil_rate_m3s  = nanv;
TS.cumulative_exfiltration_m3 = nanv;

% River exchange if saved
TS.mean_q_river_mm_h = nanv;
TS.river_exchange_rate_m3s = nanv;
TS.cumulative_river_exchange_m3 = nanv;

% V-catchment contrasts
TS.center_mean_zwt_m = nanv;
TS.side_mean_zwt_m   = nanv;
TS.center_mean_depth_m = nanv;
TS.side_mean_depth_m   = nanv;
TS.center_mean_q_exf_mm_h = nanv;
TS.side_mean_q_exf_mm_h   = nanv;

TS.upstream_mean_zwt_m = nanv;
TS.midstream_mean_zwt_m = nanv;
TS.downstream_mean_zwt_m = nanv;

% Balance diagnostics
TS.total_storage_m3 = nanv;
TS.total_storage_change_m3 = nanv;
TS.mass_residual_m3 = nanv;
TS.mass_residual_mm = nanv;

%% ------------------------------------------------------------------------
% 8) OUTLET HYDROGRAPH
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

    t_sec = time_min(:) * 60;
    TS.cumulative_outlet_m3 = cumtrapz(t_sec, q_out);
else
    warning('outlet_states.outlet_hydrograph not found. Outlet volume diagnostics will be NaN.');
end

%% ------------------------------------------------------------------------
% 9) PREALLOCATE STATIC DIAGNOSTIC MAPS
% -------------------------------------------------------------------------
MaxDepth_m       = -inf(nRows, nCols);
MaxQexf_mm_h     = -inf(nRows, nCols);
MinZwt_m         = inf(nRows, nCols);
MaxSatThick_m    = -inf(nRows, nCols);
FinalZwt_m       = nan(nRows, nCols);
FinalHt_m        = nan(nRows, nCols);
FinalSatThick_m  = nan(nRows, nCols);

TimeToWet_h      = nan(nRows, nCols);
TimeToExfil_h    = nan(nRows, nCols);
TimeToSurfaceSat_h = nan(nRows, nCols);

EverWetMask      = false(nRows, nCols);
EverExfilMask    = false(nRows, nCols);
EverSurfaceSatMask = false(nRows, nCols);

% Store profiles for snapshots
ProfileStore = struct();
ProfileStore.times_h = time_h(snapshotIdx);
ProfileStore.zwt_centerline = {};
ProfileStore.depth_centerline = {};
ProfileStore.qexf_centerline = {};
ProfileStore.cross_zwt_mid = {};
ProfileStore.cross_depth_mid = {};
ProfileStore.cross_qexf_mid = {};

%% ------------------------------------------------------------------------
% 10) MAIN LOOP OVER SAVED MAPS
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

    % ---------------------------------------------------------------------
    % Read available maps
    % ---------------------------------------------------------------------
    d_mm = get_map_by_candidates(Maps, ...
        {'Hydro.d','Hydro.depth','Hydro.depths'}, localID);

    h_t_m = get_map_by_candidates(Maps, ...
        {'Hydro.h_t','Hydro.GW_head','Hydro.water_table_elevation','Hydro.groundwater_head'}, localID);

    zwt_m = get_map_by_candidates(Maps, ...
        {'Hydro.zwt','Hydro.zwt_m','Hydro.water_table_depth','Hydro.depth_to_water_table'}, localID);

    GW_depth_m = get_map_by_candidates(Maps, ...
        {'Hydro.GW_depth','Hydro.gw_depth','Hydro.saturated_thickness','Hydro.sat_thick'}, localID);

    q_exf_raw = get_map_by_candidates(Maps, ...
        {'Hydro.q_exf','Hydro.exfiltration','Hydro.exfiltration_rate','Hydro.q_exfiltration'}, localID);

    q_river_raw = get_map_by_candidates(Maps, ...
        {'Hydro.q_river','Hydro.river_exchange','Hydro.qriv'}, localID);

    d_mm = clean_map(d_mm, diagMask);

    if all(isnan(d_mm(:)))
        error('Could not find Maps.Hydro.d in stored maps.');
    end

    % ---------------------------------------------------------------------
    % Build groundwater state from whatever was saved
    % ---------------------------------------------------------------------
    hasHead = ~all(isnan(h_t_m(:)));
    hasZwt  = ~all(isnan(zwt_m(:)));
    hasGWDepth = ~all(isnan(GW_depth_m(:)));

    if hasHead
        h_t_m = clean_map(h_t_m, diagMask);
        zwt_current_m = DEM - h_t_m;
        zwt_current_m = max(zwt_current_m, 0);
        zwt_current_m = min(zwt_current_m, SoilDepth);

        sat_thick_m = h_t_m - z_bed;
        sat_thick_m = max(sat_thick_m, 0);
        sat_thick_m = min(sat_thick_m, SoilDepth);

    elseif hasZwt
        zwt_current_m = clean_map(zwt_m, diagMask);
        zwt_current_m = max(zwt_current_m, 0);
        zwt_current_m = min(zwt_current_m, SoilDepth);

        h_t_m = DEM - zwt_current_m;
        sat_thick_m = SoilDepth - zwt_current_m;
        sat_thick_m = max(sat_thick_m, 0);
        sat_thick_m = min(sat_thick_m, SoilDepth);
        hasHead = true;

    elseif hasGWDepth
        sat_thick_m = clean_map(GW_depth_m, diagMask);
        sat_thick_m = max(sat_thick_m, 0);
        sat_thick_m = min(sat_thick_m, SoilDepth);

        zwt_current_m = SoilDepth - sat_thick_m;
        zwt_current_m = max(zwt_current_m, 0);
        zwt_current_m = min(zwt_current_m, SoilDepth);

        h_t_m = z_bed + sat_thick_m;
        hasHead = true;

    else
        zwt_current_m = nan(nRows, nCols);
        sat_thick_m = nan(nRows, nCols);
        h_t_m = nan(nRows, nCols);
    end

    hasGWState = hasHead || hasZwt || hasGWDepth;

    % ---------------------------------------------------------------------
    % Exfiltration fields
    % ---------------------------------------------------------------------
    hasQexf = ~all(isnan(q_exf_raw(:)));
    if hasQexf
        [q_exf_mm_h, q_exf_mps] = convert_flux_to_mm_h_and_mps(q_exf_raw, diagMask);
    else
        q_exf_mm_h = nan(nRows, nCols);
        q_exf_mps  = nan(nRows, nCols);
    end

    hasQriver = ~all(isnan(q_river_raw(:)));
    if hasQriver
        [q_river_mm_h, q_river_mps] = convert_flux_to_mm_h_and_mps(q_river_raw, diagMask);
    else
        q_river_mm_h = nan(nRows, nCols);
        q_river_mps  = nan(nRows, nCols);
    end

    % Masks
    wetMask = d_mm > Opt.depth_threshold_mm;
    if hasGWState
        surfaceSatMask = zwt_current_m <= Opt.surface_sat_zwt_threshold_m;
    else
        surfaceSatMask = false(nRows, nCols);
    end

    if hasQexf
        exfilMask = q_exf_mm_h > Opt.exfil_threshold_mm_h;
    else
        exfilMask = wetMask & surfaceSatMask;
    end

    wetMask = wetMask & diagMask;
    surfaceSatMask = surfaceSatMask & diagMask;
    exfilMask = exfilMask & diagMask;

    % ---------------------------------------------------------------------
    % Time-series statistics
    % ---------------------------------------------------------------------
    TS.mean_depth_m(it) = area_weighted_mean(d_mm / 1000, cellArea, diagMask);
    TS.max_depth_m(it)  = max_valid(d_mm(diagMask)) / 1000;
    TS.wet_area_pct(it) = 100 * nansum(cellArea(wetMask)) / areaTotal_m2;
    TS.surface_storage_m3(it) = nansum((d_mm(diagMask) / 1000) .* cellArea(diagMask));

    if hasGWState
        TS.mean_h_t_m(it) = area_weighted_mean(h_t_m, cellArea, diagMask);
        TS.mean_zwt_m(it) = area_weighted_mean(zwt_current_m, cellArea, diagMask);
        TS.min_zwt_m(it)  = min_valid(zwt_current_m(diagMask));
        TS.p05_zwt_m(it)  = percentile_valid(zwt_current_m(diagMask), 5);

        TS.mean_saturated_thickness_m(it) = area_weighted_mean(sat_thick_m, cellArea, diagMask);
        TS.max_saturated_thickness_m(it)  = max_valid(sat_thick_m(diagMask));
        TS.surface_saturated_area_pct(it) = 100 * nansum(cellArea(surfaceSatMask)) / areaTotal_m2;

        GW_storage_m3 = nansum(sat_thick_m(diagMask) .* Sy(diagMask) .* cellArea(diagMask));
        TS.mean_gw_storage_m3(it) = GW_storage_m3;
        TS.gw_storage_change_m3(it) = GW_storage_m3 - GW_storage_initial_m3;

        TS.center_mean_zwt_m(it) = area_weighted_mean(zwt_current_m, cellArea, centerMask);
        TS.side_mean_zwt_m(it)   = area_weighted_mean(zwt_current_m, cellArea, sideMask);

        TS.upstream_mean_zwt_m(it)   = area_weighted_mean(zwt_current_m, cellArea, upstreamMask);
        TS.midstream_mean_zwt_m(it)  = area_weighted_mean(zwt_current_m, cellArea, midstreamMask);
        TS.downstream_mean_zwt_m(it) = area_weighted_mean(zwt_current_m, cellArea, downstreamMask);
    end

    if hasQexf
        qpos_mps = max(q_exf_mps, 0);
        TS.mean_q_exf_mm_h(it) = area_weighted_mean(q_exf_mm_h, cellArea, diagMask);
        TS.max_q_exf_mm_h(it)  = max_valid(q_exf_mm_h(diagMask));
        TS.exfil_area_pct(it)  = 100 * nansum(cellArea(exfilMask)) / areaTotal_m2;
        TS.exfil_rate_m3s(it)  = nansum(qpos_mps(diagMask) .* cellArea(diagMask));

        TS.center_mean_q_exf_mm_h(it) = area_weighted_mean(q_exf_mm_h, cellArea, centerMask);
        TS.side_mean_q_exf_mm_h(it)   = area_weighted_mean(q_exf_mm_h, cellArea, sideMask);
    end

    if hasQriver
        TS.mean_q_river_mm_h(it) = area_weighted_mean(q_river_mm_h, cellArea, diagMask);
        TS.river_exchange_rate_m3s(it) = nansum(q_river_mps(diagMask) .* cellArea(diagMask));
    end

    TS.center_mean_depth_m(it) = area_weighted_mean(d_mm / 1000, cellArea, centerMask);
    TS.side_mean_depth_m(it)   = area_weighted_mean(d_mm / 1000, cellArea, sideMask);

    if hasGWState
        TS.total_storage_m3(it) = TS.surface_storage_m3(it) + TS.mean_gw_storage_m3(it);
    else
        TS.total_storage_m3(it) = TS.surface_storage_m3(it);
    end

    % ---------------------------------------------------------------------
    % Static maps
    % ---------------------------------------------------------------------
    MaxDepth_m(diagMask) = max(MaxDepth_m(diagMask), d_mm(diagMask) / 1000);

    if hasQexf
        MaxQexf_mm_h(diagMask) = max(MaxQexf_mm_h(diagMask), q_exf_mm_h(diagMask));
    end

    if hasGWState
        MinZwt_m(diagMask) = min(MinZwt_m(diagMask), zwt_current_m(diagMask));
        MaxSatThick_m(diagMask) = max(MaxSatThick_m(diagMask), sat_thick_m(diagMask));

        if it == nTimes
            FinalZwt_m = zwt_current_m;
            FinalHt_m = h_t_m;
            FinalSatThick_m = sat_thick_m;
        end
    end

    newlyWet = wetMask & isnan(TimeToWet_h);
    TimeToWet_h(newlyWet) = time_h(it);

    newlyExfil = exfilMask & isnan(TimeToExfil_h);
    TimeToExfil_h(newlyExfil) = time_h(it);

    newlySurfSat = surfaceSatMask & isnan(TimeToSurfaceSat_h);
    TimeToSurfaceSat_h(newlySurfSat) = time_h(it);

    EverWetMask = EverWetMask | wetMask;
    EverExfilMask = EverExfilMask | exfilMask;
    EverSurfaceSatMask = EverSurfaceSatMask | surfaceSatMask;

    % ---------------------------------------------------------------------
    % Snapshot figures / rasters / profiles
    % ---------------------------------------------------------------------
    if ismember(it, snapshotIdx)
        make_snapshot_figure( ...
            Dirs.Snaps, Opt.show_figures, time_h(it), ...
            DEM, d_mm, h_t_m, zwt_current_m, sat_thick_m, q_exf_mm_h, ...
            wetMask, exfilMask, surfaceSatMask, diagMask, hasGWState, hasQexf);

        if Opt.write_geotiff
            write_snapshot_rasters( ...
                Dirs.Maps, W.DEM_raster, time_h(it), ...
                d_mm/1000, h_t_m, zwt_current_m, sat_thick_m, q_exf_mm_h, ...
                double(wetMask), double(exfilMask), double(surfaceSatMask), ...
                diagMask, hasGWState, hasQexf);
        end

        % Profiles
        midRow = round(nRows/2);
        centerCol = round(midCol);

        ProfileStore.zwt_centerline{end+1} = zwt_current_m(:, centerCol); %#ok<AGROW>
        ProfileStore.depth_centerline{end+1} = d_mm(:, centerCol) / 1000; %#ok<AGROW>
        ProfileStore.qexf_centerline{end+1} = q_exf_mm_h(:, centerCol); %#ok<AGROW>

        ProfileStore.cross_zwt_mid{end+1} = zwt_current_m(midRow, :); %#ok<AGROW>
        ProfileStore.cross_depth_mid{end+1} = d_mm(midRow, :) / 1000; %#ok<AGROW>
        ProfileStore.cross_qexf_mid{end+1} = q_exf_mm_h(midRow, :); %#ok<AGROW>
    end
end

%% ------------------------------------------------------------------------
% 11) CUMULATIVE EXFILTRATION / RIVER EXCHANGE AND MASS BALANCE
% -------------------------------------------------------------------------
t_sec = TS.time_min(:) * 60;

if any(isfinite(TS.exfil_rate_m3s))
    TS.cumulative_exfiltration_m3 = cumtrapz(t_sec, fillmissing_linear(TS.exfil_rate_m3s));
end

if any(isfinite(TS.river_exchange_rate_m3s))
    TS.cumulative_river_exchange_m3 = cumtrapz(t_sec, fillmissing_linear(TS.river_exchange_rate_m3s));
end

storage0 = TS.total_storage_m3(1);
TS.total_storage_change_m3 = TS.total_storage_m3 - storage0;

% No rainfall in this benchmark.
% If no external groundwater boundary exists, approximately:
%   cumulative_outlet ≈ - total_storage_change
% so residual = -dS - Vout.
TS.mass_residual_m3 = -TS.total_storage_change_m3 - TS.cumulative_outlet_m3;
TS.mass_residual_mm = TS.mass_residual_m3 / areaTotal_m2 * 1000;

%% ------------------------------------------------------------------------
% 12) EVENT SUMMARY
% -------------------------------------------------------------------------
Event = table();
Event.metric = strings(0,1);
Event.value  = zeros(0,1);
Event.units  = strings(0,1);
Event.note   = strings(0,1);

Event = add_event_metric(Event, 'max_outlet_Q', max_valid(TS.outlet_Q_m3s), 'm3/s', ...
    'Maximum outlet discharge generated without rainfall');

Event = add_event_metric(Event, 'final_outlet_Q', TS.outlet_Q_m3s(end), 'm3/s', ...
    'Final outlet discharge');

Event = add_event_metric(Event, 'max_wet_area', max_valid(TS.wet_area_pct), '%', ...
    'Maximum surface wet area');

Event = add_event_metric(Event, 'max_surface_saturated_area', max_valid(TS.surface_saturated_area_pct), '%', ...
    'Maximum area with water table near land surface');

Event = add_event_metric(Event, 'max_exfil_area', max_valid(TS.exfil_area_pct), '%', ...
    'Maximum area with positive exfiltration above threshold');

Event = add_event_metric(Event, 'time_first_wet_10pct', first_time_threshold(TS.time_h, TS.wet_area_pct, 10), 'h', ...
    'First time wet area exceeds 10%');

Event = add_event_metric(Event, 'time_first_exfil_10pct', first_time_threshold(TS.time_h, TS.exfil_area_pct, 10), 'h', ...
    'First time exfiltration area exceeds 10%');

Event = add_event_metric(Event, 'time_first_surface_sat_10pct', first_time_threshold(TS.time_h, TS.surface_saturated_area_pct, 10), 'h', ...
    'First time surface-saturated area exceeds 10%');

Event = add_event_metric(Event, 'final_mean_zwt', TS.mean_zwt_m(end), 'm', ...
    'Final mean depth to water table');

Event = add_event_metric(Event, 'final_min_zwt', TS.min_zwt_m(end), 'm', ...
    'Final minimum depth to water table');

Event = add_event_metric(Event, 'final_GW_storage_change', TS.gw_storage_change_m3(end), 'm3', ...
    'Final groundwater storage change relative to initial');

Event = add_event_metric(Event, 'final_surface_storage', TS.surface_storage_m3(end), 'm3', ...
    'Final surface water storage');

Event = add_event_metric(Event, 'final_cumulative_outlet', TS.cumulative_outlet_m3(end), 'm3', ...
    'Cumulative outlet volume');

Event = add_event_metric(Event, 'final_mass_residual', TS.mass_residual_mm(end), 'mm', ...
    'No-rainfall residual: -dStorage - outlet volume');

%% ------------------------------------------------------------------------
% 13) EXPORT TABLES
% -------------------------------------------------------------------------
writetable(TS, fullfile(Dirs.Tables, 'Case04_GW_TimeSeries.csv'));
writetable(Event, fullfile(Dirs.Tables, 'Case04_GW_EventSummary.csv'));

%% ------------------------------------------------------------------------
% 14) CLEAN STATIC MAPS AND EXPORT RASTERS
% -------------------------------------------------------------------------
MaxDepth_m(~diagMask) = nan;
MaxQexf_mm_h(~diagMask) = nan;
MinZwt_m(~diagMask) = nan;
MaxSatThick_m(~diagMask) = nan;
FinalZwt_m(~diagMask) = nan;
FinalHt_m(~diagMask) = nan;
FinalSatThick_m(~diagMask) = nan;
TimeToWet_h(~diagMask) = nan;
TimeToExfil_h(~diagMask) = nan;
TimeToSurfaceSat_h(~diagMask) = nan;

MaxDepth_m(isinf(MaxDepth_m)) = nan;
MaxQexf_mm_h(isinf(MaxQexf_mm_h)) = nan;
MinZwt_m(isinf(MinZwt_m)) = nan;
MaxSatThick_m(isinf(MaxSatThick_m)) = nan;

if Opt.write_geotiff
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_MaxDepth_m.tif'), W.DEM_raster, MaxDepth_m);
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_TimeToWet_h.tif'), W.DEM_raster, TimeToWet_h);
    write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_EverWetMask.tif'), W.DEM_raster, double(EverWetMask));

    if any(isfinite(MinZwt_m(:)))
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_MinZWT_m.tif'), W.DEM_raster, MinZwt_m);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_MaxSaturatedThickness_m.tif'), W.DEM_raster, MaxSatThick_m);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_FinalZWT_m.tif'), W.DEM_raster, FinalZwt_m);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_FinalGWHead_m.tif'), W.DEM_raster, FinalHt_m);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_FinalSaturatedThickness_m.tif'), W.DEM_raster, FinalSatThick_m);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_TimeToSurfaceSaturation_h.tif'), W.DEM_raster, TimeToSurfaceSat_h);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_EverSurfaceSaturatedMask.tif'), W.DEM_raster, double(EverSurfaceSatMask));
    end

    if any(isfinite(MaxQexf_mm_h(:)))
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_MaxExfiltration_mm_h.tif'), W.DEM_raster, MaxQexf_mm_h);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_TimeToExfiltration_h.tif'), W.DEM_raster, TimeToExfil_h);
        write_geotiff_safe(fullfile(Dirs.Maps, 'Case04_EverExfiltrationMask.tif'), W.DEM_raster, double(EverExfilMask));
    end
end

%% ------------------------------------------------------------------------
% 15) FIGURES
% -------------------------------------------------------------------------
make_timeseries_figure(Dirs.Figures, TS, Opt);
make_volume_figure(Dirs.Figures, TS, Opt);
make_profile_figure(Dirs.Profiles, ProfileStore, Opt);
make_event_summary_text_file(Dirs.Tables, Event);

%% ------------------------------------------------------------------------
% 16) RETURN RESULTS
% -------------------------------------------------------------------------
Results = struct();
Results.Options = Opt;
Results.TimeSeries = TS;
Results.EventSummary = Event;
Results.StaticMaps.MaxDepth_m = MaxDepth_m;
Results.StaticMaps.MaxQexf_mm_h = MaxQexf_mm_h;
Results.StaticMaps.MinZwt_m = MinZwt_m;
Results.StaticMaps.MaxSatThick_m = MaxSatThick_m;
Results.StaticMaps.TimeToWet_h = TimeToWet_h;
Results.StaticMaps.TimeToExfil_h = TimeToExfil_h;
Results.StaticMaps.TimeToSurfaceSat_h = TimeToSurfaceSat_h;
Results.StaticMaps.EverWetMask = EverWetMask;
Results.StaticMaps.EverExfilMask = EverExfilMask;
Results.StaticMaps.EverSurfaceSatMask = EverSurfaceSatMask;
Results.OutputDirs = Dirs;

save(fullfile(Dirs.Root, 'Case04_GW_Advanced_Analysis_Results.mat'), 'Results', '-v7.3');

fprintf('\n============================================================\n');
fprintf('Advanced Case 4 groundwater analysis complete.\n');
fprintf('Tables:   %s\n', Dirs.Tables);
fprintf('Figures:  %s\n', Dirs.Figures);
fprintf('Profiles: %s\n', Dirs.Profiles);
fprintf('Rasters:  %s\n', Dirs.Maps);
fprintf('============================================================\n\n');

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function S = set_default(S, fieldName, value)
if ~isfield(S, fieldName) || isempty(S.(fieldName))
    S.(fieldName) = value;
end
end

function ensure_folder(p)
if ~isfolder(p)
    mkdir(p);
end
end

function x = gather_if_needed(x)
try
    if isa(x, 'gpuArray')
        x = gather(x);
    end
catch
end
end

function val = get_field_or(S, fieldName, defaultVal)
if isfield(S, fieldName)
    val = S.(fieldName);
else
    val = defaultVal;
end
end

function A = to_grid(x, nRows, nCols)
x = double(gather_if_needed(x));
if isscalar(x)
    A = x + zeros(nRows, nCols);
elseif isequal(size(x), [nRows, nCols])
    A = x;
else
    A = mean_valid(x(:)) + zeros(nRows, nCols);
end
end

function v = read_table_var(T, preferredName, defaultVal)
v = defaultVal;
if isempty(T)
    return
end

names = T.Properties.VariableNames;
idx = find(strcmpi(names, preferredName), 1, 'first');

if isempty(idx)
    idx = find(contains(lower(names), lower(preferredName)), 1, 'first');
end

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

function val = infer_scalar_soil_depth(W)
val = 1.0;
if isfield(W, 'Soil_Properties') && isfield(W.Soil_Properties, 'Soil_Depth')
    x = double(gather_if_needed(W.Soil_Properties.Soil_Depth));
    val = median(x(isfinite(x)));
end
end

function cellArea = infer_cell_area_grid(W, nRows, nCols)
if isfield(W.Wshed_Properties, 'cell_area')
    cellArea = W.Wshed_Properties.cell_area;
elseif isfield(W.Wshed_Properties, 'Resolution')
    cellArea = W.Wshed_Properties.Resolution.^2;
else
    cellArea = 1;
end
cellArea = to_grid(cellArea, nRows, nCols);
end

function t_min = get_time_minutes(running_control, flags)
t_min = [];

if ~isfield(running_control, 'time_records')
    return
end

t = running_control.time_records(:);
t = gather_if_needed(t);

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
candidates = { ...
    fullfile(tempDir, sprintf('save_map_hydro_%d.mat', storeID)), ...
    fullfile(tempDir, sprintf('save_map_hydro_%d', storeID))};

mapFile = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        mapFile = candidates{i};
        return
    end
end

error('Could not find map chunk %d in:\n%s', storeID, tempDir);
end

function nTotal = count_total_saved_frames(tempDir)
nTotal = 0;
store = 1;

while true
    try
        f = find_map_file(tempDir, store);
    catch
        break
    end

    Tmp = load(f, 'Maps');
    Maps = Tmp.Maps;

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
            if localID <= size(val,3)
                Z = double(val(:,:,localID));
                return
            end
        elseif ismatrix(val)
            Z = double(val);
            return
        end
    end
end

if isfield(Maps, 'Hydro') && isfield(Maps.Hydro, 'd')
    d = Maps.Hydro.d;
    Z = nan(size(d,1), size(d,2));
else
    Z = nan;
end
end

function [tf, val] = get_nested_field(S, pathStr)
tf = false;
val = [];

parts = strsplit(pathStr, '.');
val = S;

for i = 1:numel(parts)
    p = parts{i};
    if isstruct(val) && isfield(val, p)
        val = val.(p);
    else
        return
    end
end

tf = true;
end

function Z = clean_map(Z, validMask)
if isscalar(Z) && isnan(Z)
    Z = nan(size(validMask));
end
Z = double(gather_if_needed(Z));
Z(~validMask) = nan;
Z(~isfinite(Z)) = nan;
end

function [q_mm_h, q_mps] = convert_flux_to_mm_h_and_mps(q_raw, validMask)
q = clean_map(q_raw, validMask);
q_abs_max = max_valid(abs(q(validMask)));

if isnan(q_abs_max)
    q_mm_h = q;
    q_mps = q;
    return
end

% Heuristic:
% - Boussinesq q_exf is usually m/s, around 1e-10 to 1e-4.
% - If values are already > 0.05 in magnitude, assume mm/h.
if q_abs_max < 0.05
    q_mps = q;
    q_mm_h = q * 1000 * 3600;
else
    q_mm_h = q;
    q_mps = q / 1000 / 3600;
end
end

function m = area_weighted_mean(Z, A, validMask)
idx = validMask & isfinite(Z) & isfinite(A);
if ~any(idx(:))
    m = nan;
else
    m = nansum(Z(idx) .* A(idx)) / nansum(A(idx));
end
end

function p = percentile_valid(x, q)
x = x(isfinite(x));
if isempty(x)
    p = nan;
else
    p = prctile(x, q);
end
end

function m = mean_valid(x)
x = x(isfinite(x));
if isempty(x)
    m = nan;
else
    m = mean(x);
end
end

function m = max_valid(x)
x = x(isfinite(x));
if isempty(x)
    m = nan;
else
    m = max(x);
end
end

function m = min_valid(x)
x = x(isfinite(x));
if isempty(x)
    m = nan;
else
    m = min(x);
end
end

function y2 = fillmissing_linear(y)
y2 = y;
if all(~isfinite(y2))
    y2(:) = 0;
    return
end
idx = isfinite(y2);
x = (1:numel(y2))';
y2(~idx) = interp1(x(idx), y2(idx), x(~idx), 'linear', 'extrap');
end

function t0 = first_time_threshold(t, y, threshold)
idx = find(isfinite(y) & y >= threshold, 1, 'first');
if isempty(idx)
    t0 = nan;
else
    t0 = t(idx);
end
end

function T = add_event_metric(T, metric, value, units, note)
newRow = table(string(metric), value, string(units), string(note), ...
    'VariableNames', {'metric','value','units','note'});
T = [T; newRow];
end

function make_snapshot_figure(outDir, showFigures, time_h, DEM, d_mm, h_t_m, zwt_m, sat_thick_m, q_exf_mm_h, wetMask, exfilMask, surfaceSatMask, validMask, hasGWState, hasQexf)

if showFigures
    vis = 'on';
else
    vis = 'off';
end

fig = figure('Color','w','Units','pixels','Position',[50 50 1800 1050], ...
    'Visible', vis);

tiledlayout(2,4,'Padding','compact','TileSpacing','compact');

nexttile
plot_map(DEM, validMask, 'DEM [m]');
colorbar

nexttile
plot_map(d_mm/1000, validMask, sprintf('Surface depth [m]\nt = %.2f h', time_h));
colorbar

nexttile
if hasGWState
    plot_map(zwt_m, validMask, 'Depth to WT z_{wt} [m]');
else
    plot_map(nan(size(DEM)), validMask, 'z_{wt} not stored');
end
colorbar

nexttile
if hasGWState
    plot_map(sat_thick_m, validMask, 'Saturated thickness [m]');
else
    plot_map(nan(size(DEM)), validMask, 'GW depth not stored');
end
colorbar

nexttile
if hasGWState
    plot_map(h_t_m, validMask, 'Groundwater head h_t [m]');
else
    plot_map(nan(size(DEM)), validMask, 'h_t not stored');
end
colorbar

nexttile
if hasQexf
    plot_map(q_exf_mm_h, validMask, 'q_{exf} [mm/h]');
else
    plot_map(nan(size(DEM)), validMask, 'q_{exf} not stored');
end
colorbar

nexttile
plot_map(double(surfaceSatMask), validMask, 'Surface-saturated mask');
colorbar
caxis([0 1])

nexttile
plot_map(double(wetMask | exfilMask), validMask, 'Wet or exfiltration mask');
colorbar
caxis([0 1])

sgtitle(sprintf('Case 4 V-tilted groundwater diagnostics at t = %.2f h', time_h), ...
    'Interpreter','none');

outPng = fullfile(outDir, sprintf('Case04_snapshot_t_%07.3f_h.png', time_h));
exportgraphics(fig, outPng, 'Resolution', 250);
close(fig);

end

function plot_map(Z, validMask, ttl)
Z = double(Z);
Z(~validMask) = nan;
imagesc(Z);
axis image tight
set(gca,'YDir','normal');
title(ttl, 'Interpreter','tex');
xlabel('Column');
ylabel('Row');
end

function write_snapshot_rasters(outDir, DEM_raster, time_h, depth_m, h_t_m, zwt_m, sat_thick_m, q_exf_mm_h, wetMask, exfilMask, surfaceSatMask, validMask, hasGWState, hasQexf)

tag = sprintf('t_%07.3f_h', time_h);
tag = strrep(tag, '.', 'p');

depth_m(~validMask) = nan;
wetMask(~validMask) = nan;
exfilMask(~validMask) = nan;
surfaceSatMask(~validMask) = nan;

write_geotiff_safe(fullfile(outDir, ['Case04_Depth_' tag '.tif']), DEM_raster, depth_m);
write_geotiff_safe(fullfile(outDir, ['Case04_WetMask_' tag '.tif']), DEM_raster, wetMask);
write_geotiff_safe(fullfile(outDir, ['Case04_ExfilMask_' tag '.tif']), DEM_raster, exfilMask);
write_geotiff_safe(fullfile(outDir, ['Case04_SurfaceSatMask_' tag '.tif']), DEM_raster, surfaceSatMask);

if hasGWState
    h_t_m(~validMask) = nan;
    zwt_m(~validMask) = nan;
    sat_thick_m(~validMask) = nan;

    write_geotiff_safe(fullfile(outDir, ['Case04_GWHead_' tag '.tif']), DEM_raster, h_t_m);
    write_geotiff_safe(fullfile(outDir, ['Case04_ZWT_' tag '.tif']), DEM_raster, zwt_m);
    write_geotiff_safe(fullfile(outDir, ['Case04_SatThickness_' tag '.tif']), DEM_raster, sat_thick_m);
end

if hasQexf
    q_exf_mm_h(~validMask) = nan;
    write_geotiff_safe(fullfile(outDir, ['Case04_Qexf_mm_h_' tag '.tif']), DEM_raster, q_exf_mm_h);
end

end

function write_geotiff_safe(outPath, DEM_raster, Z)
try
    if isfield(DEM_raster, 'georef') && isfield(DEM_raster.georef, 'GeoKeyDirectoryTag')
        geotiffwrite(outPath, single(Z), DEM_raster.georef.SpatialRef, ...
            'GeoKeyDirectoryTag', DEM_raster.georef.GeoKeyDirectoryTag);
    else
        geotiffwrite(outPath, single(Z), DEM_raster.georef.SpatialRef);
    end
catch ME
    warning('Could not write GeoTIFF:\n%s\nReason: %s', outPath, ME.message);
end
end

function make_timeseries_figure(outDir, TS, Opt)

fig = figure('Color','w','Units','pixels','Position',[60 60 1800 1050]);

tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

nexttile
plot(TS.time_h, TS.outlet_Q_m3s, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Outlet Q [m^3/s]');
title('Surface outlet response without rainfall');
grid on

nexttile
plot(TS.time_h, TS.mean_depth_m, 'LineWidth', 2); hold on
plot(TS.time_h, TS.max_depth_m, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Surface depth [m]');
legend('Mean','Max','Location','best');
title('Surface-water depth generated by exfiltration');
grid on

nexttile
plot(TS.time_h, TS.mean_zwt_m, 'LineWidth', 2); hold on
plot(TS.time_h, TS.min_zwt_m, 'LineWidth', 2);
plot(TS.time_h, TS.p05_zwt_m, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('z_{wt} [m]');
legend('Mean','Min','P05','Location','best');
title('Water-table depth');
grid on

nexttile
plot(TS.time_h, TS.surface_saturated_area_pct, 'LineWidth', 2); hold on
plot(TS.time_h, TS.exfil_area_pct, 'LineWidth', 2);
plot(TS.time_h, TS.wet_area_pct, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Area [%]');
legend('Surface saturated','Exfiltration','Wet surface','Location','best');
title('Expansion of saturation/exfiltration/wet area');
grid on

nexttile
plot(TS.time_h, TS.mean_q_exf_mm_h, 'LineWidth', 2); hold on
plot(TS.time_h, TS.max_q_exf_mm_h, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('q_{exf} [mm/h]');
legend('Mean','Max','Location','best');
title('Groundwater exfiltration flux if stored');
grid on

nexttile
plot(TS.time_h, TS.center_mean_zwt_m, 'LineWidth', 2); hold on
plot(TS.time_h, TS.side_mean_zwt_m, 'LineWidth', 2);
plot(TS.time_h, TS.upstream_mean_zwt_m, '--', 'LineWidth', 1.5);
plot(TS.time_h, TS.downstream_mean_zwt_m, '--', 'LineWidth', 1.5);
xlabel('Time [h]');
ylabel('z_{wt} [m]');
legend('Center valley','Side hillslopes','Upstream','Downstream','Location','best');
title('Spatial groundwater redistribution');
grid on

sgtitle('Case 4 — V-tilted 2D groundwater propagation diagnostics');

exportgraphics(fig, fullfile(outDir, 'Case04_GW_TimeSeries.png'), 'Resolution', 250);
saveas(fig, fullfile(outDir, 'Case04_GW_TimeSeries.fig'));
close(fig);

end

function make_volume_figure(outDir, TS, Opt)

fig = figure('Color','w','Units','pixels','Position',[80 80 1700 900]);

tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile
plot(TS.time_h, TS.cumulative_outlet_m3, 'LineWidth', 2); hold on
plot(TS.time_h, TS.cumulative_exfiltration_m3, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Cumulative volume [m^3]');
legend('Outlet','Exfiltration to surface','Location','best');
title('Cumulative volumes');
grid on

nexttile
plot(TS.time_h, TS.surface_storage_m3, 'LineWidth', 2); hold on
plot(TS.time_h, TS.mean_gw_storage_m3, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Storage [m^3]');
legend('Surface storage','GW storage','Location','best');
title('Water storage components');
grid on

nexttile
plot(TS.time_h, TS.gw_storage_change_m3, 'LineWidth', 2); hold on
plot(TS.time_h, TS.total_storage_change_m3, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('\Delta storage [m^3]');
legend('\Delta GW','\Delta total','Location','best');
title('Storage depletion / redistribution');
grid on

nexttile
plot(TS.time_h, TS.mass_residual_mm, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Residual [mm]');
title('No-rainfall balance: -\DeltaS - V_{out}');
grid on

sgtitle('Case 4 — Volume diagnostics');

exportgraphics(fig, fullfile(outDir, 'Case04_GW_VolumeDiagnostics.png'), 'Resolution', 250);
saveas(fig, fullfile(outDir, 'Case04_GW_VolumeDiagnostics.fig'));
close(fig);

end

function make_profile_figure(outDir, ProfileStore, Opt)

if isempty(ProfileStore.times_h)
    return
end

n = numel(ProfileStore.times_h);

% Centerline profiles along flow direction
fig = figure('Color','w','Units','pixels','Position',[80 80 1700 900]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile
hold on
for i = 1:n
    y = ProfileStore.zwt_centerline{i};
    plot(y, 'LineWidth', 1.5);
end
xlabel('Row index along flow direction');
ylabel('z_{wt} [m]');
title('Centerline depth to water table');
grid on

nexttile
hold on
for i = 1:n
    y = ProfileStore.depth_centerline{i};
    plot(y, 'LineWidth', 1.5);
end
xlabel('Row index along flow direction');
ylabel('Depth [m]');
title('Centerline surface depth');
grid on

nexttile
hold on
for i = 1:n
    y = ProfileStore.qexf_centerline{i};
    plot(y, 'LineWidth', 1.5);
end
xlabel('Row index along flow direction');
ylabel('q_{exf} [mm/h]');
title('Centerline exfiltration');
grid on

legend(compose('%.2f h', ProfileStore.times_h), 'Location','bestoutside');
sgtitle('Case 4 — Longitudinal center-valley profiles');

exportgraphics(fig, fullfile(outDir, 'Case04_GW_LongitudinalProfiles.png'), 'Resolution', 250);
saveas(fig, fullfile(outDir, 'Case04_GW_LongitudinalProfiles.fig'));
close(fig);

% Cross-valley profiles at mid-row
fig = figure('Color','w','Units','pixels','Position',[80 80 1700 900]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile
hold on
for i = 1:n
    y = ProfileStore.cross_zwt_mid{i};
    plot(y, 'LineWidth', 1.5);
end
xlabel('Column index across valley');
ylabel('z_{wt} [m]');
title('Cross-valley depth to water table');
grid on

nexttile
hold on
for i = 1:n
    y = ProfileStore.cross_depth_mid{i};
    plot(y, 'LineWidth', 1.5);
end
xlabel('Column index across valley');
ylabel('Depth [m]');
title('Cross-valley surface depth');
grid on

nexttile
hold on
for i = 1:n
    y = ProfileStore.cross_qexf_mid{i};
    plot(y, 'LineWidth', 1.5);
end
xlabel('Column index across valley');
ylabel('q_{exf} [mm/h]');
title('Cross-valley exfiltration');
grid on

legend(compose('%.2f h', ProfileStore.times_h), 'Location','bestoutside');
sgtitle('Case 4 — Cross-valley profiles at mid-row');

exportgraphics(fig, fullfile(outDir, 'Case04_GW_CrossValleyProfiles.png'), 'Resolution', 250);
saveas(fig, fullfile(outDir, 'Case04_GW_CrossValleyProfiles.fig'));
close(fig);

end

function make_event_summary_text_file(outDir, Event)
fid = fopen(fullfile(outDir, 'Case04_GW_EventSummary.txt'), 'w');

fprintf(fid, 'Case 4 V-Tilted Groundwater Propagation Event Summary\n');
fprintf(fid, '============================================================\n\n');

for i = 1:height(Event)
    fprintf(fid, '%-38s : %12.5g %-8s | %s\n', ...
        Event.metric(i), Event.value(i), Event.units(i), Event.note(i));
end

fclose(fid);
end
