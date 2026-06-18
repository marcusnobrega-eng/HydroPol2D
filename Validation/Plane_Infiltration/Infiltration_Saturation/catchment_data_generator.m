%% ================================================================
%  SYNTHETIC CATCHMENTS FOR INFILTRATION + GROUNDWATER TEST SUITE
%  HydroPol2D testing
%
%  Purpose:
%   - Generate one raster folder per infiltration/saturation/GW test case
%   - Uniform tilted plane geometry for Cases 0–3
%   - V-tilted convergent catchment geometry for Case 4
%   - Analytical no-infiltration kinematic-wave hydrograph where rainfall > 0
%
%  Root folder:
%   /oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Plane_Infiltration/Infiltration_Saturation
%
%  Folder structure generated:
%
%   Infiltration_Saturation/
%       Input_Data_Sheets/
%       Case0_no_infiltration/
%           Static/
%       Case1A_supply_limited_deep_soil/
%           Static/
%       Case1B_capacity_limited_deep_soil/
%           Static/
%       Case2_shallow_soil_saturation_excess/
%           Static/
%       Case3_shallow_groundwater_saturation_excess/
%           Static/
%       Case4_vtilted_groundwater_propagation/
%           Static/
%
%  IMPORTANT STORAGE CONVENTIONS
%  ------------------------------------------------------------------------
%  1) In the infiltration / recharge benchmark cases, I_t is interpreted as
%     above-residual vadose storage:
%
%       total UZ capacity       = zwt * (theta_sat - theta_r)
%       initial UZ storage      = zwt * (theta_i   - theta_r)
%       initial remaining store = zwt * (theta_sat - theta_i)
%
%  2) In the 2D groundwater propagation case, Sy is intentionally treated as
%     an aquifer drainable porosity / specific yield. It is NOT forced to be
%     theta_sat - theta_r. This lets the Boussinesq groundwater propagation
%     test use a physically distinct storage coefficient.
%
%  IMPORTANT:
%   S0 = 0.01 means 1 percent slope.
%   S0 = 0.0001 would mean 0.01 percent slope.
% ================================================================

clear; clc; close all;

%% ================= ROOT OUTPUT FOLDER =================

rootDir = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Plane_Infiltration/Infiltration_Saturation';

if ~exist(rootDir, 'dir')
    mkdir(rootDir);
end

inputSheetDir = fullfile(rootDir, 'Input_Data_Sheets');

if ~exist(inputSheetDir, 'dir')
    mkdir(inputSheetDir);
end

%% ================= COMMON GRID PARAMETERS =================

dx = 20;             % grid resolution [m]

Lx = 2000;           % flow length [m]
Ly = 1000;           % plane width [m]

S0 = 0.01;           % default bed slope [-], 0.01 = 1 percent
n_manning = 0.015;   % default Manning roughness [s/m^(1/3)]

z_outlet = 0;        % outlet / valley reference elevation [m]

dt_benchmark_sec = 10 * 60;   % analytical benchmark time step [s]

%% ================= COMMON LAND-SURFACE PROPERTIES =================

LULC_val   = 1;
Albedo_val = 0.5;
LAI_val    = 0.5;

%% ================= BASE SOIL TABLE =================
% These are the default soil hydraulic parameters to be used by HydroPol2D.
% The SOIL raster stores the soil index.
%
% Case-specific overrides are applied inside the generation loop.

SOIL_table = table();

SOIL_table.Soil_type = { ...
    'left'; ...
    'right'; ...
    'center'};

SOIL_table.Index = [1; 2; 3];

SOIL_table.ksat_mm_h = [0.5; 0.3; 0.3];
SOIL_table.n_vg = [1.2; 1.09; 1.09];
SOIL_table.alpha_vg_1_m = [2.0; 0.8; 0.8];

SOIL_table.theta_sat = [0.385; 0.385; 0.385];
SOIL_table.theta_r   = [0.068; 0.068; 0.068];
SOIL_table.theta_i   = [0.0997; 0.0997; 0.0997];

% Default for local recharge-to-water-table infiltration benchmarks.
% Case 4 overrides this with aquifer Sy = 0.05.
SOIL_table.Sy = SOIL_table.theta_sat - SOIL_table.theta_r;

% Default groundwater saturated conductivity.
% Case 4 overrides this with a much larger value to produce visible
% Boussinesq drainage/exfiltration over a synthetic benchmark timescale.
SOIL_table.ksat_gw_mm_h = [6; 6; 6];

% Default soil depth in table.
% Per-case soil depth is defined below and written into each case's DTB raster.
SOIL_table.Soil_Depth_m = ones(3,1);

% Infiltration module defaults / recommended values
Ltop_m = 0.05;
dh_max_m = 1.0;
l_vg = 0.5;

%% ================= CASE DEFINITIONS =================
% DTB_m:
%   depth to bedrock / soil depth used in the raster.
%
% zwt0_m:
%   initial water table depth below land surface.
%
% Generated groundwater table raster:
%   GW_table = DEM - zwt0_m
%
% Equivalent old HydroPol2D GW_depth value:
%   GW_depth = DTB_m - zwt0_m
%
% because:
%   GW_table = DEM - DTB + GW_depth
%
% geometry_mode:
%   'plane'     -> regular downstream tilted plane
%   'v_tilted'  -> downstream tilted valley with cross-slope convergence

Cases = struct([]);

% ------------------------------------------------
% Case 0: no-infiltration hydraulic control
% ------------------------------------------------
Cases(1).folder = 'Case0_no_infiltration';
Cases(1).description = 'No-infiltration hydraulic control. Used to verify the tilted-plane routing benchmark.';
Cases(1).geometry_mode = 'plane';
Cases(1).S_long = S0;
Cases(1).S_cross = 0.0;
Cases(1).flag_infiltration = 0;
Cases(1).flag_groundwater_modeling = 0;
Cases(1).flag_baseflow = 0;
Cases(1).rain_mm_h = 20;
Cases(1).rain_duration_min = 9 * 60;
Cases(1).sim_duration_min  = 9 * 60;
Cases(1).soil_mode = 'uniform';
Cases(1).uniform_soil_index = 1;
Cases(1).DTB_m = 5.0;
Cases(1).zwt0_m = 5.0;
Cases(1).theta_i = 0.0997;
Cases(1).Sy_override = NaN;
Cases(1).ksat_gw_override_mm_h = NaN;
Cases(1).expected_note = 'Outlet should approach the no-infiltration analytical Qss.';

% ------------------------------------------------
% Case 1A: supply-limited infiltration, deep soil
% ------------------------------------------------
Cases(2).folder = 'Case1A_supply_limited_deep_soil';
Cases(2).description = 'Supply-limited infiltration. Rainfall is below expected infiltration capacity; deep soil prevents saturation.';
Cases(2).geometry_mode = 'plane';
Cases(2).S_long = S0;
Cases(2).S_cross = 0.0;
Cases(2).flag_infiltration = 1;
Cases(2).flag_groundwater_modeling = 0;
Cases(2).flag_baseflow = 0;
Cases(2).rain_mm_h = 5;
Cases(2).rain_duration_min = 9 * 60;
Cases(2).sim_duration_min  = 9 * 60;
Cases(2).soil_mode = 'uniform';
Cases(2).uniform_soil_index = 1;
Cases(2).DTB_m = 5.0;
Cases(2).zwt0_m = 5.0;
Cases(2).theta_i = 0.0997;
Cases(2).Sy_override = NaN;
Cases(2).ksat_gw_override_mm_h = NaN;
Cases(2).expected_note = 'Mean f should approach rainfall, approximately 5 mm/h. Outlet should remain close to zero after transient.';

% ------------------------------------------------
% Case 1B: capacity-limited infiltration, deep soil
% ------------------------------------------------
Cases(3).folder = 'Case1B_capacity_limited_deep_soil';
Cases(3).description = 'Capacity-limited infiltration. Rainfall exceeds expected infiltration capacity; deep soil prevents saturation.';
Cases(3).geometry_mode = 'plane';
Cases(3).S_long = S0;
Cases(3).S_cross = 0.0;
Cases(3).flag_infiltration = 1;
Cases(3).flag_groundwater_modeling = 0;
Cases(3).flag_baseflow = 0;
Cases(3).rain_mm_h = 20;
Cases(3).rain_duration_min = 9 * 60;
Cases(3).sim_duration_min  = 9 * 60;
Cases(3).soil_mode = 'uniform';
Cases(3).uniform_soil_index = 1;
Cases(3).DTB_m = 5.0;
Cases(3).zwt0_m = 5.0;
Cases(3).theta_i = 0.0997;
Cases(3).Sy_override = NaN;
Cases(3).ksat_gw_override_mm_h = NaN;
Cases(3).expected_note = 'Mean f should approach about 10.5 mm/h if Ltop=0.05 m and dh_max=1 m. Outlet should approach rainfall-excess flow.';

% ------------------------------------------------
% Case 2: shallow soil saturation excess
% ------------------------------------------------
Cases(4).folder = 'Case2_shallow_soil_saturation_excess';
Cases(4).description = 'Shallow soil saturation-excess test. Soil storage should fill and infiltration should eventually collapse.';
Cases(4).geometry_mode = 'plane';
Cases(4).S_long = S0;
Cases(4).S_cross = 0.0;
Cases(4).flag_infiltration = 1;
Cases(4).flag_groundwater_modeling = 0;
Cases(4).flag_baseflow = 0;
Cases(4).rain_mm_h = 20;
Cases(4).rain_duration_min = 12 * 60;
Cases(4).sim_duration_min  = 12 * 60;
Cases(4).soil_mode = 'uniform';
Cases(4).uniform_soil_index = 1;
Cases(4).DTB_m = 0.30;
Cases(4).zwt0_m = 0.30;
Cases(4).theta_i = 0.0997;
Cases(4).Sy_override = NaN;
Cases(4).ksat_gw_override_mm_h = NaN;
Cases(4).expected_note = 'Mean f should be near capacity early, then approach zero after the shallow soil storage fills.';

% ------------------------------------------------
% Case 3: near-surface groundwater saturation excess
% ------------------------------------------------
Cases(5).folder = 'Case3_shallow_groundwater_saturation_excess';
Cases(5).description = ['Near-surface groundwater saturation-excess test. ', ...
    'Thin soil/aquifer column, water table starts only 5 cm below land surface, ', ...
    'and rainfall is applied for 48 h.'];
Cases(5).geometry_mode = 'plane';
Cases(5).S_long = S0;
Cases(5).S_cross = 0.0;
Cases(5).flag_infiltration = 1;
Cases(5).flag_groundwater_modeling = 1;
Cases(5).flag_baseflow = 0;
Cases(5).rain_mm_h = 20;
Cases(5).rain_duration_min = 48 * 60;
Cases(5).sim_duration_min  = 48 * 60;
Cases(5).soil_mode = 'uniform';
Cases(5).uniform_soil_index = 1;
Cases(5).DTB_m = 0.30;
Cases(5).zwt0_m = 0.05;
Cases(5).theta_i = 0.0997;
Cases(5).Sy_override = NaN;                 % local closure uses theta_sat-theta_r in code
Cases(5).ksat_gw_override_mm_h = NaN;
Cases(5).expected_note = ['Very rapid saturation-excess onset is expected. ', ...
    'Initial remaining storage is only about 14.3 mm.'];

% ------------------------------------------------
% Case 4: V-tilted 2D groundwater propagation / exfiltration
% ------------------------------------------------
Cases(6).folder = 'Case4_vtilted_groundwater_propagation';
Cases(6).description = ['V-tilted convergent catchment for testing 2D Boussinesq ', ...
    'groundwater propagation without rainfall. The initial water table is ', ...
    '0.25 m below land surface everywhere. Lateral groundwater flow should ', ...
    'converge toward the central valley and downstream outlet, potentially ', ...
    'producing localized exfiltration/surface water.'];
Cases(6).geometry_mode = 'v_tilted';

% Longitudinal/downstream slope and cross-valley side slope.
% The cross-slope creates a V-shaped valley along the centerline.
Cases(6).S_long  = 0.005;   % 0.5% downstream slope
Cases(6).S_cross = 0.050;   % 5% side slope toward center valley

Cases(6).flag_infiltration = 0;
Cases(6).flag_groundwater_modeling = 1;
Cases(6).flag_baseflow = 1;

% No rainfall. We run long enough for groundwater redistribution.
Cases(6).rain_mm_h = 0;
Cases(6).rain_duration_min = 0;
Cases(6).sim_duration_min  = 72 * 60;

Cases(6).soil_mode = 'uniform';
Cases(6).uniform_soil_index = 1;

% Aquifer / soil geometry.
Cases(6).DTB_m  = 1.0;     % soil/aquifer thickness [m]
Cases(6).zwt0_m = 0.25;    % water table depth below surface [m]
Cases(6).theta_i = 0.0997;

% Case-specific aquifer storage and conductivity.
% These are intentionally larger/faster than natural hillslope values so
% that groundwater propagation is visible over a 72 h benchmark.
Cases(6).Sy_override = 0.05;                  % aquifer drainable porosity [-]
Cases(6).ksat_gw_override_mm_h = 10000;       % 10 m/h = 2.78e-3 m/s

Cases(6).expected_note = ['No rainfall input. Any surface water should come from ', ...
    'groundwater convergence/exfiltration through the Boussinesq module. ', ...
    'Expected response: water-table redistribution toward the central valley, ', ...
    'localized q_exf > 0, then shallow overland routing to the outlet.'];

%% ================= DOMAIN =================

nx = round(Ly / dx);
ny = round(Lx / dx);

[X, Y] = meshgrid(1:nx, 1:ny);

x = (X - 1) * dx;   % width direction [m]
y = (Y - 1) * dx;   % length direction [m]

%% ================= GEOREFERENCE =================

x0 = 0;
y0 = Lx;

R = maprefcells( ...
    [x0, x0 + nx * dx], ...
    [y0 - ny * dx, y0], ...
    [ny, nx], ...
    'ColumnsStartFrom', 'north');

epsgCode = 3857;

%% ================= MASTER CASE SUMMARY TABLE =================

Master = table();

%% ================= GENERATE ALL CASES =================

for icase = 1:numel(Cases)

    Cfg = Cases(icase);

    caseRootDir = fullfile(rootDir, Cfg.folder);
    staticDir   = fullfile(caseRootDir, 'Static');

    if ~exist(caseRootDir, 'dir')
        mkdir(caseRootDir);
    end

    if ~exist(staticDir, 'dir')
        mkdir(staticDir);
    end

    fprintf('\n============================================================\n');
    fprintf('Generating %s\n', Cfg.folder);
    fprintf('%s\n', Cfg.description);
    fprintf('Case root folder:\n%s\n', caseRootDir);
    fprintf('Static raster folder:\n%s\n', staticDir);
    fprintf('============================================================\n');

    %% ---------- Case-specific DEM ----------

    DEM = build_case_dem(Cfg, x, y, Ly, z_outlet);

    %% ---------- Uniform fields ----------

    LULC = LULC_val * ones(ny, nx);
    DTB  = Cfg.DTB_m * ones(ny, nx);

    Albedo = Albedo_val * ones(ny, nx);
    LAI    = LAI_val    * ones(ny, nx);

    Initial_SM = Cfg.theta_i * ones(ny, nx);

    % Kept for reproducibility in MAT file, but not written as GeoTIFF
    % because your Static folder convention does not include Manning_n.tif.
    Manning_n = n_manning * ones(ny, nx);

    %% ---------- Soil map ----------

    switch lower(Cfg.soil_mode)

        case 'uniform'

            SOIL = Cfg.uniform_soil_index * ones(ny, nx);

        case 'three_zone'

            SOIL = zeros(ny, nx);

            c1 = floor(nx / 3);
            c2 = floor(2 * nx / 3);

            left_cols   = 1:c1;
            center_cols = (c1 + 1):c2;
            right_cols  = (c2 + 1):nx;

            SOIL(:, left_cols)   = 1;   % left soil
            SOIL(:, center_cols) = 3;   % center soil
            SOIL(:, right_cols)  = 2;   % right soil

        otherwise

            error('Unknown soil_mode: %s', Cfg.soil_mode);

    end

    %% ---------- Groundwater table ----------

    % Initial water table elevation [m].
    % zwt0_m is depth below land surface.
    GW_table = DEM - Cfg.zwt0_m;

    % Equivalent old-style GW_depth value:
    % GW_table = DEM - DTB + GW_depth
    GW_depth_equivalent = Cfg.DTB_m - Cfg.zwt0_m;

    if GW_depth_equivalent < -1e-9
        warning('Case %s has zwt0_m deeper than DTB_m. Water table is below bedrock in this setup.', Cfg.folder);
    end

    %% ---------- Analytical no-infiltration benchmark ----------

    rain_mm_h = Cfg.rain_mm_h;
    rain_duration_min = Cfg.rain_duration_min;
    sim_duration_min = Cfg.sim_duration_min;

    sim_duration_sec  = sim_duration_min  * 60;
    t_sec = (0:dt_benchmark_sec:sim_duration_sec)';

    if rain_mm_h > 0

        r = rain_mm_h / 1000 / 3600;     % [m/s]
        rain_duration_sec = rain_duration_min * 60;

        alpha_kw = sqrt(max(Cfg.S_long, 1e-12)) / n_manning;

        % Time to equilibrium / time of concentration [s]
        t_eq_sec = (Lx / (alpha_kw * r^(2/3)))^(3/5);

        % Steady unit-width discharge [m2/s]
        q_eq = r * Lx;

        % Total steady no-infiltration discharge [m3/s]
        Q_eq_no_infiltration = q_eq * Ly;

        Q_analytic = zeros(size(t_sec));

        for ii = 1:numel(t_sec)

            tt = t_sec(ii);

            if tt <= rain_duration_sec

                if tt < t_eq_sec
                    q_unit = alpha_kw * (r * tt)^(5/3);
                else
                    q_unit = q_eq;
                end

            else

                % Recession intentionally left as NaN.
                q_unit = NaN;

            end

            Q_analytic(ii) = q_unit * Ly;

        end

    else

        % No rainfall benchmark. The no-rain hydraulic input is zero.
        t_eq_sec = NaN;
        Q_eq_no_infiltration = 0;
        Q_analytic = zeros(size(t_sec));

    end

    Benchmark = table();
    Benchmark.time_sec = t_sec;
    Benchmark.time_min = t_sec / 60;
    Benchmark.Q_no_infiltration_m3s = Q_analytic;

    writetable(Benchmark, fullfile(staticDir, 'Analytical_KinematicWave_Hydrograph.csv'));

    %% ---------- Per-case soil table ----------

    SOIL_table_case = SOIL_table;

    % Overwrite Soil_Depth_m in the table to match the generated DTB raster.
    SOIL_table_case.Soil_Depth_m(:) = Cfg.DTB_m;

    % Overwrite theta_i for this case.
    SOIL_table_case.theta_i(:) = Cfg.theta_i;

    % Optional case-specific aquifer Sy override.
    if isfield(Cfg, 'Sy_override') && isfinite(Cfg.Sy_override)
        SOIL_table_case.Sy(:) = Cfg.Sy_override;
    end

    % Optional case-specific groundwater conductivity override.
    if isfield(Cfg, 'ksat_gw_override_mm_h') && isfinite(Cfg.ksat_gw_override_mm_h)
        SOIL_table_case.ksat_gw_mm_h(:) = Cfg.ksat_gw_override_mm_h;
    end

    % Add infiltration model parameters to the per-case table.
    SOIL_table_case.Ltop_m   = Ltop_m   * ones(height(SOIL_table_case), 1);
    SOIL_table_case.dh_max_m = dh_max_m * ones(height(SOIL_table_case), 1);
    SOIL_table_case.l_vg     = l_vg     * ones(height(SOIL_table_case), 1);

    %% ---------- Expected infiltration / storage metrics ----------

    if strcmpi(Cfg.soil_mode, 'uniform')

        sidx = Cfg.uniform_soil_index;

        soil_row = find(SOIL_table_case.Index == sidx, 1, 'first');

        ksat_mm_h = SOIL_table_case.ksat_mm_h(soil_row);
        theta_s = SOIL_table_case.theta_sat(soil_row);
        theta_r = SOIL_table_case.theta_r(soil_row);
        theta_i = Cfg.theta_i;
        Sy_case = SOIL_table_case.Sy(soil_row);
        ksat_gw_case = SOIL_table_case.ksat_gw_mm_h(soil_row);

        Cmax_mm_h = ksat_mm_h * (dh_max_m / Ltop_m + 1);

        storage_capacity_total_above_residual_mm = ...
            Cfg.zwt0_m * (theta_s - theta_r) * 1000;

        initial_storage_above_residual_mm = ...
            Cfg.zwt0_m * max(theta_i - theta_r, 0) * 1000;

        remaining_storage_initial_mm = ...
            Cfg.zwt0_m * (theta_s - theta_i) * 1000;

        aquifer_saturated_thickness_initial_m = ...
            max(Cfg.DTB_m - Cfg.zwt0_m, 0);

        aquifer_drainable_storage_initial_mm = ...
            aquifer_saturated_thickness_initial_m * Sy_case * 1000;

        if Cfg.flag_infiltration == 1 && Cfg.rain_mm_h > 0
            expected_f_early_mm_h = min(Cfg.rain_mm_h, Cmax_mm_h);
        else
            expected_f_early_mm_h = 0;
        end

        if expected_f_early_mm_h > 0
            expected_fill_time_h = remaining_storage_initial_mm / expected_f_early_mm_h;
        else
            expected_fill_time_h = NaN;
        end

        expected_excess_mm_h = max(Cfg.rain_mm_h - expected_f_early_mm_h, 0);

        expected_Q_after_infiltration_m3s = ...
            expected_excess_mm_h / 1000 / 3600 * Lx * Ly;

    else

        Cmax_mm_h = NaN;
        Sy_case = NaN;
        ksat_gw_case = NaN;
        storage_capacity_total_above_residual_mm = NaN;
        initial_storage_above_residual_mm = NaN;
        remaining_storage_initial_mm = NaN;
        aquifer_saturated_thickness_initial_m = NaN;
        aquifer_drainable_storage_initial_mm = NaN;
        expected_f_early_mm_h = NaN;
        expected_fill_time_h = NaN;
        expected_excess_mm_h = NaN;
        expected_Q_after_infiltration_m3s = NaN;

    end

    %% ---------- Case config table ----------

    Case_Config = table();

    Case_Config.CaseFolder = string(Cfg.folder);
    Case_Config.Description = string(Cfg.description);

    Case_Config.geometry_mode = string(Cfg.geometry_mode);
    Case_Config.S_long = Cfg.S_long;
    Case_Config.S_cross = Cfg.S_cross;

    Case_Config.flag_infiltration = Cfg.flag_infiltration;
    Case_Config.flag_groundwater_modeling = Cfg.flag_groundwater_modeling;
    Case_Config.flag_baseflow = Cfg.flag_baseflow;

    Case_Config.Lx_m = Lx;
    Case_Config.Ly_m = Ly;
    Case_Config.dx_m = dx;
    Case_Config.nx = nx;
    Case_Config.ny = ny;
    Case_Config.Area_m2 = Lx * Ly;

    Case_Config.S0_default = S0;
    Case_Config.n_manning = n_manning;

    Case_Config.rain_mm_h = Cfg.rain_mm_h;
    Case_Config.rain_duration_min = Cfg.rain_duration_min;
    Case_Config.sim_duration_min = Cfg.sim_duration_min;

    Case_Config.DTB_m = Cfg.DTB_m;
    Case_Config.zwt0_m = Cfg.zwt0_m;
    Case_Config.GW_depth_equivalent_m = GW_depth_equivalent;

    Case_Config.theta_i = Cfg.theta_i;
    Case_Config.theta_r = SOIL_table_case.theta_r(1);
    Case_Config.theta_sat = SOIL_table_case.theta_sat(1);
    Case_Config.Sy_case = Sy_case;
    Case_Config.ksat_gw_case_mm_h = ksat_gw_case;

    Case_Config.Ltop_m = Ltop_m;
    Case_Config.dh_max_m = dh_max_m;
    Case_Config.l_vg = l_vg;

    Case_Config.Qss_no_infiltration_m3s = Q_eq_no_infiltration;
    Case_Config.tc_no_infiltration_min = t_eq_sec / 60;

    Case_Config.expected_Cmax_uniform_mm_h = Cmax_mm_h;
    Case_Config.expected_f_early_mm_h = expected_f_early_mm_h;
    Case_Config.expected_excess_mm_h = expected_excess_mm_h;
    Case_Config.expected_Q_after_infiltration_m3s = expected_Q_after_infiltration_m3s;

    Case_Config.storage_capacity_total_above_residual_mm = storage_capacity_total_above_residual_mm;
    Case_Config.initial_storage_above_residual_mm = initial_storage_above_residual_mm;
    Case_Config.remaining_storage_initial_mm = remaining_storage_initial_mm;
    Case_Config.expected_fill_time_h = expected_fill_time_h;

    Case_Config.aquifer_saturated_thickness_initial_m = aquifer_saturated_thickness_initial_m;
    Case_Config.aquifer_drainable_storage_initial_mm = aquifer_drainable_storage_initial_mm;

    Case_Config.ExpectedNote = string(Cfg.expected_note);

    writetable(Case_Config, fullfile(caseRootDir, 'Case_Config.csv'));
    writetable(SOIL_table_case, fullfile(caseRootDir, 'SOIL_table.csv'));

    %% ---------- Safety check ----------

    rasterNames = {'DEM','LULC','SOIL','DTB','Albedo','LAI','Initial_SM','GW_table'};

    for kk = 1:numel(rasterNames)

        A = eval(rasterNames{kk});

        if any(isnan(A(:)))
            error('%s contains NaN values in case %s.', rasterNames{kk}, Cfg.folder);
        end

        if any(isinf(A(:)))
            error('%s contains Inf values in case %s.', rasterNames{kk}, Cfg.folder);
        end

    end

    %% ---------- Write GeoTIFF rasters into Static folder ----------

    geotiffwrite(fullfile(staticDir, 'DEM.tif'),        single(DEM),        R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'LULC.tif'),       single(LULC),       R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'SOIL.tif'),       single(SOIL),       R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'DTB.tif'),        single(DTB),        R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'Albedo.tif'),     single(Albedo),     R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'LAI.tif'),        single(LAI),        R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'Initial_SM.tif'), single(Initial_SM), R, 'CoordRefSysCode', epsgCode);
    geotiffwrite(fullfile(staticDir, 'GW_table.tif'),   single(GW_table),   R, 'CoordRefSysCode', epsgCode);

    %% ---------- Save MAT file for reproducibility ----------

    save(fullfile(caseRootDir, 'Case_Rasters_And_Config.mat'), ...
        'Cfg', 'Case_Config', 'SOIL_table_case', ...
        'DEM', 'LULC', 'SOIL', 'DTB', 'Albedo', 'LAI', ...
        'Initial_SM', 'GW_table', 'Manning_n', ...
        'Benchmark', 'R', 'epsgCode', ...
        'x', 'y', 'X', 'Y');

    %% ---------- Quick overview plot ----------

    fig = figure('Color','w','Position',[100 100 1600 900]);

    subplot(2,3,1)
    imagesc(DEM)
    axis image off
    colorbar
    title('DEM [m]')

    subplot(2,3,2)
    imagesc(SOIL)
    axis image off
    colorbar
    title('SOIL index')

    subplot(2,3,3)
    imagesc(DTB)
    axis image off
    colorbar
    title('DTB / Soil depth [m]')

    subplot(2,3,4)
    imagesc(GW_table)
    axis image off
    colorbar
    title('Initial GW table elevation [m]')

    subplot(2,3,5)
    imagesc(DEM - GW_table)
    axis image off
    colorbar
    title('Initial z_{wt} [m]')

    subplot(2,3,6)
    plot(Benchmark.time_min, Benchmark.Q_no_infiltration_m3s, 'LineWidth', 2)
    grid on
    xlabel('Time [min]')
    ylabel('Q [m^3/s]')
    title(sprintf('No-rain/no-inf benchmark Q=%.3f m^3/s', Q_eq_no_infiltration))

    sgtitle(strrep(Cfg.folder, '_', '\_'));

    exportgraphics(fig, fullfile(staticDir, 'overview_plane_with_flow.png'), 'Resolution', 300);
    exportgraphics(fig, fullfile(staticDir, 'overview_plane.png'), 'Resolution', 300);

    close(fig);

    %% ---------- Append to master table ----------

    Master = [Master; Case_Config]; %#ok<AGROW>

    %% ---------- Console summary ----------

    fprintf('Geometry mode                      = %s\n', Cfg.geometry_mode);
    fprintf('Longitudinal slope                 = %.5f\n', Cfg.S_long);
    fprintf('Cross slope                        = %.5f\n', Cfg.S_cross);
    fprintf('Rainfall                           = %.3f mm/h\n', Cfg.rain_mm_h);
    fprintf('Rain duration                      = %.3f min\n', Cfg.rain_duration_min);
    fprintf('Simulation duration                = %.3f min\n', Cfg.sim_duration_min);
    fprintf('DTB                                = %.3f m\n', Cfg.DTB_m);
    fprintf('Initial water table depth          = %.3f m below surface\n', Cfg.zwt0_m);
    fprintf('Initial saturated thickness        = %.3f m\n', aquifer_saturated_thickness_initial_m);
    fprintf('Sy case                            = %.5f\n', Sy_case);
    fprintf('Ksat_gw case                       = %.3f mm/h\n', ksat_gw_case);
    fprintf('Aquifer drainable storage          = %.3f mm\n', aquifer_drainable_storage_initial_mm);
    fprintf('No-infiltration Qss                = %.6f m3/s\n', Q_eq_no_infiltration);

    if strcmpi(Cfg.soil_mode, 'uniform') && Cfg.rain_mm_h > 0
        fprintf('Expected Cmax                      = %.6f mm/h\n', Cmax_mm_h);
        fprintf('Expected early f                   = %.6f mm/h\n', expected_f_early_mm_h);
        fprintf('Expected excess                    = %.6f mm/h\n', expected_excess_mm_h);
        fprintf('Expected Q after infiltration      = %.6f m3/s\n', expected_Q_after_infiltration_m3s);
        fprintf('Initial remaining storage          = %.6f mm\n', remaining_storage_initial_mm);
        fprintf('Expected fill time                 = %.6f h\n', expected_fill_time_h);
    end

    fprintf('Finished writing case root folder:\n%s\n', caseRootDir);
    fprintf('Finished writing Static raster folder:\n%s\n', staticDir);

end

%% ================= WRITE MASTER FILES =================

writetable(Master, fullfile(inputSheetDir, 'Plane_Infiltration_Case_Master_Config.csv'));
writetable(SOIL_table, fullfile(inputSheetDir, 'Base_SOIL_table.csv'));

save(fullfile(inputSheetDir, 'Plane_Infiltration_Case_Master_Config.mat'), ...
    'Cases', 'Master', 'SOIL_table', ...
    'rootDir', 'inputSheetDir', ...
    'dx', 'Lx', 'Ly', 'S0', 'n_manning', ...
    'Ltop_m', 'dh_max_m', 'l_vg');

fprintf('\n============================================================\n');
fprintf('ALL SYNTHETIC CASES GENERATED SUCCESSFULLY\n');
fprintf('Root folder:\n%s\n', rootDir);
fprintf('Master config:\n%s\n', fullfile(inputSheetDir, 'Plane_Infiltration_Case_Master_Config.csv'));
fprintf('Base soil table:\n%s\n', fullfile(inputSheetDir, 'Base_SOIL_table.csv'));
fprintf('============================================================\n');

%% ========================================================================
% LOCAL HELPER FUNCTION
% ========================================================================

function DEM = build_case_dem(Cfg, x, y, Ly, z_outlet)

    switch lower(Cfg.geometry_mode)

        case 'plane'

            % Simple longitudinal tilted plane.
            DEM = z_outlet + Cfg.S_long .* y;

        case 'v_tilted'

            % V-shaped convergent valley with downstream slope.
            %
            % x is the cross-width coordinate.
            % The valley centerline is at x = Ly/2.
            % Elevation increases away from the valley centerline by:
            %
            %   S_cross * abs(x - Ly/2)
            %
            % and increases upstream/downstream by:
            %
            %   S_long * y
            %
            % This produces lateral groundwater gradients toward the valley.
            valley_x = Ly / 2;
            cross_distance = abs(x - valley_x);

            DEM = z_outlet + Cfg.S_long .* y + Cfg.S_cross .* cross_distance;

        otherwise

            error('Unknown geometry_mode: %s', Cfg.geometry_mode);

    end

end
