%% ========================================================================
% input_data_bypass_script.m
% ========================================================================
% HydroPol2D | Standalone MATLAB definition of ALL inputs normally handled
% inside input_data_script through spreadsheets.
%
% PURPOSE
% ------------------------------------------------------------------------
% This script is the USER-EDITABLE MATLAB alternative to the spreadsheet-
% based input workflow used by HydroPol2D.
%
% When bypass mode is active later in input_data_script, this file will be
% executed and must define a single struct:
%
%     InputData_Bypass
%
% The revised input_data_script will then read this struct and translate its
% contents into the SAME legacy variables traditionally built from Excel,
% such as:
%
%   flags
%   running_control
%   Courant_Parameters
%   Design_Storm_Parameters
%   Human_Instability
%   LULC_name, LULC_index, lulc_parameters, n_lulc, imp_index
%   SOIL_name, SOIL_index, soil_parameters, n_soil
%   Rainfall_Parameters
%   Input_Rainfall
%   Input_Transpiration
%   Input_Evaporation
%   Inflow_Parameters
%   Stage_Parameters
%   Reservoir_Data
%   Snow_Properties
%
% IMPORTANT PHILOSOPHY
% ------------------------------------------------------------------------
% 1) This file is the CLEAN USER INTERFACE.
% 2) It is NOT the final internal HydroPol2D format.
% 3) The updated input_data_script will later translate these sections into
%    the old legacy variables used by the rest of the model.
%
% In other words:
%
%   YOU edit InputData_Bypass here.
%   input_data_script later converts it into the classic HydroPol2D inputs.
%
% UNITS
% ------------------------------------------------------------------------
% Unless otherwise noted, keep the same units expected by your original
% model and spreadsheets:
%
% - time_step_model, min_time_step, max_time_step, routing times: minutes
% - rainfall intensity: mm/h
% - map times for rainfall/ETP schedules: minutes from simulation start
% - inflow discharge: m^3/s
% - stage hydrograph: m
% - projected coordinates x/y: map CRS units (typically meters)
%
% HOW TO USE THIS FILE
% ------------------------------------------------------------------------
% 1) Copy this file into the HydroPol2D working folder.
% 2) Rename it exactly to:
%       input_data_bypass_script.m
% 3) In HydroPol2D_V115, set:
%       use_inputdata_bypass = 1;
% 4) Edit the sections below.
%
% SAFE WORKFLOW RECOMMENDATION
% ------------------------------------------------------------------------
% Start by bypassing only one section (for example, inflow hydrographs),
% test the model, and only then migrate the remaining inputs.
%
% CONSISTENCY WITH UPDATED input_data_script
% ------------------------------------------------------------------------
% This version was written to match the updated input_data_script that
% expects, in bypass mode:
%
%   InputData_Bypass.general
%   InputData_Bypass.flags
%   InputData_Bypass.LULC
%   InputData_Bypass.SOIL
%
% and conditionally:
%
%   InputData_Bypass.Human_Instability
%   InputData_Bypass.Snow_Properties
%   InputData_Bypass.Rainfall_Parameters
%   InputData_Bypass.Input_Rainfall
%   InputData_Bypass.Input_Transpiration
%   InputData_Bypass.Input_Evaporation
%   InputData_Bypass.Inflow
%   InputData_Bypass.Stage
%   InputData_Bypass.Reservoir_Data
%
% IMPORTANT
% ------------------------------------------------------------------------
% - date_begin AND date_end are required here.
% - record_time_spatial_ETP must exist in general.
% - Design storm fields must live in InputData_Bypass.general with names:
%       RP, Rainfall_Duration, K, a, b, c, dt_design
% - If flag_resample = 1, make sure GIS_data.resolution_resample is defined
%   somewhere upstream in preprocessing, because input_data_script still
%   uses it later for inlet coordinate conversion and area calculations.
% ========================================================================

%% Initialize master struct
% ------------------------------------------------------------------------
% This is the ONLY variable that the bypass script must leave behind.
% The updated input_data_script will later read its fields.
% -------------------------------------------------------------------------
InputData_Bypass = struct();

%% ------------------------------------------------------------------------
% Check that InputPaths exists when using the bypass workflow
% -------------------------------------------------------------------------
if ~exist('InputPaths','var') || ~isstruct(InputPaths)
    error(['InputPaths was not found in the workspace. ' ...
           'Call InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides) ' ...
           'before running input_data_bypass_script.m.']);
end

%% ========================================================================
% SECTION 1 — GENERAL MODEL CONTROLS
% ========================================================================
% This section will later populate legacy variables such as:
%   time_step_model
%   running_control.*
%   Courant_Parameters.*
%   date_begin, date_end
%   slope_outlet
%   GIS_data.alfa_1, alfa_2, beta_1, beta_2
%   LULC_Parameters.River_Manning
%   River_K_coeff
%   ADD, min_Bt, Bmin, Bmax
%   GIS_data.min_area, tau, K_value, sl, slope_DTM
%   topo_path
%   Design_Storm_Parameters.*
%
% RULE
% ------------------------------------------------------------------------
% In standalone bypass mode, these fields should be treated as the source
% of truth for the input_data layer.
% -------------------------------------------------------------------------
InputData_Bypass.general = struct();

% -------------------------------------------------------------------------
% Time control
% -------------------------------------------------------------------------
% IMPORTANT:
% These values are in MINUTES in bypass mode.
%
% Example:
%   5 seconds = 5/60 minutes
%   1 minute  = 1
%   60 min    = 60
%
% Choose values consistent with your numerical setup.
InputData_Bypass.general.time_step_model      = 5/60;   % [min] Base model time step = 5 s
InputData_Bypass.general.min_time_step        = 0.01;      % [sec] Example minimum adaptive step
InputData_Bypass.general.max_time_step        = 3;  % [sec] Maximum allowed adaptive step
InputData_Bypass.general.time_step_increments = 0.001;  % [sec] Legacy increment used internally
InputData_Bypass.general.time_step_change     = 0.0001;   % [-] Legacy adaptive change factor

% -------------------------------------------------------------------------
% Courant / stability controls
% -------------------------------------------------------------------------
InputData_Bypass.general.alfa_min = 0.50;   % [-]
InputData_Bypass.general.alfa_max = 0.50;   % [-]

% -------------------------------------------------------------------------
% Simulation start/end dates
% -------------------------------------------------------------------------
% REQUIRED by updated input_data_script.
InputData_Bypass.general.date_begin = datetime(2025,5,1,0,0,0);
InputData_Bypass.general.date_end   = datetime(2025,5,1,3,0,0);

% Optional direct override of routing time [min]
% If omitted, input_data_script computes:
%   routing_time = minutes(date_end - date_begin)
% InputData_Bypass.general.routing_time = 1440;

% -------------------------------------------------------------------------
% Raster forcing time steps
% -------------------------------------------------------------------------
% These define the spacing between consecutive raster maps in minutes.
% Examples:
%   IMERG half-hourly   -> 30
%   MSWEP 3-hourly      -> 180
%   Daily ETP maps      -> 1440
InputData_Bypass.general.dt_rainfall_maps_min      = 5;   % [min]
InputData_Bypass.general.dt_transpiration_maps_min = 1440;  % [min]
InputData_Bypass.general.dt_evaporation_maps_min   = 1440;  % [min]

% -------------------------------------------------------------------------
% Representative rainfall filename (REQUIRED for fast manifest generation)
% -------------------------------------------------------------------------
% This must match the naming pattern used in your rainfall folder.
% It does NOT need to match the simulation period.
%
% Example formats:
%   IMERG_2025_05_01_00_00.tif
%   IMERG_2025-05-01-00-00.tif
%
% IMPORTANT:
% All rainfall rasters must follow the SAME naming convention.
InputData_Bypass.general.rainfall_filename_example = ...
    'IMERG_30min_mmhr_India_2010_03_12_00_00.tif';

% -------------------------------------------------------------------------
% Outlet / map saving / visualization controls
% -------------------------------------------------------------------------
InputData_Bypass.general.slope_outlet                  = 0.01; % [m/m]
InputData_Bypass.general.n_outlets_data                = 50;     % legacy default usually 1
InputData_Bypass.general.record_time_maps              = 10;    % [min]
InputData_Bypass.general.record_time_hydrographs       = 10;     % [min]
InputData_Bypass.general.Pol_min                       = 0.001; % plotting threshold
InputData_Bypass.general.depth_wse                     = 0.01;  % [m]
InputData_Bypass.general.flag_wse                      = 0;     % 0/1
InputData_Bypass.general.record_time_spatial_rainfall  = 10;    % [min]
InputData_Bypass.general.time_save_ETP                 = 1440;  % [min]
InputData_Bypass.general.record_time_spatial_ETP       = 1440;  % [min] REQUIRED
InputData_Bypass.general.Krs_ETP                       = 0.13;  % legacy constant
InputData_Bypass.general.albedo                        = 0.25;  % [-]

% -------------------------------------------------------------------------
% Routing parameterization
% -------------------------------------------------------------------------
InputData_Bypass.general.alfa_1  = 1.0;
InputData_Bypass.general.alfa_2  = 1.0;
InputData_Bypass.general.beta_1  = 1.0;
InputData_Bypass.general.beta_2  = 1.0;
InputData_Bypass.general.Manning = 0.035; % channel Manning n

% Optional river parameter
% InputData_Bypass.general.River_K_coeff = NaN;

% -------------------------------------------------------------------------
% Water-quality scalar controls
% -------------------------------------------------------------------------
InputData_Bypass.general.ADD    = 5.0;
InputData_Bypass.general.min_Bt = 0.0;
InputData_Bypass.general.Bmin   = 0.0;
InputData_Bypass.general.Bmax   = 1.0;

% -------------------------------------------------------------------------
% DEM smoothing / geomorphic preprocessing controls
% -------------------------------------------------------------------------
InputData_Bypass.general.min_area  = 0.5;      % km2
InputData_Bypass.general.tau       = 0.2;     % between 0 and 1
InputData_Bypass.general.K_value   = 10;      % between 0 and 20
InputData_Bypass.general.sl        = 0.001;   % m/m
InputData_Bypass.general.resolution_resample   = 30; % m
InputData_Bypass.general.slope_DTM = 0.05;    % 

% -------------------------------------------------------------------------
% Tool path used later by input_data_script / preprocessing
% -------------------------------------------------------------------------
InputData_Bypass.general.topo_path = InputPaths.topo_path;


% -------------------------------------------------------------------------
% Internal ETP meteorological forcing file
% -------------------------------------------------------------------------
% Used when:
%   flags.flag_ETP = 1
%   flags.flag_input_ETP_map = 0
%
% The path is resolved upstream by input_paths_bypass.m so that this
% bypass-data script does not need to hardcode machine- or case-specific
% file locations.
if isfield(InputPaths,'ETP_input_spreadsheet') && ~isempty(InputPaths.ETP_input_spreadsheet)
    InputData_Bypass.general.etp_input_spreadsheet = InputPaths.ETP_input_spreadsheet;
else
    InputData_Bypass.general.etp_input_spreadsheet = ...
        fullfile('Forcing','Evapotranspiration','ETP_input_data.xlsx');
end

% -------------------------------------------------------------------------
% Design storm parameters
% -------------------------------------------------------------------------
% These are used when:
%   flags.flag_alternated_blocks = 1
%   or
%   flags.flag_huff = 1
%
% IMPORTANT:
% The updated input_data_script reads these directly from "general".
InputData_Bypass.general.RP                = 10;    % [years]
InputData_Bypass.general.Rainfall_Duration = 180;   % [min]
InputData_Bypass.general.K                 = 1000;  % IDF coefficient
InputData_Bypass.general.a                 = 0.20;  % IDF parameter
InputData_Bypass.general.b                 = 10;    % IDF parameter
InputData_Bypass.general.c                 = 0.80;  % IDF parameter
InputData_Bypass.general.dt_design         = 5;     % [min]

%% ========================================================================
% SECTION 2 — MODEL FLAGS
% ========================================================================
% This section will later populate the legacy struct:
%   flags
%
% IMPORTANT
% ------------------------------------------------------------------------
% These names must match the flag names expected by HydroPol2D.
% Use 0 or 1 unless a flag explicitly supports additional coded values.
% -------------------------------------------------------------------------
InputData_Bypass.flags = struct();

% -------------------------------------------------------------------------
% Forcing and boundary flags
% -------------------------------------------------------------------------
InputData_Bypass.flags.flag_rainfall                     = 1;
InputData_Bypass.flags.flag_spatial_rainfall             = 0;
InputData_Bypass.flags.flag_ETP                          = 0;
InputData_Bypass.flags.flag_input_rainfall_map           = 0;
InputData_Bypass.flags.flag_rainfall_multiple_runs       = 0;
InputData_Bypass.flags.flag_data_source                  = 0;
InputData_Bypass.flags.flag_inflow                       = 0;
InputData_Bypass.flags.flag_satellite_rainfall           = 0;
InputData_Bypass.flags.flag_alternated_blocks            = 0;
InputData_Bypass.flags.flag_huff                         = 0;
InputData_Bypass.flags.flag_stage_hydrograph             = 0;
InputData_Bypass.flags.flag_input_ETP_map                = 0;

% -------------------------------------------------------------------------
% Numerical / hydrologic core flags
% -------------------------------------------------------------------------
InputData_Bypass.flags.flag_timestep                     = 1;
InputData_Bypass.flags.flag_infiltration                 = 0;
InputData_Bypass.flags.flag_critical                     = 0;
InputData_Bypass.flags.flag_D8                           = 0;
InputData_Bypass.flags.flag_CA                           = 0;
InputData_Bypass.flags.flag_inertial                     = 1;
InputData_Bypass.flags.flag_waterbalance                 = 0;
InputData_Bypass.flags.flag_waterquality                 = 0;
InputData_Bypass.flags.flag_reservoir                    = 0;
InputData_Bypass.flags.flag_wq_model                     = 0;
InputData_Bypass.flags.flag_groundwater_modeling         = 0;
InputData_Bypass.flags.flag_real_time_satellite_rainfall = 0;
InputData_Bypass.flags.flag_dam_break                    = 0;
InputData_Bypass.flags.flag_human_instability            = 0; % can be 0, 1, or 3
InputData_Bypass.flags.flag_boundary                     = 0;
InputData_Bypass.flags.flag_numerical_scheme             = 1;
InputData_Bypass.flags.flag_outlet_type                  = 1;
InputData_Bypass.flags.flag_adaptive_timestepping        = 1;
InputData_Bypass.flags.flag_neglect_infiltration_river   = 0;
InputData_Bypass.flags.flag_subgrid                      = 0;
InputData_Bypass.flags.flag_spatial_albedo               = 0;
InputData_Bypass.flags.flag_river_rasters                = 0;
InputData_Bypass.flags.flag_baseflow                     = 0;
InputData_Bypass.flags.flag_kinematic                    = 0;
InputData_Bypass.flags.flag_diffusive                    = 0;
InputData_Bypass.flags.flag_DTM                          = 0;
InputData_Bypass.flags.flag_abstraction                  = 0;
InputData_Bypass.flags.flag_overbanks                    = 0;
InputData_Bypass.flags.flag_snow_modeling                = 0;
InputData_Bypass.flags.flag_WQ_Rasters                   = 0;
InputData_Bypass.flags.flag_GPU                          = 0;
InputData_Bypass.flags.flag_single                       = 0;
InputData_Bypass.flags.flag_warmup                       = 0;
InputData_Bypass.flags.flag_initial_buildup              = 0;
InputData_Bypass.flags.flag_resample                     = 0;
InputData_Bypass.flags.flag_smoothening                  = 0;
InputData_Bypass.flags.flag_trunk                        = 0;
InputData_Bypass.flags.flag_fill_DEM                     = 0;
InputData_Bypass.flags.flag_smooth_cells                 = 0;
InputData_Bypass.flags.flag_reduce_DEM                   = 0;
InputData_Bypass.flags.flag_export_maps                  = 1;
InputData_Bypass.flags.flag_river_heigth_compensation    = 0;
InputData_Bypass.flags.flag_dashboard                    = 0;
InputData_Bypass.flags.flag_elapsed_time                 = 1;
InputData_Bypass.flags.flag_obs_gauges                   = 0;

%% ========================================================================
% SECTION 3 — HUMAN INSTABILITY PARAMETERS
% ========================================================================
% This section will later populate the legacy variable:
%   Human_Instability
%
% Use ONLY if flags.flag_human_instability is active.
%
% CASE A — flags.flag_human_instability = 1
% ------------------------------------------------------------------------
% Simple person instability parameterization.
%
% CASE B — flags.flag_human_instability = 3
% ------------------------------------------------------------------------
% Detailed demographic groups. If you use that option, define ALL required
% fields expected by the legacy code.
% -------------------------------------------------------------------------
InputData_Bypass.Human_Instability = struct();
InputData_Bypass.Human_Instability.mu            = 0.6;
InputData_Bypass.Human_Instability.Cd            = 1.0;
InputData_Bypass.Human_Instability.ro_person     = 985;   % [kg/m^3]
InputData_Bypass.Human_Instability.weight_person = 75;    % [kg]
InputData_Bypass.Human_Instability.height_person = 1.75;  % [m]
InputData_Bypass.Human_Instability.width1_person = 0.35;  % [m]
InputData_Bypass.Human_Instability.width2_person = 0.20;  % [m]
InputData_Bypass.Human_Instability.ro_water      = 1000;  % [kg/m^3]
InputData_Bypass.Human_Instability.gravity       = 9.81;  % [m/s^2]

%% ========================================================================
% SECTION 4 — LULC CLASS TABLE
% ========================================================================
% PURPOSE
%   Defines all land-cover classes and their hydraulic / hydrologic
%   parameters directly in MATLAB, replacing LULC_parameters.xlsx when
%   bypass mode is active.
%
% HOW THIS MAPS TO THE ORIGINAL MODEL
%   Later, input_data_script will translate this block into the same legacy
%   variables used by HydroPol2D:
%
%       LULC_name
%       LULC_index
%       lulc_parameters
%       n_lulc
%       imp_index
%
% COLUMN MEANING (must preserve this order)
%   1) LC                  : class name
%   2) Index               : raster code used in the LULC map
%   3) roughness           : Manning n [m^(-1/3) s]
%   4) h_0_mm              : interception / storage-related parameter [mm]
%   5) d_0_mm              : depression storage / related parameter [mm]
%   6) C1                  : empirical parameter 1
%   7) C2                  : empirical parameter 2
%   8) C3                  : empirical parameter 3
%   9) C4                  : empirical parameter 4
%  10) index_impervious    : impervious class ID
% ========================================================================

LULC_table = table();

LULC_table.LC = { ...
    'left'; ...
    'center'; ...
    'right'};

LULC_table.Index = [1;2;3];

LULC_table.roughness = [ ...
    0.015; ...
    1e-4; ...
    1e-4];

LULC_table.h_0_mm = zeros(3,1);
LULC_table.d_0_mm = zeros(3,1);

LULC_table.C1 = [10;10;10];
LULC_table.C2 = [0.2;0.20;0.20];
LULC_table.C3 = [800;800;800];
LULC_table.C4 = [1.20;1.20;1.20];

% IMPORTANT:
% The updated input_data_script extracts:
%   imp_index = parameters(1,end)
% so the legacy behavior expects the impervious ID to be stored in the LAST
% column, and in the FIRST retained row.
%
% To remain consistent with that legacy extraction:
%   - place the impervious class ID in row 1, last column
%   - other rows may be NaN
%
% If you later improve input_data_script, this can be made cleaner.
LULC_table.index_impervious = [ ...
    4; ...
    NaN; ...
    NaN; ...
    ];

InputData_Bypass.LULC.table = LULC_table;

%% ========================================================================
% SECTION 5 — SOIL CLASS TABLE
% ========================================================================
% PURPOSE
%   Defines all soil classes and hydraulic parameters directly in MATLAB,
%   replacing SOIL_parameters.xlsx when bypass mode is active.
%
% HOW THIS MAPS TO THE ORIGINAL MODEL
%   Later, input_data_script will translate this block into the same legacy
%   variables used by HydroPol2D:
%
%       SOIL_name
%       SOIL_index
%       soil_parameters
%       n_soil
%
% COLUMN MEANING (must preserve this order)
%   1) Soil_type      : class name
%   2) Index          : raster code used in the SOIL map
%   3) ksat_mm_h      : saturated hydraulic conductivity [mm/h]
%   4) n_vg           : van Genuchten n [-]
%   5) alpha_vg_1_m   : van Genuchten alpha [1/m]
%   6) theta_sat      : saturated volumetric water content [cm3/cm3]
%   7) theta_r        : residual volumetric water content [cm3/cm3]
%   8) theta_i        : initial volumetric water content [cm3/cm3]
%   9) Sy             : specific yield [-]
%  10) ksat_gw_mm_h   : groundwater saturated conductivity [mm/h]
%  11) Soil_Depth_m   : effective soil depth / DTB [m]
% ========================================================================

SOIL_table = table();

SOIL_table.Soil_type = { ...
    'left'; ...
    'right'; ...
    'center';};

SOIL_table.Index = [1;2;3];

SOIL_table.ksat_mm_h = [0.3;0.3;0.3];
SOIL_table.n_vg = [1.09;1.09;1.09];
SOIL_table.alpha_vg_1_m = [0.8;0.8;0.8];
SOIL_table.theta_sat = [0.385;0.385;0.385];
SOIL_table.theta_r   = [0.068;0.068;0.068];
SOIL_table.theta_i   = [0.0997;0.0997;0.0997];
SOIL_table.Sy        = [0.0997;0.0997;0.0997];
SOIL_table.ksat_gw_mm_h = [6;6;6];
SOIL_table.Soil_Depth_m = ones(3,1);

InputData_Bypass.SOIL.table = SOIL_table;

%% ========================================================================
% SECTION 6 — SNOW PARAMETERS
% ========================================================================
% PURPOSE
%   Defines snow-model parameters directly in MATLAB.
%
% HOW THIS MAPS TO THE ORIGINAL MODEL
%   Later, input_data_script / preprocessing can read this section and
%   populate the legacy Snow_Properties structure expected by the model.
%
% IMPORTANT
%   The raster-like state variables below depend on DEM_raster already
%   existing in the workspace when this script is executed.
% ========================================================================

Snow_Properties = struct();

Snow_Properties.alpha           = 0.8;
Snow_Properties.epsilon         = 0.98;
Snow_Properties.C_e             = 0.001;
Snow_Properties.DDF             = 2;      % [mm/°C/day]
Snow_Properties.T_thresh        = 0;      % [°C]
Snow_Properties.rho_snow_init   = 100;    % [kg/m^3]
Snow_Properties.rho_max         = 400;    % [kg/m^3]
Snow_Properties.k_t             = 0.1;
Snow_Properties.k_swe           = 0.001;
Snow_Properties.k_D             = 0.02;
Snow_Properties.snow_fraction_a = 0.2;

InputData_Bypass.Snow_Properties = Snow_Properties;

%% ========================================================================
% SECTION 7 — LUMPED RAINFALL TIME SERIES
% ========================================================================
% This section will later populate the legacy struct:
%   Rainfall_Parameters
%
% Use this section when:
%   flags.flag_rainfall = 1
% and NONE of the following are active:
%   flag_input_rainfall_map
%   flag_spatial_rainfall
%   flag_satellite_rainfall
%   flag_real_time_satellite_rainfall
%   flag_alternated_blocks
%   flag_huff
% ========================================================================

Rainfall_Parameters = struct();

% -------------------------------------------------------------------------
% Lumped / spatially invariant rainfall forcing
% -------------------------------------------------------------------------
% Used when:
%   flags.flag_rainfall = 1
%   flags.flag_spatial_rainfall = 0
%   flags.flag_input_rainfall_map = 0
%   flags.flag_satellite_rainfall = 0
%   flags.flag_real_time_satellite_rainfall = 0
%   flags.flag_alternated_blocks = 0
%   flags.flag_huff = 0
%
% In bypass mode, the preferred source is a time-series file resolved by
% input_paths_bypass.m:
%
%   InputPaths.Rainfall_Timeseries_File
%
% Expected file format:
%   Column 1 -> time in minutes from simulation start
%   Column 2 -> rainfall intensity in mm/h
%
% Accepted file types:
%   .xlsx, .xls, .csv
%
% If the file does not exist, the script falls back to the simple example
% series below.
% -------------------------------------------------------------------------
rainfall_file = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Stanford/Forcing/Rainfall/Rainfall_Intensity_Data.xlsx';

if isfield(InputPaths,'Rainfall_Timeseries_File') && ~isempty(InputPaths.Rainfall_Timeseries_File)
    rainfall_file = InputPaths.Rainfall_Timeseries_File;
end

if ~isempty(rainfall_file) && exist(rainfall_file,'file') == 2

    T_rain = readtable(rainfall_file);

    if width(T_rain) < 2
        error(['Rainfall time-series file must contain at least 2 columns: ' ...
               '[time_min, intensity_mm_h]. File: %s'], rainfall_file);
    end

    Rainfall_Parameters.time_rainfall      = table2array(T_rain(:,1));
    Rainfall_Parameters.intensity_rainfall = table2array(T_rain(:,2));

    Rainfall_Parameters.time_rainfall      = double(Rainfall_Parameters.time_rainfall(:));
    Rainfall_Parameters.intensity_rainfall = double(Rainfall_Parameters.intensity_rainfall(:));

    valid = isfinite(Rainfall_Parameters.time_rainfall) & ...
            isfinite(Rainfall_Parameters.intensity_rainfall) & ...
            (Rainfall_Parameters.intensity_rainfall >= 0);

    Rainfall_Parameters.time_rainfall      = Rainfall_Parameters.time_rainfall(valid);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(valid);

    if numel(Rainfall_Parameters.time_rainfall) < 2
        error('Rainfall time-series file must contain at least 2 valid rows.');
    end

    dt_rain = Rainfall_Parameters.time_rainfall(2) - Rainfall_Parameters.time_rainfall(1);

    if dt_rain <= 0
        error('Rainfall time-series file must have a positive time step.');
    end

    if any(abs(diff(Rainfall_Parameters.time_rainfall) - dt_rain) > 1e-9)
        error('Rainfall time-series file must have a constant time step.');
    end

    Rainfall_Parameters.time_step_rainfall = dt_rain;
    Rainfall_Parameters.rainfall_duration  = Rainfall_Parameters.time_rainfall(end);
    Rainfall_Parameters.n_obs_rainfall     = numel(Rainfall_Parameters.time_rainfall);

else
    % ------------------------------------------------------------
    % Fallback example (used only if no file is found)
    % ------------------------------------------------------------
    Rainfall_Parameters.time_rainfall      = [0; 240]; % [min]
    Rainfall_Parameters.intensity_rainfall = [100;100];      % [mm/h]
    Rainfall_Parameters.time_step_rainfall = 240;                      % [min]
    Rainfall_Parameters.rainfall_duration  = 240;                     % [min]
    Rainfall_Parameters.n_obs_rainfall     = numel(Rainfall_Parameters.time_rainfall);
end

InputData_Bypass.Rainfall_Parameters = Rainfall_Parameters;

%% ========================================================================
% SECTION 8 — SPATIAL RAINFALL MAP SCHEDULE
% ========================================================================
% FAST RAINFALL FILE DISCOVERY WITHOUT SCANNING THE WHOLE ARCHIVE
%
% STRATEGY
% ------------------------------------------------------------------------
% Instead of scanning a folder containing many years of rasters, we:
%
%   1) read the rainfall folder path from InputPaths
%   2) look at one sample filename pattern
%   3) use:
%        - date_begin
%        - date_end
%        - dt_rainfall_maps_min
%      to generate the exact sequence of expected timestamps
%   4) build the expected filenames directly
%   5) keep only the files that exist
%
% This is much faster than scanning very large folders on network storage.
%
% REQUIREMENT
% ------------------------------------------------------------------------
% All rainfall rasters must follow ONE consistent filename format.
%
% Example supported formats:
%   IMERG_2026_01_01_00_00.tif
%   IMERG_2026-01-01-00-00.tif
%   Rain_2026_01_01_00_00.tif
%
% The code below infers the filename template from one example file stored
% in InputData_Bypass.general.rainfall_filename_example.
% -------------------------------------------------------------------------

Input_Rainfall = struct();

rain_folder = InputPaths.Rainfall_Rasters_Folder;
date_begin_sim = InputData_Bypass.general.date_begin;
date_end_sim   = InputData_Bypass.general.date_end;
dt_rain_min    = InputData_Bypass.general.dt_rainfall_maps_min;

% -------------------------------------------------------------------------
% USER-DEFINED EXAMPLE FILENAME
% -------------------------------------------------------------------------
% This should match the naming style of the rainfall rasters in the folder.
% Example:
%   'IMERG_2026_01_01_00_00.tif'
%
% You only need ONE representative example. The code will detect the
% timestamp token inside it and generate the full manifest automatically.
if isfield(InputData_Bypass.general,'rainfall_filename_example') && ...
        ~isempty(InputData_Bypass.general.rainfall_filename_example)

    rainfall_filename_example = char(InputData_Bypass.general.rainfall_filename_example);

else
    error(['Please define InputData_Bypass.general.rainfall_filename_example ' ...
           'with one representative rainfall filename.']);
end

% -------------------------------------------------------------------------
% Build the expected simulation timestamps
% -------------------------------------------------------------------------
rain_times_expected = (date_begin_sim:minutes(dt_rain_min):date_end_sim)';
n_expected = numel(rain_times_expected);

if n_expected < 1
    error('No expected rainfall timestamps were generated. Check date_begin/date_end.');
end

% -------------------------------------------------------------------------
% Infer filename template from the example
% -------------------------------------------------------------------------
rain_template = infer_rainfall_filename_template(rainfall_filename_example);

% -------------------------------------------------------------------------
% Preallocate manifest
% -------------------------------------------------------------------------
rain_files = cell(n_expected,1);
keep_file  = false(n_expected,1);

for i = 1:n_expected
    this_time = rain_times_expected(i);

    this_name = build_rainfall_filename_from_template(rain_template, this_time);
    this_path = fullfile(rain_folder, this_name);

    if exist(this_path,'file') == 2
        rain_files{i} = this_path;
        keep_file(i)  = true;
    end
end

% Keep only files that exist
rain_files = rain_files(keep_file);
rain_times = rain_times_expected(keep_file);

n_rain = numel(rain_files);

if n_rain > 0
    Input_Rainfall.time = minutes(rain_times - rain_times(1));
    Input_Rainfall.labels_Directory = cellfun(@char, rain_files, 'UniformOutput', false);
    Input_Rainfall.num_obs_maps = n_rain;

    fprintf('[Rainfall manifest] %d file(s) found between %s and %s.\n', ...
        n_rain, char(string(date_begin_sim)), char(string(date_end_sim)));

    fprintf('[Rainfall manifest] First file: %s\n', Input_Rainfall.labels_Directory{1});
    fprintf('[Rainfall manifest] Last file : %s\n', Input_Rainfall.labels_Directory{end});

else
    Input_Rainfall.time = [];
    Input_Rainfall.labels_Directory = {};
    Input_Rainfall.num_obs_maps = 0;

    warning(['No rainfall rasters were found for the requested simulation period.\n' ...
             'Folder: %s\nExample pattern: %s'], ...
             rain_folder, rainfall_filename_example);
end

InputData_Bypass.Input_Rainfall = Input_Rainfall;

%% ========================================================================
% SECTION 9 — ETP / TRANSPIRATION / EVAPORATION MAP SCHEDULES
% ========================================================================
% Uses raster files discovered and sorted in input_paths_bypass_script.m

Input_Transpiration = struct();
Input_Evaporation   = struct();

tr_files = InputPaths.Transpiration_Raster_Files;
ev_files = InputPaths.Evaporation_Raster_Files;

n_tr = numel(tr_files);
n_ev = numel(ev_files);

if n_tr > 0
    dt_tr = InputData_Bypass.general.dt_transpiration_maps_min;

    Input_Transpiration.time = (0:n_tr-1)' * dt_tr;
    Input_Transpiration.labels_Directory = cellfun(@char, tr_files, 'UniformOutput', false);
    Input_Transpiration.num_obs_maps = n_tr;
else
    Input_Transpiration.time = [];
    Input_Transpiration.labels_Directory = {};
    Input_Transpiration.num_obs_maps = 0;
end

if n_ev > 0
    dt_ev = InputData_Bypass.general.dt_evaporation_maps_min;

    Input_Evaporation.time = (0:n_ev-1)' * dt_ev;
    Input_Evaporation.labels_Directory = cellfun(@char, ev_files, 'UniformOutput', false);
    Input_Evaporation.num_obs_maps = n_ev;
else
    Input_Evaporation.time = [];
    Input_Evaporation.labels_Directory = {};
    Input_Evaporation.num_obs_maps = 0;
end

InputData_Bypass.Input_Transpiration = Input_Transpiration;
InputData_Bypass.Input_Evaporation   = Input_Evaporation;

%% ========================================================================
% SECTION 10 — INFLOW HYDROGRAPHS
% ========================================================================
% Reads inflow forcing from CSV defined in InputPaths.Inflow_Hydrograph_CSV
%
% REQUIRED CSV FORMAT
% -------------------------------------------------------------------------
% One row per time step.
%
% Required first column:
%   time_min
%
% Then, for each inflow gauge i:
%   Q_i     -> discharge [m^3/s]
%   x_i     -> easting [m]
%   y_i     -> northing [m]
%
% Example header:
%   time_min,Q_1,x_1,y_1,Q_2,x_2,y_2
% -------------------------------------------------------------------------

Inflow = struct();

if isfield(InputPaths,'Inflow_Hydrograph_CSV') && isfile(InputPaths.Inflow_Hydrograph_CSV)

    T_in = readtable(InputPaths.Inflow_Hydrograph_CSV);

    if ~ismember('time_min', T_in.Properties.VariableNames)
        error('Inflow hydrograph CSV must contain a column named "time_min".');
    end

    Inflow.time_inflow = T_in.time_min;

    q_cols = startsWith(T_in.Properties.VariableNames, 'Q_');
    x_cols = startsWith(T_in.Properties.VariableNames, 'x_');
    y_cols = startsWith(T_in.Properties.VariableNames, 'y_');

    Q_names = T_in.Properties.VariableNames(q_cols);
    X_names = T_in.Properties.VariableNames(x_cols);
    Y_names = T_in.Properties.VariableNames(y_cols);

    nQ = numel(Q_names);
    nX = numel(X_names);
    nY = numel(Y_names);

    if ~(nQ == nX && nQ == nY)
        error('Inflow CSV must have matching numbers of Q_i, x_i, and y_i columns.');
    end

    if nQ == 0
        error('Inflow CSV contains no inflow gauges. Add at least Q_1, x_1, y_1.');
    end

    Inflow.inflow_discharge = table2array(T_in(:, Q_names));

    Inflow.x_coords = cell(nQ,1);
    Inflow.y_coords = cell(nQ,1);

    for i = 1:nQ
        xi = table2array(T_in(1, X_names(i)));
        yi = table2array(T_in(1, Y_names(i)));

        Inflow.x_coords{i} = xi;
        Inflow.y_coords{i} = yi;
    end

else
    Inflow.time_inflow = [];
    Inflow.inflow_discharge = [];
    Inflow.x_coords = {};
    Inflow.y_coords = {};
end

InputData_Bypass.Inflow = Inflow;

%% ========================================================================
% SECTION 11 — STAGE HYDROGRAPHS
% ========================================================================
% Reads stage forcing from CSV defined in InputPaths.Stage_Hydrograph_CSV
%
% REQUIRED CSV FORMAT
% -------------------------------------------------------------------------
% One row per time step.
%
% Required first column:
%   time_min
%
% Then, for each stage gauge i:
%   H_i     -> stage [m]
%   x_i     -> easting [m]
%   y_i     -> northing [m]
%
% Example header:
%   time_min,H_1,x_1,y_1,H_2,x_2,y_2
% -------------------------------------------------------------------------

Stage = struct();

if isfield(InputPaths,'Stage_Hydrograph_CSV') && isfile(InputPaths.Stage_Hydrograph_CSV)

    T_st = readtable(InputPaths.Stage_Hydrograph_CSV);

    if ~ismember('time_min', T_st.Properties.VariableNames)
        error('Stage hydrograph CSV must contain a column named "time_min".');
    end

    Stage.time_stage = T_st.time_min;

    h_cols = startsWith(T_st.Properties.VariableNames, 'H_');
    x_cols = startsWith(T_st.Properties.VariableNames, 'x_');
    y_cols = startsWith(T_st.Properties.VariableNames, 'y_');

    H_names = T_st.Properties.VariableNames(h_cols);
    X_names = T_st.Properties.VariableNames(x_cols);
    Y_names = T_st.Properties.VariableNames(y_cols);

    nH = numel(H_names);
    nX = numel(X_names);
    nY = numel(Y_names);

    if ~(nH == nX && nH == nY)
        error('Stage CSV must have matching numbers of H_i, x_i, and y_i columns.');
    end

    if nH == 0
        error('Stage CSV contains no stage gauges. Add at least H_1, x_1, y_1.');
    end

    Stage.stage = table2array(T_st(:, H_names));

    Stage.x_coords = cell(nH,1);
    Stage.y_coords = cell(nH,1);

    for i = 1:nH
        xi = table2array(T_st(1, X_names(i)));
        yi = table2array(T_st(1, Y_names(i)));

        Stage.x_coords{i} = xi;
        Stage.y_coords{i} = yi;
    end

else
    Stage.time_stage = [];
    Stage.stage = [];
    Stage.x_coords = {};
    Stage.y_coords = {};
end

InputData_Bypass.Stage = Stage;

%% ========================================================================
% SECTION 12 — RESERVOIR / CONTROL STRUCTURE DATA
% ========================================================================
% PURPOSE
% -------------------------------------------------------------------------
% This section defines hydraulic control structures such as:
%   • reservoirs
%   • detention basins
%   • outlet structures
%   • spillways
%   • orifices
%
% It will later populate:
%   Reservoir_Data
%
% IMPORTANT:
% All coordinate fields must use the SAME projected CRS and units as the
% DEM.
% ========================================================================

Reservoir_Data = struct();

% Unique identifier
Reservoir_Data.index = [1];

% Upstream control point
Reservoir_Data.x_us = [];   % [m]
Reservoir_Data.y_us = [];  % [m]

% Primary outlet
Reservoir_Data.k1 = [];
Reservoir_Data.h1 = [];          % [m]
Reservoir_Data.k2 = [];
Reservoir_Data.x_ds1 = [];  % [m]
Reservoir_Data.y_ds1 = []; % [m]

% Secondary outlet / spillway
Reservoir_Data.k3 = [];
Reservoir_Data.h2 = [];          % [m]
Reservoir_Data.k4 = [];
Reservoir_Data.x_ds2 = [];  % [m]
Reservoir_Data.y_ds2 = []; % [m]

% % Upstream control point
% Reservoir_Data.x_us = [500100.0];   % [m]
% Reservoir_Data.y_us = [4099900.0];  % [m]
% 
% % Primary outlet
% Reservoir_Data.k1 = [1.0];
% Reservoir_Data.h1 = [0.5];          % [m]
% Reservoir_Data.k2 = [1.0];
% Reservoir_Data.x_ds1 = [500150.0];  % [m]
% Reservoir_Data.y_ds1 = [4099850.0]; % [m]
% 
% % Secondary outlet / spillway
% Reservoir_Data.k3 = [1.0];
% Reservoir_Data.h2 = [1.0];          % [m]
% Reservoir_Data.k4 = [1.0];
% Reservoir_Data.x_ds2 = [500200.0];  % [m]
% Reservoir_Data.y_ds2 = [4099800.0]; % [m]

InputData_Bypass.Reservoir_Data = Reservoir_Data;

%% ========================================================================
% SECTION 13 — OBSERVED GAUGES
% ========================================================================
% Reads observed gauges from the CSV defined in:
%   InputPaths.Observed_Gauges_CSV
%
% REQUIRED CSV FORMAT
% -------------------------------------------------------------------------
% Required columns:
%   Gauge
%   Easting_m
%   Northing_m
%   Label_Name
%
% Example:
%   Gauge,Easting_m,Northing_m,Label_Name
%   1,500000,4100000,Gauge_A
%   2,500500,4099500,Gauge_B
%   3,501000,4099000,Gauge_C
% ========================================================================

Obs = table();

if isfield(InputPaths,'Observed_Gauges_CSV') && isfile(InputPaths.Observed_Gauges_CSV)

    T_obs = readtable(InputPaths.Observed_Gauges_CSV);

    requiredObsCols = {'Gauge','Easting_m','Northing_m','Label_Name'};
    missingObsCols = requiredObsCols(~ismember(requiredObsCols, T_obs.Properties.VariableNames));

    if ~isempty(missingObsCols)
        error('Observed gauges CSV is missing required column(s): %s', strjoin(missingObsCols, ', '));
    end

    Obs.("Gauge") = T_obs.Gauge;
    Obs.("Easting (m)") = T_obs.Easting_m;
    Obs.("Northing (m)") = T_obs.Northing_m;
    Obs.("Label Name") = string(T_obs.Label_Name);

else
    Obs.("Gauge") = zeros(0,1);
    Obs.("Easting (m)") = zeros(0,1);
    Obs.("Northing (m)") = zeros(0,1);
    Obs.("Label Name") = strings(0,1);
end

InputData_Bypass.Obs = Obs;
InputData_Bypass.Obs = Obs;

%% ========================================================================
% FINAL NOTES FOR USERS
% ========================================================================
% 1) This script intentionally defines EVERYTHING in one place.
% 2) You may keep unused sections present; later logic can ignore them based
%    on the flags.
% 3) The updated input_data_script should validate consistency, for example:
%      - do not activate inflow and stage simultaneously unless supported
%      - if flag_input_rainfall_map = 1, Input_Rainfall must be valid
%      - if flag_input_ETP_map = 1, both Input_Transpiration and
%        Input_Evaporation must be valid
% 4) The most important thing is to keep IDs and units consistent with your
%    rasters and legacy HydroPol2D conventions.
% ========================================================================


function dt = parse_datetime_from_filename_local(fname)
%PARSE_DATETIME_FROM_FILENAME_LOCAL
% Attempts to extract a datetime from filenames like:
%   rain_2025_06_01_00_00.tif
%   GPM_180min_2025_06_01_03_00.tif
%   2025-06-01-03-00.tif
%   2025_06_01.tif
%
% Returns NaT if parsing fails.

    dt = NaT;

    if ~(ischar(fname) || isstring(fname))
        return;
    end

    fname = char(fname);

    % Pattern 1: YYYY_MM_DD_HH_mm or YYYY-MM-DD-HH-mm
    token = regexp(fname, ...
        '(\d{4})[_-](\d{2})[_-](\d{2})[_-](\d{2})[_-](\d{2})', ...
        'tokens', 'once');

    if isempty(token)
        % Pattern 2: YYYY_MM_DD or YYYY-MM-DD
        token = regexp(fname, ...
            '(\d{4})[_-](\d{2})[_-](\d{2})', ...
            'tokens', 'once');

        if isempty(token)
            return;
        else
            yyyy = str2double(token{1});
            mm   = str2double(token{2});
            dd   = str2double(token{3});
            HH   = 0;
            MM   = 0;
        end
    else
        yyyy = str2double(token{1});
        mm   = str2double(token{2});
        dd   = str2double(token{3});
        HH   = str2double(token{4});
        MM   = str2double(token{5});
    end

    try
        dt = datetime(yyyy, mm, dd, HH, MM, 0);
    catch
        dt = NaT;
    end
end

function tpl = infer_rainfall_filename_template(example_name)
%INFER_RAINFALL_FILENAME_TEMPLATE
% Detects the datetime token format inside one example filename and stores
% the pieces needed to rebuild filenames for arbitrary datetimes.

    tpl = struct();
    tpl.example_name = char(example_name);

    % Pattern 1: YYYY_MM_DD_HH_mm
    expr1 = '(\d{4}_\d{2}_\d{2}_\d{2}_\d{2})';
    token1 = regexp(example_name, expr1, 'match', 'once');

    if ~isempty(token1)
        split_parts = regexp(example_name, expr1, 'split', 'once');
        tpl.prefix = split_parts{1};
        tpl.suffix = split_parts{2};
        tpl.mode   = 'underscore_datetime';
        return;
    end

    % Pattern 2: YYYY-MM-DD-HH-mm
    expr2 = '(\d{4}-\d{2}-\d{2}-\d{2}-\d{2})';
    token2 = regexp(example_name, expr2, 'match', 'once');

    if ~isempty(token2)
        split_parts = regexp(example_name, expr2, 'split', 'once');
        tpl.prefix = split_parts{1};
        tpl.suffix = split_parts{2};
        tpl.mode   = 'dash_datetime';
        return;
    end

    % Pattern 3: YYYY_MM_DD
    expr3 = '(\d{4}_\d{2}_\d{2})';
    token3 = regexp(example_name, expr3, 'match', 'once');

    if ~isempty(token3)
        split_parts = regexp(example_name, expr3, 'split', 'once');
        tpl.prefix = split_parts{1};
        tpl.suffix = split_parts{2};
        tpl.mode   = 'underscore_date';
        return;
    end

    % Pattern 4: YYYY-MM-DD
    expr4 = '(\d{4}-\d{2}-\d{2})';
    token4 = regexp(example_name, expr4, 'match', 'once');

    if ~isempty(token4)
        split_parts = regexp(example_name, expr4, 'split', 'once');
        tpl.prefix = split_parts{1};
        tpl.suffix = split_parts{2};
        tpl.mode   = 'dash_date';
        return;
    end

    error(['Could not infer rainfall filename template from example:\n  %s\n' ...
           'Supported formats include YYYY_MM_DD_HH_mm and YYYY-MM-DD-HH-mm.'], ...
           example_name);
end

function fname = build_rainfall_filename_from_template(tpl, dt)
%BUILD_RAINFALL_FILENAME_FROM_TEMPLATE
% Reconstructs a rainfall filename from the inferred template and datetime.

    switch tpl.mode
        case 'underscore_datetime'
            token = sprintf('%04d_%02d_%02d_%02d_%02d', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt));

        case 'dash_datetime'
            token = sprintf('%04d-%02d-%02d-%02d-%02d', ...
                year(dt), month(dt), day(dt), hour(dt), minute(dt));

        case 'underscore_date'
            token = sprintf('%04d_%02d_%02d', ...
                year(dt), month(dt), day(dt));

        case 'dash_date'
            token = sprintf('%04d-%02d-%02d', ...
                year(dt), month(dt), day(dt));

        otherwise
            error('Unsupported rainfall filename template mode: %s', tpl.mode);
    end

    fname = [tpl.prefix token tpl.suffix];
end