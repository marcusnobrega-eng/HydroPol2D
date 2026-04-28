%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%
% HydroPol2D input-data script with two modes:
%
%   use_inputdata_bypass = 0  -> ORIGINAL Excel/spreadsheet workflow
%   use_inputdata_bypass = 1  -> STANDALONE MATLAB bypass workflow
%
% IMPORTANT
% -------------------------------------------------------------------------
% 1) Excel mode is preserved as the legacy path.
% 2) In bypass mode, this script does NOT read the auxiliary spreadsheet
%    files handled here (LULC_parameters.xlsx, SOIL_parameters.xlsx,
%    Rainfall_Intensity_Data.xlsx, Inflow_Hydrograph.xlsx,
%    Stage_Hydrograph.xlsx, BC_Control_structure.xlsx, human_risk.xlsx).
% 3) In bypass mode, the file input_data_bypass_script.m must define the
%    full struct InputData_Bypass with every field needed by the selected
%    model options.
% 4) Preprocessing still reads General_Data / Flags / raster paths before
%    calling this script. Therefore, this bypass is standalone for the
%    INPUT_DATA layer, not for the whole preprocessing pipeline.
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Detect mode and load bypass file if requested
% -------------------------------------------------------------------------
if ~exist('use_inputdata_bypass','var') || isempty(use_inputdata_bypass)
    use_inputdata_bypass = 0;
end

if use_inputdata_bypass == 1
    if ~exist('input_data_bypass_script_path','var') || isempty(input_data_bypass_script_path)
        error(['use_inputdata_bypass = 1, but input_data_bypass_script_path was not ' ...
               'provided by the main wrapper.']);
    end

    if exist(input_data_bypass_script_path,'file') ~= 2
        error(['use_inputdata_bypass = 1, but the bypass input script was not found:\n  %s'], ...
              input_data_bypass_script_path);
    end

    clear InputData_Bypass
    run(input_data_bypass_script_path);

    if ~exist('InputData_Bypass','var') || ~isstruct(InputData_Bypass)
        error('The bypass script must create a struct named InputData_Bypass.');
    end

    fprintf('input_data_script: MATLAB standalone bypass mode is ACTIVE.\n');

    %% ====================================================================
    % STANDALONE BYPASS MODE
    % ====================================================================
    requiredTop = {'general','flags','LULC','SOIL'};
    for ireq = 1:numel(requiredTop)
        if ~isfield(InputData_Bypass, requiredTop{ireq})
            error('InputData_Bypass must contain the top-level field "%s" in bypass mode.', requiredTop{ireq});
        end
    end

    G     = InputData_Bypass.general;
    flags = normalize_flags_struct(InputData_Bypass.flags);

    % ---------------- Non-raster ETP forcing file ----------------
    % Used later in preprocessing when:
    %   flags.flag_ETP = 1
    %   flags.flag_input_ETP_map = 0
    %
    % This can be either .xlsx or .csv and should follow the same layout as
    % the legacy ETP_input_data.xlsx file.
    ETP_input_spreadsheet = get_optional_field(G,'etp_input_spreadsheet','');

    requiredFlags = { ...
        'flag_rainfall','flag_spatial_rainfall','flag_ETP','flag_input_rainfall_map', ...
        'flag_rainfall_multiple_runs','flag_data_source','flag_inflow','flag_satellite_rainfall', ...
        'flag_alternated_blocks','flag_huff','flag_stage_hydrograph','flag_input_ETP_map', ...
        'flag_timestep','flag_infiltration','flag_critical','flag_D8','flag_CA','flag_inertial', ...
        'flag_waterbalance','flag_waterquality','flag_reservoir','flag_wq_model','flag_groundwater_modeling', ...
        'flag_real_time_satellite_rainfall','flag_dam_break','flag_human_instability','flag_boundary', ...
        'flag_numerical_scheme','flag_outlet_type','flag_adaptive_timestepping','flag_neglect_infiltration_river', ...
        'flag_subgrid','flag_spatial_albedo','flag_river_rasters','flag_baseflow','flag_kinematic','flag_diffusive', ...
        'flag_DTM','flag_abstraction','flag_overbanks','flag_snow_modeling','flag_WQ_Rasters', ...
        'flag_GPU','flag_single','flag_warmup','flag_initial_buildup', ...
        'flag_resample','flag_smoothening','flag_trunk','flag_fill_DEM','flag_smooth_cells','flag_reduce_DEM', ...
        'flag_export_maps','flag_river_heigth_compensation','flag_dashboard','flag_elapsed_time','flag_obs_gauges' ...
    };
    missing = requiredFlags(~isfield(flags, requiredFlags));
    if ~isempty(missing)
        error("Missing required flags in InputData_Bypass.flags: %s", strjoin(missing, ", "));
    end

    % ---------------- Running Control ----------------
    time_step_model = require_field(G,'time_step_model');
    running_control.min_time_step = require_field(G,'min_time_step');
    running_control.max_time_step = require_field(G,'max_time_step');
    running_control.time_step_increments = get_optional_field(G,'time_step_increments',0.001);
    running_control.time_step_change = require_field(G,'time_step_change');

    Courant_Parameters.alfa_min = require_field(G,'alfa_min');
    Courant_Parameters.alfa_max = require_field(G,'alfa_max');

    date_begin = cast_to_datetime_local(require_field(G,'date_begin'));
    date_end   = cast_to_datetime_local(require_field(G,'date_end'));

    running_control.routing_time = round(minutes(date_end - date_begin), 4);
    if isfield(G,'routing_time')
        running_control.routing_time = G.routing_time;
    end
    if running_control.routing_time < 0
        error('Please make sure Date End is later than Date Begin.')
    end

    % ---------------- Little Constraints ----------------
    if flags.flag_infiltration == 0
        flags.flag_ETP = 0;
        warning('Since no infiltration occurs, HydroPol2D is assuming no evapotranspiration is happening. Please check accordingly.')
        pause(0.5)
    end

    running_control.volume_error = get_optional_field(G,'volume_error',5);
    running_control.factor_reduction = get_optional_field(G,'factor_reduction',2);

    if flags.flag_input_rainfall_map + flags.flag_satellite_rainfall + flags.flag_real_time_satellite_rainfall > 1
        error('Please choose only one type of spatial rainfall data.')
    end

    if flags.flag_inertial == 1 && flags.flag_CA == 1
        error('Please, add either diffusive or inertial flag.')
    end

    % ---------------- Watershed / outlet ----------------
    outlet_type = get_optional_field(G,'outlet_type',1);

    slope_outlet   = require_field(G,'slope_outlet');
    n_outlets_data = get_optional_field(G,'n_outlets_data',1);

    % ---------------- Maps and Plots Control ----------------
    running_control.record_time_maps             = require_field(G,'record_time_maps');
    running_control.record_time_hydrographs      = require_field(G,'record_time_hydrographs');
    Pol_min                                      = require_field(G,'Pol_min');
    depth_wse                                    = require_field(G,'depth_wse');
    flags.flag_wse                               = require_field(G,'flag_wse');
    running_control.record_time_spatial_rainfall = require_field(G,'record_time_spatial_rainfall');
    time_save_ETP                                = require_field(G,'time_save_ETP');
    running_control.record_time_spatial_etp      = require_field(G,'record_time_spatial_ETP');
    Krs_ETP                                      = get_optional_field(G,'Krs_ETP',0.13);
    albedo                                       = get_optional_field(G,'albedo',0.25);

    if flags.flag_input_rainfall_map == 1
        if running_control.record_time_maps < running_control.record_time_spatial_rainfall
            error('If the input rainfall maps are provided, It is not possible to save maps with fewer temporal resolution than rainfall')
        end
    end

    % ---------------- Routing Parameters ----------------
    CA_States.depth_tolerance = get_optional_field(G,'CA_depth_tolerance',1);

    GIS_data.alfa_1 = require_field(G,'alfa_1');
    GIS_data.alfa_2 = require_field(G,'alfa_2');
    GIS_data.beta_1 = require_field(G,'beta_1');
    GIS_data.beta_2 = require_field(G,'beta_2');
    LULC_Parameters.River_Manning = require_field(G,'Manning');
    River_K_coeff = get_optional_field(G,'River_K_coeff',nan);

    % ---------------- Water Quality Inputs ----------------
    ADD    = require_field(G,'ADD');
    min_Bt = require_field(G,'min_Bt');
    Bmin   = require_field(G,'Bmin');
    Bmax   = require_field(G,'Bmax');

    % ---------------- DEM Smoothing, Imposemin, Resample ----------------
    GIS_data.min_area  = require_field(G,'min_area');
    GIS_data.tau       = require_field(G,'tau');
    GIS_data.K_value   = require_field(G,'K_value');
    GIS_data.sl        = require_field(G,'sl');
    GIS_data.slope_DTM = require_field(G,'slope_DTM');

    % ---------------- TopoToolBox Folder ----------------
    topo_path = require_field(G,'topo_path');

    % ---------------- Human Instability ----------------
    if flags.flag_human_instability == 1 || flags.flag_human_instability == 3
        if ~isfield(InputData_Bypass,'Human_Instability')
            error('Bypass mode requires InputData_Bypass.Human_Instability when flag_human_instability is active.');
        end
        Human_Instability = InputData_Bypass.Human_Instability;
    else
        Human_Instability = [];
    end

    % ---------------- Snow Properties ----------------
    if flags.flag_snow_modeling == 1
        if ~isfield(InputData_Bypass,'Snow_Properties')
            error('Bypass mode requires InputData_Bypass.Snow_Properties when flag_snow_modeling = 1.');
        end
        Snow_Properties = InputData_Bypass.Snow_Properties;
    end

    % ---------------- Design Storms ----------------
    % In the finalized bypass script, these live under InputData_Bypass.general
    Design_Storm_Parameters.RP                = require_field(G,'RP');
    Design_Storm_Parameters.Rainfall_Duration = require_field(G,'Rainfall_Duration');
    Design_Storm_Parameters.K                 = require_field(G,'K');
    Design_Storm_Parameters.a                 = require_field(G,'a');
    Design_Storm_Parameters.b                 = require_field(G,'b');
    Design_Storm_Parameters.c                 = require_field(G,'c');
    Design_Storm_Parameters.time_step         = require_field(G,'dt_design');

    if flags.flag_huff == 1 && flags.flag_alternated_blocks == 1
        error('Please, enter either Alternated Blocks or Huff hyetograph.')
    end

    % ---------------- Input Rainfall Maps ----------------
    if flags.flag_input_rainfall_map == 1 || flags.flag_satellite_rainfall == 1 || flags.flag_real_time_satellite_rainfall == 1

        if flags.flag_input_rainfall_map == 1
            if ~isfield(InputData_Bypass,'Input_Rainfall')
                error('Bypass mode requires InputData_Bypass.Input_Rainfall when flag_input_rainfall_map = 1.');
            end
            Input_Rainfall = normalize_map_input(InputData_Bypass.Input_Rainfall,'Input_Rainfall');
        else
            Input_Rainfall = struct();
        end

        flags.flag_spatial_rainfall = 1;
        flags.flag_rainfall = 1;
    end

    % ---------------- Input transpiration and evaporation Maps ----------------
    if flags.flag_input_ETP_map == 1
        if ~isfield(InputData_Bypass,'Input_Transpiration') || ~isfield(InputData_Bypass,'Input_Evaporation')
            error('Bypass mode requires InputData_Bypass.Input_Transpiration and InputData_Bypass.Input_Evaporation when flag_input_ETP_map = 1.');
        end

        % In standalone mode, respect the user-provided schedules
        Input_Transpiration = normalize_map_input(InputData_Bypass.Input_Transpiration,'Input_Transpiration');
        Input_Evaporation   = normalize_map_input(InputData_Bypass.Input_Evaporation,'Input_Evaporation');
    end

    % ---------------- Satellite / real-time rainfall generated schedules ----------------
    if flags.flag_satellite_rainfall == 1
        time_maps = 60; % Only valid for 1h maps
        Input_Rainfall.time = (0:60:running_control.routing_time); % time in minutes
        Input_Rainfall.num_obs_maps = sum(~isnan(Input_Rainfall.time));
        Input_Rainfall.time = Input_Rainfall.time(1:Input_Rainfall.num_obs_maps);

        flags.flag_spatial_rainfall = 1;
        flags.flag_rainfall = 1;
        flags.flag_real_time_satellite_rainfall = 0;
        flags.flag_input_rainfall_map = 0;

    elseif flags.flag_real_time_satellite_rainfall == 1
        time_maps = 60; % Only valid for 1h maps
        Input_Rainfall.time = (0:60:360); % time in minutes
        Input_Rainfall.num_obs_maps = sum(~isnan(Input_Rainfall.time));
        Input_Rainfall.time = Input_Rainfall.time(1:Input_Rainfall.num_obs_maps);

        flags.flag_spatial_rainfall = 1;
        flags.flag_rainfall = 1;
        flags.flag_satellite_rainfall = 0;
    end

    % ---------------- LULC DATA ----------------
    [LULC_name, lulc_parameters, n_lulc, imp_index, LULC_index] = unpack_class_table(InputData_Bypass.LULC,'LULC');

    % ---------------- SOIL DATA ----------------
    [SOIL_name, soil_parameters, n_soil, ~, SOIL_index] = unpack_class_table(InputData_Bypass.SOIL,'SOIL');

    % ---------------- Rainfall ----------------
    if flags.flag_rainfall == 1 && flags.flag_alternated_blocks ~= 1 && flags.flag_huff ~= 1 && ...
            flags.flag_input_rainfall_map ~= 1 && flags.flag_spatial_rainfall ~= 1 && ...
            flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1

        if ~isfield(InputData_Bypass,'Rainfall_Parameters')
            error('Bypass mode requires InputData_Bypass.Rainfall_Parameters for lumped rainfall mode.');
        end
        Rainfall_Parameters = normalize_rainfall_parameters(InputData_Bypass.Rainfall_Parameters);

    elseif flags.flag_alternated_blocks == 1 && flags.flag_input_rainfall_map ~= 1 && ...
            flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && ...
            flags.flag_satellite_rainfall ~= 1

        [Rainfall_Parameters.time_rainfall,Rainfall_Parameters.intensity_rainfall,~,idf] = ...
            alternated_blocks(Design_Storm_Parameters.Rainfall_Duration(end), Design_Storm_Parameters.time_step, ...
            Design_Storm_Parameters.K, Design_Storm_Parameters.a, Design_Storm_Parameters.b, ...
            Design_Storm_Parameters.c, Design_Storm_Parameters.RP, 1, Paths.FigPDF);

        Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
        Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
        Rainfall_Parameters.time_step_rainfall = Design_Storm_Parameters.time_step;
        Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end);
        running_control.routing_time = max(running_control.routing_time, Rainfall_Parameters.rainfall_duration);

    elseif flags.flag_huff == 1 && flags.flag_input_rainfall_map ~= 1 && ...
            flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && ...
            flags.flag_satellite_rainfall ~= 1

        [Rainfall_Parameters.time_rainfall,Rainfall_Parameters.intensity_rainfall,~] = ...
            Huff_Curves(Design_Storm_Parameters.Rainfall_Duration(end), Design_Storm_Parameters.time_step, ...
            Design_Storm_Parameters.K, Design_Storm_Parameters.a, Design_Storm_Parameters.b, ...
            Design_Storm_Parameters.c, Design_Storm_Parameters.RP, 1, Paths.FigPDF);

        Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
        Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
        Rainfall_Parameters.time_step_rainfall = Design_Storm_Parameters.time_step;
        Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end);
        running_control.routing_time = max(running_control.routing_time, Rainfall_Parameters.rainfall_duration);
    end

    % ---------------- Load TopoToolBox Tools ----------------
    addpath(genpath(char(topo_path)));

    % ---------------- Inflow Hydrograph ----------------
    if flags.flag_inflow == 1
        if ~isfield(InputData_Bypass,'Inflow')
            error('Bypass mode requires InputData_Bypass.Inflow when flag_inflow = 1.');
        end

        [Inflow_Parameters, Wshed_Properties] = build_inflow_from_bypass(InputData_Bypass.Inflow, flags, GIS_data, Wshed_Properties, DEM_raster);

        if isempty(Wshed_Properties.n_inlets) || all(Wshed_Properties.n_inlets == 0)
            error('Please, insert the inlet coordinate(s) in the bypass inflow data.')
        end
    else
        Inflow_Parameters = init_empty_inflow();
        Wshed_Properties.n_inlets = 0;
    end

    if flags.flag_inflow == 1 && flags.flag_stage_hydrograph == 1
        error('There is no way to enter both boundary conditions (inflow and stage-hydrograph). Choose one or none of them.')
    end

    % ---------------- Stage Hydrograph ----------------
    if flags.flag_stage_hydrograph == 1
        if ~isfield(InputData_Bypass,'Stage')
            error('Bypass mode requires InputData_Bypass.Stage when flag_stage_hydrograph = 1.');
        end

        [Stage_Parameters, Wshed_Properties] = build_stage_from_bypass(InputData_Bypass.Stage, flags, GIS_data, Wshed_Properties, DEM_raster);
    else
        Stage_Parameters = init_empty_stage();
        Wshed_Properties.n_inlets_stage = 0;
    end

    % Entering a boundary condition for rainfall if design storms are used
    if flags.flag_huff == 1 || flags.flag_alternated_blocks == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_real_time_satellite_rainfall == 0 && flags.flag_satellite_rainfall == 0
        flags.flag_satellite_rainfall = 0;
        flags.flag_real_time_satellite_rainfall = 0;
        flags.flag_rainfall = 1;
        flags.flag_spatial_rainfall = 0;
    end

    % ---------------- Inflow hydrograph rate ----------------
    if flags.flag_inflow == 1
        if flags.flag_resample
            Area = GIS_data.resolution_resample^2 * Wshed_Properties.n_inlets; % m2
            Inflow_Parameters.inflow_discharge = Inflow_Parameters.inflow_discharge(:,1:Inflow_Parameters.n_stream_gauges);
            Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_discharge ./ Area * 1000 * 3600; % mm/h
        else
            Area = Wshed_Properties.Resolution^2 * Wshed_Properties.n_inlets; % m2
            Inflow_Parameters.inflow_discharge = Inflow_Parameters.inflow_discharge(:,1:Inflow_Parameters.n_stream_gauges);
            Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_discharge ./ Area * 1000 * 3600; % mm/h
        end
    else
        Inflow_Parameters.inflow_hydrograph_rate = [];
    end

    % ---------------- Reservoir Data ----------------
    if flags.flag_reservoir == 1
        if ~isfield(InputData_Bypass,'Reservoir_Data')
            error('Bypass mode requires InputData_Bypass.Reservoir_Data when flag_reservoir = 1.');
        end
        Reservoir_Data = normalize_reservoir_data(InputData_Bypass.Reservoir_Data);
    end

    % ---------------- Observed Gauges ----------------
    if flags.flag_obs_gauges == 1
        if ~isfield(InputData_Bypass,'Obs')
            error('Bypass mode requires InputData_Bypass.Obs when flag_obs_gauges = 1.');
        end
        Obs = InputData_Bypass.Obs;
    end

    clear input_data input_table

else
    %% ====================================================================
    % ORIGINAL EXCEL / SPREADSHEET WORKFLOW (preserved)
    % ====================================================================

%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%
% GD = readcell(model_folder,'Sheet','General_Data');  % raw grid (names anywhere)

% ---------------- Running Control (name-based) ----------------
time_step_model = xlnum(GD,'time_step_model')/60;   % sheet is seconds -> minutes
running_control.min_time_step = xlnum(GD,'min_time_step');
running_control.max_time_step = xlnum(GD,'max_time_step');

running_control.time_step_increments = 0.001;       % keep as constant (your choice)
running_control.time_step_change = xlnum(GD,'time_step_change');

Courant_Parameters.alfa_min = xlnum(GD,'alfa_min');
Courant_Parameters.alfa_max = xlnum(GD,'alfa_max');

date_begin = xldatetime(GD,'date_begin');
date_end   = xldatetime(GD,'date_end');

running_control.routing_time = round(minutes(date_end - date_begin), 4);  % minutes
if running_control.routing_time < 0
    error('Please make sure Date End is later than Date Begin.')
end

% ==========================================================
% FLAGS (name-based read from 'Flags' sheet)
% ==========================================================
FlagsGrid = readcell(model_folder,'Sheet','Flags');

flags = read_flags_sheet(FlagsGrid);

requiredFlags = { ...
    'flag_rainfall','flag_spatial_rainfall','flag_ETP','flag_input_rainfall_map', ...
    'flag_rainfall_multiple_runs','flag_data_source','flag_inflow','flag_satellite_rainfall', ...
    'flag_alternated_blocks','flag_huff','flag_stage_hydrograph','flag_input_ETP_map', ...
    'flag_timestep','flag_infiltration','flag_critical','flag_D8','flag_CA','flag_inertial', ...
    'flag_waterbalance','flag_waterquality','flag_reservoir','flag_wq_model','flag_groundwater_modeling', ...
    'flag_real_time_satellite_rainfall','flag_dam_break','flag_human_instability','flag_boundary', ...
    'flag_numerical_scheme','flag_outlet_type','flag_adaptive_timestepping','flag_neglect_infiltration_river', ...
    'flag_subgrid','flag_spatial_albedo','flag_river_rasters','flag_baseflow','flag_kinematic','flag_diffusive', ...
    'flag_DTM','flag_abstraction','flag_overbanks','flag_snow_modeling','flag_WQ_Rasters', ...
    'flag_GPU','flag_single','flag_warmup','flag_initial_buildup', ...
    'flag_resample','flag_smoothening','flag_trunk','flag_fill_DEM','flag_smooth_cells','flag_reduce_DEM', ...
    'flag_export_maps','flag_river_heigth_compensation','flag_dashboard','flag_elapsed_time','flag_obs_gauges' ...
};

missing = requiredFlags(~isfield(flags, requiredFlags));
if ~isempty(missing)
    error("Missing flags in Flags sheet: %s", strjoin(missing, ", "));
end

% Little Constraints
if flags.flag_infiltration == 0
    flags.flag_ETP = 0;
    warning('Since no infiltration occurs, HydroPol2D is assuming no evapotranspiration is happening. Please check accordingly.')
    pause(0.5)
end

running_control.volume_error = 5; % m3 [not needed and not considerd as input anymore]
running_control.factor_reduction = 2; % [not considered as input anymore]

if flags.flag_input_rainfall_map + flags.flag_satellite_rainfall + flags.flag_real_time_satellite_rainfall > 1
    error('Please choose only one type of spatial rainfall data.')
end

if flags.flag_inertial == 1 && flags.flag_CA == 1
    error('Please, add either diffusive or inertial flag.')
end

% ---------------- Watershed / outlet ----------------
outlet_type = single(flags.flag_outlet_type);


slope_outlet   = xlnum(GD,'slope_outlet');
n_outlets_data = 1; % Fixed in one

% ---------------- Maps and Plots Control ----------------
running_control.record_time_maps             = xlnum(GD,'record_time_maps');
running_control.record_time_hydrographs      = xlnum(GD,'record_time_hydrographs');
Pol_min                                      = xlnum(GD,'Pol_min');
depth_wse                                    = xlnum(GD,'depth_wse');
flags.flag_wse                               = xlnum(GD,'flag_wse');
running_control.record_time_spatial_rainfall = xlnum(GD,'record_time_spatial_rainfall');
time_save_ETP                                = xlnum(GD,'time_save_ETP');
running_control.record_time_spatial_etp      = xlnum(GD,'record_time_spatial_ETP');

Krs_ETP = 0.13; % Defaults
albedo  = 0.25; % Defaults

if flags.flag_input_rainfall_map == 1
    if running_control.record_time_maps < running_control.record_time_spatial_rainfall
        error('If the input rainfall maps are provided, It is not possible to save maps with fewer temporal resolution than rainfall')
    end
end

% ---------------- Routing Parameters ----------------
CA_States.depth_tolerance = 1; % Default 1 mm

GIS_data.alfa_1 = xlnum(GD,'alfa_1');
GIS_data.alfa_2 = xlnum(GD,'alfa_2');
GIS_data.beta_1 = xlnum(GD,'beta_1');
GIS_data.beta_2 = xlnum(GD,'beta_2');
LULC_Parameters.River_Manning = xlnum(GD,'Manning');
River_K_coeff = nan; % Deactivated

% ---------------- Water Quality Inputs ----------------
ADD    = xlnum(GD,'ADD');
min_Bt = xlnum(GD,'min_Bt');
Bmin   = xlnum(GD,'Bmin');
Bmax   = xlnum(GD,'Bmax');

% ---------------- DEM Smoothing, Imposemin, Resample ----------------
GIS_data.min_area  = xlnum(GD,'min_area');
GIS_data.tau       = xlnum(GD,'tau');
GIS_data.K_value   = xlnum(GD,'K_value');
GIS_data.sl        = xlnum(GD,'sl');
GIS_data.slope_DTM = xlnum(GD,'slope_DTM');

% ---------------- TopoToolBox Folder ----------------
topo_path = xlget(GD,'topo_path');

% Human Instability
if flags.flag_human_instability == 1
    human_table = readtable('human_risk.xlsx');
    Human_Instability.mu = table2array(human_table(3,2)); % friction coefficient
    Human_Instability.Cd = table2array(human_table(4,2));
    Human_Instability.ro_person = table2array(human_table(5,2)); % kg/m3
    Human_Instability.weight_person = table2array(human_table(6,2)); % kg
    Human_Instability.height_person = table2array(human_table(7,2)); % meters
    Human_Instability.width1_person = table2array(human_table(8,2)); % meters
    Human_Instability.width2_person = table2array(human_table(9,2)); % meters
    Human_Instability.ro_water = table2array(human_table(1,2)); % kg/m3
    Human_Instability.gravity = table2array(human_table(2,2)); % m/s2
elseif flags.flag_human_instability == 3
    human_table = readtable('human_risk.xlsx');
    Human_Instability.ro_water = table2array(human_table(1,10)); % kg/m3
    Human_Instability.gravity = table2array(human_table(2,10)); % m/s2
    Human_Instability.mu = table2array(human_table(3,10)); % friction coefficient
    Human_Instability.Cd = table2array(human_table(4,10));
    Human_Instability.m_c_m = table2array(human_table(5,10));
    Human_Instability.y_c_m = table2array(human_table(6,10));
    Human_Instability.m_t_m = table2array(human_table(7,10));
    Human_Instability.y_t_m = table2array(human_table(8,10));
    Human_Instability.m_a_m = table2array(human_table(9,10));
    Human_Instability.y_a_m = table2array(human_table(10,10));
    Human_Instability.m_o_m = table2array(human_table(11,10));
    Human_Instability.y_o_m = table2array(human_table(12,10));
    Human_Instability.m_c_f = table2array(human_table(13,10));
    Human_Instability.y_c_f = table2array(human_table(14,10));
    Human_Instability.m_t_f = table2array(human_table(15,10));
    Human_Instability.y_t_f = table2array(human_table(16,10));
    Human_Instability.m_a_f = table2array(human_table(17,10));
    Human_Instability.y_a_f = table2array(human_table(18,10));
    Human_Instability.m_o_f = table2array(human_table(19,10));
    Human_Instability.y_o_f = table2array(human_table(20,10));
    Human_Instability.w_c_m = sqrt((Human_Instability.m_c_m*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_c_m));
    Human_Instability.d_c_m = Human_Instability.w_c_m/2;
    Human_Instability.w_t_m = sqrt((Human_Instability.m_t_m*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_t_m));
    Human_Instability.d_t_m = Human_Instability.w_t_m/2;
    Human_Instability.w_a_m = sqrt((Human_Instability.m_a_m*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_a_m));
    Human_Instability.d_a_m = Human_Instability.w_a_m/2;
    Human_Instability.w_o_m = sqrt((Human_Instability.m_o_m*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_o_m));
    Human_Instability.d_o_m = Human_Instability.w_o_m/2;
    Human_Instability.w_c_f = sqrt((Human_Instability.m_c_f*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_c_f));
    Human_Instability.d_c_f = Human_Instability.w_c_f/2;
    Human_Instability.w_t_f = sqrt((Human_Instability.m_t_f*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_t_f));
    Human_Instability.d_t_f = Human_Instability.w_t_f/2;
    Human_Instability.w_a_f = sqrt((Human_Instability.m_a_f*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_a_f));
    Human_Instability.d_a_f = Human_Instability.w_a_f/2;
    Human_Instability.w_o_f = sqrt((Human_Instability.m_o_f*16)/(3*pi()*Human_Instability.ro_water*Human_Instability.y_o_f));
    Human_Instability.d_o_f = Human_Instability.w_o_f/2;
    temp = cell2mat({Human_Instability.m_c_m/Human_Instability.y_c_m, Human_Instability.m_t_m/Human_Instability.y_t_m, Human_Instability.m_a_m/Human_Instability.y_a_m, Human_Instability.m_o_m/Human_Instability.y_o_m,...
        Human_Instability.m_c_f/Human_Instability.y_c_f, Human_Instability.m_t_f/Human_Instability.y_t_f, Human_Instability.m_a_f/Human_Instability.y_a_f, Human_Instability.m_o_f/Human_Instability.y_o_f});
    unique_vals = unique(temp);
    list={'_cm','_tm','_am','_om','_cf','_tf','_af','_of'};
    names={'Child Boy','Teen Male','Adult Man','Elderly Man',...
        'Child Girl','Teen Female','Adult Women','Elderly Women'};
    names={[strcat(num2str(round(Human_Instability.m_c_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_c_m,1)),'m') newline names{1}],...
        [strcat(num2str(round(Human_Instability.m_t_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_t_m,1)),'m') newline names{2}],...
        [strcat(num2str(round(Human_Instability.m_a_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_a_m,1)),'m') newline names{3}],...
        [strcat(num2str(round(Human_Instability.m_o_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_o_m,1)),'m') newline names{4}],...
        [strcat(num2str(round(Human_Instability.m_c_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_c_f,1)),'m') newline names{5}],...
        [strcat(num2str(round(Human_Instability.m_t_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_t_f,1)),'m') newline names{6}],...
        [strcat(num2str(round(Human_Instability.m_a_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_a_f,1)),'m') newline names{7}],...
        [strcat(num2str(round(Human_Instability.m_o_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_o_f,1)),'m') newline names{8}]};
    for i = 1:length(unique_vals)
        Human_Instability.order(i) = find(temp==unique_vals(i),1,'first');
    end
    for i = 1:8
        Human_Instability.list{i} = list{Human_Instability.order(i)};
        Human_Instability.names{i} = names{Human_Instability.order(i)};
    end
else
    Human_Instability = [];
end

% ---------------- Design Storms (from General_Data by label) ----------------
Design_Storm_Parameters.RP                = xlnum(GD,'RP');
Design_Storm_Parameters.Rainfall_Duration = xlnum(GD,'Rainfall Duration');
Design_Storm_Parameters.K                 = xlnum(GD,'K');
Design_Storm_Parameters.a                 = xlnum(GD,'a');
Design_Storm_Parameters.b                 = xlnum(GD,'b');
Design_Storm_Parameters.c                 = xlnum(GD,'c');
Design_Storm_Parameters.time_step         = xlnum(GD,'dt_design');

if flags.flag_huff == 1 && flags.flag_alternated_blocks == 1
    error('Please, enter either Alternated Blocks or Huff hyetograph.')
end

% Input Rainfall Maps
if flags.flag_input_rainfall_map == 1 || flags.flag_satellite_rainfall == 1 || flags.flag_real_time_satellite_rainfall == 1

    T_rain = xlblock_2col( ...
        GD, ...
        'Sattelite or Radar Rainfall', ...
        'Time [min]', ...
        'Raster Directory with values in mm/h' ...
    );

    Input_Rainfall.time = T_rain.time;
    Input_Rainfall.num_obs_maps = numel(Input_Rainfall.time);
    Input_Rainfall.labels_Directory = cellstr(T_rain.directory);

    flags.flag_spatial_rainfall = 1;
    flags.flag_rainfall = 1;
end

% Input transpiration and evaporation Maps
if flags.flag_input_ETP_map == 1

    T_tr = xlblock_2col( ...
        GD, ...
        'Sattelite transpiration', ...
        'Time [day]', ...
        'Raster Directory with values in mm/day' ...
    );

    T_ev = xlblock_2col( ...
        GD, ...
        'Sattelite Evaporation', ...
        'Time [day]', ...
        'Raster Directory with values in mm/day' ...
    );

    Input_Transpiration.time = T_tr.time;
    Input_Transpiration.num_obs_maps = numel(Input_Transpiration.time);
    Input_Transpiration.labels_Directory = cellstr(T_tr.directory);

    Input_Evaporation.time = T_ev.time;
    Input_Evaporation.num_obs_maps = numel(Input_Evaporation.time);
    Input_Evaporation.labels_Directory = cellstr(T_ev.directory);
end

if flags.flag_satellite_rainfall == 1
    time_maps = 60; % Only valid for 1h maps
    Input_Rainfall.time = (0:60:running_control.routing_time); % time in minutes
    Input_Rainfall.num_obs_maps = sum(~isnan(Input_Rainfall.time));
    Input_Rainfall.time = Input_Rainfall.time(1:Input_Rainfall.num_obs_maps);

    flags.flag_spatial_rainfall = 1;
    flags.flag_rainfall = 1;
    flags.flag_real_time_satellite_rainfall = 0;
    flags.flag_input_rainfall_map = 0;
elseif flags.flag_real_time_satellite_rainfall == 1
    time_maps = 60; % Only valid for 1h maps
    Input_Rainfall.time = (0:60:360); % time in minutes, set up for 6 hours until save a register
    Input_Rainfall.num_obs_maps = sum(~isnan(Input_Rainfall.time));
    Input_Rainfall.time = Input_Rainfall.time(1:Input_Rainfall.num_obs_maps);

    flags.flag_spatial_rainfall = 1;
    flags.flag_rainfall = 1;
    flags.flag_satellite_rainfall = 0;
end

if flags.flag_input_ETP_map == 1
    time_maps = 1440; % Only valid for 1day maps
    Input_Transpiration.time = (0:time_maps:running_control.routing_time); % time in minutes
    Input_Transpiration.num_obs_maps = sum(~isnan(Input_Transpiration.time));
    Input_Transpiration.time = Input_Transpiration.time(1:Input_Transpiration.num_obs_maps);

    Input_Evaporation.time = (0:time_maps:running_control.routing_time); % time in minutes
    Input_Evaporation.num_obs_maps = sum(~isnan(Input_Evaporation.time));
    Input_Evaporation.time = Input_Evaporation.time(1:Input_Evaporation.num_obs_maps);
end

%%%%%%%%%%%%%% LULC DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable('LULC_parameters.xlsx');
input_data = table2array(input_table(:,2:end)); % numbers
LULC_name = input_table(:,1);
lulc_parameters = input_data(:,2:end);
n_lulc = sum(lulc_parameters(:,1)>=0); % Number of LULC
lulc_parameters = lulc_parameters(1:n_lulc,:);
imp_index = lulc_parameters(1,end);
LULC_index = input_data(:,1);

%%%%%%%%%%%%%% SOIL DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable('SOIL_parameters.xlsx');
input_data = table2array(input_table(:,2:end)); % numbers
soil_parameters = input_data(:,2:end); % Number of Soil Types
n_soil = sum(soil_parameters(:,1)>=0);
soil_parameters = soil_parameters(1:n_soil,:);
SOIL_index = input_data(:,1);
SOIL_name = input_table(:,1);

% Rainfall
if flags.flag_rainfall == 1 && flags.flag_alternated_blocks ~= 1 && flags.flag_huff ~= 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
    if flags.flag_rainfall_multiple_runs == 0
        input_table = readtable('Rainfall_Intensity_Data.xlsx');
        input_data = table2array(input_table);
        Rainfall_Parameters.intensity_rainfall = input_data(:,2); % mm/h
        Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
        Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
        Rainfall_Parameters.time_rainfall = input_data(:,1); % min
        Rainfall_Parameters.time_step_rainfall = input_data(2,1) - input_data(1,1); % min
        Rainfall_Parameters.rainfall_duration = input_data(2,1) - input_data(1,1); % min
    else
        path = Paths.FigPDF;
        input_table = readall(spreadsheetDatastore('rainfalls_cc.xlsx'));
        for i = 1:length(input_table.Properties.VariableNames)
            if sum(input_table.Properties.VariableNames{i}(1:4) == 'time') > 1
                if sum(find(ismember(ls(path), input_table.Properties.VariableNames(i))))>1
                    continue
                else
                    name_time = input_table.Properties.VariableNames{i};
                    mkdir(strcat(path,'\',name_time));
                    Rainfall_Parameters.name_time = i;
                    break
                end
            end
        end
        input_table.(input_table.Properties.VariableNames{name_time}) = datetime(input_table.(input_table.Properties.VariableNames{name_time}),'InputFormat','yyyy-MM-dd HH:mm:ss');
        column_time = find(strcmpi(input_table.Properties.VariableNames, input_table.Properties.VariableNames{name_time}));
        column_rain = column_time + 1;
        temp_rain = table2array(input_table(:,column_rain));
        temp_time = table2array(input_table(:,column_time));
        Rainfall_Parameters.intensity_rainfall = temp_rain(~isnat(temp_time));
        Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
        Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
        Rainfall_Parameters.time_rainfall = input_table.(input_table.Properties.VariableNames{name_time})(~isnat(temp_time));
        Rainfall_Parameters.time_step_rainfall = minutes(Rainfall_Parameters.time_rainfall(2) - Rainfall_Parameters.time_rainfall(1));
        Rainfall_Parameters.rainfall_duration = minutes(Rainfall_Parameters.time_rainfall(2) - Rainfall_Parameters.time_rainfall(1));
        date_begin = Rainfall_Parameters.time_rainfall(1);
        date_end = Rainfall_Parameters.time_rainfall(end);
        Rainfall_Parameters.time_rainfall = minutes(Rainfall_Parameters.time_rainfall(:)-Rainfall_Parameters.time_rainfall(1));
        running_control.routing_time = (Rainfall_Parameters.n_obs_rainfall-1)*Rainfall_Parameters.time_step_rainfall;
    end
elseif flags.flag_alternated_blocks == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
    [Rainfall_Parameters.time_rainfall,Rainfall_Parameters.intensity_rainfall,~,idf] = ...
        alternated_blocks(Design_Storm_Parameters.Rainfall_Duration(end),Design_Storm_Parameters.time_step,Design_Storm_Parameters.K,Design_Storm_Parameters.a,Design_Storm_Parameters.b,Design_Storm_Parameters.c,Design_Storm_Parameters.RP,1,Paths.FigPDF);
    Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
    Rainfall_Parameters.time_step_rainfall = Design_Storm_Parameters.time_step; % min
    Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end); % min
    running_control.routing_time = max(running_control.routing_time,Rainfall_Parameters.rainfall_duration);
elseif flags.flag_huff == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
    [Rainfall_Parameters.time_rainfall,Rainfall_Parameters.intensity_rainfall,~] = ...
        Huff_Curves(Design_Storm_Parameters.Rainfall_Duration(end), Design_Storm_Parameters.time_step, Design_Storm_Parameters.K,Design_Storm_Parameters.a,Design_Storm_Parameters.b,Design_Storm_Parameters.c,Design_Storm_Parameters.RP,1,Paths.FigPDF);
    Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
    Rainfall_Parameters.time_step_rainfall = Design_Storm_Parameters.time_step; % min
    Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end); % min
    running_control.routing_time = max(running_control.routing_time,Rainfall_Parameters.rainfall_duration);
end

% Load TopoToolBox Tools
addpath(genpath(char(topo_path)));

%%% ---- Inflow Hydrograph ---- %%%
input_table = readtable('Inflow_Hydrograph.xlsx');
input_table_labels = input_table(1:2,:); %#ok<NASGU>
input_table_values = input_table(3:end,:);
input_table_values = convertvars(input_table_values, 1:size(input_table_values,2), 'string');
input_table_values = convertvars(input_table_values, 1:size(input_table_values,2), 'double');

n_inflows_max = 50;
for i = 1:n_inflows_max
    Inflow_Parameters.inflow_discharge(:,i) = input_table_values{:, (i-1)*5 + 2};
end
Inflow_Parameters.n_stream_gauges = sum(sum(Inflow_Parameters.inflow_discharge(1,:)>=0));

Inflow_Parameters.easting_inlet_cells = nan(n_inflows_max,n_inflows_max);
Inflow_Parameters.northing_inlet_cells = nan(n_inflows_max,n_inflows_max);

for i = 1:Inflow_Parameters.n_stream_gauges
    if flags.flag_resample
        nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
        x_coordinates = round((-DEM_raster.georef.SpatialRef.XWorldLimits(1) + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/GIS_data.resolution_resample);
        y_coordinates = round((DEM_raster.georef.SpatialRef.YWorldLimits(2) - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/GIS_data.resolution_resample);
        points = [x_coordinates,y_coordinates];
        if size(points,1) ~= nonNanCount
            warning('Some of the inlet points are located at the same cell. Therefore, we are only considering one of them.')
        end
        effective_inlets = unique(points,'rows','stable');
        Inflow_Parameters.easting_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,1);
        Inflow_Parameters.northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2);
        Wshed_Properties.n_inlets(:,i) = size(effective_inlets,1);

        Inflow_Parameters.easting_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,1);
        northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2); %#ok<NASGU>
        Inflow_Parameters.n_inlets(:,i) = size(effective_inlets,1);
    else
        nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
        x_coordinates = round((-GIS_data.xulcorner  + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/Wshed_Properties.Resolution);
        y_coordinates = round((GIS_data.yulcorner - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/Wshed_Properties.Resolution);
        points = [x_coordinates,y_coordinates];
        if size(points,1) ~= nonNanCount
            warning('Some of the inlet points are located at the same cell. Therefore, we are only considering one of them.')
        end
        effective_inlets = unique(points,'rows','stable');
        Inflow_Parameters.easting_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,1);
        Inflow_Parameters.northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2);
        Wshed_Properties.n_inlets(:,i) = size(effective_inlets,1);
    end
end

time_inflow = double(table2array(input_table_values(:,1)));
Qall        = double(Inflow_Parameters.inflow_discharge(:,1:Inflow_Parameters.n_stream_gauges));

valid_rows = isfinite(time_inflow) & any(isfinite(Qall), 2);

time_inflow = time_inflow(valid_rows);
Qall        = Qall(valid_rows,:);

Qall(~isfinite(Qall)) = NaN;

if numel(time_inflow) < 2
    error('Inflow_Hydrograph.xlsx must contain at least 2 valid inflow rows.');
end

dt_inflow = time_inflow(2) - time_inflow(1);
if dt_inflow <= 0
    error('Inflow forcing time step must be positive.');
end

if any(abs(diff(time_inflow) - dt_inflow) > 1e-9)
    error('Inflow forcing time step must be constant.');
end

Inflow_Parameters.time_inflow      = time_inflow(:);
Inflow_Parameters.inflow_discharge = Qall;
Inflow_Parameters.n_stream_obs     = numel(time_inflow);
Inflow_Parameters.time_step_inflow = dt_inflow;

if flags.flag_inflow == 1
    if Wshed_Properties.n_inlets == 0
        error('Please, insert the inlet coordinate(s) in the Inflow_Hydrograph.xlsx file')
    end
end

if flags.flag_inflow == 1 && flags.flag_stage_hydrograph == 1
    error('There is no way to enter both boundary conditions (inflow and stage-hydrograph). Choose one or none of them.')
end

%%% ---- Stage Hydrograph ---- %%%
input_table = readtable('Stage_Hydrograph.xlsx');
input_table_labels = input_table(1:2,:); %#ok<NASGU>
input_table_values = input_table(3:end,:);
input_table_values = convertvars(input_table_values, 1:size(input_table_values,2), 'string');
input_table_values = convertvars(input_table_values, 1:size(input_table_values,2), 'double');

n_stage_max = 50;
for i = 1:n_stage_max
    Stage_Parameters.stage(:,i) = table2array(input_table(3:end,(i-1)*5 + 2));
end
Stage_Parameters.n_stage_gauges = sum(sum(Stage_Parameters.stage(1,:)>=0));

Stage_Parameters.easting_inlet_cells = nan(n_stage_max,n_stage_max);
Stage_Parameters.northing_inlet_cells = nan(n_stage_max,n_stage_max);

for i = 1:Stage_Parameters.n_stage_gauges
    if flags.flag_resample
        nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
        x_coordinates = round((-DEM_raster.georef.SpatialRef.XWorldLimits(1) + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/GIS_data.resolution_resample);
        y_coordinates = round((DEM_raster.georef.SpatialRef.YWorldLimits(2) - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/GIS_data.resolution_resample);
        points = [x_coordinates,y_coordinates];
        if size(points,1) ~= nonNanCount
            warning('Some of the inlet points are located at the same cell. Therefore, we are only considering one of them.')
        end
        effective_inlets = unique(points,'rows','stable');
        Stage_Parameters.easting_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,1);
        Stage_Parameters.northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2);
        Wshed_Properties.n_inlets_stage(:,i) = size(effective_inlets,1);

        Stage_Parameters.easting_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,1);
        Stage_Parameters.northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2);
        Wshed_Properties.n_inlets_stage(:,i) = size(effective_inlets,1);
    else
        nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
        x_coordinates = round((-GIS_data.xulcorner  + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/Wshed_Properties.Resolution);
        y_coordinates = round((GIS_data.yulcorner - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/Wshed_Properties.Resolution);
        points = [x_coordinates,y_coordinates];
        if size(points,1) ~= nonNanCount
            warning('Some of the inlet points are located at the same cell. Therefore, we are only considering one of them.')
        end
        effective_inlets = unique(points,'rows','stable');
        Stage_Parameters.easting_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,1);
        Stage_Parameters.northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2);
        Wshed_Properties.n_inlets_stage(:,i) = size(effective_inlets,1);
    end
end

time_stage = double(table2array(input_table_values(:,1)));
Hall       = double(Stage_Parameters.stage(:,1:Stage_Parameters.n_stage_gauges));

valid_rows = isfinite(time_stage) & any(isfinite(Hall), 2);

time_stage = time_stage(valid_rows);
Hall       = Hall(valid_rows,:);

Hall(~isfinite(Hall)) = NaN;

if numel(time_stage) < 2
    error('Stage_Hydrograph.xlsx must contain at least 2 valid stage rows.');
end

dt_stage = time_stage(2) - time_stage(1);
if dt_stage <= 0
    error('Stage forcing time step must be positive.');
end

if any(abs(diff(time_stage) - dt_stage) > 1e-3)
    error('Stage forcing time step must be constant.');
end

Stage_Parameters.time_stage      = time_stage(:);
Stage_Parameters.stage           = Hall;
Stage_Parameters.n_stage_obs     = numel(time_stage);
Stage_Parameters.time_step_stage = dt_stage;

% Entering a boundary condition for rainfall if design storms are used
if flags.flag_huff == 1 || flags.flag_alternated_blocks == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_real_time_satellite_rainfall == 0 && flags.flag_satellite_rainfall == 0
    flags.flag_satellite_rainfall = 0;
    flags.flag_real_time_satellite_rainfall = 0;
    flags.flag_rainfall = 1;
    flags.flag_spatial_rainfall = 0;
end

if flags.flag_resample
    Area = GIS_data.resolution_resample^2*Wshed_Properties.n_inlets; % m2
    Inflow_Parameters.inflow_discharge = Inflow_Parameters.inflow_discharge(:,1:Inflow_Parameters.n_stream_gauges);
    Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_discharge./Area*1000*3600; % mm/h
else
    Area = Wshed_Properties.Resolution^2*Wshed_Properties.n_inlets; % m2
    Inflow_Parameters.inflow_discharge = Inflow_Parameters.inflow_discharge(:,1:Inflow_Parameters.n_stream_gauges);
    Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_discharge./Area*1000*3600; % mm/h
end
clear input_data input_table

% Reservoir Data
if flags.flag_reservoir == 1
    input_table = readtable('BC_Control_structure.xlsx');
    Reservoir_Data.index = table2array(input_table(:,1));
    Reservoir_Data.x_us = table2array(input_table(:,2));
    Reservoir_Data.y_us = table2array(input_table(:,3));
    Reservoir_Data.k1 = table2array(input_table(:,4));
    Reservoir_Data.h1 = table2array(input_table(:,5));
    Reservoir_Data.k2 = table2array(input_table(:,6));
    Reservoir_Data.x_ds1 = table2array(input_table(:,7));
    Reservoir_Data.y_ds1 = table2array(input_table(:,8));
    Reservoir_Data.k3 = table2array(input_table(:,9));
    Reservoir_Data.h2 = table2array(input_table(:,10));
    Reservoir_Data.k4 = table2array(input_table(:,11));
    Reservoir_Data.x_ds2 = table2array(input_table(:,12));
    Reservoir_Data.y_ds2 = table2array(input_table(:,13));
end

clear input_table
end

%% ========================================================================
% LOCAL HELPER FUNCTIONS
% ========================================================================

function v = require_field(S, fieldName)
    if ~isfield(S, fieldName)
        error('Missing required field "%s" in bypass input.', fieldName);
    end
    v = S.(fieldName);
end

function v = get_optional_field(S, fieldName, defaultValue)
    if isfield(S, fieldName)
        v = S.(fieldName);
    else
        v = defaultValue;
    end
end

function flags = normalize_flags_struct(flagsIn)
    flags = struct();
    fn = fieldnames(flagsIn);
    for i = 1:numel(fn)
        field = matlab.lang.makeValidName(fn{i});
        flags.(field) = double(flagsIn.(fn{i}));
    end
    if ~isfield(flags,'flag_warmup')
        flags.flag_warmup = 0;
    end
end

function M = normalize_map_input(M, label)
    if ~isfield(M,'time') || ~isfield(M,'labels_Directory')
        error('%s must contain fields time and labels_Directory.', label);
    end

    M.time = double(M.time(:));

    if isstring(M.labels_Directory)
        M.labels_Directory = cellstr(M.labels_Directory);
    elseif ischar(M.labels_Directory)
        M.labels_Directory = {M.labels_Directory};
    end

    if numel(M.labels_Directory) ~= numel(M.time)
        error('%s.labels_Directory must have the same length as %s.time.', label, label);
    end

    valid = isfinite(M.time);
    M.time = M.time(valid);
    M.labels_Directory = M.labels_Directory(valid);

    if numel(M.time) < 2
        error('%s must contain at least 2 valid time entries.', label);
    end

    dt_map = M.time(2) - M.time(1);
    if dt_map <= 0
        error('%s time step must be positive.', label);
    end

    if any(abs(diff(M.time) - dt_map) > 1e-9)
        error('%s forcing time step must be constant.', label);
    end

    M.num_obs_maps = numel(M.time);
end

function [class_names, parameters, n_classes, imp_index, class_index] = unpack_class_table(S, label)
    imp_index = [];
    if istable(S)
        T = S;
        class_names = T(:,1);
        input_data = table2array(T(:,2:end));
    elseif isstruct(S) && isfield(S,'table')
        T = S.table;
        class_names = T(:,1);
        input_data = table2array(T(:,2:end));
    else
        class_names = require_field(S,'names');
        class_index = require_field(S,'index');
        parameters = require_field(S,'parameters');
        n_classes = size(parameters,1);
        if isfield(S,'impervious_index')
            imp_index = S.impervious_index;
        end
        return
    end

    class_index = input_data(:,1);
    parameters  = input_data(:,2:end);
    n_classes   = sum(parameters(:,1) >= 0);
    parameters  = parameters(1:n_classes,:);

    if strcmpi(label,'LULC')
        imp_index = parameters(1,end);
    end
end

function RP = normalize_rainfall_parameters(RP)
    required = {'time_rainfall','intensity_rainfall'};
    for i = 1:numel(required)
        if ~isfield(RP, required{i})
            error('InputData_Bypass.Rainfall_Parameters must contain field "%s".', required{i});
        end
    end

    RP.time_rainfall      = double(RP.time_rainfall(:));
    RP.intensity_rainfall = double(RP.intensity_rainfall(:));

    if numel(RP.time_rainfall) ~= numel(RP.intensity_rainfall)
        error('Rainfall_Parameters.time_rainfall and intensity_rainfall must have the same length.');
    end

    valid = isfinite(RP.time_rainfall) & isfinite(RP.intensity_rainfall) & (RP.intensity_rainfall >= 0);
    RP.time_rainfall      = RP.time_rainfall(valid);
    RP.intensity_rainfall = RP.intensity_rainfall(valid);

    if numel(RP.time_rainfall) < 2
        error('Rainfall_Parameters must contain at least 2 valid records.');
    end

    RP.time_step_rainfall = RP.time_rainfall(2) - RP.time_rainfall(1);

    if RP.time_step_rainfall <= 0
        error('Rainfall time step must be positive.');
    end

    if any(abs(diff(RP.time_rainfall) - RP.time_step_rainfall) > 1e-9)
        error('Rainfall forcing time step must be constant.');
    end

    RP.n_obs_rainfall    = numel(RP.intensity_rainfall);
    RP.rainfall_duration = RP.time_step_rainfall;   % keep legacy meaning
    RP.rainfall_end_time = RP.time_rainfall(end);   % explicit full forcing span
end

function [Inflow_Parameters, Wshed_Properties] = build_inflow_from_bypass(B, flags, GIS_data, Wshed_Properties, DEM_raster)
    required = {'time_inflow','inflow_discharge','x_coords','y_coords'};
    for i = 1:numel(required)
        if ~isfield(B, required{i})
            error('InputData_Bypass.Inflow must contain field "%s".', required{i});
        end
    end

    time_inflow = double(B.time_inflow(:));
    inflow_discharge = double(B.inflow_discharge);

    if size(inflow_discharge,1) ~= numel(time_inflow)
        error('InputData_Bypass.Inflow.inflow_discharge must have the same number of rows as time_inflow.');
    end

    valid_rows = isfinite(time_inflow) & any(isfinite(inflow_discharge), 2);
    time_inflow = time_inflow(valid_rows);
    inflow_discharge = inflow_discharge(valid_rows,:);

    inflow_discharge(~isfinite(inflow_discharge)) = NaN;

    if numel(time_inflow) < 2
        error('InputData_Bypass.Inflow must contain at least 2 valid records.');
    end

    dt_inflow = time_inflow(2) - time_inflow(1);
    if dt_inflow <= 0
        error('InputData_Bypass.Inflow time step must be positive.');
    end

    if any(abs(diff(time_inflow) - dt_inflow) > 1e-9)
        error('InputData_Bypass.Inflow forcing time step must be constant.');
    end

    x_coords = normalize_coord_input(B.x_coords, size(inflow_discharge,2), 'Inflow.x_coords');
    y_coords = normalize_coord_input(B.y_coords, size(inflow_discharge,2), 'Inflow.y_coords');

    n_stream_gauges = size(inflow_discharge,2);
    max_pts = max(cellfun(@numel, x_coords));
    max_pts = max(1, max_pts);

    Inflow_Parameters.inflow_discharge = inflow_discharge;
    Inflow_Parameters.n_stream_gauges = n_stream_gauges;
    Inflow_Parameters.easting_inlet_cells = nan(max_pts, n_stream_gauges);
    Inflow_Parameters.northing_inlet_cells = nan(max_pts, n_stream_gauges);
    Wshed_Properties.n_inlets = zeros(1,n_stream_gauges);

    for g = 1:n_stream_gauges
        [xc, yc, n_eff] = build_cells_from_xy(x_coords{g}, y_coords{g}, flags, GIS_data, Wshed_Properties, DEM_raster);
        Inflow_Parameters.easting_inlet_cells(1:n_eff,g)  = xc(:);
        Inflow_Parameters.northing_inlet_cells(1:n_eff,g) = yc(:);
        Wshed_Properties.n_inlets(:,g) = n_eff;
    end

    Inflow_Parameters.n_stream_obs = numel(time_inflow);
    Inflow_Parameters.time_inflow = time_inflow;
    Inflow_Parameters.time_step_inflow = dt_inflow;
end

function [Stage_Parameters, Wshed_Properties] = build_stage_from_bypass(B, flags, GIS_data, Wshed_Properties, DEM_raster)
    required = {'time_stage','stage','x_coords','y_coords'};
    for i = 1:numel(required)
        if ~isfield(B, required{i})
            error('InputData_Bypass.Stage must contain field "%s".', required{i});
        end
    end

    time_stage = double(B.time_stage(:));
    stage_vals = double(B.stage);

    if size(stage_vals,1) ~= numel(time_stage)
        error('InputData_Bypass.Stage.stage must have the same number of rows as time_stage.');
    end

    valid_rows = isfinite(time_stage) & any(isfinite(stage_vals), 2);
    time_stage = time_stage(valid_rows);
    stage_vals = stage_vals(valid_rows,:);

    stage_vals(~isfinite(stage_vals)) = NaN;

    if numel(time_stage) < 2
        error('InputData_Bypass.Stage must contain at least 2 valid records.');
    end

    dt_stage = time_stage(2) - time_stage(1);
    if dt_stage <= 0
        error('InputData_Bypass.Stage time step must be positive.');
    end

    if any(abs(diff(time_stage) - dt_stage) > 1e-3)
        error('InputData_Bypass.Stage forcing time step must be constant.');
    end

    x_coords = normalize_coord_input(B.x_coords, size(stage_vals,2), 'Stage.x_coords');
    y_coords = normalize_coord_input(B.y_coords, size(stage_vals,2), 'Stage.y_coords');

    n_stage_gauges = size(stage_vals,2);
    max_pts = max(cellfun(@numel, x_coords));
    max_pts = max(1, max_pts);

    Stage_Parameters.stage = stage_vals;
    Stage_Parameters.n_stage_gauges = n_stage_gauges;
    Stage_Parameters.easting_inlet_cells = nan(max_pts, n_stage_gauges);
    Stage_Parameters.northing_inlet_cells = nan(max_pts, n_stage_gauges);
    Wshed_Properties.n_inlets_stage = zeros(1,n_stage_gauges);

    for g = 1:n_stage_gauges
        [xc, yc, n_eff] = build_cells_from_xy(x_coords{g}, y_coords{g}, flags, GIS_data, Wshed_Properties, DEM_raster);
        Stage_Parameters.easting_inlet_cells(1:n_eff,g)  = xc(:);
        Stage_Parameters.northing_inlet_cells(1:n_eff,g) = yc(:);
        Wshed_Properties.n_inlets_stage(:,g) = n_eff;
    end

    Stage_Parameters.n_stage_obs = numel(time_stage);
    Stage_Parameters.time_stage = time_stage;
    Stage_Parameters.time_step_stage = dt_stage;
end

function [xc, yc, n_eff] = build_cells_from_xy(x_in, y_in, flags, GIS_data, Wshed_Properties, DEM_raster)
    if flags.flag_resample
        x_coordinates = round((-DEM_raster.georef.SpatialRef.XWorldLimits(1) + x_in) / GIS_data.resolution_resample);
        y_coordinates = round(( DEM_raster.georef.SpatialRef.YWorldLimits(2) - y_in) / GIS_data.resolution_resample);
    else
        x_coordinates = round((-GIS_data.xulcorner + x_in) / Wshed_Properties.Resolution);
        y_coordinates = round(( GIS_data.yulcorner - y_in) / Wshed_Properties.Resolution);
    end

    points = [x_coordinates(:), y_coordinates(:)];
    effective_inlets = unique(points,'rows','stable');

    xc = effective_inlets(:,1);
    yc = effective_inlets(:,2);
    n_eff = size(effective_inlets,1);
end

function C = normalize_coord_input(v, nExpected, label)
    if iscell(v)
        C = v;
    elseif isnumeric(v)
        if size(v,2) == nExpected
            C = cell(1,nExpected);
            for i = 1:nExpected
                tmp = v(:,i);
                C{i} = tmp(~isnan(tmp));
            end
        elseif size(v,1) == nExpected
            C = cell(1,nExpected);
            for i = 1:nExpected
                tmp = v(i,:);
                C{i} = tmp(~isnan(tmp));
            end
        else
            error('%s must be a cell array or numeric array compatible with the number of gauges.', label);
        end
    else
        error('%s must be a cell array or numeric array.', label);
    end

    if numel(C) ~= nExpected
        error('%s must contain exactly %d entries.', label, nExpected);
    end

    for i = 1:numel(C)
        C{i} = C{i}(:);
    end
end

function R = normalize_reservoir_data(R)
    required = {'index','x_us','y_us','k1','h1','k2','x_ds1','y_ds1','k3','h2','k4','x_ds2','y_ds2'};
    for i = 1:numel(required)
        if ~isfield(R, required{i})
            error('InputData_Bypass.Reservoir_Data must contain field "%s".', required{i});
        end
    end
end

function S = init_empty_inflow()
    S = struct('inflow_discharge',[],'n_stream_gauges',0,'easting_inlet_cells',[], ...
               'northing_inlet_cells',[],'n_stream_obs',0,'time_inflow',[],'time_step_inflow',[]);
end

function S = init_empty_stage()
    S = struct('stage',[],'n_stage_gauges',0,'easting_inlet_cells',[], ...
               'northing_inlet_cells',[],'n_stage_obs',0,'time_stage',[],'time_step_stage',[]);
end

function v = xlget(GD, key)
    S = strings(size(GD));
    for r = 1:size(GD,1)
        for c = 1:size(GD,2)
            if ischar(GD{r,c}) || isstring(GD{r,c})
                S(r,c) = string(GD{r,c});
            end
        end
    end

    [rr, cc] = find(strcmpi(strtrim(S), key), 1, 'first');
    if isempty(rr)
        error("General_Data: key '%s' not found.", key);
    end
    if cc == size(GD,2)
        error("General_Data: key '%s' found in last column; no value to the right.", key);
    end

    v = GD{rr, cc+1};
end

function x = xlnum(GD, key)
    v = xlget(GD, key);
    if isempty(v) || (isstring(v) && strlength(v)==0)
        x = NaN;
        return
    end
    x = double(v);
end

function t = xldatetime(GD, key)
    v = xlget(GD, key);

    if isdatetime(v)
        t = v;
        return
    end

    if isnumeric(v)
        t = datetime(v + datenum('30-Dec-1899'),'ConvertFrom','datenum');
        return
    end

    t = datetime(string(v));
end

function t = cast_to_datetime_local(v)
    if isdatetime(v)
        t = v;
    elseif isnumeric(v)
        t = datetime(v + datenum('30-Dec-1899'),'ConvertFrom','datenum');
    else
        t = datetime(string(v));
    end
end

function flags = read_flags_sheet(FlagsGrid)
    flags = struct();

    nR = size(FlagsGrid,1);
    nC = size(FlagsGrid,2);

    for r = 1:nR
        for c = 1:nC
            keyRaw = FlagsGrid{r,c};

            if iscell(keyRaw) && numel(keyRaw) >= 1
                keyRaw = keyRaw{1};
            end

            try
                keyStr = string(keyRaw);
            catch
                continue
            end

            if numel(keyStr) == 0
                continue
            end

            keyStr = keyStr(1);
            keyStr = strtrim(keyStr);

            if strlength(keyStr) == 0
                continue
            end

            if ~startsWith(lower(keyStr), "flag_")
                continue
            end

            if c == nC
                valRaw = [];
            else
                valRaw = FlagsGrid{r,c+1};
            end

            if iscell(valRaw) && numel(valRaw) >= 1
                valRaw = valRaw{1};
            end

            val = NaN;

            if isempty(valRaw)
                val = NaN;
            elseif ismissing(string(valRaw))
                val = NaN;
            elseif islogical(valRaw)
                val = double(valRaw);
            elseif isnumeric(valRaw)
                if isscalar(valRaw)
                    val = double(valRaw);
                else
                    val = double(valRaw(1));
                end
            else
                s = strtrim(string(valRaw));
                if strlength(s) == 0 || ismissing(s)
                    val = NaN;
                else
                    sLow = lower(s);
                    if any(sLow == ["true","t","yes","y","on"])
                        val = 1;
                    elseif any(sLow == ["false","f","no","n","off"])
                        val = 0;
                    else
                        tmp = str2double(sLow);
                        if ~isnan(tmp)
                            val = tmp;
                        else
                            val = NaN;
                        end
                    end
                end
            end

            fieldName = matlab.lang.makeValidName(char(keyStr));
            flags.(fieldName) = val;
        end
    end

    if ~isfield(flags,'flag_warmup')
        flags.flag_warmup = 0;
    end
end

function T = xlblock_2col(GD, headerText, col1Header, col2Header)
    S = strings(size(GD));
    for r = 1:size(GD,1)
        for c = 1:size(GD,2)
            v = GD{r,c};
            if ischar(v) || isstring(v)
                S(r,c) = string(v);
            end
        end
    end

    [r0,~] = find(strcmpi(strtrim(S), headerText), 1, 'first');
    if isempty(r0), error("General_Data: header '%s' not found.", headerText); end

    rh = r0 + 1;
    c1 = find(strcmpi(strtrim(S(rh,:)), col1Header), 1, 'first');
    c2 = find(strcmpi(strtrim(S(rh,:)), col2Header), 1, 'first');
    if isempty(c1) || isempty(c2)
        error("General_Data: could not find column headers '%s' and '%s' under '%s'.", col1Header, col2Header, headerText);
    end

    rd = rh + 1;
    times = [];
    paths = strings(0,1);

    r = rd;
    while r <= size(GD,1)
        vtime = GD{r,c1};
        vpath = GD{r,c2};

        if isempty(vtime) || (isstring(vtime) && strlength(vtime)==0)
            break
        end

        times(end+1,1) = double(vtime); %#ok<AGROW>

        if ismissing(string(vpath)) || strlength(string(vpath))==0
            paths(end+1,1) = "";
        else
            paths(end+1,1) = string(vpath);
        end

        r = r + 1;
    end

    T = table(times, paths, 'VariableNames', {'time','directory'});
end