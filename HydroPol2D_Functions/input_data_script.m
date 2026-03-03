%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%
GD = readcell(model_folder,'Sheet','General_Data');  % raw grid (names anywhere)

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

% --- Optional: sanity checks (fail fast if a required flag is missing) ---
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

% Volume error
running_control.volume_error = 5; % m3 [not needed and not considerd as input anymore]
running_control.factor_reduction = 2; % [not considered as input anymore]

if flags.flag_input_rainfall_map + flags.flag_satellite_rainfall + flags.flag_real_time_satellite_rainfall > 1
    error('Please choose only one type of spatial rainfall data.')
end

if flags.flag_inertial == 1 && flags.flag_CA == 1
    error('Please, add either diffusive or inertial flag.')
end


% ---------------- Watershed / outlet ----------------
if flags.flag_outlet_type == 1
    outlet_type = 1;
end

slope_outlet   = xlnum(GD,'slope_outlet');

% This key does NOT exist in your current General_Data -> add it or delete this line:
n_outlets_data = 1; % Fixed in one



% Lateral_Groundwater_Flux = input_table_watershed_inputs(3); % m3/s/m

% ---------------- Maps and Plots Control ----------------
running_control.record_time_maps            = xlnum(GD,'record_time_maps');
running_control.record_time_hydrographs     = xlnum(GD,'record_time_hydrographs');
Pol_min                                    = xlnum(GD,'Pol_min');
depth_wse                                   = xlnum(GD,'depth_wse');
flags.flag_wse                               = xlnum(GD,'flag_wse');
running_control.record_time_spatial_rainfall = xlnum(GD,'record_time_spatial_rainfall');
time_save_ETP                                = xlnum(GD,'time_save_ETP');
running_control.record_time_spatial_etp      = xlnum(GD,'record_time_spatial_ETP');

% These keys do NOT exist in your current General_Data -> add them or remove from code:
Krs_ETP = 0.13; % Defaults
albedo  = 0.25; % Defaults


if flags.flag_input_rainfall_map == 1
    if running_control.record_time_maps < running_control.record_time_spatial_rainfall
        error('If the input rainfall maps are provided, It is not possible to save maps with fewer temporal resolution than rainfall')
    end
end

% ---------------- Routing Parameters ----------------
CA_States.depth_tolerance = 1; % Default 1 mm


% River Heigth and Width (Deactivated. Now it is spatially distributed)
% ---------------- River Height and Width ----------------
GIS_data.alfa_1 = xlnum(GD,'alfa_1');
GIS_data.alfa_2 = xlnum(GD,'alfa_2');
GIS_data.beta_1 = xlnum(GD,'beta_1');
GIS_data.beta_2 = xlnum(GD,'beta_2');
LULC_Parameters.River_Manning = xlnum(GD,'Manning');

% This key does NOT exist in your current General_Data -> add it or remove:
River_K_coeff = nan; % Deactivated

% GIS_data.xulcorner = input_table_abstraction(5);
% GIS_data.yulcorner = input_table_abstraction(6);

% ---------------- Water Quality Inputs ----------------
ADD    = xlnum(GD,'ADD');
min_Bt = xlnum(GD,'min_Bt');
Bmin   = xlnum(GD,'Bmin');
Bmax   = xlnum(GD,'Bmax');

% DEM Smoothing, Imposemin & Resample
% ---------------- DEM Smoothing, Imposemin, Resample ----------------
GIS_data.min_area            = xlnum(GD,'min_area');
GIS_data.tau                 = xlnum(GD,'tau');
GIS_data.K_value             = xlnum(GD,'K_value');
GIS_data.sl                  = xlnum(GD,'sl');
GIS_data.resolution_resample = xlnum(GD,'resolution_resample');
GIS_data.slope_DTM           = xlnum(GD,'slope_DTM');

% ---------------- TopoToolbox Folder ----------------
topo_path = xlget(GD,'topo_path');   % <-- rename 'Path' -> 'topo_path' in Excel

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
    unique = unique(temp);
    list={'_cm','_tm','_am','_om','_cf','_tf','_af','_of'};
    names={'Child Boy','Teen Male','Adult Man','Elderly Man',...
        'Child Girl','Teen Female','Adult Women','Elderly Women'};
    names={[strcat(num2str(round(Human_Instability.m_c_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_c_m,1)),'m') newline names{1}],
        [strcat(num2str(round(Human_Instability.m_t_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_t_m,1)),'m') newline names{2}],
        [strcat(num2str(round(Human_Instability.m_a_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_a_m,1)),'m') newline names{3}],
        [strcat(num2str(round(Human_Instability.m_o_m,1)),'Kg', '| ', num2str(round(Human_Instability.y_o_m,1)),'m') newline names{4}],
        [strcat(num2str(round(Human_Instability.m_c_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_c_f,1)),'m') newline names{5}],
        [strcat(num2str(round(Human_Instability.m_t_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_t_f,1)),'m') newline names{6}],
        [strcat(num2str(round(Human_Instability.m_a_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_a_f,1)),'m') newline names{7}],
        [strcat(num2str(round(Human_Instability.m_o_f,1)),'Kg', '| ', num2str(round(Human_Instability.y_o_f,1)),'m') newline names{8}]
        };
    for i = 1:length(unique)
        Human_Instability.order(i) = find(temp==unique(i));
    end
    for i = 1:8
        Human_Instability.list{i} = list{Human_Instability.order(i)};
        Human_Instability.names{i} = names{Human_Instability.order(i)};
    end
else
    Human_Instability = [];
end

% ---------------- Design Storms (from General_Data by label) ----------------
Design_Storm_Parameters.RP              = xlnum(GD,'RP');                % years
Design_Storm_Parameters.Rainfall_Duration = xlnum(GD,'Rainfall Duration'); % minutes (see note below)
Design_Storm_Parameters.K               = xlnum(GD,'K');
Design_Storm_Parameters.a               = xlnum(GD,'a');
Design_Storm_Parameters.b               = xlnum(GD,'b');
Design_Storm_Parameters.c               = xlnum(GD,'c');
Design_Storm_Parameters.time_step       = xlnum(GD,'dt_design');                % minutes

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

    Input_Rainfall.time = T_rain.time;                     % minutes
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

    Input_Transpiration.time = T_tr.time;                       % days
    Input_Transpiration.num_obs_maps = numel(Input_Transpiration.time);
    Input_Transpiration.labels_Directory = cellstr(T_tr.directory);

    Input_Evaporation.time = T_ev.time;                         % days
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
    % flags.flag_input_rainfall_map = 0;
end

if flags.flag_input_ETP_map == 1
    time_maps = 1440; % Only valid for 1day maps
    Input_Transpiration.time = (0:time_maps:running_control.routing_time); % time in minutes
    Input_Transpiration.num_obs_maps = sum(~isnan(Input_Transpiration.time));
    Input_Transpiration.time = Input_Transpiration.time(1:Input_Transpiration.num_obs_maps);
    
    Input_Evaporation.time = (0:time_maps:running_control.routing_time); % time in minutes
    Input_Evaporation.num_obs_maps = sum(~isnan(Input_Evaporation.time));
    Input_Evaporation.time = Input_Evaporation.time(1:Input_Evaporation.num_obs_maps);

    % flags.flag_spatial_rainfall = 1;
    % flags.flag_rainfall = 1;
    % flags.flag_real_time_satellite_rainfall = 0;
    % flags.flag_input_rainfall_map = 0;
end

%%%%%%%%%%%%%% LULC DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable('LULC_parameters.xlsx');
input_data = table2array(input_table(:,2:end)); % numbers
lulc_parameters = input_data(:,2:end);
n_lulc = sum(lulc_parameters(:,1)>=0); % Number of LULC
lulc_parameters = lulc_parameters(1:n_lulc,:);
imp_index = lulc_parameters(1,end);
% LULC Index
LULC_index = input_data(:,1);

%%%%%%%%%%%%%% SOIL DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable('SOIL_parameters.xlsx');
input_data = table2array(input_table(:,2:end)); % numbers
soil_parameters = input_data(:,2:end); % Number of Soil Types
n_soil = sum(soil_parameters(:,1)>=0);
soil_parameters = soil_parameters(1:n_soil,:);
% SOIL Index
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
        path = what;
        path = strcat(path.path,'\','Modeling_Results');
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
        % change the date format
        input_table.(input_table.Properties.VariableNames{name_time})=datetime(input_table.(input_table.Properties.VariableNames{name_time}),'InputFormat','yyyy-MM-dd HH:mm:ss');
        column_time = find(strcmpi(input_table.Properties.VariableNames, input_table.Properties.VariableNames{name_time}));
        column_rain = column_time + 1;
        % Gathering data
        temp_rain = table2array(input_table(:,column_rain));
        temp_time = table2array(input_table(:,column_time));
        Rainfall_Parameters.intensity_rainfall = temp_rain(~isnat(temp_time));
        Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
        Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
        Rainfall_Parameters.time_rainfall = input_table.(input_table.Properties.VariableNames{name_time})(~isnat(temp_time)); % min
        Rainfall_Parameters.time_step_rainfall = minutes(Rainfall_Parameters.time_rainfall(2) - Rainfall_Parameters.time_rainfall(1)); % min
        Rainfall_Parameters.rainfall_duration = minutes(Rainfall_Parameters.time_rainfall(2) - Rainfall_Parameters.time_rainfall(1)); % min
        date_begin = Rainfall_Parameters.time_rainfall(1);
        date_end = Rainfall_Parameters.time_rainfall(end);
        Rainfall_Parameters.time_rainfall = minutes(Rainfall_Parameters.time_rainfall(:)-Rainfall_Parameters.time_rainfall(1));

        % --- Routing Time --- %
        % sum(~isnat(temp_time))>sum(~isnan(temp_rain))
        running_control.routing_time = (Rainfall_Parameters.n_obs_rainfall-1)*Rainfall_Parameters.time_step_rainfall; % Minutes
    end
elseif flags.flag_alternated_blocks == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
    % Run Alternated Blocks Model
    [Rainfall_Parameters.time_rainfall,Rainfall_Parameters.intensity_rainfall,~,idf] = alternated_blocks(Design_Storm_Parameters.Rainfall_Duration(end),Design_Storm_Parameters.time_step,Design_Storm_Parameters.K,Design_Storm_Parameters.a,Design_Storm_Parameters.b,Design_Storm_Parameters.c,Design_Storm_Parameters.RP,1);
    Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
    Rainfall_Parameters.time_step_rainfall = Design_Storm_Parameters.time_step; % min
    Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end); % min
    % Routing Time >= Rainfall Duration
    running_control.routing_time = max(running_control.routing_time,Rainfall_Parameters.rainfall_duration);
elseif flags.flag_huff == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_spatial_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1
    % Run Huff Model
    [Rainfall_Parameters.time_rainfall,Rainfall_Parameters.intensity_rainfall,~] = Huff_Curves(Design_Storm_Parameters.Rainfall_Duration(end), Design_Storm_Parameters.time_step, Design_Storm_Parameters.K,Design_Storm_Parameters.a,Design_Storm_Parameters.b,Design_Storm_Parameters.c,Design_Storm_Parameters.RP,1);
    Rainfall_Parameters.n_obs_rainfall = sum(Rainfall_Parameters.intensity_rainfall >= 0);
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall(1:Rainfall_Parameters.n_obs_rainfall);
    Rainfall_Parameters.time_step_rainfall = Design_Storm_Parameters.time_step; % min
    Rainfall_Parameters.rainfall_duration = Rainfall_Parameters.time_rainfall(end); % min
    % Routing Time >= Rainfall Duration
    running_control.routing_time = max(running_control.routing_time,Rainfall_Parameters.rainfall_duration);
end

% Load TopoToolBox Tools
addpath(genpath(char(topo_path)));

%%% ---- Inflow Hydrograph ---- %%%
    input_table = readtable('Inflow_Hydrograph.xlsx');
    input_table_labels = input_table(1:2,:);
    input_table_values = (input_table(3:end,:));
    % --- Convert Everything to String
    input_table_values = convertvars(input_table_values,[1:1:size(input_table_values,2)],'string');
    % --- Converting Everything to Double
    input_table_values = convertvars(input_table_values,[1:1:size(input_table_values,2)],'double');

    n_inflows_max = 50;
    for i = 1:n_inflows_max
        Inflow_Parameters.inflow_discharge(:,i) = table2array(input_table(3:end,(i-1)*5 + 2));
    end
    Inflow_Parameters.n_stream_gauges = sum(sum(Inflow_Parameters.inflow_discharge(1,:)>=0)); % Number of stream gauges

    %Prealocating array for relative coordinates of inlets
    Inflow_Parameters.easting_inlet_cells = nan(n_inflows_max,n_inflows_max);
    Inflow_Parameters.northing_inlet_cells = nan(n_inflows_max,n_inflows_max);

    for i = 1:Inflow_Parameters.n_stream_gauges
        if flags.flag_resample
            % Do we have to include Resolution/2 in the calculations?
            nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
            x_coordinates = round((-DEM_raster.georef.SpatialRef.XWorldLimits(1) + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/GIS_data.resolution_resample); % Check R/2
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
            northing_inlet_cells(1:size(effective_inlets,1),i) = effective_inlets(:,2);
            Inflow_Parameters.n_inlets(:,i) = size(effective_inlets,1);
        else
            % Do we have to include Resolution/2 in the calculations?
            nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
            x_coordinates = round((-GIS_data.xulcorner  + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/Wshed_Properties.Resolution);
            y_coordinates = round((GIS_data.yulcorner - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/Wshed_Properties.Resolution); % Check R/2
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
    Inflow_Parameters.n_stream_obs = sum(sum(Inflow_Parameters.inflow_discharge(:,1)>=0)); % Number of stream data per gauge
    Inflow_Parameters.time_inflow = table2array(input_table(3:Inflow_Parameters.n_stream_obs,1)); % min
    Inflow_Parameters.time_step_inflow = (table2array(input_table(4,1)) - table2array(input_table(3,1))); % min

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
    input_table_labels = input_table(1:2,:);
    input_table_values = (input_table(3:end,:));
    % --- Convert Everything to String
    input_table_values = convertvars(input_table_values,[1:1:size(input_table_values,2)],'string');
    % --- Converting Everything to Double
    input_table_values = convertvars(input_table_values,[1:1:size(input_table_values,2)],'double');

    n_stage_max = 50;
    for i = 1:n_stage_max
        Stage_Parameters.stage(:,i) = table2array(input_table(3:end,(i-1)*5 + 2));
    end
    Stage_Parameters.n_stage_gauges = sum(sum(Stage_Parameters.stage(1,:)>=0)); % Number of stream gauges

    %Prealocating array for relative coordinates of inlets
    Stage_Parameters.easting_inlet_cells = nan(n_stage_max,n_stage_max);
    Stage_Parameters.northing_inlet_cells = nan(n_stage_max,n_stage_max);

    for i = 1:Stage_Parameters.n_stage_gauges
        if flags.flag_resample
            % Do we have to include Resolution/2 in the calculations?
            nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
            x_coordinates = round((-DEM_raster.georef.SpatialRef.XWorldLimits(1) + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/GIS_data.resolution_resample); % Check R/2
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
            % Do we have to include Resolution/2 in the calculations?
            nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
            x_coordinates = round((-GIS_data.xulcorner  + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/Wshed_Properties.Resolution);
            y_coordinates = round((GIS_data.yulcorner - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/Wshed_Properties.Resolution); % Check R/2
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
    Stage_Parameters.n_stage_obs = sum(sum(Stage_Parameters.stage(:,1)>=0)); % Number of stream data per gauge
    Stage_Parameters.time_stage = table2array(input_table(3:Stage_Parameters.n_stage_obs,1)); % min
    Stage_Parameters.time_step_stage = (table2array(input_table(4,1)) - table2array(input_table(3,1))); % min
%  Entering a boundary condition for rainfall if design storms are used
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


function v = xlget(GD, key)
%XLGET Find a label (key) anywhere in a readcell() grid and return the cell to the right.
%   v = xlget(GD,'time_step_model') returns the value in the cell right of the key.
%
% Notes:
% - Case-insensitive
% - Errors if key not found
% - If multiple matches exist, returns the first match (top-to-bottom, left-to-right)

    S = strings(size(GD));
    for r = 1:size(GD,1)
        for c = 1:size(GD,2)
            if ischar(GD{r,c}) || isstring(GD{r,c})
                S(r,c) = string(GD{r,c});
            end
        end
    end

    [rr, cc] = find(strcmpi(S, key), 1, 'first');
    if isempty(rr)
        error("General_Data: key '%s' not found.", key);
    end

    if cc == size(GD,2)
        error("General_Data: key '%s' found in last column; no value to the right.", key);
    end

    v = GD{rr, cc+1};
end


function x = xlnum(GD, key)
%XLNUM numeric version of xlget
    v = xlget(GD, key);
    if isempty(v) || (isstring(v) && strlength(v)==0)
        x = NaN;
        return
    end
    x = double(v);
end


function t = xldatetime(GD, key)
%XLDATETIME datetime version of xlget (handles Excel serials or datetimes/strings)
    v = xlget(GD, key);

    if isdatetime(v)
        t = v;
        return
    end

    if isnumeric(v)
        % Excel serial date number -> MATLAB datetime (Excel origin: 30-Dec-1899)
        t = datetime(v + datenum('30-Dec-1899'),'ConvertFrom','datenum');
        return
    end

    % string parse fallback
    t = datetime(string(v));
end

function flags = read_flags_sheet(FlagsGrid)
%READ_FLAGS_SHEET Builds flags struct by scanning for any cell 'flag_*'
% and reading the value in the cell to the right.
%
% Supports:
% - Flags sheet is a readcell() grid
% - flag names anywhere in the grid (case-insensitive)
% - value is assumed to be in the cell immediately to the right
% - values can be numeric, logical, char, string, empty, missing

    flags = struct();

    nR = size(FlagsGrid,1);
    nC = size(FlagsGrid,2);

    for r = 1:nR
        for c = 1:nC

            % ---- read potential key cell safely ----
            keyRaw = FlagsGrid{r,c};

            % unwrap cell-in-cell
            if iscell(keyRaw) && numel(keyRaw) >= 1
                keyRaw = keyRaw{1};
            end

            % convert to string for checking (safe)
            try
                keyStr = string(keyRaw);
            catch
                continue
            end

            if numel(keyStr) == 0
                continue
            end
            keyStr = keyStr(1);                 % scalar
            keyStr = strtrim(keyStr);

            if strlength(keyStr) == 0
                continue
            end

            % only accept keys like "flag_*"
            if ~startsWith(lower(keyStr), "flag_")
                continue
            end

            % ---- value is the cell to the right ----
            if c == nC
                % flag in last column: no value cell
                valRaw = [];
            else
                valRaw = FlagsGrid{r,c+1};
            end

            % unwrap cell-in-cell
            if iscell(valRaw) && numel(valRaw) >= 1
                valRaw = valRaw{1};
            end

            % ---- normalize value to numeric (0/1) or NaN ----
            val = NaN;

            if isempty(valRaw)
                val = NaN;

            elseif ismissing(string(valRaw))  % handles <missing>, NaT, etc
                val = NaN;

            elseif islogical(valRaw)
                val = double(valRaw);

            elseif isnumeric(valRaw)
                if isscalar(valRaw)
                    val = double(valRaw);
                else
                    % if someone pasted a vector by mistake, take first
                    val = double(valRaw(1));
                end

            else
                % char/string like '1', '0', 'true', 'false', 'yes', 'no'
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

            % store field (make valid)
            fieldName = matlab.lang.makeValidName(char(keyStr));
            flags.(fieldName) = val;

        end
    end

    % Optional default (prevents crash if user forgot it)
    if ~isfield(flags,'flag_warmup')
        flags.flag_warmup = 0;
    end
end

function T = xlblock_2col(GD, headerText, col1Header, col2Header)
%XLBLOCK_2COL Finds a 2-column block under a header and returns it as a table.
%   headerText: the title above the block (e.g., 'Sattelite or Radar Rainfall')
%   col1Header/col2Header: the header texts in the first row of the block.

    % make a string grid for searching
    S = strings(size(GD));
    for r = 1:size(GD,1)
        for c = 1:size(GD,2)
            v = GD{r,c};
            if ischar(v) || isstring(v)
                S(r,c) = string(v);
            end
        end
    end

    [r0,c0] = find(strcmpi(strtrim(S), headerText), 1, 'first');
    if isempty(r0), error("General_Data: header '%s' not found.", headerText); end

    % find the header row immediately below
    rh = r0 + 1;

    % find column indices for the two columns (same row)
    c1 = find(strcmpi(strtrim(S(rh,:)), col1Header), 1, 'first');
    c2 = find(strcmpi(strtrim(S(rh,:)), col2Header), 1, 'first');
    if isempty(c1) || isempty(c2)
        error("General_Data: could not find column headers '%s' and '%s' under '%s'.", col1Header, col2Header, headerText);
    end

    % data start row
    rd = rh + 1;

    % read down until first column is empty
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