%%%%%%%%%%%%%% INPUT DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable(model_folder);

% Running Control
input_table_running_control = table2array(input_table(:,2));
time_step_model = input_table_running_control(1)/60; % Dividing to convert to min
running_control.min_time_step = input_table_running_control(2);
running_control.max_time_step = input_table_running_control(3);
running_control.time_step_increments = input_table_running_control(4);
running_control.time_step_change = input_table_running_control(5);
Courant_Parameters.alfa_max = input_table_running_control(6);
Courant_Parameters.alfa_min = input_table_running_control(7);
Courant_Parameters.v_threshold = input_table_running_control(8);
Courant_Parameters.slope_alfa = input_table_running_control(9);
date_begin = input_table_running_control(10);
date_end = input_table_running_control(11);
% --- Routing Time --- %
running_control.routing_time = round((date_end - date_begin)*1440,4); % Minutes
date_begin = datetime(datestr(date_begin+datenum('30-Dec-1899')));
date_end = datetime(datestr(date_end+datenum('30-Dec-1899')));
if running_control.routing_time < 0
    error('Please make sure Date End is later than Date Begin.')
end

% Flags
input_table_flags = table2array(input_table(:,5));
flags.flag_rainfall = input_table_flags(1);
flags.flag_abstraction = input_table_flags(2);
flags.flag_inflow = input_table_flags(3);
flags.flag_waterbalance = input_table_flags(4);
flags.flag_waterquality = input_table_flags(5);
flags.flag_timestep = input_table_flags(6);
flags.flag_warmup = input_table_flags(7);
flags.flag_initial_buildup = input_table_flags(8);
flags.flag_wq_model = input_table_flags(9);
flags.flag_infiltration = input_table_flags(10);
flags.flag_critical = input_table_flags(11);
flags.flag_spatial_rainfall	 = input_table_flags(12);
flags.flag_D8	= input_table_flags(13);
flags.flag_diffusive = input_table_flags(14);
flags.flag_resample = input_table_flags(15);
flags.flag_smoothening = input_table_flags(16);
flags.flag_trunk = input_table_flags(17);
flags.flag_export_maps = input_table_flags(18);
flags.flag_fill_DEM = input_table_flags(19);
flags.flag_smooth_cells = input_table_flags(20);
flags.flag_reduce_DEM = input_table_flags(21);
flags.flag_ETP = input_table_flags(22);
flags.flag_obs_gauges = input_table_flags(23);
flags.flag_GPU = input_table_flags(24);
flags.flag_single = input_table_flags(25);
flags.flag_reservoir = input_table_flags(26);
flags.flag_human_instability = input_table_flags(27);
flags.flag_input_rainfall_map = input_table_flags(28);
flags.flag_satellite_rainfall = input_table_flags(29);
flags.flag_real_time_satellite_rainfall = input_table_flags(30);
flags.flag_unit_for_forecasting = input_table_flags(31);
flags.flag_dam_break = input_table_flags(32);
flags.flag_groundwater_modeling  = input_table_flags(33);
flags.flag_river_heigth_compensate = input_table_flags(34);
flags.flag_rainfall_multiple_runs = input_table_flags(35);
flags.flag_inertial = input_table_flags(36);
% Matricial Variables
% input_table_matricial = table2array(input_table(:,9));
% running_control.time_step_matrices = input_table_matricial(1);
% factor_cells = input_table_matricial(2);

if flags.flag_input_rainfall_map + flags.flag_satellite_rainfall + flags.flag_real_time_satellite_rainfall > 1
    error('Please choose only one type of spatial rainfall data.')
end

if flags.flag_inertial == 1 && flags.flag_diffusive == 1
    error('Please, add either diffusive or inertial flag.')
end

% Watershed Inputs and Cuts
input_table_watershed_inputs = table2array(input_table(:,12));
outlet_type = input_table_watershed_inputs(1);
slope_outlet = input_table_watershed_inputs(2);
n_outlets_data = input_table_watershed_inputs(3);
Lateral_Groundwater_Flux = input_table_watershed_inputs(4); % m3/s/m

% Maps and Plots Control
input_table_map_plots = table2array(input_table(:,15));
running_control.record_time_maps = input_table_map_plots(1);
running_control.record_time_hydrographs = input_table_map_plots(2);
Pol_min = input_table_map_plots(3);
depth_wse = input_table_map_plots(4);
flags.flag_wse = input_table_map_plots(5);
flags.flag_elapsed_time = input_table_map_plots(6);
running_control.record_time_spatial_rainfall = input_table_map_plots(7);
time_save_ETP = input_table_map_plots(8);

% Routing Parameters
routing_parameters = table2array(input_table(:,18));
CA_States.depth_tolerance = routing_parameters(1);


% River Heigth and Width
input_table_river = table2array(input_table(:,21));
GIS_data.alfa_1 = input_table_river(1);
GIS_data.alfa_2 = input_table_river(2);
GIS_data.beta_1 = input_table_river(3);
GIS_data.beta_2 = input_table_river(4);

% GIS_data.xulcorner = input_table_abstraction(5);
% GIS_data.yulcorner = input_table_abstraction(6);

% Water Quality Inputs
input_table_WQ_parameter = table2array(input_table(:,24));
ADD = input_table_WQ_parameter(1);
min_Bt = input_table_WQ_parameter(2);
Bmin = input_table_WQ_parameter(3);
Bmax = input_table_WQ_parameter(4);

% DEM Smoothing, Imposemin & Resample
input_table_DEM = table2array(input_table(:,28));
GIS_data.min_area = input_table_DEM(1);
GIS_data.tau = input_table_DEM(2);
GIS_data.K_value = input_table_DEM(3);
GIS_data.sl = input_table_DEM(4);
GIS_data.resolution_resample = input_table_DEM(5);


% TopoToolbox Folder
topo_path = table2cell(input_table(1,31));

% Human Instability
if flags.flag_human_instability == 1
    human_table = readtable('human_risk.xlsx');
    Human_Instability.mu = table2array(human_table(4,2)); % friction coefficient
    Human_Instability.Cd = table2array(human_table(3,2));
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

% Design Storms
input_table_design = table2array((input_table(1:9,43)));
flags.flag_alternated_blocks = input_table_design(1);
flags.flag_huff = input_table_design(2);
Design_Storm_Parameters.RP = input_table_design(3); % year
Design_Storm_Parameters.Rainfall_Duration = input_table_design(4); % min
Design_Storm_Parameters.K = input_table_design(5); % K
Design_Storm_Parameters.a = input_table_design(6); % a
Design_Storm_Parameters.b = input_table_design(7); % b
Design_Storm_Parameters.c = input_table_design(8); % c
Design_Storm_Parameters.time_step = input_table_design(9); % min

if flags.flag_huff == 1 && flags.flag_alternated_blocks == 1
    error('Please, enter either Alternated Blocks or Huff hyetograph.')
end

% Input Rainfall Maps
if flags.flag_input_rainfall_map == 1
    input_table_rainfall = ((input_table(2:end,45:46)));
    Input_Rainfall.time = table2array(input_table_rainfall(:,1)); % % min
    Input_Rainfall.num_obs_maps = sum(~isnan(Input_Rainfall.time));
    Input_Rainfall.time = Input_Rainfall.time(1:Input_Rainfall.num_obs_maps);
    input_table_dir = input_table_rainfall(:,2);

    for i = 1:Input_Rainfall.num_obs_maps
        Input_Rainfall.labels_Directory{i,:} = input_table_dir{i,:};
    end
    flags.flag_spatial_rainfall = 1;
    flags.flag_rainfall = 1;
    flags.flag_real_time_satellite_rainfall = 0;
    flags.flag_satellite_rainfall = 0;
    % flags.flag_input_rainfall_map = 0;
else
    Input_Rainfall = [];
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

% Inflow
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
easting_inlet_cells = nan(n_inflows_max,n_inflows_max);
northing_inlet_cells = nan(n_inflows_max,n_inflows_max);

for i = 1:Inflow_Parameters.n_stream_gauges
    if flags.flag_resample
        % Do we have to include Resolution/2 in the calculations?
        nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
        easting_inlet_cells(1:nonNanCount,i) = round((-GIS_data.xulcorner + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/GIS_data.resolution_resample);
        northing_inlet_cells(1:nonNanCount,i) = round((DEM_raster.georef.SpatialRef.YWorldLimits(2) - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/GIS_data.resolution_resample);
        Wshed_Properties.n_inlets(:,i) = sum(easting_inlet_cells(1:nonNanCount,i)>=0);
    else
        % Do we have to include Resolution/2 in the calculations?
        nonNanCount = sum(~ismissing(input_table_values(1:end,(i-1)*5 + 3)));
        easting_inlet_cells(1:nonNanCount,i) = round((-GIS_data.xulcorner + table2array(input_table_values(1:nonNanCount,(i-1)*5 + 3)))/Wshed_Properties.Resolution);
        northing_inlet_cells(1:nonNanCount,i) = round((GIS_data.yulcorner - table2array(input_table_values(1:nonNanCount,(i-1)*5 + 4)))/Wshed_Properties.Resolution);
        Wshed_Properties.n_inlets(:,i) = sum(sum(easting_inlet_cells(1:nonNanCount,i)>=0));
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
    input_table = readtable('Reservoir_Data_ori.xlsx');
    Reservoir_Data.index = table2array(input_table(:,1));
    Reservoir_Data.x = table2array(input_table(:,2));
    Reservoir_Data.y = table2array(input_table(:,3));
    Reservoir_Data.Lef = table2array(input_table(:,4));
	Reservoir_Data.Hv = table2array(input_table(:,5));
	Reservoir_Data.Cd = table2array(input_table(:,6));
	Reservoir_Data.Aculv = table2array(input_table(:,7));
    Reservoir_Data.Ho = table2array(input_table(:,8));
    Reservoir_Data.Co = table2array(input_table(:,9));
    Reservoir_Data.Dir = table2array(input_table(:,10));
    Reservoir_Data.Qreg = table2array(input_table(:,11));
    % Calculating Kv = Cd * Lef and orifice parameters
    Reservoir_Data.Kv = Reservoir_Data.Cd.*Reservoir_Data.Lef;
	Reservoir_Data.Ko = Reservoir_Data.Co.*sqrt(2*9.81*Reservoir_Data.Aculv);
end

clear input_table