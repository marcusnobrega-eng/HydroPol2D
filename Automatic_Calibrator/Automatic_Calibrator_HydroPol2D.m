%% HydroPol2D Solver

% We have 9 variables in total, to calibrate
% Let nl be the number of LULC and ns the number of soil classifications
% Therefore, we have 9 * nl + 9 * ns decision variables
% This variable collects all these data following:

% x = [n1, n2 ... n_{nl} ...
%      h0_1, h_0,2 ... h_0_{nl}
%      C1_1, C1_2 ... C1_{nl}
%      C2_1, C2_2 ... C2_{nl}
%      C3_1, C3_2 ... C3_{nl}
%      C4_1, C4_2 ... C4_{nl}
%      ksat1, ksat2 ... ksat{ns}
%      dtheta1, dtheta2 ... dtheta{ns}
%      ps1, psi2 ... psi{ns}]^T

%% Load HydroPol2D Preprocessing Data
clear all; clc;
model_folder = 'General_Data_HydroPol2D.xlsx';
input_table = readtable(model_folder); 

% Load Model Functions
HydroPol2D_tools = char(table2cell(input_table(9,31)));
addpath(genpath(char(HydroPol2D_tools)));

% Pre-Processing
HydroPol2D_preprocessing
save('HydroPol2D_preprocessing_input_data.mat')


%% Defining Minimum and Maximum Constraints
clear all
clc
load('HydroPol2D_preprocessing_input_data.mat') % Loading only important data

% Automatic Calibrator Spreadsheet
calibrator_sheet = 'HydroPol2D_Automatic_Calibrator_Data_3.xlsx';

% Reading Parameters
parameters_data = readtable(calibrator_sheet,'Sheet','Parameters');
n_LULC = table2array(parameters_data(1,2));
n_SOIL = table2array(parameters_data(3,2));

% LULC-Based Parameters
if n_LULC == 1
    parameters_data.Var9 = str2double(parameters_data.Var9);
end
lulc_based_parameters = table2array(parameters_data(1:(1 + n_LULC- 1),5:16));
% ---- Parameter Ranges ----- %
n_range = lulc_based_parameters(:,1:2);
h_0_range = lulc_based_parameters(:,3:4);
C_1_range = lulc_based_parameters(:,5:6);
C_2_range = lulc_based_parameters(:,7:8);
C3_range = lulc_based_parameters(:,9:10);
C4_range = lulc_based_parameters(:,11:12);

% SOIL-Based Parameters
soil_based_parameters = table2array(parameters_data(1:(1 + n_SOIL- 1),19:24));
ksat_range = soil_based_parameters(:,1:2);
dtheta_range = soil_based_parameters(:,3:4);
psi_range = soil_based_parameters(:,5:6);


%% Load Initial Maps and Input Data
addpath(genpath(char(topo_path)));

% Reading Initial Maps and Input Data
initial_maps_data = readtable(calibrator_sheet,'Sheet','Initial_Maps');
flags.flag_calibrate_wq =  table2array(initial_maps_data(1,5));
n_events =  table2array(initial_maps_data(1,6));
warmup_depths_file = initial_maps_data(1:(1 + n_events - 1),2);
warmup_I_0_file = initial_maps_data(1:(1 + n_events - 1),3);
number_of_gauges = table2array(initial_maps_data(1:(1 + n_events - 1),4));

% Deciding which events will be used
% gauges_used = ones(1,size(number_of_gauges))
gauges_used = [1 0 1]; % ones for the events considered in the objective functions

% Reading Input Data
observed_data = readtable(calibrator_sheet,'Sheet','Observed_Data');

for i = 1:n_events
    % Reading other input
    time_observed(:,i) = table2array(observed_data(4:end,(i-1)*22 + 1));
    time_step_rainfall(:,i) = table2array(observed_data(5,(i-1)*22 + 1)) - table2array(observed_data(4,(i-1)*22 + 1));
    intensity_rainfall_observed(:,i) = table2array(observed_data(4:end,(i-1)*22 + 2));
    for j = 1:number_of_gauges(i)
        observed_flow(:,j,i) = table2array(observed_data(4:end,(i-1)*2 + (3 + (j-1)*2) +(i-1)*20)); % page is the event
        pollutant_concentration(:,j,i) = table2array(observed_data(4:end,(i-1)*2+ (4 + (j-1)*2) + +(i-1)*20)); % page is the event
        x_observed(:,j,i) = table2array(observed_data(1,(i-1)*2+ (4 + (j-1)*2) + +(i-1)*20));
        y_observed(:,j,i) = table2array(observed_data(2,(i-1)*2+ (4 + (j-1)*2) + +(i-1)*20));        
    end
    ADD_events(1,i) = table2array(observed_data(2,(i-1)*22 + 2));
end

% Observed Gauges
% Raster Extent
xulcorner = DEM_raster.refmat(3,1); % Up Left Corner
yulcorner = DEM_raster.refmat(3,2);

% --- Converting coordinates to local coordinates in pixels
easting_observed_cells = round((-xulcorner + x_observed)/Wshed_Properties.Resolution);
northing_observed_cells = round((yulcorner - y_observed)/Wshed_Properties.Resolution);


for i = 1:n_events
% Extent Problem   
    if flags.flag_resample ~= 1
        % Extent Issues
        zzz = GRIDobj(string(warmup_depths_file{i,:}));
        zzz = (resample(zzz,DEM_raster));
        warmup_depths(:,:,i) = double(zzz.Z)*1000; % mm 
        yyy = GRIDobj(string(warmup_I_0_file{i,:}));
        yyy = resample(yyy,DEM_raster);
        warmup_I_0(:,:,i) = double(yyy.Z); % mm       
    else
        resolution = resolution_resample; % m
        % Resample
        zzz = (GRIDobj(string(warmup_depths_file{i,:})));
        zzz = resample(zzz,resolution);
        warmup_depths(:,:,i) = double(zzz.Z);

        zzz = (GRIDobj(string(warmup_I_0_file{i,:})));
        zzz = resample(zzz,resolution);
        warmup_I_0(:,:,i) = double(zzz.Z);        
    end
end

for ii = 1:n_events
    % Finding the size of delta_p_obs
    n_observations_max = sum(~isnan(time_observed(:,ii)));
    steps_max_events(ii) = ceil(time_observed(n_observations_max)/time_step_model);
end

steps_max = max(steps_max_events); 
delta_p_obs = NaN(n_events,steps_max);

for ii = 1:n_events
    % Calculation of delta_p from Rainfall Data

    n_observations = sum(~isnan(time_observed(:,ii)));


    % Simulation Time    
%     routing_time = time_observed(n_observations,ii);   
    routing_time = steps_max_events(ii)*time_step_model;
    steps = routing_time/time_step_model; % number of calculation steps

    % Rainfall
    intensity_rainfall = intensity_rainfall_observed(1:n_observations,ii); % mm/h
    intensity_rainfall = intensity_rainfall*flags.flag_rainfall;
    % Conversion of rainfall into the time-step of calculations for concentrated rainfall
    if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~=1  % Only for concentrated rainfall
        z2 = 0;
        intensity_rainfall_length = length(intensity_rainfall) - 1; % Number of intervals        
        intensity_discretized = zeros(1,steps_max); % Preallocating
        for i =1:steps
            time = i*time_step_model; % min
            z1 = find(time_observed(:,ii) <= time,1,'last');
            if z1 > z2
                intensity_discretized(i) = intensity_rainfall(z1);
            else
                intensity_discretized(i) = intensity_rainfall(z2);
            end
            z2 = z1;
        end
    end
    [~,delta_p_obs(ii,1:steps),~] = accumulated_incremental(steps,intensity_discretized,time_step_model);
    clear intensity_discretized
end


% [Obj_fun] = Optimization_Function(x,delta_p_obs,warmup_depths,warmup_I_0,n_events,time_observed,time_step_rainfall,intensity_rainfall_observed,observed_flow,pollutant_concentration,ADD_events,time_obs , alfa_albedo_input , alfa_max , alfa_min , alfa_save , avgtemp_stations , B_t , Bmax , Bmin , C , C_3 , C_4 , Cd , cell_area , climatologic_spatial_duration , col_outlet , coordinate_x , coordinate_y , coordinates_stations , d , d_0 , d_p , d_t , date_begin , delta_p , delta_p_agg , dem , DEM_etp , DEM_raster , depth_tolerance , dmax_final , drainage_area , easting_obs_gauges , elevation , elevation_down_t , elevation_left_t , elevation_right_t , elevation_up_t , EMC_outlet , ETP , ETP_save , factor_cells , flags.flag_critical , flags.flag_D8 , flags.flag_ETP , flags.flag_infiltration , flags.flag_inflow , flags.flag_obs_gauges , flags.flag_rainfall , flags.flag_spatial_rainfall , flags.flag_timestep , flags.flag_waterbalance , flags.flag_waterquality , flags.flag_wq_model , flow_acceleration , flow_tolerance , flows_cells , G_stations , gravity , h_0 , h_0_fulldomain , I_0 , I_p , I_t , I_tot_end , I_tot_end_cell, idx_nan , idx_nan_5 , idx_not_nan , idx_outlet , inflow , inflow_cells , inflow_volume , k , k_out , Krs , ksat , ksat_fulldomain , last_record_hydrograph , last_record_maps , lat , mass_lost , mass_outlet , max_time_step , maxtemp_stations , min_Bt , min_time_step , mintemp_stations , mu , n_inlets , n_stream_gauges , northing_obs_gauges , num_obs_gauges , nx_max , ny_max , Out_Conc , outet_hydrograph , outflow_volume , outlet_index , outlet_index_fulldomain , outlet_runoff_volume , outlet_type , P_conc , psi , psi_fulldomain , rainfall_matrix , rainfall_matrix_full_domain , relative_vol_error , Resolution , risk , ro_water , roughness , roughness_fulldomain , routing_time , row_outlet , slope_alfa , slope_outlet , spatial_domain , t , t_previous , teta_i , teta_i_fulldomain , teta_sat , teta_sat_fulldomain , time_calculation_routing , time_change_matrices , time_change_records , time_deltap , time_ETP , time_hydrograph , time_record_hydrograph , time_records , time_save_previous , time_step , time_step_change , time_step_increments , time_step_model , time_step_save , tmin_wq , Tot_Washed , Tr , u2_stations , ur_stations , v_threshold , vel_down , vel_left , vel_right , vel_up , vol_outlet , weight_person, width1_person, width2_person);


%% Genetic Algorithm Optimization Problem
pop_size = 100;
generations = 10;
stall_generations = 100;
tolerance = 1e-6;
stall_limit = 720*60; % If no improvement is found in this time in seconds, the algorithm stops.

% Number of variables
nvars = 6*n_LULC + 3*n_SOIL;

% Minimum and Maximum Ranges
% Reading Parameters
% parameters_data = readtable('HydroPol2D_Automatic_Calibrator_Data.xlsx','Sheet','Parameters');
minmax_parameters = table2array(parameters_data(3:(3 + n_LULC-1),5:16));
for i = 1:6
    min_parameters_LULC(:,i) = minmax_parameters(1:end,(i-1)*2 + 1);
    max_parameters_LULC(:,i) = minmax_parameters(1:end,(i-1)*2 + 2);
end

minmax_parameters = table2array(parameters_data(3:(3 + n_SOIL - 1),19:24));
for i = 1:3
    min_parameters_SOIL(:,i) = minmax_parameters(1:end,(i-1)*2 + 1);
    max_parameters_SOIL(:,i) = minmax_parameters(1:end,(i-1)*2 + 2);
end

lb = [min_parameters_LULC(:); min_parameters_SOIL(:)]; % Lower bound
ub = [max_parameters_LULC(:) ; max_parameters_SOIL(:)]; % Upper bound

A = []; % Inequality constraint given by A*x <= b
b = []; 
Aeq = []; % Equality constraint given by Aeq*x = beq
beq = []; % 
nonlcon = []; % Non-linear constraints of x (function of x)
intcon = []; % Integer only constraints of x
options = optimoptions('ga','PlotFcns', {@gaplotbestf,@gaplotstopping,@gaplotscores,@gaplotscorediversity}, ...
    'OutputFcn',{@ga_save_each_gen},'PopulationSize',pop_size,'Generations',generations,'StallTimeLimit',stall_limit,'MaxStallGenerations',stall_generations,'FunctionTolerance',tolerance);

% %% Boundary Condition in the Observed Flow
% for i = 1:size(observed_flow,3)
%     for j = 1:size(observed_flow,2)
%         zzz = observed_flow(:,j,i);
%         if sum(zzz) == 0
%             zzz(:) = 1e-6; % Numerical Boundary condition for the fitness functions that have divisions
%             observed_flow(:,j,i) = zzz;
%         end
%     end
% end

%% Converting All Automatic Calibration Data to a Struct
Automatic_Calibrator.time_observed = time_observed;
Automatic_Calibrator.observed_flow = observed_flow;
Automatic_Calibrator.observed_flow = observed_flow;
Automatic_Calibrator.pollutant_concentration = pollutant_concentration;
Automatic_Calibrator.warmup_I_0 = warmup_I_0;
Automatic_Calibrator.warmup_depths = warmup_depths;
Automatic_Calibrator.n_events = n_events;
Automatic_Calibrator.n_LULC = n_LULC;
Automatic_Calibrator.n_SOIL = n_SOIL;
Automatic_Calibrator.easting_observed_cells = easting_observed_cells;
Automatic_Calibrator.northing_observed_cells = northing_observed_cells;
Automatic_Calibrator.observed_flow = observed_flow;
Automatic_Calibrator.delta_p_obs = delta_p_obs;

% Ojective Function
fun = @(x)Objective_Function_HydroPol2D(x,min_soil_moisture,Automatic_Calibrator,flags.flag_calibrate_wq,Reservoir_Data,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties);
% GA problem
[x,fval,exitflag,output,population,scores] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);

save('Workspace_Optimization.mat');


% % Observed
% scatter(time_observed,pollutant_concentration)
% 
% for i = 1:size(population_gen,2)
%     plot(time_observed,Cmod_best(i,:));
%     hold on
% end