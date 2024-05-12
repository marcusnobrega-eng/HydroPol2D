% Sensitivity Analysis Code
% Developer: Marcus Nobrega
% Goal - Do sensitivity analysis from modeling results

%% Pre-Processing 2
input_table = readtable(model_folder);

% Load Model Functions
HydroPol2D_tools = char(table2cell(input_table(9,31)));
addpath(genpath(char(HydroPol2D_tools)));


%% Baseline Parameters

HydroPol2D_preprocessing

filename = "Pre_Processing_Data.mat";
save(filename);

% HydroPol2D_Main_While
load('Pre_Processing_Data.mat'); k = 1;
if flags.flag_D8 == 1
    [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA_D8(wse_slope_zeros,Distance_Matrix,min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall);
else
    [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA(min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall);
end
% Save Data - Baseline
Flow = Qmod; Baseline_Outputs(1,:,1) = Flow;
Depth = Dmod; Baseline_Outputs(1,:,2) = Depth;
if flags.flag_waterquality == 1
    Conc = Cmod; Baseline_Outputs(1,:,3) = Conc;
end
Baseline_Outputs(1,1,4) = Flooded_Area; % Flooded areas (m2)
Baseline_Outputs(1,1,5) = Risk_Area; % Areas with drag risk (m2)

% %%%%%%%%%%%%%% LULC DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable('LULC_parameters.xlsx');
input_data = table2array(input_table(:,2:end)); % numbers
lulc_parameters = input_data(:,2:end);
n_lulc = sum(lulc_parameters(:,1)>=0); % Number of LULC
lulc_parameters = lulc_parameters(1:n_lulc,:);
% imp_index = lulc_parameters(1,end);
% % LULC Index
% % LULC_index = input_data(:,1);
%
% %%%%%%%%%%%%%% SOIL DATA %%%%%%%%%%%%%%%%%%%
input_table = readtable('SOIL_parameters.xlsx');
input_data = table2array(input_table(:,2:end)); % numbers
soil_parameters = input_data(:,2:end); % Number of Soil Types
n_soil = sum(soil_parameters(:,1)>=0);
soil_parameters = soil_parameters(1:n_soil,:);
soil_parameters(:,[3 5]) = []; % BE CAREFULL HERE
% SOIL Index
% SOIL_index = input_data(:,1);
% SOIL_name = input_table(:,1);

% Baseline Parameters
% x_LULC_hydro = LULC_Properties.;
% x_LULC_WQ = ;
% x_SOIL_hydro = ;

% Pre-Processing
% Note that the model will assume the input data from the sheets as the
% baseline scenario for the pre-processing analysis


%% LULC Calculations
% Preallocate Array
data_SA_LULC = zeros(n_lulc*6,length(running_control.time_record_hydrograph),length(var_range),3); % Saving Flow, Depth, and Concentration
% - All variables, All times, All ranges, All States
m = 0; mm = 0;
for jj = 1:2 % Number of LULC-Based Parameters
    for ii = 1:n_lulc
        mm = mm + 1;
        perc = mm/(n_lulc*6)*100
        pause(0.25)
        for kk = 1:length(var_range)
            % Baseline Parameters
            load('Pre_Processing_Data.mat'); k = 1;
            m = m + 1;
            if var_range(kk,1) == 1 % No change at all, jump everything
                % Save Data
                Flow = Baseline_Outputs(1,:,1); data_SA_LULC(mm,:,kk,1) = Flow;
                Depth = Baseline_Outputs(1,:,2); data_SA_LULC(mm,:,kk,2) = Depth;
                if flags.flag_waterquality == 1
                    Conc = Baseline_Outputs(1,:,3); data_SA_LULC(mm,:,kk,3) = Conc;
                end
                data_SA_LULC(mm,1,kk,4) = Baseline_Outputs(1,1,4); % Flooded areas (m2)
                data_SA_LULC(mm,1,kk,5) = Baseline_Outputs(1,1,5); % Areas with drag risk (m2)
                
            else % Here parameters are different than baseline
                if flags.flag_waterquality ~= 1 && jj <= 2 % Only Hydrodynamics
                    if jj == 1
                        LULC_Properties.roughness(LULC_Properties.idx_lulc(:,:,ii)) = var_range(kk,1)*lulc_parameters(ii,1); % New value from Sensitivity Analysis
                    elseif jj == 2
                        LULC_Properties.h_0(LULC_Properties.idx_lulc(:,:,ii)) = var_range(kk,1)*lulc_parameters(ii,2); % New value from Sensitivity Analysis
                    end
                else
                    if jj == 3
                        LULC_Properties.C_1(LULC_Properties.idx_lulc(:,:,ii)) = var_range(kk,1)*lulc_parameters(ii,3); % New value from Sensitivity Analysis
                    elseif jj == 4
                        LULC_Properties.C_2(LULC_Properties.idx_lulc(:,:,ii)) = var_range(kk,1)*lulc_parameters(ii,4); % New value from Sensitivity Analysis
                    elseif jj == 5
                        LULC_Properties.C_3(LULC_Properties.idx_lulc(:,:,ii)) = var_range(kk,1)*lulc_parameters(ii,5); % New value from Sensitivity Analysis
                    else
                        LULC_Properties.C_4(LULC_Properties.idx_lulc(:,:,ii)) = var_range(kk,1)*lulc_parameters(ii,6); % New value from Sensitivity Analysis
                    end
                end

                % Run HydroPol2D Model and Retrieve Data
                if flags.flag_waterquality ~= 1
                    idx_nan_5 = [];
                    WQ_States = [];
                end
                if flags.flag_D8 == 1
                    [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA_D8(wse_slope_zeros,Distance_Matrix,min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall);
                else
                    [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA(min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall);
                end
                % Save Data
                Flow = Qmod; data_SA_LULC(mm,:,kk,1) = Flow;
                Depth = Dmod; data_SA_LULC(mm,:,kk,2) = Depth;
                if flags.flag_waterquality == 1
                    Conc = Cmod; data_SA_LULC(mm,:,kk,3) = Conc;
                end
                Flood = Flooded_Area; data_SA_LULC(mm,1,kk,4) = Flooded_Area;
                Risk = Risk_Area; data_SA_LULC(mm,1,kk,5) = Risk_Area;
            end
        end
    end
end


%% Calculating Metrics from Results
% Peak Flow
baseline_peakflow = max(Baseline_Outputs(1,:,1));
baseline_volume = sum(Baseline_Outputs(1,:,1)*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1))*60);
baseline_time_to_peak = (find(Baseline_Outputs(1,:,1) == max(Baseline_Outputs(1,:,1)),1,'first') - 1)*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1)); % min;
for ii = 1:size(data_SA_LULC,1) % Number of parameters
    for jj = 1:length(var_range) % Number of parameter variation
        % Hydrologic-Hydrodynamic Functions
        flow = data_SA_LULC(ii,:,jj,1); % m3/s
        peakflow_LULC(ii,jj) = max(flow); % m3/s
        volume_LULC(ii,jj) = sum(flow*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1))*60); % m3
        time_to_peak_LULC(ii,jj) = (find(flow == max(flow),1,'first') - 1)*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1)); % min
        flood_areas = data_SA_LULC(ii,1,jj,4); % m2
        risk_areas = data_SA_LULC(ii,1,jj,5); % m2
        peakflow_variance_LULC(ii,jj) = (peakflow_LULC(ii,jj) - baseline_peakflow)/baseline_peakflow*100; % percentage
        volume_variance_LULC(ii,jj) = (volume_LULC(ii,jj) - baseline_volume)/baseline_volume*100; % percentage
        time_to_peak_variance_LULC(ii,jj) = (time_to_peak_LULC(ii,jj) - baseline_time_to_peak)/baseline_time_to_peak*100; % percentage
        flood_area_variance_LULC(ii,jj) = (flood_areas - Baseline_Outputs(1,1,4))/Baseline_Outputs(1,1,4)*100; % percentage
        risk_area_variance_LULC(ii,jj) = (risk_areas - Baseline_Outputs(1,1,5))/Baseline_Outputs(1,1,5)*100; % percentage
        % WQ Functions
        if flags.flag_waterquality == 1
            % WQ Functions
        end
    end
end


%% SOIL Calculations
% Preallocate Array
data_SA_SOIL = zeros(n_soil*3,length(running_control.time_record_hydrograph),length(var_range),3); % Saving Flow, Depth, and Concentration
% - All variables, All times, All ranges, All States
m = 0; mm = 0;
for jj = 1:3 % Number of SOIL-Based Parameters
    for ii = 1:n_soil
        mm = mm + 1;
        perc = mm/(n_soil*3)*100
        pause(0.25)
        for kk = 1:length(var_range)
            load('Pre_Processing_Data.mat'); k = 1;
            if var_range(kk,1) == 1 % No change at all, jump everything
                % Save Data
                Flow = Baseline_Outputs(1,:,1); data_SA_SOIL(mm,:,kk,1) = Flow;
                Depth = Baseline_Outputs(1,:,2); data_SA_SOIL(mm,:,kk,2) = Depth;
                if flags.flag_waterquality == 1
                    Conc = Baseline_Outputs(1,:,3); data_SA_SOIL(mm,:,kk,3) = Conc;
                end
                data_SA_SOIL(mm,1,kk,4) = Baseline_Outputs(1,1,4); % Flooded areas (m2)
                data_SA_SOIL(mm,1,kk,5) = Baseline_Outputs(1,1,5); % Areas with drag risk (m2)
            else % Here parameters are different than baseline
                if jj == 1
                    Soil_Properties.ksat(Soil_Properties.idx_soil(:,:,ii)) = var_range(kk,1)*soil_parameters(ii,1); % New value from Sensitivity Analysis
                elseif jj == 2
                    Soil_Properties.psi(Soil_Properties.idx_soil(:,:,ii)) =  var_range(kk,1)*soil_parameters(ii,2); % New value from Sensitivity Analysis
                elseif jj == 3
                    Soil_Properties.teta_sat(Soil_Properties.idx_soil(:,:,ii)) = var_range(kk,1)*soil_parameters(ii,3); % New value from Sensitivity Analysis
                    Soil_Properties.teta_i(Soil_Properties.idx_soil(:,:,ii)) = 0; % Imposition
                end
                % Boundary Condition in Ksat
                Soil_Properties.ksat(LULC_Properties.idx_imp) = 0; % Impervious areas

                % Run HydroPol2D Model and Retrieve Data
                if flags.flag_D8 == 1
                    [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA_D8(wse_slope_zeros,Distance_Matrix,min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall);
                else
                    [Qmod, Cmod,Dmod,Flooded_Area,Risk_Area] = HydroPol2D_Routing_Solver_SA(min_soil_moisture,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties,idx_rivers,Lateral_Groundwater_Flux,Reservoir_Data,Input_Rainfall);
                end
                % Save Data
                Flow = Qmod; data_SA_SOIL(mm,:,kk,1) = Flow;
                Depth = Dmod; data_SA_SOIL(mm,:,kk,2) = Depth;
                if flags.flag_waterquality == 1
                    Conc = Cmod; data_SA_SOIL(mm,:,kk,3) = Conc;
                end
                Flood = Flooded_Area; data_SA_SOIL(mm,1,kk,4) = Flooded_Area;
                Risk = Risk_Area; data_SA_SOIL(mm,1,kk,5) = Risk_Area;                
            end
        end
    end
end


%% Calculating Metrics from Results - SOIL
% Peak Flow
baseline_peakflow = max(Baseline_Outputs(1,:,1));
baseline_volume = sum(Baseline_Outputs(1,:,1)*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1))*60);
baseline_time_to_peak = (find(Baseline_Outputs(1,:,1) == max(Baseline_Outputs(1,:,1)),1,'first') - 1)*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1)); % min;

for ii = 1:size(data_SA_SOIL,1) % Number of parameters
    for jj = 1:length(var_range) % Number of parameter variation
        % Hydrologic-Hydrodynamic Functions
        flow = data_SA_SOIL(ii,:,jj,1); % m3/s
        peakflow_SOIL(ii,jj) = max(flow); % m3/s
        volume_SOIL(ii,jj) = sum(flow*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1))*60); % m3
        time_to_peak_SOIL(ii,jj) = (find(flow == max(flow),1,'first') - 1)*(running_control.time_record_hydrograph(2) - running_control.time_record_hydrograph(1)); % min
        flood_areas = data_SA_SOIL(ii,1,jj,4); % m2
        risk_areas = data_SA_SOIL(ii,1,jj,5); % m2
        peakflow_variance_SOIL(ii,jj) = (peakflow_SOIL(ii,jj) - baseline_peakflow)/baseline_peakflow*100; % percentage
        volume_variance_SOIL(ii,jj) = (volume_SOIL(ii,jj) - baseline_volume)/baseline_volume*100; % percentage
        time_to_peak_variance_SOIL(ii,jj) = (time_to_peak_SOIL(ii,jj) - baseline_time_to_peak)/baseline_time_to_peak*100; % percentage
        flood_area_variance_SOIL(ii,jj) = (flood_areas - Baseline_Outputs(1,1,4))/Baseline_Outputs(1,1,4)*100; % percentage
        risk_area_variance_SOIL(ii,jj) = (risk_areas - Baseline_Outputs(1,1,5))/Baseline_Outputs(1,1,5)*100; % percentage

        % WQ Functions
        if flags.flag_waterquality == 1
            % WQ Functions
        end
    end
end


%% Plotting Results - LULC Based
close all
% Creating Modeling Results Folder
% Create the folder name
folderName = 'Modeling_Results';

% Check if the folder already exists
if ~exist(folderName, 'dir')
    % If it doesn't exist, create the folder
    mkdir(folderName);
    disp('Folder "Modeling_Results" created successfully!');
else
    disp('Data sucessfully exported in Modeling_Results Folder');
end
index = 1;
for i = 1:2 % Variables
    for j = 1:n_lulc
        if i == 1
            x_labels{index,1} = sprintf('$n_%d$',j);
        elseif i == 2
            x_labels{index,1} = sprintf('$h_{0,%d}$',j);
        elseif i == 3
            x_labels{index,1} = sprintf('$C_{1,%d}$',j);
        elseif i == 4
            x_labels{index,1} = sprintf('$C_{2,%d}$',j);
        elseif i == 5
            x_labels{index,1} = sprintf('$C_{3,%d}$',j);
        elseif i == 6
            x_labels{index,1} = sprintf('$C_{4,%d}$',j);
        end
        index = index + 1;
    end
end
var_range_plot = var_range - 1;
figure % Plotting Intensities
set(gcf,'units','inches','position',[3,3,6.5,4])
colors = linspecer(n_lulc);
index = 0;
linetypes = {'-',':','-.','--','-',':','-.','--','-'};
for jj = 1:4 % Number of Output Functions
    index = index + 1;
    subplot(4,2,index)
    title_label = {'Peak Flow','Runoff Volume','Time to Peak','Flooded Areas','Risk Areas'};
    for ii = 1:n_lulc-1 % Manning and h_0 only
        if jj == 1
            plot(var_range_plot*100,peakflow_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 2
            plot(var_range_plot*100,volume_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 3
            plot(var_range_plot*100,time_to_peak_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 4
            plot(var_range_plot*100,flood_area_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 5
            plot(var_range_plot*100,risk_area_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        end
        hold on
    end
    xlabel('Parameter Variation [\%]','Interpreter','latex','FontSize',12)
    ylabel('Output Variation [\%]','Interpreter','latex','FontSize',12)
    title(title_label(jj),'interpreter','latex','FontSize',12)
    legend(x_labels(1:(n_lulc-1)),'Location','Best','Interpreter','latex','FontSize',12)
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    grid on
    box on

    index = index + 1;
    subplot(4,2,index)
    colors = linspecer(n_lulc);
    for ii = 1:n_lulc-1  % Manning and h_0 only
        if jj == 1
            plot(var_range_plot*100,peakflow_variance_LULC(n_lulc + ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 2
            plot(var_range_plot*100,volume_variance_LULC(n_lulc + ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 3
            plot(var_range_plot*100,time_to_peak_variance_LULC(n_lulc + ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 4
            plot(var_range_plot*100,flood_area_variance_LULC(n_lulc + ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 5
            plot(var_range_plot*100,risk_area_variance_LULC(n_lulc + ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        end
        hold on
    end
    xlabel('Parameter Variation [\%]','Interpreter','latex','FontSize',12)
    ylabel('Output Variation [\%]','Interpreter','latex','FontSize',12)
    title(title_label(jj),'interpreter','latex','FontSize',12)
    legend(x_labels((n_lulc + 1):(2*n_lulc-1)),'Location','Best','Interpreter','latex','FontSize',12)
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    grid on
    box on
    hold on
end
file = 'Sensitivity_Analysis_LULC.pdf';
exportgraphics(gcf,fullfile(folderName,file),'ContentType','vector')

%% Plotting Results - Only Manning
close all
% Creating Modeling Results Folder
% Create the folder name
folderName = 'Modeling_Results';

% Check if the folder already exists
if ~exist(folderName, 'dir')
    % If it doesn't exist, create the folder
    mkdir(folderName);
    disp('Folder "Modeling_Results" created successfully!');
else
    disp('Data sucessfully exported in Modeling_Results Folder');
end
index = 1;
for i = 1:2 % Variables
    for j = 1:n_lulc
        if i == 1
            x_labels{index,1} = sprintf('$n_%d$',j);
        elseif i == 2
            x_labels{index,1} = sprintf('$h_{0,%d}$',j);
        elseif i == 3
            x_labels{index,1} = sprintf('$C_{1,%d}$',j);
        elseif i == 4
            x_labels{index,1} = sprintf('$C_{2,%d}$',j);
        elseif i == 5
            x_labels{index,1} = sprintf('$C_{3,%d}$',j);
        elseif i == 6
            x_labels{index,1} = sprintf('$C_{4,%d}$',j);
        end
        index = index + 1;
    end
end
var_range_plot = var_range - 1;
figure % Plotting Intensities
set(gcf,'units','inches','position',[3,3,6.5,4])
colors = linspecer(n_lulc);
index = 0;
linetypes = {'-',':','-.','--','-',':','-.','--','-'};
for jj = 1:3 % Number of Output Functions
    index = index + 1;
    subplot(3,1,index)
    title_label = {'Peak Flow','Runoff Volume','Time to Peak'};
    for ii = 1:n_lulc-1 % Manning and h_0 only
        if jj == 1
            plot(var_range_plot*100,peakflow_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 2
            plot(var_range_plot*100,volume_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        elseif jj == 3
            plot(var_range_plot*100,time_to_peak_variance_LULC(ii,:),'color',colors(ii,:),'linewidth',2,'linestyle',linetypes(ii));
        end
        hold on
    end
    xlabel('Parameter Variation [\%]','Interpreter','latex','FontSize',12)
    ylabel('Output Variation [\%]','Interpreter','latex','FontSize',12)
    title(title_label(jj),'interpreter','latex','FontSize',12)
    legend(x_labels(1:(n_lulc-1)),'Location','Best','Interpreter','latex','FontSize',12)
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    grid on
    box on

end
%% Plotting Results - SOIL Based
close all
% Creating Modeling Results Folder
% Create the folder name
folderName = 'Modeling_Results';

% Check if the folder already exists
if ~exist(folderName, 'dir')
    % If it doesn't exist, create the folder
    mkdir(folderName);
    disp('Folder "Modeling_Results" created successfully!');
else
    disp('Data sucessfully exported in Modeling_Results Folder');
end
index = 1;
for i = 1:3 % Variables
    for j = 1:n_soil
        if i == 1
            x_labels{index,1} = sprintf('$k_{sat,%d}$',j);
        elseif i == 2
            x_labels{index,1} = sprintf('$\\psi_{%d}$ ',j);
        elseif i == 3
            x_labels{index,1} = sprintf('$\\Delta \\theta_{%d}$',j);
        end
        index = index + 1;
    end

end

var_range_plot = var_range - 1;
figure % Plotting Intensities
set(gcf,'units','inches','position',[2.5,2.5,8,4])
colors = linspecer(n_soil*3);
index = 1;
for jj = 1:4 % Number of Output Functions
    for kk = 1:3 % Number of Decision Variables
        subplot(4,3,index)
        index = index + 1;
        % Peak Flow / Volume / Time to Peak
        title_label = {'Peak Flow','Runoff Volume','Time to Peak','Flooded Areas','Risk Areas'};
        for ii = 1:n_soil % All soil types
            if jj == 1
                plot(var_range_plot*100,peakflow_variance_SOIL((kk-1)*n_soil + ii,:),'color',colors((kk-1)*n_soil + ii,:),'linewidth',2);
                legend(x_labels((kk-1)*n_soil + 1:kk*n_soil),'Location','NorthEastOutside','Interpreter','latex','FontSize',9)
            elseif jj == 2
                plot(var_range_plot*100,volume_variance_SOIL((kk-1)*n_soil + ii,:),'color',colors((kk-1)*n_soil + ii,:),'linewidth',2);
            elseif jj == 3
                plot(var_range_plot*100,time_to_peak_variance_SOIL((kk-1)*n_soil + ii,:),'color',colors((kk-1)*n_soil + ii,:),'linewidth',2);
            elseif jj == 4
                plot(var_range_plot*100,flood_area_variance_SOIL((kk-1)*n_soil + ii,:),'color',colors((kk-1)*n_soil + ii,:),'linewidth',2);
            % elseif jj == 5
            %     plot(var_range_plot*100,risk_area_variance_SOIL((kk-1)*n_soil + ii,:),'color',colors((kk-1)*n_soil + ii,:),'linewidth',2);
            end
            hold on
        end
        xlabel('Parameter Variation [\%]','Interpreter','latex','FontSize',12)
        ylabel('Output Variation [\%]','Interpreter','latex','FontSize',12)
        title(title_label(jj),'interpreter','latex','FontSize',12)
        
        % legend(x_labels((kk-1)*n_soil + 1:kk*n_soil),'Location','Best','Interpreter','latex','FontSize',12)
        set(gca, 'TickLength', [0.015 0.01]);
        set(gca,'Tickdir','out')
        set(gca, 'FontName', 'Garamond', 'FontSize', 12)
        grid on
        box on
    end
end


file = 'Sensitivity_Analysis_SOIL.pdf';
exportgraphics(gcf,fullfile(folderName,file),'ContentType','vector')

%% Only Risk and Flood Areas for Soils
var_range_plot = var_range - 1;
figure % Plotting Intensities
set(gcf,'units','inches','position',[3,3,6.5,4])
colors = linspecer(n_soil*3);
index = 1;
subplot(2,1,1)
for i = 1:3
    plot(var_range_plot*100,flood_area_variance_SOIL(i,:),'color',colors(i,:),'linewidth',2);
    hold on
end
xlabel('Parameter Variation [\%]','Interpreter','latex','FontSize',12)
ylabel('Output Variation [\%]','Interpreter','latex','FontSize',12)
title('Flooded Areas','interpreter','latex','FontSize',12)
legend(x_labels((kk-1)*n_soil + 1:kk*n_soil),'Location','Best','Interpreter','latex','FontSize',12)
set(gca, 'TickLength', [0.015 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
grid on
box on

subplot(2,1,2)
for i = 1:3
    plot(var_range_plot*100,risk_area_variance_SOIL(i,:),'color',colors(i,:),'linewidth',2);
    hold on
end
xlabel('Parameter Variation [\%]','Interpreter','latex','FontSize',12)
ylabel('Output Variation [\%]','Interpreter','latex','FontSize',12)
title('Risk Areas','interpreter','latex','FontSize',12)
legend(x_labels((kk-1)*n_soil + 1:kk*n_soil),'Location','Best','Interpreter','latex','FontSize',12)
set(gca, 'TickLength', [0.015 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
grid on
box on
legend('$k_{sat}$','$\Delta \theta$','$\psi$','interpreter','latex')

%% Plotting Hydrographs
color_plot = linspecer(8);
line_width = [1.5 1.25 1 0.75 1 1.25 1.5];
line_type = {':','-.','--','-','--','-.',':'};
for i = 1:16
    for j = 1:size(data_SA_LULC,3)
        for k = 1:2 % Two Variables
            index_plot = mod(i,8);
            if index_plot == 0
                index_plot = 1;
            end
            subplot(1,2,k)
            dim = 1; % 1 = Flow, 2 = Depth, 3 = Conc
            variable = data_SA_LULC(i,:,j,dim);
            plot(running_control.time_record_hydrograph,variable,'color',color_plot(index_plot,:),'linewidth',line_width(j),'linestyle',line_type(j))
            hold on
        end
    end
end

save('Sensitivity_workspace');
