% Running HydroPol2D and Retrieving A few States


function [Obj_fun,Qmod,Cmod,Qmod_gauges,Cmod_gauges] = HydroPol2D_Solver_with_x(x,min_soil_moisture,Automatic_Calibrator,flag_calibrate_wq,Reservoir_Data,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties)
x = x';

% function [Obj_fun] = Objective_Function_HydroPol2D(x,min_soil_moisture,time_observed,observed_flow,pollutant_concentration,delta_p_obs,warmup_I_0,warmup_depths,n_events,n_LULC,n_SOIL,easting_observed_cells_calibration,northing_observed_cells_calibration,flag_calibrate_wq,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, steps, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties)
% --- LULC Parameter Assigning --- %
for i = 1:Automatic_Calibrator.n_LULC %
    % Only Roughness and h_0
    LULC_Properties.roughness(LULC_Properties.idx_lulc(:,:,i)) = x(i,1); % assigning values for roughness at impervious areas
    LULC_Properties.h_0(LULC_Properties.idx_lulc(:,:,i)) = x(i + Automatic_Calibrator.n_LULC,1); % Initial Abstraction
end


% Calculation of Spatial Inputs
% -- Soil Parameter Assigning ---- %
Soil_Properties.ksat = zeros(size(Elevation_Properties.elevation_cell)); Soil_Properties.psi = Soil_Properties.ksat; Soil_Properties.teta_s = Soil_Properties.ksat;
for i = 1:Automatic_Calibrator.n_SOIL
    Soil_Properties.ksat(Soil_Properties.idx_soil(:,:,i)) = x(i + 6*Automatic_Calibrator.n_LULC,1); % similar things happening below
    Soil_Properties.teta_sat(Soil_Properties.idx_soil(:,:,i)) = x(i + 6*Automatic_Calibrator.n_LULC + 1*Automatic_Calibrator.n_SOIL,1);
    Soil_Properties.psi(Soil_Properties.idx_soil(:,:,i)) = x(i + 6*Automatic_Calibrator.n_LULC + 2*Automatic_Calibrator.n_SOIL,1);
end
Soil_Properties.teta_i = zeros(size(Soil_Properties.teta_s));


%%
% Mask in Impervious Areas
Soil_Properties.ksat(LULC_Properties.idx_imp) = 0;



% Preallocation
LULC_Properties.C_1 = zeros(size(Elevation_Properties.elevation_cell));
LULC_Properties.C_2 = zeros(size(Elevation_Properties.elevation_cell));
LULC_Properties.C_3 = zeros(size(Elevation_Properties.elevation_cell));
LULC_Properties.C_4  = zeros(size(Elevation_Properties.elevation_cell));
%% Looping for n_events

for ii = 1:Automatic_Calibrator.n_events
    sprintf('Event %d',ii)


    % Warmup Files - Initial Conditions
    depths.d_t = double(Automatic_Calibrator.warmup_depths(:,:,ii));
    depths.d_t(isnan(depths.d_t)) = 0;
    depths.d_t(isnan(Elevation_Properties.elevation_cell)) = nan; depths.d_t(isinf(Elevation_Properties.elevation_cell)) = nan;
    Soil_Properties.I_0 = double(Automatic_Calibrator.warmup_I_0(:,:,ii));
    Soil_Properties.I_0(isnan(Soil_Properties.I_0)) = 0;
    Soil_Properties.I_0(isnan(Elevation_Properties.elevation_cell)) = nan; Soil_Properties.I_0(isinf(Elevation_Properties.elevation_cell)) = nan;
    Soil_Properties.I_t = Soil_Properties.I_0; Soil_Properties.I_p = Soil_Properties.I_0;


    % Calculation of B_t
    if flag_calibrate_wq == 1 %
        for i = 1:Automatic_Calibrator.n_LULC
            LULC_Properties.C_1(LULC_Properties.idx_lulc(:,:,i)) =  x(i + 2*Automatic_Calibrator.n_LULC,1);
            LULC_Properties.C_2(LULC_Properties.idx_lulc(:,:,i)) =  x(i + 3*Automatic_Calibrator.n_LULC,1);
            LULC_Properties.C_3(LULC_Properties.idx_lulc(:,:,i)) =  x(i + 4 *Automatic_Calibrator.n_LULC,1);
            LULC_Properties.C_4(LULC_Properties.idx_lulc(:,:,i)) =  x(i + 5*Automatic_Calibrator.n_LULC,1);
        end
        WQ_States.B_t = LULC_Properties.C_1.*(1 - exp(1).^(-LULC_Properties.C_2*ADD_events(ii))); % kg/ha
        WQ_States.B_t = WQ_States.B_t*Wshed_Properties.cell_areasting_observed_cells_calibrationea/10^4; % kg per cell
    end

    % Entering Observations
    n_observations = sum(~isnan(Automatic_Calibrator.time_observed(:,ii)));
    obs_discharge = Automatic_Calibrator.observed_flow(1:n_observations,:,ii);
    obs_pollutograph = Automatic_Calibrator.pollutant_concentration(1:n_observations,:,ii);

    % Observations
    easting_obs_cell = Automatic_Calibrator.easting_observed_cells(1,:,ii); % Values where observations are made
    northing_obs_cell = Automatic_Calibrator.northing_observed_cells(1,:,ii); % Values where observations are made

    time_obs = Automatic_Calibrator.time_observed(1:n_observations,ii);

    % Entering delta_p
    count = sum(~isnan(Automatic_Calibrator.delta_p_obs(ii,:)));
    delta_p = Automatic_Calibrator.delta_p_obs(ii,1:count);
    time_deltap = cumsum(ones(1,length(delta_p))*time_step_model);
    running_control.time_deltap = time_deltap;
    
    % Adding data to BC_States
    BC_States.delta_p = delta_p;

    % Simulation Time
    running_control.routing_time = time_obs(end,1);
    t = running_control.time_step_model;

    % Run HydroPol2D Routing Solver
    if flags.flag_D8 ~= 1                                  
        [Qmod(:,ii), Cmod(:,ii),Qmod_gauges(:,:,ii),Cmod_gauges(:,:,ii)] = HydroPol2D_Routing_Solver(t,time_obs,easting_obs_cell,northing_obs_cell,min_soil_moisture,Reservoir_Data,BC_States, CA_States, Courant_Parameters, date_begin, DEM_raster, depths, Elevation_Properties, flags, gauges, GIS_data, Human_Instability, Hydro_States, idx_nan, idx_nan_5, idx_outlet, Inflow_Parameters, LULC_Properties, Maps, nx_max, ny_max, outlet_index, outlet_runoff_volume, outlet_type, Rainfall_Parameters, recording_parameters, running_control, slope_outlet, Soil_Properties, t_previous, time_calculation_routing, time_step, time_step_model, tmin_wq, WQ_States, Wshed_Properties);
    else
        % Function for D8
    end
    flag_outlet = 0; % Not using the outlet, 1 using the outlet

    flag_opt_fun = 5; % Decides which objective functions is being used

    small_number = 1e-6; % This number is used in functions that have division so that the fraction does not become null
    

    if flag_outlet ~= 1
        if flags.flag_waterquality ~= 1
            if flag_opt_fun == 1 % NSE
                % Nash-Suctclife-Efficiency and Objective function
                weights = 1/size(Qmod_gauges(:,:,ii),2)*ones(1,size(Qmod_gauges(:,:,ii),2));
                size_mod = size(Qmod_gauges(:,:,ii),1);
                NSE = 1 - sum((obs_discharge(1:size_mod,:) - Qmod_gauges(:,:,ii)).^2)./(sum((obs_discharge(1:size_mod,:) - mean(Qmod_gauges(:,:,ii))).^2)); % We want to maximize it               
                if sum(isnan(NSE))
                    error('Some of the NSE are becoming NaN.')
                end
                Obj_fun(ii,1) = (-1)*sum(weights.*NSE); % Therefore, we want to minimize it
            end

            if flag_opt_fun == 2 % RMSE
                weights = 1/size(Qmod_gauges(:,:,ii),2)*ones(1,size(Qmod_gauges(:,:,ii),2));                
                n_elements = length(obs_discharge);
                RMSE = sqrt(sum((obs_discharge(1:length(Qmod_gauges(:,:,ii)),:) - Qmod_gauges(:,:,ii)).^2/n_elements));
                Obj_fun(ii,1) = sum(weights.*RMSE); % minimize RMSE
            end

            if flag_opt_fun == 3 % R2
                r2 = corrcoef(obs_discharge(1:length(Qmod_gauges(:,:,ii)),:),Qmod_gauges(:,:,ii));
                Obj_fun(ii,1) = -r2(1,2); % maximize r2
            end

            if flag_opt_fun == 4 % Peak Flow
                z1 = max(obs_discharge(1:length(Qmod_gauges(:,:,ii)),:));
                z2 = max(Qmod_gauges(:,:,ii));
                Obj_fun(ii,1) = abs(z1 - z2); % Absolute value of peak flows
            end
            if flag_opt_fun == 5 % NSE and Volume Error
                % NSE and Volume Error Objective
                % Function
                weights = 1/size(Qmod_gauges(:,:,ii),2)*ones(1,size(Qmod_gauges(:,:,ii),2));
                size_mod = size(Qmod_gauges(:,:,ii),1);
                NSE = 1 - sum((obs_discharge(1:size_mod,:) - Qmod_gauges(:,:,ii)).^2)./(sum((obs_discharge(1:size_mod,:) - mean(Qmod_gauges(:,:,ii))).^2) + small_number); % We want to maximize it               
                if sum(isnan(NSE))
                    error('Some of the NSE are becoming NaN.')
                end               
                % Volume Error
                VE = abs(sum((obs_discharge(1:size_mod,:) - Qmod_gauges(:,:,ii)))./(sum(obs_discharge(1:size_mod,:)) + small_number));
                Obj_fun(ii,1) = (-1)*sum(weights.*(1*NSE - 0.5*VE)); % Therefore, we want to minimize it
            end  
            if flag_opt_fun == 6 % RMSE and Volume Error
                % Nash-Suctclife-Efficiency and Volume Error Objective
                % Function
                weights = 1/size(Qmod_gauges(:,:,ii),2)*ones(1,size(Qmod_gauges(:,:,ii),2));                
                n_elements = length(obs_discharge);
                RMSE = sqrt(sum((obs_discharge(1:length(Qmod_gauges(:,:,ii)),:) - Qmod_gauges(:,:,ii)).^2/n_elements));
                % Volume Error
                size_mod = size(Qmod_gauges(:,:,ii),1);
                VE = abs(sum((obs_discharge(1:size_mod,:) - Qmod_gauges(:,:,ii)))./(sum(obs_discharge(1:size_mod,:)) + small_number));
                Obj_fun(ii,1) = 1*RMSE + 0.5*VE; % Therefore, we want to minimize it
            end             
        else % Modelling Water Quality
            if flag_opt_fun == 1 % NSE
                % Nash-Suctclife-Efficiency and Objective function
                NSE = 1 - sum((obs_pollutograph - Cmod_gauges(:,:,ii)).^2)/(sum((obs_pollutograph - mean(obs_pollutograph)).^2)); % We want to maximize it
                Obj_fun(ii,1) = -NSE; % Therefore, we want to minimize it
            end

            if flag_opt_fun == 2 % RMSE
                n_elements = length(obs_pollutograph);
                RMSE = sqrt(sum((obs_pollutograph - Cmod_gauges(:,:,ii)).^2/n_elements));
                Obj_fun(ii,1) = RMSE; % minimize RMSE
            end

            if flag_opt_fun == 3 % R2
                r2 = corrcoef(obs_pollutograph,Cmod_gauges(:,:,ii));
                Obj_fun(ii,1) = -r2(1,2); % maximize r2
            end

            if flag_opt_fun == 4 % Peak Flow
                z1 = max(obs_pollutograph);
                z2 = max(Cmod_gauges(:,:,ii));
                Obj_fun(ii,1) = abs(z1 - z2); % Absolute value of peak flows
            end          
        end

    else

        if flag_waterquality ~= 1
            if flag_opt_fun == 1 % NSE
                % Nash-Suctclife-Efficiency and Objective function
                NSE = 1 - sum((obs_discharge - Qmod).^2)/(sum((obs_discharge - mean(obs_discharge)).^2)); % We want to maximize it
                Obj_fun(ii,1) = -NSE; % Therefore, we want to minimize it
            end

            if flag_opt_fun == 2 % RMSE
                n_elements = length(obs_discharge);
                RMSE = sqrt(sum((obs_discharge - Qmod).^2/n_elements));
                Obj_fun(ii,1) = RMSE; % minimize RMSE
            end

            if flag_opt_fun == 3 % R2
                r2 = corrcoef(obs_discharge,Qmod);
                Obj_fun(ii,1) = -r2(1,2); % maximize r2
            end

            if flag_opt_fun == 4 % Peak Flow
                z1 = max(obs_discharge);
                z2 = max(Qmod);
                Obj_fun(ii,1) = abs(z1 - z2); % Absolute value of peak flows
            end
        else % Modelling Water Quality
            if flag_opt_fun == 1 % NSE
                % Nash-Suctclife-Efficiency and Objective function
                NSE = 1 - sum((obs_pollutograph - Cmod).^2)/(sum((obs_pollutograph - mean(obs_pollutograph)).^2)); % We want to maximize it
                Obj_fun(ii,1) = -NSE; % Therefore, we want to minimize it
            end

            if flag_opt_fun == 2 % RMSE
                n_elements = length(obs_pollutograph);
                RMSE = sqrt(sum((obs_pollutograph - Cmod).^2/n_elements));
                Obj_fun(ii,1) = RMSE; % minimize RMSE
            end

            if flag_opt_fun == 3 % R2
                r2 = corrcoef(obs_pollutograph,Cmod);
                Obj_fun(ii,1) = -r2(1,2); % maximize r2
            end

            if flag_opt_fun == 4 % Peak Flow
                z1 = max(obs_pollutograph);
                z2 = max(Cmod);
                Obj_fun(ii,1) = abs(z1 - z2); % Absolute value of peak flows
            end
        end
    end
    sprintf('End of all events')
end

%% Plot a Quick Chart
% figure(2)
% if flag_waterquality ~= 1
%     for jj = 1:n_events
%         if n_events > 1
%             col = 2;
%         else
%             col = 1;
%         end
%         subplot(n_events,col,ii);
%         plot(time_obs,Qmod,'-black','Marker','*')
%         hold on
%         plot(time_obs,obs_discharge,'--b','Marker','*')
%         xlabel('Elapsed time (min)','FontSize',12,'Interpreter','latex')
%         ylabel('Discharge ($m^3/s$)','FontSize',12,'Interpreter','latex')
%         legend('Modeled','Observed','interpreter','latex')
%     end
% else
%     for jj = 1:n_events
%         if n_events > 1
%             col = 2;
%         else
%             col = 1;
%         end
%         subplot(n_events,col,ii);
%         plot(time_obs,Cmod,'-black','Marker','*')
%         hold on
%         plot(time_obs,obs_pollutograph,'--b','Marker','*')
%         xlabel('Elapsed time (min)','FontSize',12,'Interpreter','latex')
%         ylabel('Pollutant Concentration ($mg/L$)','FontSize',12,'Interpreter','latex')
%         legend('Modeled','Observed','interpreter','latex')
%     end
% end
% pause(0.25)
% close(2)
% Final Objective Function
Obj_fun = mean(Obj_fun);

end

