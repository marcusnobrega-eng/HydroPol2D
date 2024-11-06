% Infiltration Routine
% Developer: Marcus Nobrega, Ph.D.
% Goal: Estimate effective precipitation and infiltration at cells domain
% Date: 8/15/2024


if k == 1 % First time-step
    if flags.flag_infiltration == 1
        % Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix ...
        %                  + BC_States.inflow + depths.d_0 - Hydro_States.ETP/24)./(time_step/60); % Inflow rate [mm/h]
        % Hydro_States.i_a = (depths.d_0 + BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow)./(time_step/60);  % Inflow rate [mm/h]
        Hydro_States.i_a = (depths.d_0)./(time_step/60);  % Inflow rate [mm/h]
        % Hydro_States.i_a = (depths.d_0 + BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix)./(time_step/60);  % Inflow rate [mm/h]


        % ---- Forward Explicit Green-Ampt ---- %
        C = Soil_Properties.ksat.*(1 +  ...
            ((depths.d_0 + Soil_Properties.psi).*(Soil_Properties.teta_sat -  ...
            Soil_Properties.teta_i))./Soil_Properties.I_0); % matrix form of Infiltration Capacity [mm/h]
        % ---- Green-Ampt ---- %
        % Infiltration rate
        if time_step*60 < 0.1 % If we are using high resolution time-steps
            % We can approximate the soil infiltration capacity curve 
            Hydro_States.f = min(C,Hydro_States.i_a); % Infiltration Rate            
            % Soil matrix mass balance
            Soil_Properties.I_t = max(Soil_Properties.I_0 +  ...
                                     (Hydro_States.f)*time_step/60,min_soil_moisture); % Extracting ETP from soil. Limiting I_t to the minimum of typically 5 mm        else
            inf_volume = nansum(nansum((Soil_Properties.I_t - Soil_Properties.I_0).*C_a)); % m3
            inf_volume_cell = Soil_Properties.I_t - Soil_Properties.I_0;
            depths.d_p = depths.d_p - inf_volume_cell; % Taking away infiltration prior calculating ETP and other fluxes
        else
            % We need to solve the implicit GA equation
            [Soil_Properties.I_t,Hydro_States.f] = GA_Newton_Raphson(Soil_Properties.I_0,time_step/60 ...
                                                   ,Soil_Properties.ksat,Soil_Properties.psi,Soil_Properties.teta_sat - Soil_Properties.teta_i,depths.d_0,Hydro_States.i_a);
            % Infiltrated Volume
            inf_volume = 1/1000*nansum(nansum((Soil_Properties.I_t - Soil_Properties.I_0).*C_a)); % m3            
            inf_volume_cell = Soil_Properties.I_t - Soil_Properties.I_p; % mm            
            % Taking away deep percolation
            Soil_Properties.I_t = max(Soil_Properties.I_t,min_soil_moisture);
            depths.d_p = depths.d_p - inf_volume_cell; % Taking away infiltration prior calculating ETP and other fluxes
        end
        if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1
            BC_States.delta_ETR_agg = zeros(size(Hydro_States.ETR));
            % Real evapotranspiration data
            Hydro_States.ETR  = input_evaporation + input_transpiration; % % Real evapotranspiration [mm/day] using the input data         
            % Areas with little soil moisture
            Hydro_States.idx_ETR = and(Soil_Properties.I_t == min_soil_moisture,depths.d_p == 0); % Adopting the minimum depth
            Hydro_States.ETR(isnan(Hydro_States.ETR)) = 0; % Fixing nan values to 0
            Hydro_States.ETR((idx_nan)) = nan; % Cells outside of the domain remain nan
            % Areas with little soil moisture but surface water
            idx_open_water =  depths.d_p > 0;
            Hydro_States.ETR(idx_open_water) =  min((depths.d_p(idx_open_water))/(time_step/60/24),Hydro_States.ETR(idx_open_water)); % Limiting to the available surface water                      
            BC_States.delta_ETR_agg(idx_open_water) = Hydro_States.ETR(idx_open_water)*(time_step/60/24); % mm
            idx_soil_available = and(depths.d_p == 0, Soil_Properties.I_t > min_soil_moisture); % Cells with no surface depth and with available soil moisture
            Hydro_States.ETR(idx_soil_available) =  min((Soil_Properties.I_t(idx_soil_available) - min_soil_moisture(idx_soil_available))/(time_step/60/24),Hydro_States.ETR(idx_soil_available)); % Limiting to the available soil moisture
            Soil_Properties.I_t(idx_soil_available) = Soil_Properties.I_t(idx_soil_available) - Hydro_States.ETR(idx_soil_available)*(time_step/60/24); % mm in the soil
            Hydro_States.ETP = zeros(size(Soil_Properties.I_t)); % Setting ETP to 0
            BC_States.delta_ETR_agg(idx_soil_available) = 0; % mm (already taken from the soil)
        elseif flags.flag_ETP == 1
            Hydro_States.idx_ETR = Soil_Properties.I_t == min_soil_moisture; % Adopting the minimum depth
            Hydro_States.ETR  = Hydro_States.ETP; % Real evapotranspiration [mm/h]
            Hydro_States.ETR(Hydro_States.idx_ETR) = Hydro_States.ETP(Hydro_States.idx_ETR) - ...
                                                      1/(time_step/60)*(min_soil_moisture(Hydro_States.idx_ETR) - (Soil_Properties.I_0(Hydro_States.idx_ETR) + (Hydro_States.f(Hydro_States.idx_ETR))*time_step/60));
        else
            BC_States.delta_ETR_agg = 0;
        end
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ...
            BC_States.inflow - inf_volume_cell + ...
            idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux ...
            - BC_States.delta_ETR_agg; % Effective precipitation within 1 computation time-step [mm]

        depths.d_t = depths.d_0 + depths.pef; % Adding pef from the previous depth.
        if min(min(depths.d_t)) < -1e-4
            factor_time = 1;
            catch_index = catch_index + 1;
            warning('Negative depths. Please reduce the time-step.')
        else
            catch_index = 1;
            factor_time = factor_time + 1;
            running_control.max_time_step = min(running_control.max_time_step*(1+0.1*factor_time),max_dt);            
            depths.d_t(depths.d_t < -1e-4) = 0;
        end
        Soil_Properties.T = Soil_Properties.Tr; % Beginning to track the replenishing time
    else
        depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ...
            BC_States.inflow + ...
            idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux; % Effective precipitation within 1 computation time-step [mm]
        depths.d_t = depths.d_0 + depths.pef;
        Hydro_States.f = 0*size(depths.d_0);
    end
else
    if flags.flag_infiltration == 1
        % Inflow Rate 
        % Hydro_States.i_a = (depths.d_p + BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow)./(time_step/60);  % Inflow rate [mm/h]
        Hydro_States.i_a = (depths.d_p)./(time_step/60);  % Inflow rate [mm/h]
        % Hydro_States.i_a = (depths.d_p + BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix)./(time_step/60);  % Inflow rate [mm/h]
        
        % Infiltration Capacity
        C = Soil_Properties.ksat.*(1 + ((depths.d_p + ...
            Soil_Properties.psi).*(Soil_Properties.teta_sat - Soil_Properties.teta_i))./Soil_Properties.I_p); % matrix form of Infiltration Capacity [mm/h]
        
        % Cells with top layer exceeding root zone 
        Non_C_idx = Soil_Properties.I_t > Soil_Properties.Lu*1000; % Cells that exceed the top layer infiltrated depth
        
        % Limiting infiltration in these areas
        C(Non_C_idx) = Soil_Properties.ksat(Non_C_idx); % If that condition ocurrs, reduce the capacity to the saturation
        
        % Cells where infiltration capacity is higher than inflow
        Hydro_States.idx_C = (Hydro_States.i_a <= C); % Cells where the inflow rate is smaller than the infiltration capacity
        
        % Cells where the inflow rate is larger than the infiltration capacity 
        idx_i_a = (Hydro_States.i_a > C); % Cells where the inflow rate is larger than the infiltration capacity
        
        % Recovery time
        Soil_Properties.T = Soil_Properties.T - time_step/60; % Recoverying time (hours)
        Soil_Properties.T(idx_i_a) = Soil_Properties.Tr(idx_i_a); % Saturated Areas we begin the process again
        Soil_Properties.idx_T = Soil_Properties.T < 0; % Cells where the replenishing time is reached
        
        % Infiltration rate
        if time_step*60 < 5*60 % If we are using high resolution time-steps
            % We can approximate the soil infiltration capacity curve 
            Hydro_States.f = min(C,Hydro_States.i_a); % Infiltration rate (mm/hr)
            
            % Soil matrix mass balance
            Soil_Properties.I_t = max(Soil_Properties.I_p + Hydro_States.f*time_step/60 - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
            inf_volume_cell = Soil_Properties.I_t - Soil_Properties.I_p;
            inf_volume = 1/1000*nansum(nansum((Soil_Properties.I_t - Soil_Properties.I_p).*C_a)); % m3   
            depths.d_p = depths.d_p - inf_volume_cell; % Taking away infiltration prior calculating ETP and other fluxes
        else
            % We need to solve the implicit GA equation
            [Soil_Properties.I_t,Hydro_States.f] = GA_Newton_Raphson(Soil_Properties.I_p,time_step/60 ...
                                                   ,Soil_Properties.ksat,Soil_Properties.psi,Soil_Properties.teta_sat - Soil_Properties.teta_i,depths.d_p,Hydro_States.i_a);
            % Infiltrated Volume
            inf_volume = 1/1000*nansum(nansum((Soil_Properties.I_t - Soil_Properties.I_p).*C_a)); % m3
            inf_volume_cell = Soil_Properties.I_t - Soil_Properties.I_p; % mm
            % Taking away deep percolation
            Soil_Properties.I_t = max(Soil_Properties.I_t - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
            depths.d_p = depths.d_p - inf_volume_cell; % Taking away infiltration prior calculating ETP and other fluxes
        end
            
        % Limiting evapotranspiration in soils with low moisture
        if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1
            BC_States.delta_ETR_agg = zeros(size(Hydro_States.ETR));
            % Real evapotranspiration data
            Hydro_States.ETR  = input_evaporation + input_transpiration; % % Real evapotranspiration [mm/day] using the input data         
            % Areas with little soil moisture
            Hydro_States.idx_ETR = and(Soil_Properties.I_t == min_soil_moisture,depths.d_p == 0); % Adopting the minimum depth
            Hydro_States.ETR(isnan(Hydro_States.ETR)) = 0; % Fixing nan values to 0
            Hydro_States.ETR((idx_nan)) = nan; % Cells outside of the domain remain nan
            % Areas with little soil moisture but surface water
            idx_open_water =  depths.d_p > 0;
            Hydro_States.ETR(idx_open_water) =  min((depths.d_p(idx_open_water))/(time_step/60/24),Hydro_States.ETR(idx_open_water)); % Limiting to the available surface water                      
            BC_States.delta_ETR_agg(idx_open_water) = Hydro_States.ETR(idx_open_water)*(time_step/60/24); % mm
            idx_soil_available = and(depths.d_p == 0, Soil_Properties.I_t > min_soil_moisture); % Cells with no surface depth and with available soil moisture
            Hydro_States.ETR(idx_soil_available) =  min((Soil_Properties.I_t(idx_soil_available) - min_soil_moisture(idx_soil_available))/(time_step/60/24),Hydro_States.ETR(idx_soil_available)); % Limiting to the available soil moisture
            Soil_Properties.I_t(idx_soil_available) = Soil_Properties.I_t(idx_soil_available) - Hydro_States.ETR(idx_soil_available)*(time_step/60/24); % mm in the soil
            Hydro_States.ETP = zeros(size(Soil_Properties.I_t)); % Setting ETP to 0
            BC_States.delta_ETR_agg(idx_soil_available) = 0; % mm (already taken from the soil)
        elseif flags.flag_ETP == 1
            Hydro_States.idx_ETR = Soil_Properties.I_t == min_soil_moisture; % Adopting the minimum depth
            Hydro_States.ETR  = Hydro_States.ETP; % Real evapotranspiration [mm/h]
            Hydro_States.ETR(Hydro_States.idx_ETR) = Hydro_States.ETP(Hydro_States.idx_ETR) - ...
                                                      1/(time_step/60)*(min_soil_moisture(Hydro_States.idx_ETR) - (Soil_Properties.I_0(Hydro_States.idx_ETR) + (Hydro_States.f(Hydro_States.idx_ETR))*time_step/60));
        else
            BC_States.delta_ETR_agg = 0;
        end
        % Refreshing Soil_Properties.I_t to I_0 for cases where idx_T > 0
        Soil_Properties.I_t(Soil_Properties.idx_T) = min_soil_moisture(Soil_Properties.idx_T);

        % Effective Precipitation at Domain Cells
        depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ...
                     BC_States.inflow + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux ...
                   - BC_States.delta_ETR_agg;
        % depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix - Hydro_States.f*time_step/60 + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
        depths.d_t = depths.d_p + depths.pef; %% ATTENTION HERE

        if min(min(depths.d_t)) < -1e-8
            catch_index = catch_index + 1;
            warning('Negative depths. Please reduce the time-step.')
        else
            depths.d_t(depths.d_t < 1e-6) = 0;
        end        
    else
        Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow + depths.d_p - Hydro_States.ETP/24)./(time_step/60);
        depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.ETP/24  + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
        % depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix - Hydro_States.ETP/24  + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
        depths.d_t = depths.d_p + depths.pef;
        Hydro_States.f = 0*size(depths.d_p);
    end
end

% Effective Precipitation (Available depth at the beginning of the time-step
depths.d_tot = depths.d_t;

% Correcting depths in case it reach overbanks
if flags.flag_subgrid == 1 % Maybe we have a change from inbank <-> overbank
    % Eq. 15 and 16 of A subgrid channel model for simulating river hydraulics andfloodplain inundation over large and data sparse areas
    % Inbank - Overbank
    idx = depths.d_tot/1000 > Wshed_Properties.River_Depth & depths.d_p/1000 <= Wshed_Properties.River_Depth & idx_rivers; % Cells in which there is a change from inbank to overbank
    if sum(sum(idx)) > 0
        depths.d_tot(idx) = 1000*(depths.d_tot(idx)/1000 - (depths.d_tot(idx)/1000 - Elevation_Properties.elevation_cell(idx) + (Elevation_Properties.elevation_cell(idx) - Wshed_Properties.River_Depth(idx))).*(1 - Wshed_Properties.River_Width(idx)/Wshed_Properties.Resolution)); % mm 
        C_a(idx) = Wshed_Properties.Resolution^2;
    end
    % Overbank - Inbank
    idx = depths.d_tot/1000 <= Wshed_Properties.River_Depth & depths.d_p/1000 >= Wshed_Properties.River_Depth & idx_rivers; % Cells in which there is a change from overbank to inbank
    if sum(sum(idx)) > 0
        factor = depths.d_tot(idx)/1000 - Elevation_Properties.elevation_cell(idx) + (Elevation_Properties.elevation_cell(idx) - Wshed_Properties.River_Depth(idx));
        depths.d_tot(idx) = 1000*(depths.d_tot(idx)/1000 + (Wshed_Properties.Resolution*(factor))./Wshed_Properties.River_Width(idx) ...
            - (depths.d_tot(idx)/1000 - Elevation_Properties.elevation_cell(idx) + (Elevation_Properties.elevation_cell(idx) - Wshed_Properties.River_Depth(idx)))); % mm
        C_a(idx) = Wshed_Properties.River_Width(idx)*Wshed_Properties.Resolution;
    end 
end
