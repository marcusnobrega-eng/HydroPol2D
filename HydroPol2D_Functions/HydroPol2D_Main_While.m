%% Main HydroPol2D While Loop
% Developer: Marcus Nobrega
% Date 5/21/2023
% Goal - Run the main modeling process of the model

% Loading Input File in case you want to avoid doing all preprocessing
% clear all
% load workspace_14_de_julho.mat
% clear all
% load workspace_franquinho.mat
% save('preprocessing_input.mat');
format short g

tic 
k = 1; % time-step counter
C = 0; % initial infiltration capacity
t = running_control.min_time_step; % inital time
time_step = running_control.min_time_step/60; % time-step of the model in min
update_spatial_BC;
Flooded_Area = 0; % initial flooded area
velocities.velocity_raster = 0; % initial velocity
Risk_Area = 0; % initial risk area
saver_count = 1; % starts in 1 but the next pointer should be 2, this is auto fixed when t reach the next time aggregation.
store = 1; % (meaning, luis?)
%flags.flag_inertial = 1; % Using Inertial Model
t_previous = 0;
factor_time = 0;
max_dt = running_control.max_time_step;
    current_storage = nansum(nansum(Wshed_Properties.cell_area*depths.d_t/1000)); % m3
if flags.flag_inertial == 1
    outflow_prev = outflow_bates;
end
catch_index = 1;
% ---- Plotting Results in Real-Time ---- %
% n_snaps = 10; % Number of plots. Time will be divided equally
% dt_snap = running_control.routing_time/n_snaps; time_snap = [1:1:n_snaps]*dt_snap; z2_snap = 0;
if flags.flag_dashboard == 1
    ax.flags = flags;
    ax = HydroPol2D_running_dashboard(ax,Maps, zeros(size(DEM_raster.Z)), DEM_raster, ...
        subsref({gauges, extra_parameters.gauges}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})),...
        1,1);
end
% ---- Main Loop --- %
while t <= (running_control.routing_time + running_control.min_time_step/60) % Running up to the end of the simulation
    try
    % Snapshot results Results
    % z1_snap = find(time_snap>=t,1,'first');
    % if z1_snap > z2_snap
    %    Snapshot_Results
    %    pause(0.05);
    % end
    % z2_snap = z1_snap;

    % Infiltration and Effective Precipitation Calculation
    % Show stats
    
    if tmin_wq < 0 || isnan(tmin_wq) || isinf(tmin_wq)
        error('Instability. in the Water Quality Model.')
    end
    
    if k == 1 % First time-step
        if flags.flag_infiltration == 1
            % Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix ...
            %                  + BC_States.inflow + depths.d_0 - Hydro_States.ETP/24)./(time_step/60); % Inflow rate [mm/h]
            Hydro_States.i_a = (depths.d_0)./(time_step/60); % Inflow rate [mm/h]            
            C = Soil_Properties.ksat.*(1 +  ...
                                      ((depths.d_0 + Soil_Properties.psi).*(Soil_Properties.teta_sat -  ...
                                      Soil_Properties.teta_i))./Soil_Properties.I_0); % matrix form of Infiltration Capacity [mm/h]
            Hydro_States.f = min(C,Hydro_States.i_a); % Infiltration Rate
            Soil_Properties.I_t = max(Soil_Properties.I_0 +  ...
                                     (Hydro_States.f)*time_step/60,min_soil_moisture); % Extracting ETP from soil. Limiting I_t to the minimum of typically 5 mm
            Hydro_States.idx_ETR = Soil_Properties.I_t == min_soil_moisture; % Adopting the minimum depth 
            Hydro_States.ETR  = Hydro_States.ETP; % Real evapotranspiration [mm/h]
            Hydro_States.ETR (Hydro_States.idx_ETR) = Hydro_States.ETP(Hydro_States.idx_ETR) - (min_soil_moisture(Hydro_States.idx_ETR) - (Soil_Properties.I_0(Hydro_States.idx_ETR) + (Hydro_States.f(Hydro_States.idx_ETR))*time_step/60));
            % depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ... 
            %              BC_States.inflow - Hydro_States.f*time_step/60 + ...
            %              idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux; % Effective precipitation within 1 computation time-step [mm]
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ... 
                         - Hydro_States.f*time_step/60 + ...
                         idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux; % Effective precipitation within 1 computation time-step [mm]            
            depths.d_t = depths.d_0 + depths.pef; % Adding pef from the previous depth.
            if min(min(depths.pef)) < -1e-8
                error('Negative depths. Please reduce the time-step.')
            else
                depths.d_t(depths.d_t < 1e-6) = 0;
            end
            Soil_Properties.T = Soil_Properties.Tr; % Beginning to track the replenishing time
        else
            % depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ...
            %              BC_States.inflow + ...
            %              idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux; % Effective precipitation within 1 computation time-step [mm]
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + ...
                         idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux; % Effective precipitation within 1 computation time-step [mm]            
            depths.d_t = depths.d_0 + depths.pef;
        end
    else
        if flags.flag_infiltration == 1
            % Effective precipitation - Green-Ampt(1911)
            % Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix +  ...
            %                     BC_States.inflow + depths.d_p - Hydro_States.ETP/24)./(time_step/60);  % Inflow rate [mm/h]  
            Hydro_States.i_a = (depths.d_p)./(time_step/60);  % Inflow rate [mm/h]             
            C = Soil_Properties.ksat.*(1 + ((depths.d_p + ...
                                       Soil_Properties.psi).*(Soil_Properties.teta_sat - Soil_Properties.teta_i))./Soil_Properties.I_p); % matrix form of Infiltration Capacity [mm/h]
            Non_C_idx = Soil_Properties.I_t > Soil_Properties.Lu*1000; % Cells that exceed the top layer infiltrated depth
            C(Non_C_idx) = Soil_Properties.ksat(Non_C_idx); % If that condition ocurrs, reduce the capacity to the saturation
            
            Hydro_States.idx_C = (Hydro_States.i_a <= C); % Cells where the inflow rate is smaller than the infiltration capacity 
            idx_i_a = (Hydro_States.i_a > C); % Cells where the inflow rate is larger than the infiltration capacit
            Soil_Properties.T = Soil_Properties.T - time_step/60; % Recoverying time (hours)
            Soil_Properties.T(idx_i_a) = Soil_Properties.Tr(idx_i_a); % Saturated Areas we begin the process again
            Soil_Properties.idx_T = Soil_Properties.T < 0; % Cells where the replenishing time is reached
            Hydro_States.f = min(C,Hydro_States.i_a); % Infiltration rate (mm/hr)
            Soil_Properties.I_t = max(Soil_Properties.I_p + Hydro_States.f*time_step/60 - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
            Hydro_States.idx_ETR = Soil_Properties.I_t == min_soil_moisture;
            Hydro_States.ETR  = Hydro_States.ETP;
            Hydro_States.ETR (Hydro_States.idx_ETR) = Hydro_States.ETP(Hydro_States.idx_ETR) - (min_soil_moisture(Hydro_States.idx_ETR) - (Soil_Properties.I_p(Hydro_States.idx_ETR) + Hydro_States.f(Hydro_States.idx_ETR)*time_step/60 - Soil_Properties.k_out(Hydro_States.idx_ETR)*time_step/60));
            % Refreshing Soil_Properties.I_t to I_0 for cases where idx_T > 0
            Soil_Properties.I_t(Soil_Properties.idx_T) = min_soil_moisture(Soil_Properties.idx_T);
            % depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.f*time_step/60 + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix - Hydro_States.f*time_step/60 + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_p + depths.pef; %% ATTENTION HERE
        else
            Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow + depths.d_p - Hydro_States.ETP/24)./(time_step/60);
            % depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.ETP/24  + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix - Hydro_States.ETP/24  + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;            
            depths.d_t = depths.d_p + depths.pef;
        end
    end

    % Effective Precipitation (Available depth at the beginning of the time-step
    depths.d_tot = depths.d_t;

    if flags.flag_D8 == 1 && flags.flag_diffusive == 1
        %%%% WATER DEPTHS %%%
        depths.d_left_cell = [zeros(ny,1),depths.d_tot(:,1:(nx-1))];
        depths.d_right_cell = [depths.d_tot(:,(2:(nx))) zeros(ny,1)];
        depths.d_up_cell = [zeros(1,nx) ; depths.d_tot(1:(ny-1),:)];
        depths.d_down_cell = [depths.d_tot(2:(ny),:) ; zeros(1,nx)];
        if flags.flag_D8 == 1
            % Simulate with D-8 flow direction
            % --- Adding NE, SE, SW, NW --- %
            depths.d_NE_t(2:(ny),1:(nx-1)) = depths.d_tot(1:(ny-1),2:nx); % OK
            depths.d_SE_t(1:(ny-1),1:(nx-1)) = depths.d_tot(2:ny,2:nx); % OK
            depths.d_SW_t(1:(ny-1),2:(nx)) = depths.d_tot(2:(ny),1:(nx-1)); % OK
            depths.d_NW_t(2:ny,2:nx) = depths.d_tot(1:(ny-1),1:(nx-1)); % OK
        else
            % Simulate with D-4 flow direction
            % Everything already calculated
        end
    end
    % 
    %% Flood Routing (Cellular Automata or Fully Hydrodynamic Model)
    if flags.flag_D8 == 1 % D-8
        if flags.flag_inertial ~= 1
            [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,flow_rate.qout_ne_t,flow_rate.qout_se_t,flow_rate.qout_sw_t,flow_rate.qout_nw_t,depths.d_t,CA_States.I_tot_end_cell] = ...
                CA_Routing_8D(Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical,CA_States.depth_tolerance);
        else
            % -------------------- Local Inertial Formulation ----------------%
            [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf] = ...
                Bates_Inertial_8D(Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                depths.d_tot, depths.d_p,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,Wshed_Properties.Resolution,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,CA_States.depth_tolerance,outflow_prev,idx_nan,flags.flag_critical);
        end
    else % 4-D
        % CA
        if flags.flag_inertial ~= 1
            [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell] = ...
                CA_Routing(Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical);
        else
            % -------------------- Local Inertial Formulation ----------------%
            [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf] = ...
                Bates_Inertial_4D(Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                depths.d_tot, depths.d_p,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,Wshed_Properties.Resolution,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,CA_States.depth_tolerance,outflow_prev,idx_nan,flags.flag_critical);
        end
    end


    %% Outflows become Inflows
    flow_rate.qin_left_t = [zeros(ny,1),flow_rate.qout_right_t(:,1:(nx-1))];
    flow_rate.qin_right_t = [flow_rate.qout_left_t(:,(2:(nx))) zeros(ny,1)];
    flow_rate.qin_up_t = [zeros(1,nx) ; flow_rate.qout_down_t(1:(ny-1),:)];
    flow_rate.qin_down_t = [flow_rate.qout_up_t(2:(ny),:) ; zeros(1,nx)];

    % Inflows - Inclined Directions
    if flags.flag_D8 == 1
        if flags.flag_GPU == 1
            zero_matrix = gpuArray(zeros(size(Elevation_Properties.elevation_cell)));
        else
            zero_matrix = zeros(size(Elevation_Properties.elevation_cell));
        end
        flow_rate.qin_ne_t = zero_matrix; flow_rate.qin_se_t = zero_matrix; flow_rate.qin_sw_t = zero_matrix; flow_rate.qin_nw_t = zero_matrix;
        flow_rate.qin_ne_t(2:(ny),1:(nx-1)) = flow_rate.qout_sw_t(1:(ny-1),2:(nx)); % OK
        flow_rate.qin_se_t(1:(ny-1),1:(nx-1)) = flow_rate.qout_nw_t(2:ny,2:nx); % OK
        flow_rate.qin_sw_t(1:(ny-1),2:(nx)) = flow_rate.qout_ne_t(2:(ny),1:(nx-1)); % OK
        flow_rate.qin_nw_t(2:ny,2:nx) = flow_rate.qout_se_t(1:(ny-1),1:(nx-1)); % OK
    end

    if flags.flag_D8 == 1
        flow_rate.qin_t = flow_rate.qin_left_t + flow_rate.qin_right_t + flow_rate.qin_up_t + flow_rate.qin_down_t ...
                        + flow_rate.qin_ne_t + flow_rate.qin_se_t + flow_rate.qin_sw_t + flow_rate.qin_nw_t;
    else
        flow_rate.qin_t = flow_rate.qin_left_t + flow_rate.qin_right_t + flow_rate.qin_up_t + flow_rate.qin_down_t;
    end
    
    % Including Source Terms (inflow hydrograph or water taken from the
    % domain)
    flow_rate.qin_t = flow_rate.qin_t + BC_States.inflow/(time_step/60); % mm/h
    idx3 = logical(isnan(flow_rate.qin_t) + isinf(flow_rate.qin_t)); % Mask to take out non cells
    flow_rate.qin_t(idx3) = 0; % No flow at these cells

    % Mass Balance Equation (outflow already taken)
    depths.d_t = depths.d_t + flow_rate.qin_t*time_step/60;

    % Water Balance Error
    water_balance_error_volume = abs(sum(sum(depths.d_t(depths.d_t<0))))*Wshed_Properties.Resolution^2*0.001; % m3
    water_balance_error_mm = water_balance_error_volume/Wshed_Properties.drainage_area*1000; % mm

    if water_balance_error_volume > 0.5 % We need to define better this parameter
        factor_time = 1;
        catch_index = catch_index + 1;
        error('Mass balance error too high.')
    else
        catch_index = 1;
        factor_time = factor_time + 1;
        running_control.max_time_step = min(running_control.max_time_step*(1+0.1*factor_time),max_dt);        
    end

    depths.d_t = max(depths.d_t,0);  % Taking away negative masses
    

    % ------ Water Quality Modeling ----- %
    % Water Quality Parameters for f(B(t))
    if flags.flag_waterquality == 1
        if flags.flag_D8 ~= 1
            [WQ_States.B_t,WQ_States.P_conc,Out_Conc,tmin_wq,tot_W_out,WQ_States.mass_lost,WQ_States.Tot_Washed] = build_up_wash_off(LULC_Properties.C_3,LULC_Properties.C_4,flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,WQ_States.B_t,time_step,nx,ny,Wshed_Properties.cell_area,outlet_index,idx_nan_5,flags.flag_wq_model,WQ_States.mass_lost,WQ_States.Tot_Washed,LULC_Properties.Bmin,LULC_Properties.Bmax,LULC_Properties.min_Bt);
        else
            [WQ_States.B_t,WQ_States.P_conc,Out_Conc,tmin_wq,tot_W_out,WQ_States.mass_lost,WQ_States.Tot_Washed] = build_up_wash_off_8D(LULC_Properties.C_3,LULC_Properties.C_4,flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,flow_rate.qout_ne_t,flow_rate.qout_se_t,flow_rate.qout_sw_t,flow_rate.qout_nw_t,WQ_States.B_t,time_step,nx,ny,Wshed_Properties.cell_area,outlet_index,idx_nan_5,flags.flag_wq_model,WQ_States.mass_lost,WQ_States.Tot_Washed,LULC_Properties.Bmin,LULC_Properties.Bmax,LULC_Properties.min_Bt);
        end
    end
    %%%% Checking Mass Balance
    if flags.flag_waterquality == 1
        if sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t)))) > 1.2*initial_mass  % More than 5%
            error('Brutal instability in B(t). More than 20% difference')
        end
    end

    %% Refreshing Time-step
    % running_control.pos_save = find(running_control.time_change_records < t,1,'last');
    % running_control.time_save = running_control.time_change_records(running_control.pos_save); % min
    % running_control.delta_time_save = running_control.time_save - running_control.time_save_previous;
    % running_control.time_save_previous = running_control.time_save;
    % running_control.actual_record_timestep = find(running_control.time_change_records < t,1,'last');

    running_control.pos_save = ceil((t*60)/running_control.time_step_change);
    running_control.time_save = (running_control.pos_save - 1)*running_control.time_step_change/60;
    running_control.delta_time_save = running_control.time_save - running_control.time_save_previous;
    running_control.time_save_previous = running_control.time_save;
    running_control.actual_record_timestep = ceil((t*60)/running_control.time_step_change);

    % Runoff Coefficient Calculation
    BC_States.outflow_volume  = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600*time_step*60 + BC_States.outflow_volume ;
    if flags.flag_spatial_rainfall == 1
        if flags.flag_inflow == 1
            inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area)) + ...
                         nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area;
        else
            inflow_vol = nansum(nansum(nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area));
        end
            BC_States.inflow_volume = inflow_vol + BC_States.inflow_volume; % m3        
    elseif flags.flag_spatial_rainfall ~= 1 && flags.flag_inflow == 1
        inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area)) + BC_States.delta_p_agg/1000*Wshed_Properties.drainage_area;
        BC_States.inflow_volume = inflow_vol +  BC_States.inflow_volume; % check future
    else
        inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area));
        BC_States.inflow_volume = inflow_vol + BC_States.delta_p_agg/1000*Wshed_Properties.drainage_area + BC_States.inflow_volume; % check future
    end

    % Storage Calculation
    previous_storage = current_storage;
    current_storage = nansum(nansum(Wshed_Properties.cell_area*depths.d_t/1000)); % m3
    if flags.flag_infiltration == 0
        Hydro_States.f = zeros(size(outlet_states.outlet_flow,1),size(outlet_states.outlet_flow,2));
    end    
    % dS/dt = Qin - Qout = Rain + Inflow - Outflow - infiltration
    delta_storage = current_storage - previous_storage;
    outflow_flux = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600 + ...
                    nansum(nansum(Hydro_States.f/1000/3600*Wshed_Properties.cell_area)); % m3 per sec
    flux_volumes = inflow_vol - outflow_flux*time_step*60; % dt(Qin - Qout) m3
    volume_error = flux_volumes - delta_storage;    

    % Introducing the volume error to the inflow cells
    if flags.flag_inflow == 1
        mask = logical(Wshed_Properties.inflow_cells);
        if any(any(mask))
            % depths.d_t(mask) = depths.d_t(mask) - 1000*(volume_error)/sum(sum(mask))/Wshed_Properties.cell_area;
        end
    end
    % Refreshing time-step script
    refreshing_timestep;
    
    %% Updating Boundary Conditions
    update_spatial_BC; % Updating rainfall and etr B.C.

    % Saving Plotting Values - Recording Time
    % Maps of Flood Depths, WSE and Pollutant Concentrations
    % --- Calculating EMC --- %
    if  flags.flag_automatic_calibration ~= 1
        if flags.flag_waterquality == 1
            WQ_States.mass_outlet = max(WQ_States.mass_outlet + Out_Conc*((nansum(nansum(outlet_states.outlet_flow)/1000/3600*1000)))*(time_step*60),0); % mg
            WQ_States.vol_outlet = max((nansum(nansum(outlet_states.outlet_flow))/1000/3600*1000)*(time_step*60) + WQ_States.vol_outlet,0);
        end
    end

    % Previous Time-step
    if k == 1
        t_previous = running_control.time_calculation_routing(k,1)/60;
    else
        t_previous = t;
    end
    % Current time
    t_save = t + running_control.time_calculation_routing(k,1)/60;

    %% Saving Output Maps
    save_output_maps;

    % Clearing stored values
    if flags.flag_GPU == 1
        flow_rate.qout_left_t = gpuArray(zeros(ny,nx));
        flow_rate.qout_right_t = gpuArray(zeros(ny,nx));
        flow_rate.qout_up_t = gpuArray( zeros(ny,nx));
        flow_rate.qout_down_t = gpuArray(zeros(ny,nx));
        flow_rate.qout_ne_t = gpuArray(zeros(ny,nx));
        flow_rate.qout_se_t = gpuArray(zeros(ny,nx));
        flow_rate.qout_sw_t = gpuArray(zeros(ny,nx));
        flow_rate.qout_nw_t = gpuArray(zeros(ny,nx));
        flow_rate.qin_t = gpuArray(zeros(ny,nx));
    else
        flow_rate.qout_left_t = zeros(ny,nx);
        flow_rate.qout_right_t = zeros(ny,nx);
        flow_rate.qout_up_t = zeros(ny,nx);
        flow_rate.qout_down_t = zeros(ny,nx);
        flow_rate.qout_ne_t = zeros(ny,nx);
        flow_rate.qout_se_t = zeros(ny,nx);
        flow_rate.qout_sw_t = zeros(ny,nx);
        flow_rate.qout_nw_t = zeros(ny,nx);
        flow_rate.qin_t = zeros(ny,nx);
    end

    % Previous Depths and Moisture
    depths.d_p = depths.d_t;

    Soil_Properties.I_p = Soil_Properties.I_t;
    if flags.flag_inertial == 1
        % saving previous outflows
        outflow_prev = outflow_bates; % Corrected previous outflow
    end    

    % Saving Results in time_observations - Only Valid for Calibration
    save_automatic_cabralition_outputs;

    % Refreshing Time-step
    t = running_control.time_calculation_routing(k,1)/60 + t;
    time_step_save(k,2) = running_control.time_calculation_routing(k,1);
    time_step_save(k,1) = t;
    k = k + 1;

    % Simulating Dams failure    
    if flags.flag_dam_break
        if depths.d_t(gauges.northing_obs_gauges(1),gauges.easting_obs_gauges(1))/1000 > 24 && flag_break_1 == 1
            for i_breaker = 1:57
                Elevation_Properties.elevation_cell(breakers.northing(i_breaker),breakers.easting(i_breaker)) = Elevation_Properties.elevation_cell(gauges.northing_obs_gauges(1),gauges.easting_obs_gauges(1));
            end
            flag_break_1 = 0;
        elseif depths.d_t(gauges.northing_obs_gauges(10),gauges.easting_obs_gauges(10))/1000 > 10.5 && flag_break_2 == 1
            for i_breaker =58:length(breakers.northing)
                Elevation_Properties.elevation_cell(breakers.northing(i_breaker),breakers.easting(i_breaker)) = Elevation_Properties.elevation_cell(gauges.northing_obs_gauges(10),gauges.easting_obs_gauges(10));
            end
            flag_break_2 = 0;    
        end
    end
    
    % Show Stats
    if flags.flag_waterquality == 1
        perc______duremain______tsec_______dtmm______infmmhr____CmgL_____dtmWQ___VolErrorm3 = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(Hydro_States.f)), max(max((WQ_States.P_conc))), tmin_wq,volume_error]
    else
        per__duremain__tsec___dtmm__infmmhr__Vmax___VolErrorm3 = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(Hydro_States.f)), max(max(velocities.velocity_raster)),volume_error]
    end

    catch ME % Reduce the time-step
        t = t - time_step; 
        t_previous = t;
        % time_step = time_step/2; % min
        wave_celerity = sqrt(9.81*(max(max(max(depths.d_tot/1000)),max(max(depths.d_p/1000))))); % Using d_p, which is the depth at the end of the time-step
        max_vel = max(max(velocities.velocity_raster));
        factor = 1/catch_index;
        new_timestep = factor*(min(Courant_Parameters.alfa_min*Wshed_Properties.Resolution./(max_vel+wave_celerity))); % alpha of 0.4
        % new_timestep = factor*(min(0.25*Wshed_Properties.Resolution./(wave_celerity))); % alpha of 0.4
        
        % dt_water_balance = min(min(depths.d_p/1000*Wshed_Properties.cell_area./(CA_States.I_tot_end_cell/(time_step*60)))); % sec
        dt_water_balance = new_timestep;        
        new_timestep = min(new_timestep,dt_water_balance);
        if catch_index > 1
            running_control.max_time_step = new_timestep*1.5; % Assuming a smaller max time-step to avoid large integrations
        end
        new_timestep = min(new_timestep,running_control.max_time_step); % sec
        time_step = new_timestep/60; % min
        t = t + time_step;
        depths.d_t = depths.d_p;
        Soil_Properties.I_t = Soil_Properties.I_p;
        update_spatial_BC
        % outflow_bates = outflow_prev; % Not actually necessary anymore
    end    
end

% Cloasing the dashboard
if flags.flag_dashboard == 1
    delete(ax.app)
end
% Saving the last modeled data
Maps.Hydro.d=Maps.Hydro.d(:,:,1:saver_count);
if flags.flag_infiltration == 1
    Maps.Hydro.I_t=Maps.Hydro.I_t(:,:,1:saver_count);
end
if flags.flag_rainfall > 0 && flags.flag_spatial_rainfall == 1
    Maps.Hydro.spatial_rainfall_maps=Maps.Hydro.spatial_rainfall_maps(:,:,1:saver_count);
end
if flags.flag_human_instability == 1
    Maps.Hydro.risk=Maps.Hydro.risk(:,:,1:saver_count);
elseif flags.flag_human_instability == 2
elseif flags.flag_human_instability == 3
    Maps.Hydro.risk_cm = Maps.Hydro.risk_cm(:,:,1:saver_count);
    Maps.Hydro.risk_tm = Maps.Hydro.risk_tm(:,:,1:saver_count);
    Maps.Hydro.risk_am = Maps.Hydro.risk_am(:,:,1:saver_count);
    Maps.Hydro.risk_om = Maps.Hydro.risk_om(:,:,1:saver_count);
    Maps.Hydro.risk_cf = Maps.Hydro.risk_cf(:,:,1:saver_count);
    Maps.Hydro.risk_tf = Maps.Hydro.risk_tf(:,:,1:saver_count);
    Maps.Hydro.risk_af = Maps.Hydro.risk_af(:,:,1:saver_count);
    Maps.Hydro.risk_of = Maps.Hydro.risk_of(:,:,1:saver_count);
end
if flags.flag_ETP == 1
    Maps.Hydro.ETP_save=Maps.Hydro.ETP_save(:,:,1:saver_count);
end
if flags.flag_waterquality == 1
    Maps.WQ_States.Pol_Conc_Map=Maps.WQ_States.Pol_Conc_Map(:,:,1:saver_count);
    Maps.WQ_States.Pol_mass_map=Maps.WQ_States.Pol_mass_map(:,:,1:saver_count);
end
save(strcat('Temporary_Files\save_map_hydro_',num2str(store),'.mat'),'Maps');

% Returning Variables to CPU
if flags.flag_GPU == 1
    % Structure Arrays
    BC_States = structfun(@gather, BC_States, 'UniformOutput', false);
    CA_States = structfun(@gather, CA_States, 'UniformOutput', false);
    Courant_Parameters = structfun(@gather, Courant_Parameters, 'UniformOutput', false);
    depths = structfun(@gather, depths, 'UniformOutput', false);
    Elevation_Properties = structfun(@gather, Elevation_Properties, 'UniformOutput', false);
    flags = structfun(@gather, flags, 'UniformOutput', false);
    % if flags.flag_D8 == 1
    %     wse_slope_zeros = structfun(@gather, wse_slope_zeros, 'UniformOutput', false);
    %     Distance_Matrix = structfun(@gather, Distance_Matrix, 'UniformOutput', false);
    % end
    % Gauges Label
    if flags.flag_obs_gauges == 1
        gauges = structfun(@gather, gauges, 'UniformOutput', false);
        gauges.labels_observed_string = extra_parameters.gauges.label_observed_string;
    end
    GIS_data = structfun(@gather, GIS_data, 'UniformOutput', false);
    if flags.flag_human_instability > 0
        Human_Instability = structfun(@gather, Human_Instability, 'UniformOutput', false);
    end
    Hydro_States = structfun(@gather, Hydro_States, 'UniformOutput', false);
    Inflow_Parameters = structfun(@gather, Inflow_Parameters, 'UniformOutput', false);
    LULC_Properties = structfun(@gather, LULC_Properties, 'UniformOutput', false);
    Rainfall_Parameters = structfun(@gather, Rainfall_Parameters, 'UniformOutput', false);
    recording_parameters = structfun(@gather, recording_parameters, 'UniformOutput', false);
    running_control = structfun(@gather, running_control, 'UniformOutput', false);
    Soil_Properties = structfun(@gather, Soil_Properties, 'UniformOutput', false);
    if flags.flag_waterquality == 1
        WQ_States = structfun(@gather, WQ_States, 'UniformOutput', false);
    end
    Wshed_Properties = structfun(@gather, Wshed_Properties, 'UniformOutput', false);
    outlet_states = structfun(@gather, outlet_states, 'UniformOutput', false);
    velocities = structfun(@gather, velocities, 'UniformOutput', false);
    if flags.flag_spatial_rainfall == 1
        Spatial_Rainfall_Parameters  = structfun(@gather, Spatial_Rainfall_Parameters, 'UniformOutput', false);
    end
    if flags.flag_ETP == 1
        ETP_Parameters  = structfun(@gather, ETP_Parameters, 'UniformOutput', false);
    end

    % Double Arrays
    elevation = gather(elevation);
    idx_nan = gather(idx_nan);
    if flags.flag_waterquality == 1
        idx_nan_5 = gather(idx_nan_5);
    end
    idx_outlet = gather(idx_outlet);
    k = gather(k);
    nx = gather(nx);
    if flags.flag_waterquality == 1
        Out_Conc = gather(Out_Conc);
    end
    outlet_index = gather(outlet_index);
    outlet_runoff_volume = gather(outlet_runoff_volume);
    outlet_type = gather(outlet_type);
    slope_outlet = gather(slope_outlet);
    t = gather(t);
    t_previous = gather(t_previous);
    time_calculation_routing = gather(time_calculation_routing);
    time_step = gather(time_step);
    time_step_model = gather(time_step_model);
    tmin_wq = gather(tmin_wq);
    C = gather(C);
    min_soil_moisture = gather(min_soil_moisture);
end


function [M_left,M_right,M_up,M_down,M_NE,M_SE,M_SW,M_NW] = forward_matrices(M,nx,ny,flag_D8)
    % Generates Forward Matrices of a given Matix
    M_left = [zeros(ny,1),M(:,1:(nx-1))];
    M_right = [M(:,(2:(nx))) zeros(ny,1)];
    M_up = [zeros(1,nx) ; M(1:(ny-1),:)];
    M_down = [M(2:(ny),:) ; zeros(1,nx)];
    if flag_D8 == 1
        % Simulate with D-8 flow direction
        % --- Adding NE, SE, SW, NW --- %
        M_NE(2:(ny),1:(nx-1)) = M(1:(ny-1),2:nx); % OK
        M_SE(1:(ny-1),1:(nx-1)) = M(2:ny,2:nx); % OK
        M_SW(1:(ny-1),2:(nx)) = M(2:(ny),1:(nx-1)); % OK
        M_NW(2:ny,2:nx) = M(1:(ny-1),1:(nx-1)); % OK
    end
end
