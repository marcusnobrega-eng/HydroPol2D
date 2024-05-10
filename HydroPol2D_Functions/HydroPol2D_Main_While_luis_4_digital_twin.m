%% Main HydroPol2D While Loop
% Developers: Marcus Nobrega & Luis Castillo
% Date 9/21/2023
% Goal - Run the main modeling process of the model
format short g

k = 1;
C = 0;
t = running_control.time_step_model;
Flooded_Area = 0;
Risk_Area = 0;
elevation = DEM_raster.Z;

if flags.flag_real_time_satellite_rainfall == 1
    current_time = datetime('now');
    [lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
    time_zone = timezone(mean(lon));
    t2 = running_control.time_step_model;
    rain_flag=0;
    end_flag=0;
    d_db = zeros(ny_max,nx_max,length(running_control.time_store));
    r_db = d_db; i_db = d_db; v_db = d_db;
    dates = strings(length(running_control.time_store),1);
    if flags.flag_unit_for_forecasting == 1
        flag_loader = 1;
        hour_to_save_previous = -1;
        save_counter = 7;
        gauges_data_all_forecast=table();
        app = HydroPol2D_Monitor_forecast();
        DCPs_all = [];
        gauges_all = table();
        temp_gauges = table();
    else
        save_counter = 0;
        app = HydroPol2D_Monitor();
    end
    ax_d = app.UIAxes; ax_r = app.UIAxes2; ax_date = app.DateTimeTextArea; ax_iter = app.iter;
    ax_i = app.UIAxes_3; ax_v = app.UIAxes_2; ax_system_output = app.SystemOutput;
    ax_list = app.gauges_list;
    % Calling the dashboard
    if flags.flag_unit_for_forecasting == 1
        [d_db,r_db,i_db,v_db,S_p,A,RA,dates,t_d,system_output,save_counter]=HydroPol2D_real_time_dashboard(ax_d,ax_r,ax_i,ax_v,ax_date,ax_iter,[],[],[],[],d_db,r_db,i_db,v_db,1,[],[],[],DEM_raster,register_data_2,dates,[],"Initializing the system...",ax_system_output,0,gauges,ax_list,flags.flag_unit_for_forecasting,save_counter);
    else
        [d_db,r_db,i_db,v_db,S_p,A,RA,dates,t_d,system_output,~]=HydroPol2D_real_time_dashboard(ax_d,ax_r,ax_i,ax_v,ax_date,ax_iter,Maps.Hydro.d(:,:,1),Maps.Hydro.spatial_rainfall_maps(:,:,1),Maps.Hydro.I_t(:,:,1),Maps.Hydro.I_t(:,:,1),d_db,r_db,i_db,v_db,1,[],[],[],DEM_raster,register_data_2,dates,[],"Initializing the system...",ax_system_output,0,gauges,ax_list,0,save_counter);
    end    
end


while t <= (running_control.routing_time + running_control.min_time_step/60)
    % Infiltration and Available Depth
    % Show stats
    if flags.flag_waterquality == 1
        perc____duremain______tsec_______dtmm______infmmhr____CmgL_____dtmWQ = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(C)), max(max((WQ_States.P_conc))), tmin_wq]
    else
        perc____duremain______tsec_______dtmm______infmmhr____ = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(C))]
    end
    if tmin_wq < 0 || isnan(tmin_wq) || isinf(tmin_wq)
        error('Instability')
    end
    % Inflows
    if flags.flag_inflow == 1
        BC_States.inflow = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2)); % This is to solve spatially, don't delete
        for i = 1:Inflow_Parameters.n_stream_gauges
            BC_States.inflow = BC_States.inflow + BC_States.delta_inflow_agg(i)*Wshed_Properties.inflow_cells(:,:,i);
        end
    end

    if k == 1
        if flags.flag_infiltration == 1
            Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow + depths.d_0 - Hydro_States.ETP/24)./(time_step/60);
            C = Soil_Properties.ksat.*(1 + ((depths.d_0 + Soil_Properties.psi).*(Soil_Properties.teta_sat - Soil_Properties.teta_i))./Soil_Properties.I_0); % matrix form
            Hydro_States.f = min(C,Hydro_States.i_a);
            Soil_Properties.I_t = max(Soil_Properties.I_0 + (Hydro_States.f)*time_step/60,min_soil_moisture); % Extracting ETP from soil
            Hydro_States.idx_ETR = Soil_Properties.I_t == min_soil_moisture;
            Hydro_States.ETR  = Hydro_States.ETP;
            Hydro_States.ETR (Hydro_States.idx_ETR) = Hydro_States.ETP(Hydro_States.idx_ETR) - (min_soil_moisture(Hydro_States.idx_ETR) - (Soil_Properties.I_0(Hydro_States.idx_ETR) + (Hydro_States.f(Hydro_States.idx_ETR))*time_step/60));
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.f*time_step/60+ idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_0 + depths.pef;
            Soil_Properties.T = Soil_Properties.Tr; % Beginning to track the replenishing time
        else
            %Hydro_States.i_a = (BC_States.delta_p_agg*rainfall_matrix + BC_States.delta_inflow_agg*Wshed_Properties.inflow_cells + d_0)./(time_step/60);
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow+ idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_0 + depths.pef;
        end
    else
        if flags.flag_infiltration == 1
            % Effective precipitation - Green-Ampt(1911)
            Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow + depths.d_p - Hydro_States.ETP/24)./(time_step/60);
            C = Soil_Properties.ksat.*(1 + ((depths.d_p + Soil_Properties.psi).*(Soil_Properties.teta_sat - Soil_Properties.teta_i))./Soil_Properties.I_p); % matrix form
            Non_C_idx = Soil_Properties.I_t > Soil_Properties.Lu*1000;
            C(Non_C_idx) = Soil_Properties.ksat(Non_C_idx);

            Hydro_States.idx_C = (Hydro_States.i_a <= C); % Values where Hydro_States.i_a is below C
            idx_i_a = (Hydro_States.i_a > C); % Values where Hydro_States.i_a is larger than C
            Soil_Properties.T = Soil_Properties.T - time_step/60/60; % Recoverying time (hours)
            Soil_Properties.T(idx_i_a) = Soil_Properties.Tr(idx_i_a); % Saturated Areas we begin the process again
            Soil_Properties.idx_T = Soil_Properties.T < 0; % Cells where the replenishing time is reached
            Hydro_States.f = min(C,Hydro_States.i_a); % Infiltration rate (mm/hr)
            Soil_Properties.I_t = max(Soil_Properties.I_p + Hydro_States.f*time_step/60 - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
            Hydro_States.idx_ETR = Soil_Properties.I_t == min_soil_moisture;
            Hydro_States.ETR  = Hydro_States.ETP;
            Hydro_States.ETR (Hydro_States.idx_ETR) = Hydro_States.ETP(Hydro_States.idx_ETR) - (min_soil_moisture(Hydro_States.idx_ETR) - (Soil_Properties.I_p(Hydro_States.idx_ETR) + Hydro_States.f(Hydro_States.idx_ETR)*time_step/60 - Soil_Properties.k_out(Hydro_States.idx_ETR)*time_step/60));
            % Refreshing Soil_Properties.I_t to I_0 for cases where idx_T > 0
            Soil_Properties.I_t(Soil_Properties.idx_T) = min_soil_moisture(Soil_Properties.idx_T);
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.f*time_step/60 + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_p + depths.pef; %% ATTENTION HERE
        else
            Hydro_States.i_a = (BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow + depths.d_p - Hydro_States.ETP/24)./(time_step/60);
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.ETP/24 + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_p + depths.pef;
        end
    end
    depths.d_tot = depths.d_t;
    if flags.flag_D8 == 1
        %%%% WATER DEPTHS %%%
        depths.d_left_cell = [zeros(ny_max,1),depths.d_tot(:,1:(nx_max-1))];
        depths.d_right_cell = [depths.d_tot(:,(2:(nx_max))) zeros(ny_max,1)];
        depths.d_up_cell = [zeros(1,nx_max) ; depths.d_tot(1:(ny_max-1),:)];
        depths.d_down_cell = [depths.d_tot(2:(ny_max),:) ; zeros(1,nx_max)];
        if flags.flag_D8 == 1
            % Simulate with D-8 flow direction
            % --- Adding NE, SE, SW, NW --- %
            depths.d_NE_t(2:(ny_max),1:(nx_max-1)) = depths.d_tot(1:(ny_max-1),2:nx_max); % OK
            depths.d_SE_t(1:(ny_max-1),1:(nx_max-1)) = depths.d_tot(2:ny_max,2:nx_max); % OK
            depths.d_SW_t(1:(ny_max-1),2:(nx_max)) = depths.d_tot(2:(ny_max),1:(nx_max-1)); % OK
            depths.d_NW_t(2:ny_max,2:nx_max) = depths.d_tot(1:(ny_max-1),1:(nx_max-1)); % OK
        else
            % Simulate with D-4 flow direction
            % Everything already calculated
        end
    end

    if flags.flag_D8 == 1 % D-8 C-A Model
%         [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,outflow_rate.qout_ne_t,outflow_rate.qout_se_t,outflow_rate.qout_sw_t,outflow_rate.qout_nw_t,depths.d_t,CA_States.I_tot_end_cell,~] = CA_Routing_8D(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.p,flags.flag_reservoir,Elevation_Properties.elevation_cell,depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical,wse_slope_zeros,Distance_Matrix);
        [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,outflow_rate.qout_ne_t,outflow_rate.qout_se_t,outflow_rate.qout_sw_t,outflow_rate.qout_nw_t,depths.d_t,CA_States.I_tot_end_cell,~] = CA_Routing_8D_not_optimized(Elevation_Properties.elevation_cell,Elevation_Properties.elevation_left_t,Elevation_Properties.elevation_right_t,Elevation_Properties.elevation_up_t,Elevation_Properties.elevation_down_t,depths.d_tot,depths.d_left_cell,depths.d_right_cell,depths.d_up_cell,depths.d_down_cell,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical,Elevation_Properties.elevation_NE_t, Elevation_Properties.elevation_SE_t, Elevation_Properties.elevation_SW_t, Elevation_Properties.elevation_NW_t, depths.d_NE_t, depths.d_SE_t, depths.d_SW_t, depths.d_NW_t,wse_slope_zeros,Distance_Matrix);
    else % 4-D CA Model
        % [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,CA_States.I_cell] = CA_Routing(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.p,flags.flag_reservoir,Elevation_Properties.elevation_cell,Elevation_Properties.elevation_left_t,Elevation_Properties.elevation_right_t,Elevation_Properties.elevation_up_t,Elevation_Properties.elevation_down_t,depths.d_tot,depths.d_left_cell,depths.d_right_cell,depths.d_up_cell,depths.d_down_cell,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical);
        [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,~] = CA_Routing(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.p,flags.flag_reservoir,Elevation_Properties.elevation_cell,depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical);
    end
    %     qout_t = qout_left + qout_right + qout_up + qout_down + qout_ne + qout_se + qout_sw + qout_ne;
    % Inflows - Von-Neuman
    outflow_rate.qin_left_t = [zeros(ny_max,1),outflow_rate.qout_right_t(:,1:(nx_max-1))];
    outflow_rate.qin_right_t = [outflow_rate.qout_left_t(:,(2:(nx_max))) zeros(ny_max,1)];
    outflow_rate.qin_up_t = [zeros(1,nx_max) ; outflow_rate.qout_down_t(1:(ny_max-1),:)];
    outflow_rate.qin_down_t = [outflow_rate.qout_up_t(2:(ny_max),:) ; zeros(1,nx_max)];

    % Inflows - Inclined Directions
    if flags.flag_D8 == 1
        if flags.flag_GPU == 1
            zero_matrix = gpuArray(zeros(size(Elevation_Properties.elevation_cell)));
        else
            zero_matrix = zeros(size(Elevation_Properties.elevation_cell));
        end
        outflow_rate.qin_ne_t = zero_matrix; outflow_rate.qin_se_t = zero_matrix; outflow_rate.qin_sw_t = zero_matrix; outflow_rate.qin_nw_t = zero_matrix;
        outflow_rate.qin_ne_t(2:(ny_max),1:(nx_max-1)) = outflow_rate.qout_sw_t(1:(ny_max-1),2:(nx_max)); % OK
        outflow_rate.qin_se_t(1:(ny_max-1),1:(nx_max-1)) = outflow_rate.qout_nw_t(2:ny_max,2:nx_max); % OK
        outflow_rate.qin_sw_t(1:(ny_max-1),2:(nx_max)) = outflow_rate.qout_ne_t(2:(ny_max),1:(nx_max-1)); % OK
        outflow_rate.qin_nw_t(2:ny_max,2:nx_max) = outflow_rate.qout_se_t(1:(ny_max-1),1:(nx_max-1)); % OK
    end
    if flags.flag_D8 == 1
        outflow_rate.qin_t = outflow_rate.qin_left_t + outflow_rate.qin_right_t + outflow_rate.qin_up_t + outflow_rate.qin_down_t + outflow_rate.qin_ne_t + outflow_rate.qin_se_t + outflow_rate.qin_sw_t + outflow_rate.qin_nw_t;
        %         outflow_rate.qin_t = outflow_rate.qin_left_t + outflow_rate.qin_right_t + outflow_rate.qin_up_t + outflow_rate.qin_down_t;
    else
        outflow_rate.qin_t = outflow_rate.qin_left_t + outflow_rate.qin_right_t + outflow_rate.qin_up_t + outflow_rate.qin_down_t;
    end
    idx3 = logical(isnan(outflow_rate.qin_t) + isinf(outflow_rate.qin_t));
    outflow_rate.qin_t(idx3) = 0;

    % Water Quality Parameters for f(B(t))
    if flags.flag_waterquality == 1
        if flags.flag_D8 ~= 1
            [WQ_States.B_t,WQ_States.P_conc,Out_Conc,tmin_wq,tot_W_out,WQ_States.mass_lost,WQ_States.Tot_Washed] = build_up_wash_off(LULC_Properties.C_3,LULC_Properties.C_4,outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,WQ_States.B_t,time_step,nx_max,ny_max,Wshed_Properties.cell_area,outlet_index,idx_nan_5,flags.flag_wq_model,WQ_States.mass_lost,WQ_States.Tot_Washed,LULC_Properties.Bmin,LULC_Properties.Bmax,LULC_Properties.min_Bt);
        else
            [WQ_States.B_t,WQ_States.P_conc,Out_Conc,tmin_wq,tot_W_out,WQ_States.mass_lost,WQ_States.Tot_Washed] = build_up_wash_off_8D(LULC_Properties.C_3,LULC_Properties.C_4,outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,outflow_rate.qout_ne_t,outflow_rate.qout_se_t,outflow_rate.qout_sw_t,outflow_rate.qout_nw_t,WQ_States.B_t,time_step,nx_max,ny_max,Wshed_Properties.cell_area,outlet_index,idx_nan_5,flags.flag_wq_model,WQ_States.mass_lost,WQ_States.Tot_Washed,LULC_Properties.Bmin,LULC_Properties.Bmax,LULC_Properties.min_Bt);
        end
    end
    %%%% Checking Mass Balance
    if flags.flag_waterquality == 1
        if sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t)))) > 1.2*initial_mass  % More than 5%
            error('Brutal instability in B(t). More than 20% difference')
        end
    end
    % New Time-step Calculation
    running_control.pos_save = find(running_control.time_change_records < t,1,'last');
    running_control.time_save = running_control.time_change_records(running_control.pos_save); % min
    running_control.delta_time_save = running_control.time_save - running_control.time_save_previous;
    running_control.time_save_previous = running_control.time_save;
    running_control.actual_record_timestep = find(running_control.time_change_records < t,1,'last');
    if running_control.delta_time_save > 0 || k == 1 % First time-step
        if flags.flag_timestep == 0
            %%% SOLUTION FOR COURANT METHOD %%%
            % depths.d_left_cell = [zeros(ny_max,1),depths.d_t(:,1:(nx_max-1))];
            % depths.d_right_cell = [depths.d_t(:,(2:(nx_max))) zeros(ny_max,1)];
            % depths.d_up_cell = [zeros(1,nx_max) ; depths.d_t(1:(ny_max-1),:)];%
            % depths.d_down_cell = [depths.d_t(2:(ny_max),:) ; zeros(1,nx_max)];
            % if flags.flag_D8 == 1
            %     depths.d_NE_cell = zero_matrix; depths.d_SE_cell = zero_matrix; depths.d_SW_cell = zero_matrix; depths.d_NW_cell = zero_matrix;
            %     depths.d_NE_cell(2:(ny_max),1:(nx_max-1)) = depths.d_t(1:(ny_max-1),2:nx_max); % OK
            %     depths.d_SE_cell(1:(ny_max-1),1:(nx_max-1)) = depths.d_t(2:ny_max,2:nx_max); % OK
            %     depths.d_SW_cell(1:(ny_max-1),2:(nx_max)) = depths.d_t(2:(ny_max),1:(nx_max-1)); % OK
            %     depths.d_NW_cell(2:ny_max,2:nx_max) = depths.d_t(1:(ny_max-1),1:(nx_max-1)); % OK
            % end
            % Considering Velocity as Celerity
            %             velocities.vel_left = (I_cell(:,:,1)./(0.5/1000.*(depths.d_t + d_left_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_left_cell)/2/1000); % m/s
            %             velocities.vel_right = (I_cell(:,:,2)./(0.5/1000.*(depths.d_t + d_right_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_right_cell)/2/1000);
            %             velocities.vel_up = (I_cell(:,:,3)./(0.5/1000.*(depths.d_t + d_up_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_up_cell)/2/1000); % m/s;;
            %             velocities.vel_down = (I_cell(:,:,4)./(0.5/1000.*(depths.d_t + d_down_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_down_cell)/2/1000); % m/s;

            z = depths.d_t;
            if flags.flag_GPU == 1
                Courant_Parameters.dt_threshold = gpuArray(1); % mm
            else
                Courant_Parameters.dt_threshold = (1); % mm
            end
            z(depths.d_t < Courant_Parameters.dt_threshold) = 1e12;
            velocities.vel_left = (outflow_rate.qout_left_t/1000/3600)*Wshed_Properties.Resolution^2./(Wshed_Properties.Resolution*z/1000); % m/s
            velocities.vel_right = (outflow_rate.qout_right_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s
            velocities.vel_up = (outflow_rate.qout_up_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s
            velocities.vel_down = (outflow_rate.qout_down_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s

            if flags.flag_D8 == 1
                %                 velocities.vel_ne = (CA_States.I_cell(:,:,6)./(0.5/1000.*(depths.d_t + d_NE_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_NE_cell)/2/1000); % m/s
                %                 velocities.vel_se = (CA_States.I_cell(:,:,7)./(0.5/1000.*(depths.d_t + d_SE_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_SE_cell)/2/1000);
                %                 velocities.vel_sw = (CA_States.I_cell(:,:,8)./(0.5/1000.*(depths.d_t + d_SW_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_SW_cell)/2/1000); % m/s;;
                %                 velocities.vel_nw = (CA_States.I_cell(:,:,9)./(0.5/1000.*(depths.d_t + d_NW_cell).*Resolution))/(time_step*60) + sqrt(9.81*(depths.d_t + d_NW_cell)/2/1000); % m/s;
                velocities.vel_ne = (outflow_rate.qout_ne_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s
                velocities.vel_se = (outflow_rate.qout_se_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s
                velocities.vel_sw = (outflow_rate.qout_sw_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s
                velocities.vel_nw = (outflow_rate.qout_nw_t/1000/3600)*Wshed_Properties.Resolution./(z/1000); % m/s
            end
            %%%%%%%%%%%%%% Find the Maximum Velocity
            velocities.max_velocity_left = max(max(velocities.vel_left));
            velocities.max_velocity_right = max(max(velocities.vel_right));
            velocities.max_velocity_up = max(max(velocities.vel_up));
            velocities.max_velocity_down = max(max(velocities.vel_down));
            if flags.flag_D8 == 1
                velocities.max_velocity_ne = max(max(velocities.vel_ne));
                velocities.max_velocity_se = max(max(velocities.vel_se));
                velocities.max_velocity_sw = max(max(velocities.vel_sw));
                velocities.max_velocity_nw = max(max(velocities.vel_nw));
            end
            % - Velocit Raster - %
            %%% Maximum of All of Them %%%
            velocities.velocity_raster = max(velocities.vel_left,velocities.vel_right);
            velocities.velocity_raster = max(velocities.velocity_raster,velocities.vel_up);
            velocities.velocity_raster = max(velocities.velocity_raster,velocities.vel_down);
            if flags.flag_D8 == 1
                velocities.velocity_raster = max(velocities.velocity_raster,velocities.vel_ne);
                velocities.velocity_raster = max(velocities.velocity_raster,velocities.vel_se);
                velocities.velocity_raster = max(velocities.velocity_raster,velocities.vel_sw);
                velocities.velocity_raster = max(velocities.velocity_raster,velocities.vel_nw);
                velocities.velocity_vector = [velocities.max_velocity_left, velocities.max_velocity_right, velocities.max_velocity_up, velocities.max_velocity_down, velocities.max_velocity_ne, velocities.max_velocity_se, velocities.max_velocity_sw, velocities.max_velocity_nw];
            else
                velocities.velocity_vector = [velocities.max_velocity_left, velocities.max_velocity_right, velocities.max_velocity_up, velocities.max_velocity_down];
            end

            velocities.max_velocity = max(velocities.velocity_vector);
            if flags.flag_D8 == 1
                Courant_Parameters.factor_grid = sqrt(1/2);
            else
                Courant_Parameters.factor_grid = 1;
            end
            if flags.flag_GPU == 1
                Courant_Parameters.factor_grid = gpuArray(Courant_Parameters.factor_grid);
            end
            if velocities.max_velocity > 0
                new_timestep = (Courant_Parameters.factor_grid*Wshed_Properties.Resolution/velocities.max_velocity); % seconds
                Courant_Parameters.time_step_factor = max(Courant_Parameters.alfa_max - Courant_Parameters.slope_alfa*(max(velocities.max_velocity - Courant_Parameters.v_threshold,0)),Courant_Parameters.alfa_min);
                Courant_Parameters.alfa_save(running_control.pos_save,1) = Courant_Parameters.time_step_factor;
                new_timestep = new_timestep*Courant_Parameters.alfa_save(running_control.pos_save,1);
                new_timestep = double(min(new_timestep,running_control.max_time_step));
            elseif velocities.max_velocity < 0
                error('Model instability. Velocities are becoming negative.')
            elseif isnan(velocities.max_velocity)
                new_timestep = running_control.max_time_step;
            end

            if velocities.max_velocity == 0
                new_timestep = running_control.max_time_step;
            end

            % ---- Calculation of Stability --- %
            if flags.flag_human_instability == 1
                Human_Instability.F_person = Human_Instability.weight_person*Human_Instability.gravity; % N
                % Buyoance
                Human_Instability.F_buoy = (Human_Instability.width1_person*Human_Instability.width2_person)*depths.d_t/1000*Human_Instability.ro_water*Human_Instability.gravity; % N
                % Available Friction
                Human_Instability.available_friction = Human_Instability.mu*(max(Human_Instability.F_person - Human_Instability.F_buoy,0));
                % Hydrodynamic Force
                Human_Instability.hydro_force = max(1/2*(Human_Instability.Cd*Human_Instability.width1_person*depths.d_t/1000.*velocities.velocity_raster.^2),0);
                % Risk Factor
                Human_Instability.risk_t = min(Human_Instability.hydro_force./Human_Instability.available_friction,1);
                Human_Instability.risk_t(idx_nan) = nan;
            end
        else             % % % % % % % %  Solution for the Stable Method - Time-step refreshment
            % Water Slopes Calculation (THERE IS A MISTAKE HERE)
            error('This method is currenly not working, please choose the courant method')
        end
        if flags.flag_waterquality == 1
            % Adding Water Quality Minimum Time-step
            Courant_Parameters.alfa_wq = 1; % Reduction factor of water quality
            new_timestep = min(Courant_Parameters.alfa_wq*tmin_wq,new_timestep);
        end
        % Rounding time-step to min_timestep or max_timestep with the required
        % precision. This is not very interesting, maybe we should delete
        % it
        running_control.time_calculation_routing(k,1) = new_timestep;
        running_control.time_calculation_routing(k,1) = max(running_control.time_step_increments*floor(running_control.time_calculation_routing(k,1)/running_control.time_step_increments),running_control.min_time_step);
        running_control.time_calculation_routing(k,1) = min(running_control.time_calculation_routing(k,1),running_control.max_time_step);
        time_step = running_control.time_calculation_routing(k,1)/60; % new time-step for the next run
        if running_control.time_calculation_routing(k,1) == running_control.min_time_step % Check if we reached the minimum time-step
            unstable = 1; % If this is 1, it means we possibly had an unstable period at least
        end
    elseif  k == 1
        running_control.time_calculation_routing(k,1) = time_step*60; % sec
    else
        running_control.time_calculation_routing(k,1) = running_control.time_calculation_routing(k-1,1);
    end
    velocities.max_velocity_previous = velocities.max_velocity; % Assigning right velocities
    % Agregating Inflows to the New Time-step
    if flags.flag_inflow > 0
        for z = 1:Inflow_Parameters.n_stream_gauges
            z1 = find(BC_States.time_deltainflow(z,:) > t_previous,1,'first'); % begin of the time-step
            z2 = find(BC_States.time_deltainflow(z,:) <= t,1,'last'); % end of the time-step
            if isempty(z1)
                z1 = 1;
            end
            if isempty(z2) || z2 < z1
                z2 = z1;
            end
            if time_step >= time_step_model
                BC_States.delta_inflow_agg(z,1) = mean(BC_States.delta_inflow(z,z1:z2))/(time_step_model*60)*time_step*60;
            else
                BC_States.delta_inflow_agg(z,1) = BC_States.delta_inflow(z,z1)/(time_step_model*60)*time_step*60;
            end
        end
    end
    % Agregating Precipitation to the New Time-step
    if flags.flag_rainfall > 0
        if flags.flag_spatial_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
            z1 = find(running_control.time_deltap > t_previous,1,'first'); % begin of the time-step
            z2 = find(running_control.time_deltap <= t,1,'last'); % end of the time-step
            if z2 < z1
                z2 = z1;
            end
            if time_step >= time_step_model
                BC_States.delta_p_agg = mean(BC_States.delta_p(1,z1:z2))/(time_step_model*60)*time_step*60;
            else
                BC_States.delta_p_agg = BC_States.delta_p(1,z1)/(time_step_model*60)*time_step*60;
            end
            if isnan(BC_States.delta_p_agg)
                ttt = 1;
            end
        elseif flags.flag_spatial_rainfall == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
            % Spatial Rainfall
            % Code for Spatial-Varying Rainfall
            % Times of full rainfall dataset
            z1 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg <= t_previous,1,'last');
            z2 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg <= t,1,'last');

            % Times for aggregated rainfall dataset
            zz1 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t_previous,1,'last');
            zz2 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t,1,'last');

            if Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2)  == Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) % No difference between them, same time-step
                Rainfall_Parameters.index_aggregation = 1; % We save only 1 map
            elseif zz1 == zz2 && z2 ~= z1 % We save more than 1 map
                Rainfall_Parameters.index_aggregation = Rainfall_Parameters.index_aggregation + (z2-z1);
            elseif zz2 > zz1
                Rainfall_Parameters.index_aggregation = 1;
            end

            if z1 ~= z2 || z2 == length(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg)
                Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1); % Coordinates (easting) of each rain gauge
                Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2); % Coordinates (northing) of each rain gauge
                Spatial_Rainfall_Parameters.x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,2)]'; % Pixel easting coordinates
                Spatial_Rainfall_Parameters.y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,1)]'; % Pixel northing coordinates

                % Spatial Rainfall
                if z2 == length(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg)
                    spatial_rainfall = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2));

                else
                    rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(z2,1:Spatial_Rainfall_Parameters.n_raingauges)'; % Values of rainfall at t for each rain gauge
                    % idx_rainfall = logical(isnan(rainfall) | rainfall == 0);
                    idx_rainfall = logical(isnan(rainfall));
                    rainfall(idx_rainfall) = []; % Taking out nans
                    Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = []; % Taking out nans
                    Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = []; % Taking out nans
                    if isempty(Spatial_Rainfall_Parameters.x_coordinate) % No rainfall
                        spatial_rainfall = zeros(size(Elevation_Properties.elevation_cell));
                        spatial_rainfall(idx_nan) = nan;
                    else
                        [spatial_rainfall] = Rainfall_Interpolator(Spatial_Rainfall_Parameters.x_coordinate,Spatial_Rainfall_Parameters.y_coordinate,rainfall,Spatial_Rainfall_Parameters.x_grid,Spatial_Rainfall_Parameters.y_grid); % Interpolated Values
                        spatial_rainfall(idx_nan) = nan;
                    end

                    if nansum(nansum(spatial_rainfall)) > 0
                        sprintf('Rainfall interpolated')
                    end
                end
                if zz2 > zz1 % Saving Maps
                    Maps.Hydro.spatial_rainfall_maps(:,:,zz2) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                    zzz = Maps.Hydro.spatial_rainfall_maps(:,:,zz2);
                    BC_States.average_spatial_rainfall(zz2,1) = mean(zzz(zzz>=0));
                end                
                BC_States.delta_p_agg = spatial_rainfall/3600*time_step*60; % Matrix of delta P for each pixel
                % Saving Maps for Spatial Aggregation
                rainfall_spatial_aggregation(:,:,Rainfall_Parameters.index_aggregation) = gather(spatial_rainfall); % Saving high resolution map
            end
        elseif flags.flag_input_rainfall_map == 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
            % Input Rainfall Maps
            % Times of full rainfall dataset
            z1_input = find(Input_Rainfall.time <= t_previous,1,'last');
            z2_input = find(Input_Rainfall.time <= t,1,'last');
            if z2_input > z1_input
                % Read Geotiff
                input_rainfall = GRIDobj(Input_Rainfall.labels_Directory{z2_input}{1});
                if flags.flag_resample == 1
                    raster_resample = DEM_raster;
                    if input_rainfall.cellsize ~= GIS_data.resolution_resample
                        if sum(input_rainfall.refmat(:) == DEM_raster.refmat(:)) ~= 6 % we have 6 information in refmat
                            input_rainfall = resample(input_rainfall,DEM_raster,'nearest');
                            MASK = DEM_raster; MASK.Z = ~isnan(MASK.Z); 
                            input_rainfall = clip(input_rainfall,MASK);
                        else
                            % Resample other two rastersMaps.Hydro
                            input_rainfall = resample(input_rainfall,DEM_raster,'nearest');
                        end   
                    end
                end
                if flags.flag_single == 1
                    input_rainfall = single(input_rainfall.Z); % Only the values
                    if flags.flag_GPU == 1
                        input_rainfall = gpuArray(input_rainfall); % Only the values
                    end
                else               
                    input_rainfall = (input_rainfall.Z); % Only the values
                    if flags.flag_GPU == 1
                        input_rainfall = gpuArray(input_rainfall); % Only the values
                    end                         
                end
                input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;flags.flag_real_time_satellite_rainfall == 1
                Maps.Hydro.spatial_rainfall_maps(:,:,z2_input) = input_rainfall;
                BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
                BC_States.average_spatial_rainfall(z2_input,1) = mean(input_rainfall(input_rainfall>=0));
            end
        elseif flags.flag_input_rainfall_map == 0 && flags.flag_satellite_rainfall == 1 && flags.flag_real_time_satellite_rainfall ~= 1
            % Satellite Rainfall
            if flags.flag_satellite_rainfall == 1
                % Input Rainfall Maps
                % Times of full rainfall dataset
                z1_input = find(Input_Rainfall.time <= t_previous,1,'last');
                z2_input = find(Input_Rainfall.time <= t,1,'last');
                if z2_input > z1_input
                    product = 'PDIRNow1hourly';
                    [rainfall_raster, register,register_data,~] = Satellite_rainfall_processing([],[],register,product,date_begin,date_end,flags.flag_satellite_rainfall,flags.flag_real_time_satellite_rainfall,DEM_raster);
                    input_rainfall = rainfall_raster.Z;
                    input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;
                    Maps.Hydro.spatial_rainfall_maps(:,:,z2_input) = input_rainfall;
                    BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
                    BC_States.average_spatial_rainfall(z2_input,1) = mean(input_rainfall(input_rainfall>=0)); 
                end
            end
        elseif flags.flag_input_rainfall_map == 0 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall == 1
            if flags.flag_real_time_satellite_rainfall == 1
                % Input Rainfall Maps
                % Times of full rainfall dataset
                z1_input = find(Input_Rainfall.time <= t_previous,1,'last');
                z2_input = find(Input_Rainfall.time <= t,1,'last');
                if z2_input ~= z1_input
                    if flags.flag_unit_for_forecasting == 1
                        if flag_loader == 1
                            %Check the last state update from the main system
                            trier = 1;
                            while trier ==1
                                load('current_run_time','hour_to_save'); 
                                if (hour_to_save == 0 || hour_to_save == 6 || hour_to_save == 12 || hour_to_save == 18) && hour_to_save_previous ~= hour_to_save    
                                    %Load the last state of the system from the near-real time simulation
                                    app_temp = app; t_d_temp=t_d;
                                    load('Last_system_update');
                                    app = app_temp; t_d = t_d_temp;
                                    % t_d.Period = 0.3;
                                    % Cumulating the DCPs data 
                                    DCPs_all = [DCPs_all;DCP_data];
                                    for i = 1:gauges.num_obs_gauges
                                        if flags.flag_GPU
                                            gauges_all = horzcat(gauges_all, table(squeeze(d_db(gauges.northing_obs_gauges(i),gauges.easting_obs_gauges(i),:)),'VariableNames',{extra_parameters.gauges.labels_observed_string{i}{1}}));
                                        else
                                            gauges_all = horzcat(gauges_all, table(squeeze(d_db(gauges.northing_obs_gauges(i),gauges.easting_obs_gauges(i),:)),'VariableNames',{gauges.labels_observed_string{i}{1}}));
                                        end
                                    end
                                    % Adding the time for the observed
                                    % gauged data
                                    gauges_all = horzcat(gauges_all, table(datetime(dates(:)),'VariableNames',{'Time'}));
                                    % updating the t_d values
                                    ax_d = app.UIAxes; ax_r = app.UIAxes2; ax_date = app.DateTimeTextArea; ax_iter = app.iter;
                                    ax_i = app.UIAxes_3; ax_v = app.UIAxes_2; ax_system_output = app.SystemOutput;
                                    ax_list = app.gauges_list;
                                    app.gauges_list.Items = char_vector_cell;
                                    gauges_data_all_forecast = [gauges_data_all_forecast;gauges_data_all];
                                    % Recovering the last observation from
                                    % the gauges
                                    app.time_zone = time_zone;
                                    app.gauges_data = gauges_data_all(gauges_data_all.Time >= gauges_data_all.Time(end)-hours(24),char_vector_cell);
                                    app.gauges_time = gauges_data_all(gauges_data_all.Time >= gauges_data_all.Time(end)-hours(24),'Time');
                                    register = 0;
                                    first_register = dates(end);
                                    flag_loader = 0;
                                    flags.flag_unit_for_forecasting = 1;
                                    trier = 0;
                                    d_db = cat(3,zeros(ny_max, nx_max, 17),d_db); r_db = cat(3,zeros(ny_max, nx_max, 17),r_db);
                                    i_db = cat(3,zeros(ny_max, nx_max, 17),i_db); v_db = cat(3,zeros(ny_max, nx_max, 17),v_db);
                                    dates = [strings(17,1);dates];
                                    %updates databases in the dashboard
                                    t_d.TimerFcn{29} = "Initializing the system for forecasting...";
                                    t_d.TimerFcn{11} = i_db;
                                    t_d.TimerFcn{15} = v_db;
                                    t_d.TimerFcn{3} = d_db;
                                    t_d.TimerFcn{7} = r_db;
                                    t_d.TimerFcn{25} = dates;
                                    t_d.TimerFcn{31} = hour_to_save;
                                    t_d.TimerFcn{34} = dates(end);

                                    if hour_to_save == 0
                                        % to choose the first cicle
                                        cycler = 1;
                                    else
                                        % to choose one of thre three left
                                        % cicles (2, 3 or 4)
                                        cycler = hour_to_save/6 +1; 
                                    end
                                    hour_to_save_previous = hour_to_save;
                                else
                                    % Waiting for the main system update
                                    pause(300)
                                end
                            end
                        end
                        [rainfall_raster, register, register_data] = NWP_rainfall_processing(t_d,ax_system_output,register,cycler,DEM_raster,dates,first_register);
                        register_data_2=register_data;
                        % here is check if a new main system is available
                        % if positive, the system ignore the current
                        % forecast and jumps into the newest
                        load('current_run_time','hour_to_save');
                        if (hour_to_save == 0 || hour_to_save == 6 || hour_to_save == 12 || hour_to_save == 18) && hour_to_save_previous ~= hour_to_save
                            % Here, even if the system didnt finish the
                            % forecast, it will start running the most recent
                            % update
                            flag_loader = 1;
                            % saving gauges observation befored exported
                            for i = 1:gauges.num_obs_gauges
                                if flags.flag_GPU
                                    temp_gauges = horzcat(temp_gauges, table(squeeze(d_db(gauges.northing_obs_gauges(i),gauges.easting_obs_gauges(i),:)),'VariableNames',{extra_parameters.gauges.labels_observed_string{i}{1}}));
                                else
                                    temp_gauges = horzcat(temp_gauges, table(squeeze(d_db(gauges.northing_obs_gauges(i),gauges.easting_obs_gauges(i),:)),'VariableNames',{gauges.labels_observed_string{i}{1}}));
                                end
                            end
                            % Adding the time for the observed
                            temp_gauges = horzcat(temp_gauges, table(squeeze(datetime(dates(:))),'VariableNames',{'Time'}));
                            % gauged data
                            gauges_all = [gauges_all; temp_gauges];
                            gauges_all = table();
                            temp_gauges= table();
                            % saving the data from the forecasted data
                            writetable(gauges_data_all_forecast,'Forecast_gauges_summary.txt','Delimiter',',')
                            writetable(DCPs_all,'DCPs_summary.txt','Delimiter',',')
                            save(strcat('data_',num2str(hour_to_save_previous')),'i_db','v_db','d_db','r_db','dates');
                        end
                        if register == 121
                            % the images for that cycle are ended, new
                            % cycle will be load
                            flag_loader = 1;
                            gauges_all = table();
                            temp_gauges= table();
                            % saving the data from the forecasted data
                            writetable(gauges_data_all_forecast,'Forecast_gauges_summary.txt','Delimiter',',')
                            writetable(DCPs_all,'DCPs_summary.txt','Delimiter',',')
                            save(strcat('data_',num2str(hour_to_save_previous')),'i_db','v_db','d_db','r_db','dates');
                        end
                    else
                        product = 'PDIRNow1hourly';
                        [rainfall_raster, register,register_data,register_data_2] = Satellite_rainfall_processing(t_d,ax_system_output,register,product,date_begin,date_end,flags.flag_satellite_rainfall,flags.flag_real_time_satellite_rainfall,DEM_raster);
                        temp_data = register_data{1};
                        hour_to_save = str2double(temp_data(16:17));
                        % This save the current workspace according to 
                        % certain time of the day (4 saves per day) because
                        % due to the NWP update frequency
                        if hour_to_save == 0 || hour_to_save == 6 || hour_to_save == 12 || hour_to_save == 18
                            % if dates(end) > datetime("15-Nov-2023 00:00:00") 
                            %     pause(7200)
                            % end
                            save('current_run_time','hour_to_save');
                            vars = whos;
                            exclude = {'app','t_d','ax_d','ax_date','ax_i','ax_iter','ax_list','ax_r','ax_system_output','ax_v'}; 
                            saveVars = {};
                            for i = 1:length(vars)
                                if ~ismember(vars(i).name, exclude)
                                    saveVars{end+1} = vars(i).name;
                                end
                            end
                            save('Last_system_update.mat',saveVars{:},'-v7.3');
                        end
                    end
                    input_rainfall = rainfall_raster.Z;
                    input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;
                    if end_flag ==1
                        Maps.Hydro.spatial_rainfall_maps(:,:,end) = input_rainfall;
                    else
                        Maps.Hydro.spatial_rainfall_maps(:,:,z2_input) = input_rainfall;
                    end
                    BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
                    BC_States.average_spatial_rainfall(z2_input,1) = mean(input_rainfall(input_rainfall>=0));
                    rain_flag=1;
                end
            end 
        end
    end

    % Aggregating ETP for next time-step

    if flags.flag_ETP == 1
        z1 = find(ETP_Parameters.climatologic_spatial_duration <= t_previous,1,'last');
        z2 = find(ETP_Parameters.climatologic_spatial_duration <= t,1,'last');
        if isempty(z1) && isempty(z2)
            % No data, No ETP.
            Hydro_States.ETP = zeros(size(Elevation_Properties.elevation_cell));
        else
            if ~isempty(z1) && z1 == z2 && z2 == length(ETP_Parameters.climatologic_spatial_duration)
                % Outside of ETP data maximum duration
                Hydro_States.ETP = zeros(size(Elevation_Properties.elevation_cell));
                Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
                ETR_save(:,:,z2) = Hydro_States.ETR ;
            elseif  (isempty(z1) && z2 > 0) || z2 > z1 && z2 < length(ETP_Parameters.climatologic_spatial_duration)
                if flags.flag_GPU == 1 || flags.flag_single == 1
                    day_of_year = day(extra_parameters.ETP.time_ETP(z2,1),'dayofyear');
                else
                    day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
                end
                [Hydro_States.ETP] = ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2),Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations,ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations,ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);
                if nansum(nansum(Hydro_States.ETP)) == 0
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,z2-1);
                    Hydro_States.ETR  = ETR_save(:,:,z2-1);
                    if nansum(nansum(Hydro_States.ETP)) == 0
                        warning('No ETP and ETR data. Assuming it equals 0')
                        Hydro_States.ETP = zeros(size(DEM));
                        Hydro_States.ETP(idx_nan) = nan;
                        Hydro_States.ETR  = zeros(size(DEM));
                        Hydro_States.ETR (idx_nan) = nan;
                    end
                end
                Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
                ETR_save(:,:,z2) = Hydro_States.ETR ;
            elseif z1 == 1 && z2 == 1 && k == 1
                % First data. We assume a constant ETP using
                % 1st data
                if flags.flag_GPU == 1 || flags.flag_single == 1
                    day_of_year = day(extra_parameters.ETP.time_ETP(z2,1),'dayofyear');
                else
                    day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
                end
                [Hydro_States.ETP] = ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2),Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations,ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations,ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);
                if nansum(nansum(Hydro_States.ETP)) == 0
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,z2-1);
                    Hydro_States.ETR  = ETR_save(:,:,z2-1);
                end
                Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
                ETR_save(:,:,z2) = Hydro_States.ETR ;
            end
        end
    end

    % Previous Time-step
    if k == 1
        t_previous = running_control.time_calculation_routing(k,1)/60;
        %             t_previous_date = date_begin + running_control.time_calculation_routing(k,1)/60/60/24; % Datetime
    else
        t_previous = t;
        %             t_previous_date = t_previous_date + t/60/24; % Datetime
    end

    % Inflows and Depth Refreshments
    depths.d_t = depths.d_t + outflow_rate.qin_t*time_step/60;

    % Clearing stored values
    if flags.flag_GPU == 1
        outflow_rate.qout_left_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_right_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_up_t =gpuArray( zeros(ny_max,nx_max));
        outflow_rate.qout_down_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_ne_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_se_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_sw_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_nw_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qin_t = gpuArray(zeros(ny_max,nx_max));
    else
        outflow_rate.qout_left_t = zeros(ny_max,nx_max);
        outflow_rate.qout_right_t = zeros(ny_max,nx_max);
        outflow_rate.qout_up_t = zeros(ny_max,nx_max);
        outflow_rate.qout_down_t = zeros(ny_max,nx_max);
        outflow_rate.qout_ne_t = zeros(ny_max,nx_max);
        outflow_rate.qout_se_t = zeros(ny_max,nx_max);
        outflow_rate.qout_sw_t = zeros(ny_max,nx_max);
        outflow_rate.qout_nw_t = zeros(ny_max,nx_max);
        outflow_rate.qin_t = zeros(ny_max,nx_max);
    end

    % Saving Plotting Values - Recording Time
    % Maps of Flood Depths, WSE and Pollutant Concentrations
    % --- Calculating EMC --- %
    if  flags.flag_automatic_calibration ~= 1
        if flags.flag_waterquality == 1
            WQ_States.mass_outlet = max(WQ_States.mass_outlet + Out_Conc*((nansum(nansum(outlet_states.outlet_flow)/1000/3600*1000)))*(time_step*60),0); % mg
            WQ_States.vol_outlet = max((nansum(nansum(outlet_states.outlet_flow))/1000/3600*1000)*(time_step*60) + WQ_States.vol_outlet,0);
        end
    end
    % Current time
    t_save = t + running_control.time_calculation_routing(k,1)/60;

    if flags.flag_real_time_satellite_rainfall == 1
        t_save_2 = t2 + running_control.time_calculation_routing(k,1)/60;
    end

    % Maps with generally coarser resolution
    recording_parameters.actual_record_state = find(running_control.time_records < t_save,1,'last');
    recording_parameters.delta_record = recording_parameters.actual_record_state - recording_parameters.last_record_maps;
    recording_parameters.last_record_maps = recording_parameters.actual_record_state;


    if flags.flag_automatic_calibration ~= 1 && flags.flag_real_time_satellite_rainfall == 1
        if k == 1 % First time-step
            Maps.Hydro.d(:,:,1) = depths.d_t;
            if flags.flag_human_instability == 1
                Maps.Hydro.risk(:,:,1) = Human_Instability.risk_t;
            end
            Maps.Hydro.I_t(:,:,1) = Soil_Properties.I_t;
        elseif recording_parameters.delta_record > 0 % Saving maps
            if t_save < 360
                t_store = recording_parameters.actual_record_state;
                Maps.Hydro.d(:,:,t_store) = depths.d_t;
                if flags.flag_human_instability == 1
                    Maps.Hydro.risk(:,:,t_store) = Human_Instability.risk_t;
                end
                Maps.Hydro.I_t(:,:,t_store) = Soil_Properties.I_t;  

                % sending info for the dashboard
                [d_db,r_db,i_db,v_db,S_p,A,RA,dates,t_d,system_output,save_counter]=HydroPol2D_real_time_dashboard(ax_d,ax_r,ax_i,ax_v,ax_date,ax_iter,Maps.Hydro.d(:,:,t_store),[],Maps.Hydro.I_t(:,:,t_store),velocities.velocity_raster,d_db,r_db,i_db,v_db,k,S_p,A,RA,DEM_raster,[],dates,t_d,[],ax_system_output,0,[],[],0,save_counter);
              
            else
                % save the last record, then register it on the host.
                t_store = recording_parameters.actual_record_state;
                Maps.Hydro.d(:,:,t_store) = depths.d_t;
                if flags.flag_human_instability == 1
                    Maps.Hydro.risk(:,:,t_store) = Human_Instability.risk_t;
                end
                Maps.Hydro.I_t(:,:,t_store) = Soil_Properties.I_t;  
                % wiping register
                t_save = t_save - 420;
                recording_parameters.actual_record_state = 0;
                recording_parameters.delta_record = 0;
                recording_parameters.last_record_maps = 0;
            end
        elseif rain_flag==1
               [d_db,r_db,i_db,v_db,S_p,A,RA,dates,t_d,system_output,save_counter]=HydroPol2D_real_time_dashboard(ax_d,ax_r,ax_i,ax_v,ax_date,ax_iter,[],Maps.Hydro.spatial_rainfall_maps(:,:,end),[],[],d_db,r_db,i_db,v_db,k,S_p,A,RA,DEM_raster,register_data_2,dates,t_d,[],ax_system_output,0,[],[],0,save_counter);
               rain_flag=0;
        end % Calls the sub
    elseif flags.flag_automatic_calibration ~= 1
        if k == 1 % First time-step
            Maps.Hydro.d(:,:,1) = depths.d_t;
            if flags.flag_human_instability == 1
                Maps.Hydro.risk(:,:,1) = Human_Instability.risk_t;
            end
            Maps.Hydro.I_t(:,:,1) = Soil_Properties.I_t;
        elseif recording_parameters.delta_record > 0 % Saving maps
            t_store = recording_parameters.actual_record_state;
            Maps.Hydro.d(:,:,t_store) = depths.d_t;
            if flags.flag_human_instability == 1
                Maps.Hydro.risk(:,:,t_store) = Human_Instability.risk_t;
            end
            Maps.Hydro.I_t(:,:,t_store) = Soil_Properties.I_t;
        end % Calls the sub
    end


    % Hydrographs and Pollutographs with same time resolution
    recording_parameters.actual_record_hydrograph = find(running_control.time_record_hydrograph <= t_save,1,'last');
    recording_parameters.delta_record_hydrograph = recording_parameters.actual_record_hydrograph - recording_parameters.last_record_hydrograph;
    recording_parameters.last_record_hydrograph = recording_parameters.actual_record_hydrograph;

    if flags.flag_automatic_calibration ~= 1
        if k == 1
            outlet_states.outlet_hydrograph(1,1) = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600; % m3/s
            running_control.time_hydrograph(1,1) = running_control.time_calculation_routing(k,1)/60;
            outlet_states.depth_outlet(1,1) = mean(depths.d_t(idx_outlet));
            % Maximum Flodded Areas
            Flooded_Area = max(sum(sum((depths.d_t > 150)))*Wshed_Properties.Resolution^2,Flooded_Area); % Areas larger than 15 cm
            % Maximum Risk Areas
            if flags.flag_human_instability == 1
                % Maximum Risk Areas
                Risk_Area = max(sum(sum((Human_Instability.risk_t == 1)))*Wshed_Properties.Resolution^2,Risk_Area);
            end
            if flags.flag_waterquality(1,1) == 1
                Maps.WQ_States.Pol_Conc_Map(:,:,1) = WQ_States.P_conc;
                Maps.WQ_States.Pol_mass_map(:,:,1) = WQ_States.B_t;
                WQ_States.EMC_outlet(1,1) = WQ_States.mass_outlet/WQ_States.vol_outlet; % mg/L
                WQ_States.mass_outlet_save(1,1) = WQ_States.mass_outlet;
                WQ_States.vol_outlet_save(1,1) = WQ_States.vol_outlet;
            end
            % Saving Data of Input Gauges
            if flags.flag_obs_gauges == 1
                for i = 1:gauges.num_obs_gauges
                    gauges.x_cell = gauges.easting_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                    gauges.wse_cell(1,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000 + Elevation_Properties.elevation_cell(gauges.y_cell,gauges.x_cell); % m
                    gauges.depth_cell(1,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000; % m
                    gauges.hydrograph_cell(1,i) = CA_States.I_tot_end_cell(gauges.y_cell,gauges.x_cell)/(running_control.time_calculation_routing(k)); % m3/s
                end
                if flags.flag_real_time_satellite_rainfall
                    % Initialize gaugues_data using app.gauges_list
                    % Extracting the names from gauges or observed points
                    char_vector_cell = {};
                    if flags.flag_GPU
                        char_num = extra_parameters.gauges.labels_observed_string; 
                    else
                        char_num = gauges.labels_observed_string;
                    end
                    for i = 1:numel(char_num)
                        if flags.flag_GPU
                            char_vector_cell = [char_vector_cell; char(extra_parameters.gauges.labels_observed_string{i})];
                            app.gauges_time = array2table(gather(double(running_control.time_hydrograph(end)/60/24)) + date_begin - hours(time_zone), 'VariableNames', {'Time'});
                        else
                            char_vector_cell = [char_vector_cell; char(gauges.labels_observed_string{i})];
                            app.gauges_time = array2table(double(running_control.time_hydrograph(end)/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'});
                        end
                    end
                    
                    if flags.flag_river_heigth_compensate
                        app.gauges_data = array2table(gauges.depth_cell.*(Wshed_Properties.Resolution./gauges.witdh'), 'VariableNames', char_vector_cell);
                    else
                        app.gauges_data = array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell);
                    end
                    app.gauges_list.Items = char_vector_cell;
                    % Putting the code to display the data as simulated
                    
                    % Making a general database for all data gathered
                    if flags.flag_GPU
                        if flags.flag_river_heigth_compensate
                            gauges_data_all = array2table(gather(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh')), 'VariableNames', char_vector_cell);
                        else
                            gauges_data_all = array2table(gather(gauges.depth_cell(end,:)), 'VariableNames', char_vector_cell);
                        end
                        gauges_data_all = [gauges_data_all, array2table(1, 'VariableNames', {'code_simu'})];
                        gauges_data_all = [gauges_data_all, array2table(gather(double(running_control.time_hydrograph(end)/60/24)) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                    else
                        if flags.flag_river_heigth_compensate
                            gauges_data_all = array2table(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh'), 'VariableNames', char_vector_cell);
                        else
                            gauges_data_all = array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell);
                        end
                        gauges_data_all = [gauges_data_all, array2table(1, 'VariableNames', {'code_simu'})];
                        gauges_data_all = [gauges_data_all, array2table(double(running_control.time_hydrograph(end)/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                    end
                    % Donwloading data for DCPs
                    if flags.flag_data_source == 1
                        DCP_data = DCP_gather_processing([],[]);
                    elseif flags.flag_data_source == 2
                        if flags.flag_GPU == 1
                            DCP_data = DCP_gather_processing(extra_parameters.gauges.labels_observed_string, date_begin);
                        else 
                            DCP_data = DCP_gather_processing(gauges.labels_observed_string, date_begin);
                        end
                    end
                    
                    DCP_data.Time(:) = DCP_data.Time(:) - hours(time_zone);
                    app.HADS_gauges = DCP_data;
                    app.HADS_geoinfo = gauges;
                    if flags.flag_GPU
                        app.HADS_geoinfo.labels_observed_string = extra_parameters.gauges.labels_observed_string;
                    end
                    % Save all gauges records into a single table

                    % Plot the shapefile with circular symbols
                    mapshow(ax_d, mappoint(gauges.x_coord_gauges', gauges.y_coord_gauges'), ...
                        'DisplayType', 'point', ...
                        'Marker', 'o', ...
                        'MarkerEdgeColor', 'k', ...
                        'MarkerSize', 10,'LineWidth', 1.1);
                    % if some simbols as alert will be implemented, this
                    % will be the right place, similar as done above
                end
            end
            if flags.flag_waterquality == 1
                WQ_States.WQ_States.outet_pollutograph(1,1) = Out_Conc; % Already averaged for all outlet cells
            end
        elseif recording_parameters.delta_record_hydrograph > 0
            t_store = recording_parameters.actual_record_hydrograph;
            outlet_states.outlet_hydrograph(t_store,1) = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600; % m3/s
            running_control.time_hydrograph(t_store,1) = t;
            outlet_states.depth_outlet(t_store,1) = mean(depths.d_t(idx_outlet));
            % Maximum Flodded Areas
            Flooded_Area = max(sum(sum((depths.d_t > 150)))*Wshed_Properties.Resolution^2,Flooded_Area); % Areas larger than 15 cm
            if flags.flag_human_instability == 1
                % Maximum Risk Areas
                Risk_Area = max(sum(sum((Human_Instability.risk_t == 1)))*Wshed_Properties.Resolution^2,Risk_Area);
            end
            if flags.flag_waterquality == 1
                Maps.WQ_States.Pol_Conc_Map(:,:,t_store) = WQ_States.P_conc;
                Maps.WQ_States.Pol_mass_map(:,:,t_store) = WQ_States.B_t;
                WQ_States.EMC_outlet(t_store,1) = WQ_States.mass_outlet/WQ_States.vol_outlet; % mg/L
                WQ_States.mass_outlet_save(t_store,1) = WQ_States.mass_outlet;
                WQ_States.vol_outlet_save(t_store,1) = WQ_States.vol_outlet;
            end
            if flags.flag_obs_gauges == 1
                % Saving Data of Input Gauges
                for i = 1:gauges.num_obs_gauges
                    gauges.x_cell = gauges.easting_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                    gauges.wse_cell(t_store,i)= depths.d_t(gauges.y_cell,gauges.x_cell)/1000 + Elevation_Properties.elevation_cell(gauges.y_cell,gauges.x_cell); % m
                    gauges.depth_cell(t_store,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000; % m
                    gauges.hydrograph_cell(t_store,i) = CA_States.I_tot_end_cell(gauges.y_cell,gauges.x_cell)/(running_control.time_calculation_routing(k)); % m3/s
                end
                if flags.flag_real_time_satellite_rainfall
                    for i = 1:gauges.num_obs_gauges
                        gauges.x_cell = gauges.easting_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                        gauges.depth_cell(end,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000; % m
                    end
                    % if t_save_2 > running_control.record_time_maps
                    % Here to gather info from de digital twin or by the
                    % forecast system
                    if flags.flag_unit_for_forecasting && hour_to_save_previous ~= -1
                        if table2array(gauges_data_all_forecast(end,'Time') - gauges_data_all_forecast(1,'Time')) > hours(168)
                            gauges_data_all_forecast(1,:)=[];
                        else
                            if flags.flag_GPU
                                if flags.flag_river_heigth_compensate
                                    temp_table = [array2table(gather(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh')), 'VariableNames', char_vector_cell),array2table(double(gather(t2)/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                else
                                    temp_table = [array2table(gather(gauges.depth_cell(end,:)), 'VariableNames', char_vector_cell),array2table(double(gather(t2)/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                end
                            else
                                if flags.flag_river_heigth_compensate
                                    temp_table = [array2table(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh'), 'VariableNames', char_vector_cell),array2table(double(t2/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                else
                                    temp_table = [array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell),array2table(double(t2/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                end
                            end
                            
                            if hour_to_save == 0 && time_zone > 0
                                temp_table = [temp_table, array2table(24-time_zone, 'VariableNames', {'code_simu'})];
                            elseif hour_to_save == 0 && time_zone < 0
                                temp_table = [temp_table, array2table(0-time_zone, 'VariableNames', {'code_simu'})];
                            else
                                temp_table = [temp_table, array2table(hour_to_save-time_zone, 'VariableNames', {'code_simu'})];
                            end
                            gauges_data_all_forecast = [gauges_data_all_forecast;temp_table];

                            % Updating data for forecast plots
                            app.gauges_forecast= gauges_data_all_forecast;
                        end

                    else
                        if table2array(app.gauges_time(end,'Time') - app.gauges_time(1,'Time')) > hours(168)
                            app.gauges_data(1,:) = []; app.gauges_time(1,:) = []; % Dropping the first line of tables when they reach 12 hours of data
                            % app.gauges_data = [app.gauges_data; array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell)];
                            % app.gauges_time = [app.gauges_time; array2table(double(t2/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                            if table2array(gauges_data_all(end,'Time')-gauges_data_all(1,'Time')) > hours(720)
                                gauges_data_all(1,:) = [];
                                % temp_table = [array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell),array2table(1, 'VariableNames', {'code_simu'})];
                                % temp_table = [temp_table,array2table(double(t2/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                % gauges_data_all = [gauges_data_all;temp_table];
                            end
                        else
                            
                            if flags.flag_GPU
                                if flags.flag_river_heigth_compensate
                                    app.gauges_data = [app.gauges_data; array2table(gather(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh')), 'VariableNames', char_vector_cell)];
                                    temp_table = [array2table(gather(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh')), 'VariableNames', char_vector_cell),array2table(1, 'VariableNames', {'code_simu'})];
                                else
                                    app.gauges_data = [app.gauges_data; array2table(gather(gauges.depth_cell(end,:)), 'VariableNames', char_vector_cell)];
                                    temp_table = [array2table(gather(gauges.depth_cell(end,:)), 'VariableNames', char_vector_cell),array2table(1, 'VariableNames', {'code_simu'})];
                                end
                                temp_table = [temp_table,array2table(gather(double(t2/60/24)) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                app.gauges_time = [app.gauges_time; array2table(gather(double(t2/60/24)) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                            else
                                if flags.flag_river_heigth_compensate
                                    app.gauges_data = [app.gauges_data; array2table(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh'), 'VariableNames', char_vector_cell)];
                                    temp_table = [array2table(gauges.depth_cell(end,:).*(Wshed_Properties.Resolution./gauges.witdh'), 'VariableNames', char_vector_cell),array2table(1, 'VariableNames', {'code_simu'})];
                                else
                                    app.gauges_data = [app.gauges_data; array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell)];
                                    temp_table = [array2table(gauges.depth_cell(end,:), 'VariableNames', char_vector_cell),array2table(1, 'VariableNames', {'code_simu'})];
                                end
                                temp_table = [temp_table,array2table(double(t2/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                                app.gauges_time = [app.gauges_time; array2table(double(t2/60/24) + date_begin - hours(time_zone), 'VariableNames', {'Time'})];
                            end
                            gauges_data_all = [gauges_data_all;temp_table];
                        end
                    end
                    % Dowloading DCPs data
                    if datetime('now') - current_time> hours(1) % check if an hour has elapsed
                        if flags.flag_data_source == 1
                            DCP_data_temp = DCP_gather_processing([],[]);
                        elseif flags.flag_data_source == 2
                            if flags.flag_GPU == 1
                                DCP_data_temp = DCP_gather_processing(extra_parameters.gauges.labels_observed_string, date_begin);
                            else
                                DCP_data_temp = DCP_gather_processing(gauges.labels_observed_string, date_begin);
                            end
                        end
                        DCP_data_temp.Time(:) = DCP_data_temp.Time(:) - hours(time_zone);

                        if DCP_data.Time(end) ~= DCP_data_temp.Time(end)
                             % Table2 has more rows, concatenate table2 to the end of table1
                             % if ~isempty(setxor(DCP_data.Properties.VariableNames, DCP_data_temp.Properties.VariableNames)) && length(DCP_data.Properties.VariableNames) < length(DCP_data_temp.Properties.VariableNames)
                             %    temp_col = setxor(DCP_data.Properties.VariableNames, DCP_data_temp.Properties.VariableNames);
                             %    DCP_data = addvars(DCP_data, NaN(size(DCP_data,1),1), 'NewVariableNames', temp_col);
                             % end
                             common_cols = intersect(DCP_data.Properties.VariableNames,DCP_data_temp.Properties.VariableNames);
                             DCP_data = outerjoin(DCP_data(:,common_cols),DCP_data_temp(:,common_cols),'keys',common_cols,'MergeKeys',true);
                             % Remove duplicates and sort the data, update
                             % current_time
                             DCP_data = unique(DCP_data);
                             DCP_data = sortrows(DCP_data,'Time');
                             current_time = datetime('now');
                        else
                            DCP_data = DCP_data_temp;
                        end
                        if DCP_data.Time(end)-DCP_data.Time(1) > hours(168)
                            DCP_data(1,:) = [];
                        end
                        app.HADS_gauges = DCP_data;
                    end                         
                end
            end

            if flags.flag_waterquality == 1
                WQ_States.outet_pollutograph(t_store,1) = Out_Conc;
                % Saving Data of Input Gauges
                if flags.flag_obs_gauges == 1
                    for i = 1:gauges.num_obs_gauges
                        gauges.x_cell = gauges.northing_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                        WQ_States.pollutograph_cell(t_store,i)= WQ_States.P_conc(gauges.y_cell,gauges.x_cell); % mg/L
                    end
                end
            end
        end % Calls the sub
    end
    % Saving Maximum Depths and Outlet Flows
    if flags.flag_automatic_calibration ~= 1
        if k == 1
            depths.dmax_final = depths.d_t;
            velocities.vmax_final = velocities.velocity_raster;
            if flags.flag_waterquality == 1
                WQ_States.max_Conc = zeros(size(Elevation_Properties.elevation_cell));
            end
        else
            depths.dmax_final = max(depths.d_t,depths.dmax_final);
            velocities.vmax_final = max(velocities.velocity_raster,velocities.vmax_final);
            if flags.flag_waterquality == 1
                WQ_States.max_Conc = max(WQ_States.max_Conc,WQ_States.P_conc);
            end
        end
    end
    outlet_runoff_volume = nansum(nansum(outlet_states.outlet_flow))*time_step/60*Wshed_Properties.cell_area/Wshed_Properties.drainage_area + outlet_runoff_volume; % mm

    % Previous Depths
    depths.d_p = depths.d_t;
    Soil_Properties.I_p = Soil_Properties.I_t;

    % Runoff Coefficient Calculation
    BC_States.outflow_volume  = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600*time_step*60 + BC_States.outflow_volume ;
    if flags.flag_spatial_rainfall == 1
        if flags.flag_inflow == 1
            BC_States.inflow_volume = sum(sum(BC_States.delta_inflow_agg'/1000.*Wshed_Properties.n_inlets*Wshed_Properties.cell_area)) + nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area + BC_States.inflow_volume; % check future
        else
            BC_States.inflow_volume =  nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area + BC_States.inflow_volume; % check future
        end
    elseif flags.flag_spatial_rainfall ~= 1 && flags.flag_inflow == 1
        BC_States.inflow_volume = sum(sum(BC_States.delta_inflow_agg'/1000.*Wshed_Properties.n_inlets*Wshed_Properties.cell_area)) + BC_States.delta_p_agg/1000*Wshed_Properties.drainage_area + BC_States.inflow_volume; % check future
    else
        BC_States.inflow_volume = BC_States.delta_p_agg/1000*Wshed_Properties.drainage_area + BC_States.inflow_volume; % check future
    end
    % Increase the counter
    if(nansum(nansum(depths.d_t))) == 0
        ttt = 1;
    end

    % Saving Results in time_observations - Only Valid for Calibration
    if flags.flag_automatic_calibration == 1
        z1_save = find(time_obs <= t,1,'last');
        if ~isempty(z1_save)
            if z1_save > z2_save
                Qmod(z1_save,1) = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600; % m3/s
                if flags.flag_waterquality == 1
                    Cmod(z1_save,1) = Out_Conc; % mg/L
                else
                    Cmod(z1_save,1) = 0; % mg/L
                end
                % Observed Gauges
                for jj = 1:length(northing_obs_cell)
                    Qmod_gauges(z1_save,jj) = CA_States.I_tot_end_cell(northing_obs_cell(jj),easting_obs_cell(jj))/(running_control.time_calculation_routing(k));
                    if flags.flag_waterquality == 1
                        Cmod_gauges(z1_save,jj) = Maps.WQ_States.Pol_Conc(northing_obs_cell(jj),easting_obs_cell(jj));
                    else
                        Cmod_gauges(z1_save,jj) = 0;
                    end
                end
            end
            z2_save = z1_save;
        else
            z1_save = 0;
            z2_save = z1_save;
        end
    end

    if isnan(depths.d_t)
        ttt = 1
    end

    % Refreshing Time-step
    if flags.flag_real_time_satellite_rainfall == 1
        if t <= 360
            t = running_control.time_calculation_routing(k,1)/60 + t;
            time_step_save(k,2) = running_control.time_calculation_routing(k,1);
            time_step_save(k,1) = t;
            if t >= (running_control.routing_time + running_control.min_time_step/60);
                t = t-360;
                % t_previous = running_control.time_calculation_routing(k,1)/60;
                % t_previous = 0;
                recording_parameters.last_record_maps = 1;
                end_flag=1;
                % time_step_save(k,2) = running_control.time_calculation_routing(k,1);
                % time_step_save(k,1) = t;
            end
        else
            t = t-360;
            recording_parameters.last_record_maps = 1;
            end_flag=1;
        end
        t2 =  running_control.time_calculation_routing(k,1)/60 + t2;
    else
        t = running_control.time_calculation_routing(k,1)/60 + t;
        time_step_save(k,2) = running_control.time_calculation_routing(k,1);
        time_step_save(k,1) = t;
    end

    k = k + 1;
end

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
        gauges.labels_observed_string = extra_parameters.gauges.labels_observed_string;
    end
    GIS_data = structfun(@gather, GIS_data, 'UniformOutput', false);
    if flags.flag_human_instability == 1
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
    nx_max = gather(nx_max);
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
