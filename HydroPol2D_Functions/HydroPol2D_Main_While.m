%% Main HydroPol2D While Loop
% Developer: Marcus Nobrega
% Date 5/21/2023
% Goal - Run the main modeling process of the model
format short g

tic 
k = 1;
C = 0;
t = running_control.time_step_model;
Flooded_Area = 0;
velocities.velocity_raster = 0;
Risk_Area = 0;
saver_count = 1; % starts in 1 but the next pointer should be 2, this is auto fixed when t reach the next time aggregation.
store = 1;

%If no dam breaking is simulated
if flags.flag_dam_break == 1
    input_table_breaker = readtable('Damns_breaks_points.xlsx'); 
    breakers.easting_absolute = table2array(input_table_breaker(:,2)); % Easting Coordinates
    breakers.northing_absolute = table2array(input_table_breaker(:,3)); % Northing Coordinates
    breakers.ID = table2array(input_table_breaker(:,1)); % ID
    % --- Converting coordinates to local coordinates in pixels
    breakers.easting= round((-GIS_data.xulcorner + breakers.easting_absolute)/Wshed_Properties.Resolution);
    breakers.northing = round((GIS_data.yulcorner - breakers.northing_absolute)/Wshed_Properties.Resolution);
    flag_break_1 = 1; flag_break_2 = 1;    
end

n_snaps = 40;
dt_snap = running_control.routing_time/n_snaps;
time_snap = [1:1:n_snaps]*dt_snap;
z2_snap = 0;
while t <= (running_control.routing_time + running_control.min_time_step/60)

    % Snap Results
    z1_snap = find(time_snap>=t,1,'first');
    if z1_snap > z2_snap
       Snapshot_Results
       pause(0.1);
    end
    z2_snap = z1_snap;

    % Infiltration and Available Depth
    % Show stats
    
    if flags.flag_waterquality == 1
        perc______duremain______tsec_______dtmm______infmmhr____CmgL_____dtmWQ = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(C)), max(max((WQ_States.P_conc))), tmin_wq]
    else
        perc______duremain______tsec_______dtmm______infmmhr____Vmax = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(C)), max(max(velocities.velocity_raster))]
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
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.f*time_step/60 + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_0 + depths.pef;
            Soil_Properties.T = Soil_Properties.Tr; % Beginning to track the replenishing time
        else
            %Hydro_States.i_a = (BC_States.delta_p_agg*rainfall_matrix + BC_States.delta_inflow_agg*Wshed_Properties.inflow_cells + d_0)./(time_step/60);
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
            depths.d_t = depths.d_0 + depths.pef;
            if depths.d_t < 0
                ttt = 1;
            end
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
            depths.pef = BC_States.delta_p_agg.*Wshed_Properties.rainfall_matrix + BC_States.inflow - Hydro_States.ETP/24  + idx_rivers*Wshed_Properties.Resolution/1000*Lateral_Groundwater_Flux;
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
        % [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,outflow_rate.qout_ne_t,outflow_rate.qout_se_t,outflow_rate.qout_sw_t,outflow_rate.qout_nw_t,depths.d_t,CA_States.I_tot_end_cell,~] = CA_Routing_8D_not_optimized(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.p,flags.flag_reservoir,Elevation_Properties.elevation_cell,depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical,wse_slope_zeros,Distance_Matrix);
        [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,outflow_rate.qout_ne_t,outflow_rate.qout_se_t,outflow_rate.qout_sw_t,outflow_rate.qout_nw_t,depths.d_t,CA_States.I_tot_end_cell] = CA_Routing_8D(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.pv,flags.flag_reservoir,Elevation_Properties.elevation_cell,depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical,Reservoir_Data.Ko,Reservoir_Data.po,CA_States.depth_tolerance);
    else % 4-D CA Model
        % [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,CA_States.I_cell] = CA_Routing(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.p,flags.flag_reservoir,Elevation_Properties.elevation_cell,depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical);
        [outflow_rate.qout_left_t,outflow_rate.qout_right_t,outflow_rate.qout_up_t,outflow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,velocity_term] = CA_Routing(Reservoir_Data.Dir,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.Kv,Reservoir_Data.p,flags.flag_reservoir,Elevation_Properties.elevation_cell,depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical,velocity_term,Reservoir_Data.Ko,Reservoir_Data.po);
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

    if min(min(outflow_rate.qin_t)) < 0
        ttt = 1;
    end

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

            flow_depth = depths.d_t; % Depth in which velocities will be calculated (mm)
            flow_depth(depths.d_t < CA_States.depth_tolerance) = 1e12;
            velocities.vel_left = (outflow_rate.qout_left_t/1000/3600)*Wshed_Properties.Resolution^2./(Wshed_Properties.Resolution*flow_depth/1000); % m/s
            velocities.vel_right = (outflow_rate.qout_right_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
            velocities.vel_up = (outflow_rate.qout_up_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
            velocities.vel_down = (outflow_rate.qout_down_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
            
            if flags.flag_D8 == 1
                velocities.vel_ne = (outflow_rate.qout_ne_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
                velocities.vel_se = (outflow_rate.qout_se_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
                velocities.vel_sw = (outflow_rate.qout_sw_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
                velocities.vel_nw = (outflow_rate.qout_nw_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
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
            % Old max velocity
            if k > 1
                old_velocity = velocities.max_velocity;
            end
                
            velocities.max_velocity = max(velocities.velocity_vector);

            % Checking changes in max velocity
            if k > 1
                if velocities.max_velocity > 25 && old_velocity > 0
                    warning('Velocities larger than 25  m per sec. Possible instability.')
                    subplot(2,1,1)
                    [row_maxvel, col_maxvel] = find(velocities.velocity_raster == max(max(velocities.velocity_raster)));
                    [X, Y] = meshgrid(1:1:size(DEM_raster.Z,2),1:1:size(DEM_raster.Z,1));
                    hold on
                    surf_plot(max(max(velocities.velocity_raster)),t,'v','m/s',velocities.velocity_raster,1,0,32,0.85,0,[0 90],X,Y)
                    subplot(2,1,2)
                    surf_plot(max(max(depths.d_t/1000)),t,'d','m',depths.d_t/1000,1,0,32,0.85,0,[0 90],X,Y)
                    pause(0.1)
                end
            end

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

            if velocities.max_velocity == 0 && k == 1
                new_timestep = running_control.time_step_model*60; % Seconds
            elseif velocities.max_velocity == 0
                new_timestep = running_control.max_time_step; % Seconds
            end

            % ---- Calculation of Stability --- %
            if flags.flag_human_instability == 1
                
            elseif flags.flag_human_instability == 2
            elseif flags.flag_human_instability == 3
                Human_Instability.risk_t_cm = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_c_m,Human_Instability.y_c_m,Human_Instability.w_c_m,Human_Instability.d_c_m);
                Human_Instability.risk_t_tm = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_t_m,Human_Instability.y_t_m,Human_Instability.w_t_m,Human_Instability.d_t_m);
                Human_Instability.risk_t_am = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_a_m,Human_Instability.y_a_m,Human_Instability.w_a_m,Human_Instability.d_a_m);
                Human_Instability.risk_t_om = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_o_m,Human_Instability.y_o_m,Human_Instability.w_o_m,Human_Instability.d_o_m);
                Human_Instability.risk_t_cf = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_c_f,Human_Instability.y_c_f,Human_Instability.w_c_f,Human_Instability.d_c_f);
                Human_Instability.risk_t_tf = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_t_f,Human_Instability.y_t_f,Human_Instability.w_t_f,Human_Instability.d_t_f);
                Human_Instability.risk_t_af = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_a_f,Human_Instability.y_a_f,Human_Instability.w_a_f,Human_Instability.d_a_f);
                Human_Instability.risk_t_of = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_o_f,Human_Instability.y_o_f,Human_Instability.w_o_f,Human_Instability.d_o_f);
            end
        else          
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
            z1 = find(BC_States.time_deltainflow > t_previous,1,'first'); % begin of the time-step
            z2 = find(BC_States.time_deltainflow <= t,1,'last'); % end of the time-step
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
            % if BC_States.delta_inflow_agg(z,1) > 0
            %     ttt = 1
            % end
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
                    % Maps.Hydro.spatial_rainfall_maps(:,:,zz2) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                    % zzz = Maps.Hydro.spatial_rainfall_maps(:,:,zz2);
                    Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = mean(rainfall_spatial_aggregation,3);
                    zzz = Maps.Hydro.spatial_rainfall_maps(:,:,saver_count);
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

            % Times for aggregated rainfall dataset
            zz1 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t_previous,1,'last');
            zz2 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t,1,'last');

            if Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2)  == Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) % No difference between them, same time-step
                Rainfall_Parameters.index_aggregation = 1; % We save only 1 map
            elseif zz1 == zz2 && z1_input ~= z2_input % We save more than 1 map
                Rainfall_Parameters.index_aggregation = Rainfall_Parameters.index_aggregation + (z2_input-z1_input);
            elseif zz2 > zz1
                Rainfall_Parameters.index_aggregation = 1;
            end            
            
            if z2_input > z1_input
              % Read Geotiff
                % input_rainfall = GRIDobj(Input_Rainfall.labels_Directory{z2_input}{1});
                try
                    [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{z2_input}{1}),'CoordinateSystemType','geographic');
                    % Reproject the coordinates from EPSG:4326 to EPSG:3857
                    rR.GeographicCRS=geocrs(4326);
                    if rR.CellExtentInLatitude ~= GIS_data.resolution_resample
                        input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,1);
                    end
                catch
                    [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{z2_input}{1}));
                    if rR.CellExtentInWorldX ~= GIS_data.resolution_resample
                        input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,0);
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
                if zz2 > zz1 % Saving Maps
                    % Maps.Hydro.spatial_rainfall_maps(:,:,zz2) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                    % zzz = Maps.Hydro.spatial_rainfall_maps(:,:,zz2);
                    Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = mean(rainfall_spatial_aggregation,3);
                    zzz = Maps.Hydro.spatial_rainfall_maps(:,:,saver_count);
                    BC_States.average_spatial_rainfall(zz2,1) = mean(zzz(zzz>=0));
                end                 
                input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;
                % Maps.Hydro.spatial_rainfall_maps(:,:,z2_input) = input_rainfall;
                BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
                % BC_States.average_spatial_rainfall(z2_input,1) = mean(input_rainfall(input_rainfall>=0));
                rainfall_spatial_aggregation(:,:,Rainfall_Parameters.index_aggregation) = gather(input_rainfall); % Saving high resolution map              
            end
        elseif flags.flag_input_rainfall_map == 0 && flags.flag_satellite_rainfall == 1 && flags.flag_real_time_satellite_rainfall ~= 1
            % Satellite Rainfall
            if flags.flag_satellite_rainfall == 1
                % Maps with generally coarser resolution
                recording_parameters.actual_record_state_rainfall = find(running_control.time_records < t,1,'last');
                recording_parameters.delta_record_rainfall = recording_parameters.actual_record_state_rainfall - recording_parameters.last_record_maps_rainfall;
                recording_parameters.last_record_maps_rainfall = recording_parameters.actual_record_state_rainfall; 

                % Input Rainfall Maps
                % Times of full rainfall dataset
                z1_input = find(Input_Rainfall.time <= t_previous,1,'last');
                z2_input = find(Input_Rainfall.time <= t,1,'last');
                if z2_input > z1_input
                    product = 'PDIRNow1hourly';
                    [rainfall_raster, register,~] = Satellite_rainfall_processing([],[],register,product,date_begin,date_end,flags.flag_satellite_rainfall,flags.flag_real_time_satellite_rainfall,DEM_raster);
                    input_rainfall = rainfall_raster.Z;
                    input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;

                    if recording_parameters.delta_record_rainfall > 0
                        record_map_indice = recording_parameters.last_record_maps_rainfall;
                        % Maps.Hydro.spatial_rainfall_maps(:,:,record_map_indice) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                        % zzz = Maps.Hydro.spatial_rainfall_maps(:,:,record_map_indice);
                        Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                        zzz = Maps.Hydro.spatial_rainfall_maps(:,:,saver_count);
                        BC_States.average_spatial_rainfall(record_map_indice,1) = mean(zzz(zzz>=0));
                        Rainfall_Parameters.index_aggregation = 0; % Returning to first rainfall
                    end 

                    Rainfall_Parameters.index_aggregation = Rainfall_Parameters.index_aggregation + (z2_input-z1_input);
                    % Map from previous timestep
                    rainfall_spatial_aggregation(:,:,Rainfall_Parameters.index_aggregation) = input_rainfall;
                     
                end
            end
        elseif flags.flag_input_rainfall_map == 0 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall == 1
            if flags.flag_real_time_satellite_rainfall == 1
                %%%%%%%%%%%%%%%%%
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
                % Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
                Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP; % Saving ETP Maps
                ETR_save(:,:,z2) = Hydro_States.ETR ;
            elseif  (isempty(z1) && z2 > 0) || z2 > z1 && z2 < length(ETP_Parameters.climatologic_spatial_duration)
                if flags.flag_GPU == 1 || flags.flag_single == 1
                    day_of_year = day(extra_parameters.ETP.time_ETP(z2,1),'dayofyear');
                else
                    day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
                end
                [Hydro_States.ETP] = ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2),Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations,ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations,ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);
                if nansum(nansum(Hydro_States.ETP)) == 0
                    if saver_count==1
                        Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,12);
                    else
                        Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,saver_count-1);
                    end
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
                % Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
                Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP; % Saving ETP Maps
                ETR_save(:,:,z2) = Hydro_States.ETR ;
            elseif z1 == 1 && z2 == 1 && k == 1
                % First data. We assume a constant ETP using
                % 1st data
                if flags.flag_GPU == 1 || flags.flag_single == 1
                    day_of_year = day(extra_parameters.ETP.time_ETP(z2,1),'dayofyear');
                end
                [Hydro_States.ETP] = ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2),Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations,ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations,ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);
                if nansum(nansum(Hydro_States.ETP)) == 0
                    if saver_count==1
                        Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,12);
                    else
                        Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,saver_count-1);
                    end
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,z2-1);
                    Hydro_States.ETR  = ETR_save(:,:,z2-1);
                end
                % Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
                Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP; % Saving ETP Maps
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
    if min(min(depths.d_t)) < 0 
        ttt = 1;
    end

    % Clearing stored values
    if flags.flag_GPU == 1
        outflow_rate.qout_left_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_right_t = gpuArray(zeros(ny_max,nx_max));
        outflow_rate.qout_up_t = gpuArray( zeros(ny_max,nx_max));
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

    % Maps with generally coarser resolution
    recording_parameters.actual_record_state = find(running_control.time_records < t_save,1,'last');
    recording_parameters.delta_record = recording_parameters.actual_record_state - recording_parameters.last_record_maps;
    recording_parameters.last_record_maps = recording_parameters.actual_record_state;

    if flags.flag_automatic_calibration ~= 1
        if k == 1 % First time-step
            Maps.Hydro.d(:,:,1) = depths.d_t;
            if flags.flag_human_instability == 1
                Maps.Hydro.risk(:,:,1) = Human_Instability.risk_t;
            elseif flags.flag_human_instability == 2
            elseif flags.flag_human_instability == 3
                Maps.Hydro.risk_cm(:,:,1) = Human_Instability.risk_t_cm;
                Maps.Hydro.risk_tm(:,:,1) = Human_Instability.risk_t_tm;
                Maps.Hydro.risk_am(:,:,1) = Human_Instability.risk_t_am;
                Maps.Hydro.risk_om(:,:,1) = Human_Instability.risk_t_om;
                Maps.Hydro.risk_cf(:,:,1) = Human_Instability.risk_t_cf;
                Maps.Hydro.risk_tf(:,:,1) = Human_Instability.risk_t_tf;
                Maps.Hydro.risk_af(:,:,1) = Human_Instability.risk_t_af;
                Maps.Hydro.risk_of(:,:,1) = Human_Instability.risk_t_of;
            end
            Maps.Hydro.I_t(:,:,1) = Soil_Properties.I_t;
        elseif recording_parameters.delta_record > 0 % Saving maps
            t_store = recording_parameters.actual_record_state;
            % Maps.Hydro.d(:,:,t_store) = depths.d_t;
            Maps.Hydro.d(:,:,saver_count) = depths.d_t;
            if flags.flag_human_instability == 1
                % Maps.Hydro.risk(:,:,t_store) = Human_Instability.risk_t;
                Maps.Hydro.risk(:,:,saver_count) = Human_Instability.risk_t;
            elseif flags.flag_human_instability == 2
            elseif flags.flag_human_instability == 3
                Maps.Hydro.risk_cm(:,:,saver_count) = Human_Instability.risk_t_cm;
                Maps.Hydro.risk_tm(:,:,saver_count) = Human_Instability.risk_t_tm;
                Maps.Hydro.risk_am(:,:,saver_count) = Human_Instability.risk_t_am;
                Maps.Hydro.risk_om(:,:,saver_count) = Human_Instability.risk_t_om;
                Maps.Hydro.risk_cf(:,:,saver_count) = Human_Instability.risk_t_cf;
                Maps.Hydro.risk_tf(:,:,saver_count) = Human_Instability.risk_t_tf;
                Maps.Hydro.risk_af(:,:,saver_count) = Human_Instability.risk_t_af;
                Maps.Hydro.risk_of(:,:,saver_count) = Human_Instability.risk_t_of;
            end
            % Maps.Hydro.I_t(:,:,t_store) = Soil_Properties.I_t;
            Maps.Hydro.I_t(:,:,saver_count) = Soil_Properties.I_t;
            if flags.flag_waterquality ~= 1
                saver_count = saver_count+1;
                if saver_count > 12
                    saver_count = 1;
                    save(strcat('Temporary_Files\save_map_hydro_',num2str(store),'.mat'),'Maps');
                    store = store + 1;
                end
            end
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
                % Maps.WQ_States.Pol_Conc_Map(:,:,t_store) = WQ_States.P_conc;
                % Maps.WQ_States.Pol_mass_map(:,:,t_store) = WQ_States.B_t;
                Maps.WQ_States.Pol_Conc_Map(:,:,saver_count) = WQ_States.P_conc;
                Maps.WQ_States.Pol_mass_map(:,:,saver_count) = WQ_States.B_t;
                saver_count = saver_count+1;
                if saver_count > 12
                    saver_count = 1;
                    save(strcat('Temporary_Files\save_map_hydro_',num2str(store),'.mat'),'Maps');
                    store = store + 1;
                end

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
        ttt = 1;
    end

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
