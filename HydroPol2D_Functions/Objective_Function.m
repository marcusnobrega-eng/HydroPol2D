%% Automatic Calibration Function
% - Developer: Marcus Nobrega

% Objective Function
function [Obj_function_value] = Objective_Function(x,observed_concentration,time_observed_concentration)
% Inputs:
% x:= decision vector collecting C3 and C4
% observed_concentration = vector of concentration in mg/L
% time_observed_concentration = vector with the respective times in minutes

% Note:
% Please, load the correct workspace below with all HydroPol2D inputs
% besides the decision variables C3 and C4. Also, please make sure that d0
% and h0 are entered below.

% h0 and d0 values
h0 = 0; % mm
d0 = 0.2; % mm

% ----- LOAD THE WORKSPACES BELOW ---- 
% --- 0.5 deg case --- %
% load workspace_73_34mmh.mat

% --- 2 deg case --- %
load workspace_2deg_62_61mmh


% Values Changed for Calibration
C_3 = double(idx_not_nan)*x(1);
C_4 = double(idx_not_nan)*x(2);

obs_index_prev = 0;
h_0 = double(idx_not_nan)*h0;
d_0 = double(idx_not_nan)*d0;

while t <= (routing_time + min_time_step/60)
    % Infiltration and Available Depth
    % Show stats
    if tmin_wq < 0 || isnan(tmin_wq) || isinf(tmin_wq)
        error('Instability')
    end
    % Inflows
    if flag_inflow == 1
        inflow = zeros(size(dem,1),size(dem,2)); % This is to solve spatially, don't delete
        for i = 1:n_stream_gauges
            inflow = inflow + delta_inflow_agg(i)*inflow_cells(:,:,i);
        end
    end

    if k == 1
        if flag_infiltration == 1
            i_a = (delta_p_agg.*rainfall_matrix + inflow + d_0 - ETP/24)./(time_step/60);
            C = ksat.*(1 + ((d_0 + psi).*(teta_sat - teta_i))./I_0); % matrix form
            f = min(C,i_a);
            I_t = I_0 + (f - ETP)*time_step/60; % Extracting ETP from soil
            pef = delta_p_agg.*rainfall_matrix + inflow - f*time_step/60;
            d_t = d_0 + pef;
            inf_m3s_t = f/1000.*cell_area/3600;
            T = Tr; % Beginning to track the replenishing time
        else
            %i_a = (delta_p_agg*rainfall_matrix + delta_inflow_agg*inflow_cells + d_0)./(time_step/60);
            pef = delta_p_agg.*rainfall_matrix + inflow;
            d_t = d_0 + pef;
        end
    else
        if flag_infiltration == 1
            % Effective precipitation - Green-Ampt(1911)
            i_a = (delta_p_agg.*rainfall_matrix + inflow + d_p - ETP/24)./(time_step/60);
            C = ksat.*(1 + ((d_p + psi).*(teta_sat - teta_i))./I_p); % matrix form
            idx_C = (i_a <= C); % Values where i_a is below C
            idx_i_a = (i_a > C); % Values where i_a is larger than C
            T = T - time_step/60/60; % Recoverying time (hours)
            T(idx_i_a) = Tr(idx_i_a); % Saturated Areas we begin the process again
            idx_T = T < 0; % Cells where the replenishing time is reached
            f = min(C,i_a); % Infiltration rate (mm/hr)
            I_t = max(I_p + f*time_step/60 - ETP/24 - k_out.*double(idx_C)*time_step/60,0);
            % Refreshing I_t to I_0 for cases where idx_T > 0
            I_t(idx_T) = I_0(idx_T);
            pef = delta_p_agg.*rainfall_matrix + inflow - f*time_step/60;
            d_t = d_p + pef; %% ATTENTION HERE
            inf_m3s_t = f/1000.*cell_area/3600;
        else
            i_a = (delta_p_agg.*rainfall_matrix + inflow + d_p - ETP/24)./(time_step/60);
            pef = delta_p_agg.*rainfall_matrix + inflow - ETP/24;
            d_t = d_p + pef;
        end
    end
    d_tot = d_t;

    %%%% WATER DEPTHS %%%
    d_left_cell = [zeros(ny_max,1),d_tot(:,1:(nx_max-1))];
    d_right_cell = [d_tot(:,(2:(nx_max))) zeros(ny_max,1)];
    d_up_cell = [zeros(1,nx_max) ; d_tot(1:(ny_max-1),:)];
    d_down_cell = [d_tot(2:(ny_max),:) ; zeros(1,nx_max)];
    if flag_D8 == 1
        % Simulate with D-8 flow direction
        % --- Adding NE, SE, SW, NW --- %
        zero_matrix = NaN(size(spatial_domain));
        elevation_NE_t = zero_matrix; elevation_SE_t = zero_matrix;
        elevation_SW_t = zero_matrix; elevation_NW_t = zero_matrix;
        d_NE_t = zero_matrix; d_SE_t = zero_matrix; d_SW_t = zero_matrix; d_NW_t = zero_matrix;

        elevation_NE_t(2:(ny_max),1:(nx_max-1)) = elevation(1:(ny_max-1),2:nx_max); % OK
        elevation_SE_t(1:(ny_max-1),1:(nx_max-1)) = elevation(2:ny_max,2:nx_max); % OK
        elevation_SW_t(1:(ny_max-1),2:(nx_max)) = elevation(2:(ny_max),1:(nx_max-1)); % OK
        elevation_NW_t(2:ny_max,2:nx_max) = elevation(1:(ny_max-1),1:(nx_max-1)); % OK

        d_NE_t(2:(ny_max),1:(nx_max-1)) = d_tot(1:(ny_max-1),2:nx_max); % OK
        d_SE_t(1:(ny_max-1),1:(nx_max-1)) = d_tot(2:ny_max,2:nx_max); % OK
        d_SW_t(1:(ny_max-1),2:(nx_max)) = d_tot(2:(ny_max),1:(nx_max-1)); % OK
        d_NW_t(2:ny_max,2:nx_max) = d_tot(1:(ny_max-1),1:(nx_max-1)); % OK
    else
        % Simulate with D-4 flow direction
        % Everything already calculated
    end

    %  Flood routing for each cell
    % Check if Nan or Inf occured in d_t
    if max(max(isinf(d_t))) == 1
        ttt = 1;
    end
    if flag_D8 == 1 % D-8 C-A Model
        [qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,qout_ne_t,qout_se_t,qout_sw_t,qout_nw_t,d_t,I_tot_end_cell,I_cell] = CA_Routing_8D(elevation_cell,elevation_left_t,elevation_right_t,elevation_up_t,elevation_down_t,d_tot,d_left_cell,d_right_cell,d_up_cell,d_down_cell,roughness,cell_area,time_step,h_0,Resolution,I_tot_end_cell,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,idx_nan,flag_critical,elevation_NE_t, elevation_SE_t, elevation_SW_t, elevation_NW_t, d_NE_t, d_SE_t, d_SW_t, d_NW_t);
    else % 4-D CA Model
        [qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,d_t,I_tot_end_cell,I_cell] = CA_Routing(elevation_cell,elevation_left_t,elevation_right_t,elevation_up_t,elevation_down_t,d_tot,d_left_cell,d_right_cell,d_up_cell,d_down_cell,roughness,cell_area,time_step,h_0,Resolution,I_tot_end_cell,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,idx_nan,flag_critical);
    end
    %     qout_t = qout_left + qout_right + qout_up + qout_down + qout_ne + qout_se + qout_sw + qout_ne;
    % Inflows - Von-Neuman
    qin_left_t = [zeros(ny_max,1),qout_right_t(:,1:(nx_max-1))];
    qin_right_t = [qout_left_t(:,(2:(nx_max))) zeros(ny_max,1)];
    qin_up_t = [zeros(1,nx_max) ; qout_down_t(1:(ny_max-1),:)];
    qin_down_t = [qout_up_t(2:(ny_max),:) ; zeros(1,nx_max)];
    % Inflows - Inclined Directions
    if flag_D8 == 1
        zero_matrix = zeros(size(spatial_domain));
        qin_ne_t = zero_matrix; qin_se_t = zero_matrix; qin_sw_t = zero_matrix; qin_nw_t = zero_matrix;
        qin_ne_t(2:(ny_max),1:(nx_max-1)) = qout_sw_t(1:(ny_max-1),2:(nx_max)); % OK
        qin_se_t(1:(ny_max-1),1:(nx_max-1)) = qout_nw_t(2:ny_max,2:nx_max); % OK
        qin_sw_t(1:(ny_max-1),2:(nx_max)) = qout_ne_t(2:(ny_max),1:(nx_max-1)); % OK
        qin_nw_t(2:ny_max,2:nx_max) = qout_se_t(1:(ny_max-1),1:(nx_max-1)); % OK
    end
    if flag_D8 == 1
        qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t + qin_ne_t + qin_se_t + qin_sw_t + qin_nw_t;
        %         qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t;
    else
        qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t;
    end
    idx = qin_t > flow_tolerance;
    idx2 = qin_t <= flow_tolerance;
    idx3 = logical(isnan(qin_t) + isinf(qin_t));
    qin_t(idx3) = 0;
    flows_cells(idx) = 1;
    flows_cells(idx2) = 0;


    % Water Quality Parameters for f(B(t))
    if flag_waterquality == 1
        if flag_D8 ~= 1
            [B_t,P_conc,Out_Conc,tmin_wq,tot_W_out,mass_lost,Tot_Washed] = build_up_wash_off(C_3,C_4,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,B_t,time_step,nx_max,ny_max,cell_area,outlet_index,idx_nan_5,flag_wq_model,mass_lost,Tot_Washed,Bmin,Bmax,min_Bt);
        else
            [B_t,P_conc,Out_Conc,tmin_wq,tot_W_out,mass_lost,Tot_Washed] = build_up_wash_off_8D(C_3,C_4,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,qout_ne_t,qout_se_t,qout_sw_t,qout_nw_t,B_t,time_step,nx_max,ny_max,cell_area,outlet_index,idx_nan_5,flag_wq_model,mass_lost,Tot_Washed,Bmin,Bmax,min_Bt);
        end
    end
    %%%% Checking Mass Balance
    if flag_waterquality == 1
        if sum(sum(B_t(~isinf(B_t)))) > 1.2*initial_mass  % More than 5%
            error('Brutal instability in B(t). More than 20% difference')
        end
    end
    % New Time-step Calculation
    pos_save = find(time_change_records < t,1,'last');
    time_save = time_change_records(pos_save); % min
    delta_time_save = time_save - time_save_previous;
    time_save_previous = time_save;
    actual_record_timestep = find(time_change_records < t,1,'last');
    if delta_time_save > 0 || k == 1 % First time-step
        if flag_timestep == 0
            %%% SOLUTION FOR COURANT METHOD %%%
            d_left_cell = [zeros(ny_max,1),d_t(:,1:(nx_max-1))];
            d_right_cell = [d_t(:,(2:(nx_max))) zeros(ny_max,1)];
            d_up_cell = [zeros(1,nx_max) ; d_t(1:(ny_max-1),:)];%
            d_down_cell = [d_t(2:(ny_max),:) ; zeros(1,nx_max)];
            if flag_D8 == 1
                d_NE_cell = zero_matrix; d_SE_cell = zero_matrix; d_SW_cell = zero_matrix; d_NW_cell = zero_matrix;
                d_NE_cell(2:(ny_max),1:(nx_max-1)) = d_t(1:(ny_max-1),2:nx_max); % OK
                d_SE_cell(1:(ny_max-1),1:(nx_max-1)) = d_t(2:ny_max,2:nx_max); % OK
                d_SW_cell(1:(ny_max-1),2:(nx_max)) = d_t(2:(ny_max),1:(nx_max-1)); % OK
                d_NW_cell(2:ny_max,2:nx_max) = d_t(1:(ny_max-1),1:(nx_max-1)); % OK
            end
            % Considering Velocity as Celerity
            %             vel_left = (I_cell(:,:,1)./(0.5/1000.*(d_t + d_left_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_left_cell)/2/1000); % m/s
            %             vel_right = (I_cell(:,:,2)./(0.5/1000.*(d_t + d_right_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_right_cell)/2/1000);
            %             vel_up = (I_cell(:,:,3)./(0.5/1000.*(d_t + d_up_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_up_cell)/2/1000); % m/s;;
            %             vel_down = (I_cell(:,:,4)./(0.5/1000.*(d_t + d_down_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_down_cell)/2/1000); % m/s;

            % Using CA Velocities
            %             vel_left = (qout_left_t/1000/3600); % m/s
            %             vel_right = (qout_right_t/1000/3600); % m/s
            %             vel_up = (qout_up_t/1000/3600); % m/s
            %             vel_down = (qout_down_t/1000/3600); % m/s

            z = d_t;
            dt_threshold = 0.05; % mm
            z(d_t < dt_threshold) = 1e12;
            vel_left = (qout_left_t/1000/3600)*Resolution^2./(Resolution*z/1000); % m/s
            vel_right = (qout_right_t/1000/3600)*Resolution./(z/1000); % m/s
            vel_up = (qout_up_t/1000/3600)*Resolution./(z/1000); % m/s
            vel_down = (qout_down_t/1000/3600)*Resolution./(z/1000); % m/s

            if flag_D8 == 1
                %                 vel_ne = (I_cell(:,:,6)./(0.5/1000.*(d_t + d_NE_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_NE_cell)/2/1000); % m/s
                %                 vel_se = (I_cell(:,:,7)./(0.5/1000.*(d_t + d_SE_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_SE_cell)/2/1000);
                %                 vel_sw = (I_cell(:,:,8)./(0.5/1000.*(d_t + d_SW_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_SW_cell)/2/1000); % m/s;;
                %                 vel_nw = (I_cell(:,:,9)./(0.5/1000.*(d_t + d_NW_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_NW_cell)/2/1000); % m/s;
                vel_ne = (qout_ne_t/1000/3600); % m/s
                vel_se = (qout_se_t/1000/3600); % m/s
                vel_sw = (qout_sw_t/1000/3600); % m/s
                vel_nw = (qout_nw_t/1000/3600); % m/s
            end
            %%%%%%%%%%%%%% Find the Maximum Velocity
            max_velocity_left = max(max(vel_left));
            max_velocity_right = max(max(vel_right));
            max_velocity_up = max(max(vel_up));
            max_velocity_down = max(max(vel_down));
            if flag_D8 == 1
                max_velocity_ne = max(max(vel_ne));
                max_velocity_se = max(max(vel_se));
                max_velocity_sw = max(max(vel_sw));
                max_velocity_nw = max(max(vel_nw));
            end
            % - Velocit Raster - %
            %%% Maximum of All of Them %%%
            velocity_raster = max(vel_left,vel_right);
            velocity_raster = max(velocity_raster,vel_up);
            velocity_raster = max(velocity_raster,vel_down);
            if flag_D8 == 1
                velocity_raster = max(velocity_raster,vel_ne);
                velocity_raster = max(velocity_raster,vel_se);
                velocity_raster = max(velocity_raster,vel_sw);
                velocity_raster = max(velocity_raster,vel_nw);
                velocity_vector = [max_velocity_left, max_velocity_right, max_velocity_up, max_velocity_down, max_velocity_ne, max_velocity_se, max_velocity_sw, max_velocity_nw];
            else
                velocity_vector = [max_velocity_left, max_velocity_right, max_velocity_up, max_velocity_down];
            end

            max_velocity = max(velocity_vector);
            if k == 1
                dvdt = 0;
            else
                dvdt = max((max_velocity - max_velocity_previous)/(time_step_change),0); % m/s2
            end
            flow_velocity(pos_save,1) = max_velocity;
            flow_acceleration(pos_save,1) =  dvdt; % m/s2
            if flag_D8 == 1
                factor_grid = sqrt(1/2);
            else
                factor_grid = 1;
            end
            if max_velocity == 0
                new_timestep = max_time_step;
            else
                new_timestep = (factor_grid*Resolution/max_velocity); % seconds
            end
            time_step_factor = max(alfa_max - slope_alfa*(max(max_velocity - v_threshold,0)),alfa_min);
            alfa_save(pos_save,1) = time_step_factor;
            new_timestep = new_timestep*alfa_save(pos_save,1);


        else             % % % % % % % %  Solution for the Stable Method - Time-step refreshment
            % Water Slopes Calculation (THERE IS A MISTAKE HERE)
            error('This method is currenly not working, please choose the courant method')
        end
        if flag_waterquality == 1
            % Adding Water Quality Minimum Time-step
            alfa_wq = 1; % Reduction factor of water quality
            new_timestep = min(alfa_wq*tmin_wq,new_timestep);
        else
            new_timestep = min(time_step_factor*new_timestep);
        end
        % Rounding time-step to min_timestep or max_timestep with the required
        % precision. This is not very interesting, maybe we should delete
        % it
        time_calculation_routing(k,1) = new_timestep;
        time_calculation_routing(k,1) = max(time_step_increments*floor(time_calculation_routing(k,1)/time_step_increments),min_time_step);
        time_calculation_routing(k,1) = min(time_calculation_routing(k,1),max_time_step);
        time_step = time_calculation_routing(k,1)/60; % new time-step for the next run
        if time_calculation_routing(k,1) == min_time_step % Check if we reached the minimum time-step
            unstable = 1; % If this is 1, it means we possibly had an unstable period at least
        end
    elseif  k == 1
        time_calculation_routing(k,1) = time_step*60; % sec
    else
        time_calculation_routing(k,1) = time_calculation_routing(k-1,1);
    end
    max_velocity_previous = max_velocity; % Assigning right velocities
    % Agregating Inflows to the New Time-step
    if flag_inflow > 0
        for z = 1:n_stream_gauges
            z1 = find(time_deltainflow(z,:) > t_previous,1,'first'); % begin of the time-step
            z2 = find(time_deltainflow(z,:) <= t,1,'last'); % end of the time-step
            if isempty(z1)
                z1 = 1;
            end
            if isempty(z2) || z2 < z1
                z2 = z1;
            end
            if time_step >= time_step_model
                delta_inflow_agg(z,1) = mean(delta_inflow(z,z1:z2))/(time_step_model*60)*time_step*60;
            else
                delta_inflow_agg(z,1) = delta_inflow(z,z1)/(time_step_model*60)*time_step*60;
            end
        end
    end
    % Agregating Precipitation to the New Time-step
    if flag_rainfall > 0
        if flag_spatial_rainfall ~= 1
            z1 = find(time_deltap > t_previous,1,'first'); % begin of the time-step
            z2 = find(time_deltap <= t,1,'last'); % end of the time-step
            if z2 < z1
                z2 = z1;
            end
            if time_step >= time_step_model
                delta_p_agg = mean(delta_p(1,z1:z2))/(time_step_model*60)*time_step*60;
            else
                delta_p_agg = delta_p(1,z1)/(time_step_model*60)*time_step*60;
            end
        else
            % Spatial Rainfall
            % Code for Spatial-Varying Rainfall
            z1 = find(rainfall_spatial_duration <= t_previous,1,'last');
            z2 = find(rainfall_spatial_duration <= t,1,'last');
            if z1 ~= z2 || z2 == length(rainfall_spatial_duration)
                x_coordinate = coordinates(1:n_raingauges,1); % Coordinates (easting) of each rain gauge
                y_coordinate = coordinates(1:n_raingauges,2); % Coordinates (northing) of each rain gauge
                x_grid = xulcorner + Resolution*[1:1:size(DEM_raster.Z,2)]'; % Pixel easting coordinates
                y_grid = yulcorner - Resolution*[1:1:size(DEM_raster.Z,1)]'; % Pixel northing coordinates

                % Spatial Rainfall
                if z2 == length(rainfall_spatial_duration)
                    spatial_rainfall = zeros(size(dem,1),size(dem,2));

                else
                    rainfall = rainfall_raingauges(z2,1:n_raingauges)'; % Values of rainfall at t for each rain gauge
                    [spatial_rainfall] = Rainfall_Interpolator(x_coordinate,y_coordinate,rainfall,x_grid,y_grid); % Interpolated Values
                end
                delta_p_agg = spatial_rainfall/3600*time_step_model*60; % Matrix of delta P for each pixel
                spatial_rainfall_maps(:,:,z2) = spatial_rainfall;
                average_spatial_rainfall(z2,1) = mean(spatial_rainfall(spatial_rainfall>=0));
            end
        end
    end


    % Inflows and Depth Refreshments
    d_t = d_t + qin_t*time_step/60;
    % Clearing stored values
    qout_left_t = zeros(ny_max,nx_max);
    qout_right_t = zeros(ny_max,nx_max);
    qout_up_t = zeros(ny_max,nx_max);
    qout_down_t = zeros(ny_max,nx_max);
    qout_ne_t = zeros(ny_max,nx_max);
    qout_se_t = zeros(ny_max,nx_max);
    qout_sw_t = zeros(ny_max,nx_max);
    qout_nw_t = zeros(ny_max,nx_max);
    qin_t = zeros(ny_max,nx_max);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Previous Depths
    d_p = d_t;
    I_p = I_t;


    % Saving Modeled Value
    obs_index = find(t >= time_observed_concentration,1,'last');  
    if isempty(obs_index)
        obs_index = 0;
    end
    if obs_index - obs_index_prev > 0
        modeled_concentration(1,obs_index) = Out_Conc;
    end
    obs_index_prev = obs_index;

    % Increase the counter
    t = time_step + t;
    k = k + 1;
end

Obj_function_value = sqrt(1/length(observed_concentration)*sum((observed_concentration' - modeled_concentration).^2));

% NSE
[NSE ~] = nashsutcliffe([time_observed_concentration',observed_concentration], [time_observed_concentration',modeled_concentration']);

% Obj_function_value = -NSE;

% R2
r2 = corrcoef(observed_concentration,modeled_concentration);

% Outputs
output_data = [time_observed_concentration',modeled_concentration',observed_concentration];
end