%% Main HydroPol2D While Loop
% Developer: Marcus Nobrega, Ph.D.
% Date 3/6/2025
% Goal - Run the main modeling process of the model

clear all
load workspace_moab_gpu.mat

tic
k = 1; % time-step counter
C = 0; % initial infiltration capacity
t = running_control.min_time_step; % inital time
time_step = running_control.min_time_step/60; % time-step of the model in min
saver_count = 1; % starts in 1 but the next pointer should be 2, this is auto fixed when t reach the next time aggregation.
update_spatial_BC;
Flooded_Area = 0; % initial flooded area
velocities.velocity_raster = 0; % initial velocity
Risk_Area = 0; % initial risk area
store = 1; % Index for saving maps
t_previous = 0;
factor_time = 0;
running_control.max_time_step = 60; % seconds
max_dt = running_control.max_time_step;

% Initial System Storage
S_c = nansum(nansum(C_a.*Hydro_States.S/1000));
S_p = nansum(nansum(Wshed_Properties.Resolution.*Wshed_Properties.River_Width.*depths.d_t/1000));
S_UZ = nansum(nansum(C_a.*Soil_Properties.I_t/1000));
S_GW = nansum(nansum(C_a.*Soil_Properties.Sy.*(BC_States.h_t - (Elevation_Properties.elevation_cell - Soil_Properties.Soil_Depth))));
S_prev = S_c + S_p + S_UZ + S_GW;

% Initial Activation
if flags.flag_inertial == 1
    outflow_prev = outflow_bates;
end
catch_index = 1;


if flags.flag_obs_gauges ~= 1
    extra_parameters.gauges = [];
end

if flags.flag_dashboard == 1
    ax.flags = flags;
    ax = HydroPol2D_running_dashboard(ax,Maps, zeros(size(DEM_raster.Z)), DEM_raster,extra_parameters.gauges,BC_States,time_step,Wshed_Properties.Resolution,1,1,C_a);
end

% #################### Main Loop (HydroPol2D)  ################ %
while t <= (running_control.routing_time + running_control.min_time_step/60) % Running up to the end of the simulation
    try
        % -------------- Hydrological Model --------------- %
        Hydrological_Model; % Runs the interception + infiltration + GW routing model

        % Preallocating cels for Cellular Automata 
        if flags.flag_D8 == 1 && flags.flag_diffusive == 1
            %%%% WATER DEPTHS %%%
            depths.d_left_cell = [zeros(ny,1),depths.d_tot(:,1:(nx-1))]; depths.d_right_cell = [depths.d_tot(:,(2:(nx))) zeros(ny,1)]; depths.d_up_cell = [zeros(1,nx) ; depths.d_tot(1:(ny-1),:)]; depths.d_down_cell = [depths.d_tot(2:(ny),:) ; zeros(1,nx)];
            if flags.flag_D8 == 1
                % Simulate with D-8 flow direction
                % --- Adding NE, SE, SW, NW --- %
                depths.d_NE_t(2:(ny),1:(nx-1)) = depths.d_tot(1:(ny-1),2:nx); depths.d_SE_t(1:(ny-1),1:(nx-1)) = depths.d_tot(2:ny,2:nx); depths.d_SW_t(1:(ny-1),2:(nx)) = depths.d_tot(2:(ny),1:(nx-1)); depths.d_NW_t(2:ny,2:nx) = depths.d_tot(1:(ny-1),1:(nx-1)); % OK
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
                % Still has to code the 8D version of the local inertial
                % formulation
                [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,flow_rate.qout_ne_t,flow_rate.qout_se_t,flow_rate.qout_sw_t,flow_rate.qout_nw_t,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf] = ...
                    Bates_Inertial_8D(Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                    flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                    depths.d_tot, depths.d_p,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,Wshed_Properties.Resolution,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,CA_States.depth_tolerance,outflow_prev,idx_nan,flags.flag_critical);
            end
        else % 4-D
            % CA
            if flags.flag_inertial ~= 1 && flags.flag_CA == 1
                [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell] = ...
                    CA_Routing(Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                    flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                    depths.d_tot,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,LULC_Properties.h_0,Wshed_Properties.Resolution,CA_States.I_tot_end_cell,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,idx_nan,flags.flag_critical);
            else
                if flags.flag_diffusive == 1
                % --------------------- Diffusive Wave Formulation -----%
                [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
                    Diffusive_Wave_Model_Implicit(flags.flag_numerical_scheme,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                    flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                    depths.d_tot, depths.d_p,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,Wshed_Properties.Resolution,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,CA_States.depth_tolerance,outflow_prev,idx_nan,flags.flag_critical,flags.flag_subgrid,Wshed_Properties.Inbank_Manning,Wshed_Properties.Overbank_Manning,Wshed_Properties.River_Width, Wshed_Properties.River_Depth,Qc,Qf,Qci,Qfi,C_a);                
                elseif flags.flag_kinematic == 1
                % --------------------- Kinematic Wave Formulation -----%
                [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
                    Kinematic_Wave_Model(flags.flag_numerical_scheme,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                    flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                    depths.d_tot, depths.d_p,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,Wshed_Properties.Resolution,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,CA_States.depth_tolerance,outflow_prev,idx_nan,flags.flag_critical,flags.flag_subgrid,Wshed_Properties.Inbank_Manning,Wshed_Properties.Overbank_Manning,Wshed_Properties.River_Width, Wshed_Properties.River_Depth,Qc,Qf,Qci,Qfi,C_a);                
                else
                % -------------------- Local Inertial Formulation ----------------%
                [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
                    Local_Inertial_Model_D4(flags.flag_numerical_scheme,Reservoir_Data.x_index,Reservoir_Data.y_index,Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3,Reservoir_Data.h2,Reservoir_Data.k4,Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index,Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index,...
                    flags.flag_reservoir,Elevation_Properties.elevation_cell,...
                    depths.d_tot, depths.d_p,LULC_Properties.roughness,Wshed_Properties.cell_area,time_step,Wshed_Properties.Resolution,outlet_index,outlet_type,slope_outlet,Wshed_Properties.row_outlet,Wshed_Properties.col_outlet,CA_States.depth_tolerance,outflow_prev,idx_nan,flags.flag_critical,flags.flag_subgrid,Wshed_Properties.Inbank_Manning,Wshed_Properties.Overbank_Manning,Wshed_Properties.River_Width, Wshed_Properties.River_Depth,Qc,Qf,Qci,Qfi,C_a);
                end
            end
        end
        
        %% Outflows become Inflows for CA
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
            flow_rate.qin_ne_t = zero_matrix; flow_rate.qin_se_t = zero_matrix; flow_rate.qin_sw_t = zero_matrix; flow_rate.qin_nw_t = zero_matrix;  flow_rate.qin_ne_t(2:(ny),1:(nx-1)) = flow_rate.qout_sw_t(1:(ny-1),2:(nx)); flow_rate.qin_se_t(1:(ny-1),1:(nx-1)) = flow_rate.qout_nw_t(2:ny,2:nx); flow_rate.qin_sw_t(1:(ny-1),2:(nx)) = flow_rate.qout_ne_t(2:(ny),1:(nx-1)); flow_rate.qin_nw_t(2:ny,2:nx) = flow_rate.qout_se_t(1:(ny-1),1:(nx-1)); % OK
        end

        if flags.flag_D8 == 1
            flow_rate.qin_t = flow_rate.qin_left_t + flow_rate.qin_right_t + flow_rate.qin_up_t + flow_rate.qin_down_t ...
                + flow_rate.qin_ne_t + flow_rate.qin_se_t + flow_rate.qin_sw_t + flow_rate.qin_nw_t;
        else
            flow_rate.qin_t = flow_rate.qin_left_t + flow_rate.qin_right_t + flow_rate.qin_up_t + flow_rate.qin_down_t;
        end

        % Including Source Terms (inflow hydrograph or water taken from the
        % domain)
        % flow_rate.qin_t = flow_rate.qin_t + BC_States.inflow/(time_step/60); % mm/h
        idx3 = logical(isnan(flow_rate.qin_t) + isinf(flow_rate.qin_t)); % Mask to take out non cells
        flow_rate.qin_t(idx3) = 0; % No flow at these cells

        % Mass Balance Equation (outflow already taken)
        if flags.flag_inertial == 1 || flags.flag_diffusive == 1 || flags.flag_kinematic == 1
            depths.d_t = depths.d_t + 0*flow_rate.qin_t*time_step/60; % All balance is inside the model
        elseif flags.flag_CA == 1
            depths.d_t = depths.d_t + flow_rate.qin_t*time_step/60;
        end

        % Water Balance Error
        water_balance_error_volume = abs(sum(sum(depths.d_t(depths.d_t<0))))*Wshed_Properties.Resolution^2*0.001; % m3
        water_balance_error_mm = water_balance_error_volume/Wshed_Properties.drainage_area*1000; % mm

        if water_balance_error_volume > running_control.volume_error % We need to define better this parameter
            factor_time = 1;
            catch_index = catch_index + 5;
            error('Mass balance error too high.')
        else
            catch_index = 1;
            factor_time = factor_time + 1;
            running_control.max_time_step = min(running_control.max_time_step*(1+0.05*factor_time),max_dt);
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

        %% Human Instability Calculations
        Human_Instability_Module

        %% Refreshing Time-step
        running_control.pos_save = ceil((t*60)/running_control.time_step_change);
        running_control.time_save = (running_control.pos_save - 1)*running_control.time_step_change/60;
        running_control.delta_time_save = running_control.time_save - running_control.time_save_previous;
        running_control.time_save_previous = running_control.time_save;
        running_control.actual_record_timestep = ceil((t*60)/running_control.time_step_change);


        %% Mass Balance Check
        mass_balance_check

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
            t_previous = running_control.time_calculation_routing/60;
        else
            t_previous = t;
        end
        % Current time
        t_save = t + running_control.time_calculation_routing/60;

        %% Saving Output Maps
        save_output_maps;

        % Clearing stored values
        if flags.flag_GPU == 1
            flow_rate.qout_left_t = gpuArray(zeros(ny,nx)); flow_rate.qout_right_t = gpuArray(zeros(ny,nx)); flow_rate.qout_up_t = gpuArray( zeros(ny,nx)); flow_rate.qout_down_t = gpuArray(zeros(ny,nx)); flow_rate.qout_ne_t = gpuArray(zeros(ny,nx)); flow_rate.qout_se_t = gpuArray(zeros(ny,nx)); flow_rate.qout_sw_t = gpuArray(zeros(ny,nx)); flow_rate.qout_nw_t = gpuArray(zeros(ny,nx)); flow_rate.qin_t = gpuArray(zeros(ny,nx));
        else
            flow_rate.qout_right_t = zeros(ny,nx); flow_rate.qout_up_t = zeros(ny,nx);flow_rate.qout_down_t = zeros(ny,nx);flow_rate.qout_ne_t = zeros(ny,nx);flow_rate.qout_se_t = zeros(ny,nx);flow_rate.qout_sw_t = zeros(ny,nx);flow_rate.qout_nw_t = zeros(ny,nx);flow_rate.qin_t = zeros(ny,nx);
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
        t = running_control.time_calculation_routing/60 + t;
        % time_step_save(k,2) = running_control.time_calculation_routing;
        % time_step_save(k,1) = t;
        k = k + 1;

        % Show Stats
        if flags.flag_waterquality == 1
            perc__duremain___tsec____dtmm___infmmhr__CmgL___dtmWQ_VolErrorm3 = [(t)/running_control.routing_time*100, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(Hydro_States.f)), max(max((WQ_States.P_conc))), tmin_wq,volume_error]
        else
            perc_t__duremain___tsec____dtmm___infmmhr__CmgL___vel__VolErrorm3 = [(t)/running_control.routing_time*100, t/60/24, (toc/((t)/running_control.routing_time) - toc)/3600,time_step*60,max(max(depths.d_t(~isinf(depths.d_t)))),max(max(Hydro_States.f)), max(max(velocities.velocity_raster)),volume_error]
        end

        % Water Quality Instability
        if tmin_wq < 0 || isnan(tmin_wq) || isinf(tmin_wq)
            error('Instability. in the Water Quality Model.')
        end

    catch ME % In case an error occurs in the model
        disp(ME.message)
        t = t - time_step; % Returning the time-step
        t_previous = t;
        wave_celerity = sqrt(9.81*(max(max(max(depths.d_tot/1000)),max(max(depths.d_p/1000))))); % Using d_p and d_tot
        max_vel = max(max(velocities.velocity_raster));
        catch_index = catch_index + 1;
        factor = 1/catch_index/running_control.factor_reduction;
        if flags.flag_adaptive_timestepping == 1
            new_timestep = factor*(min(Courant_Parameters.alfa_min*Wshed_Properties.Resolution./(wave_celerity))); % alpha of 0.4
        else
            new_timestep = factor*(min(Courant_Parameters.alfa_min*Wshed_Properties.Resolution./(max_vel+wave_celerity))); % alpha of 0.4
        end
        % dt_water_balance = min(min(depths.d_p/1000*Wshed_Properties.cell_area./(CA_States.I_tot_end_cell/(time_step*60)))); % sec
        dt_water_balance = new_timestep;
        new_timestep = min(new_timestep,dt_water_balance);
        if catch_index > 1
            running_control.max_time_step = min(new_timestep*1.5,running_control.max_time_step); % Assuming a smaller max time-step to avoid large integrations
        end        
        new_timestep = min(new_timestep,running_control.max_time_step); % sec
        if isinf(new_timestep)
           new_timestep = running_control.min_time_step;
        end
        time_step = new_timestep/60; % min
        t = t + time_step;
        depths.d_t = depths.d_p;
        Soil_Properties.I_t = Soil_Properties.I_p;
        if flags.flag_baseflow == 1
            BC_States.h_t = BC_States.h_0;
        end
        % current_storage = previous_storage;
        update_spatial_BC
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
        gauges.labels_observed_string = extra_parameters.gauges.labels_observed_string;
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


%% Save Workspace for post-processing
save('modeled_results')
