%%% ========================================================================
%  [HydroPol2D_Main_While]
%  Developer: Marcus Nobrega, Ph.D.
%  Date: [04/06/2025]
% -------------------------------------------------------------------------
%  Purpose:
%      [This function combines all internal model functions to run HydroPol2D.]
%
%  Inputs:
%      - input1: [You can enter a workspace with your model pre-processed]
%
%  Outputs:
%      - output1: [You can save your workspace with the results of your model]
%
%  Notes:
%      [To proper run the dashboard, you need Matlab 2021 or higher]
% ========================================================================

% clear all
% load workspace_CG_event.mat
% load workspace_beaver_creek_correct.mat
% clear all
% load workspace_amazon_simple.mat
% clear all
% load workspace_beaver.mat
% load example_1.mat
% load test.mat
% load test_radar.mat
% load workspace_synthetic_all_cases
% load workspace_CG_event_subgrid_functions.mat
% flags.flag_overbanks = 1;
% clear all; load pune_90m_ABM.mat;
% clear all; load pune_30m_ABM.mat
% clear all
% load pune_30m_ABM.mat
% load workspace_pune_100yr.mat
% load pune_inflow_90m.mat
% load workspace_pune_15000cms_inflow.mat;
% load pune_90m.mat;
% flags.flag_ETP = 0;
% clc
% flag.flag_ETP = 0;
% clear all; load workspace_wetbeavercreek.mat;

% clear all
% load v_tilted_workspace.mat
% load v_tilted_workspace_boundary.mat

clear all
load workspace_plane.mat
% flags.flag_subgrid = 0;
flags.flag_outlet_type = 2;
flags.flag_dashboard = 0;
running_control.report_every_percent = 0.25;
zero_matrix = zeros(ny,nx); zero_matrix(idx_nan) = nan;
tic
k = 1; % time-step counter
C = 0; % initial infiltration capacity
t = running_control.min_time_step; % inital time
progress_percent = (t / running_control.routing_time) * 100;
saver_count = 1; % starts in 1 but the next pointer should be 2, this is auto fixed when t reach the next time aggregation.
update_spatial_BC;
Flooded_Area = 0; % initial flooded area
velocities.velocity_raster = 0; % initial velocity
Risk_Area = 0; % initial risk area
store = 1; % Index for saving maps
t_previous = 0;
factor_time = 0;

running_control.max_time_step = min(running_control.record_time_hydrographs*60, running_control.max_time_step);
running_control.min_time_step = min(running_control.record_time_hydrographs*60, running_control.min_time_step);

max_dt = running_control.max_time_step;
ax.timer = minutes(double(gather(running_control.time_records(recording_parameters.actual_record_state)))) + date_begin;
ax.percentage = gather((t)/running_control.routing_time*100);

% %  DELETE ALL OF THIS
% time_step = running_control.min_time_step/60; % time-step of the model in min
% flags.flag_dashboard = 1; % DELETE
% flags.flag_subgrid = 0; % DELETE
% flags.flag_baseflow = 0; % DELETE
% flags.flag_critical = 0; % DELETE
% flags.flag_spatial_rainfall = 0; % DELETE
% flags.flag_rainfall = 1;
% flags.flag_ETP = 0; % DELETE
% SubgridTables = [];  % DELETE

% %  DELETE ALL OF THIS
% %  Assuming the average of soil depth
% avg_soil_depth = nanmean(nanmean(Soil_Properties.Soil_Depth));
% avg_soil_depth = 1;
% Soil_Properties.Soil_Depth = avg_soil_depth;
% Soil_Properties.ksat_gw = 10*Soil_Properties.ksat_gw;

% Activate this to use an average soil depth
% Soil_Properties.Soil_Depth(Soil_Properties.Soil_Depth>=0) = avg_soil_depth;

% BC_States.z0 = elevation - Soil_Properties.Soil_Depth; % Bedrock elevation [m]
% BC_States.z0(idx_nan) = nan;

% Initial Conditions for the domain
% h_depth = 0.8; % [m]
% BC_States.z0 = elevation - avg_soil_depth;
% BC_States.z0(14:16,14:16) = BC_States.z0(14:16,14:16) - 0.5;
% BC_States.h_0 = BC_States.z0 + h_depth;
% BC_States.h_t = BC_States.h_0;

% BC_States.z0(Wshed_Properties.idx_rivers) = elevation(Wshed_Properties.idx_rivers) - 0.05;
% BC_States.h_t(Wshed_Properties.idx_rivers) = BC_States.z0(Wshed_Properties.idx_rivers) + 0.05;
% BC_States.h_0 = BC_States.h_t;
%
% % Soil_Properties.ksat = 10000*Soil_Properties.ksat; % DELETE
% % Soil_Properties.ksat_gw = 10000*Soil_Properties.ksat_gw; % DELETE
% % Activate this to ensure very shallow depth in rivers
% % Soil_Properties.Soil_Depth(Wshed_Properties.idx_rivers) = 1; % m


% Close all figures
close all force;

% Find and delete all app windows (AppDesigner uifigures)
appWindows = findall(groot, 'Type', 'Figure', 'Tag', 'AppWindow');
delete(appWindows);

% Optionally clear memory/workspace
% clear all; clc;

try
    Subgrid_Properties;
catch
    Subgrid_Properties = [];
end
% Subgrid_Properties = [];
% 
roughness_squared = LULC_Properties.roughness.^2;
% 

% Initial System Storage
S_c = nansum(nansum(Wshed_Properties.Resolution^2.*Hydro_States.S/1000)); % Canopy

if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    S_p = nansum(nansum((Wshed_Properties.Resolution - Wshed_Properties.River_Width).*Wshed_Properties.Resolution.*max((depths.d_t/1000 - Wshed_Properties.River_Depth),0))) + ...
        nansum(nansum(Wshed_Properties.Resolution.*Wshed_Properties.River_Width.*depths.d_t/1000));
else
    S_p = nansum(nansum(Wshed_Properties.Resolution^2.*depths.d_t/1000));
end

S_UZ = nansum(nansum(Wshed_Properties.Resolution^2.*Soil_Properties.I_t/1000)); % UZ storage
S_GW = nansum(nansum(Wshed_Properties.Resolution^2.*Soil_Properties.Sy.*(BC_States.h_t - (Elevation_Properties.elevation_cell - Soil_Properties.Soil_Depth)))); % GW Storage
S_prev = S_c + S_p + S_UZ + S_GW;

% Initial Activation
if flags.flag_inertial == 1
    outflow_prev = outflow_bates;
end
catch_index = 1;


if flags.flag_obs_gauges ~= 1
    if isempty(extra_parameters.gauges)
        gauges = [];
    end
else
    extra_parameters = gauges;
    extra_parameters.labels_observed_string = [];
end


% if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= "WGS 84 / Pseudo-Mercator"
%     flags.flag_dashboard = 0; % Deactivating Dashboard
% end

if flags.flag_dashboard == 1
    ax.flags = flags;
    ax = HydroPol2D_running_dashboard(ax,Maps, zeros(size(DEM_raster.Z)), DEM_raster,gauges,BC_States,time_step,Wshed_Properties.Resolution,1,1,C_a);
end

%% ------------------------------------------------------------------------
% User-controlled reporting / console-print frequency
% This DOES NOT affect save_output_maps; it only affects screen printing
% -------------------------------------------------------------------------
if ~isfield(running_control,'report_every_percent') || isempty(running_control.report_every_percent)
    running_control.report_every_percent = 1;   % default = print every 1% of simulation
end

running_control.report_every_percent = max(running_control.report_every_percent, eps);
next_report_percent = running_control.report_every_percent;

%% ------------------------------------------------------------------------
% Cumulative mass-balance diagnostics
% errors = [Interception, Snow, ET/ETR, Infiltration, Groundwater] in m3
% -------------------------------------------------------------------------
cum_errors_m3 = zeros(1,5);   % cumulative signed module errors [m3]

module_names = { ...
    'Interception Module', ...
    'Snow Module', ...
    'Evaporation / Evapotranspiration Module', ...
    'Infiltration Module', ...
    'Groundwater Module'};


%% ------------------------------------------------------------------------
% Preallocate flow-rate work arrays once (reuse every iteration)
% -------------------------------------------------------------------------
if flags.flag_GPU == 1
    zero2D = gpuArray(zeros(ny,nx,'like',depths.d_t));
else
    zero2D = zeros(ny,nx,'like',depths.d_t);
end

% Outflows
flow_rate.qout_left_t  = zero2D;
flow_rate.qout_right_t = zero2D;
flow_rate.qout_up_t    = zero2D;
flow_rate.qout_down_t  = zero2D;
flow_rate.qin_t        = zero2D;


% Diagonal flows only if D8 is used
if flags.flag_D8 == 1
    flow_rate.qout_ne_t = zero2D;
    flow_rate.qout_se_t = zero2D;
    flow_rate.qout_sw_t = zero2D;
    flow_rate.qout_nw_t = zero2D;

end

%% ------------------------------------------------------------------------
% Preallocate mass-balance report history
% -------------------------------------------------------------------------
n_reports_est = ceil(100 / running_control.report_every_percent) + 10;

mass_balance_history.time_h        = zeros(n_reports_est,1);
mass_balance_history.errors_m3     = zeros(n_reports_est,5);
mass_balance_history.errors_mm_h   = zeros(n_reports_est,5);
mass_balance_history.cum_errors_m3 = zeros(n_reports_est,5);
mass_balance_history.cum_errors_mm = zeros(n_reports_est,5);
mass_balance_history.count         = 0;

%% #################### Main Loop (HydroPol2D)  ################ %

profile on
while t <= (running_control.routing_time + running_control.min_time_step/60) % Running up to the end of the simulation
    %tic
    try
        if t >= 89
            ttt = 1;
        end
        % -------------- Hydrological Model --------------- %
        % BC_States.delta_p_agg(BC_States.delta_p_agg >= 0) = 100* (time_step/60); % 1 mm/h
        % time_step = 1; % min
        
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
                
                    if flags.flag_subgrid == 1
                        % New standalone paper-style shared-face subgrid solver
                        [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t, ...
                         outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf, ...
                         Qc,Qf,Qci,Qfi,C_a,eta_t,V_t] = ...
                            Local_Inertial_Model_D4_Subgrid( ...
                            flags.flag_numerical_scheme, ...
                            Reservoir_Data.x_index,Reservoir_Data.y_index, ...
                            Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3, ...
                            Reservoir_Data.h2,Reservoir_Data.k4, ...
                            Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index, ...
                            Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index, ...
                            flags.flag_reservoir, ...
                            Elevation_Properties.elevation_cell, ...
                            depths.d_tot, depths.d_p, ...
                            LULC_Properties.roughness, roughness_squared, ...
                            Wshed_Properties.cell_area, time_step, Wshed_Properties.Resolution, ...
                            outlet_index, outlet_type, slope_outlet, ...
                            Wshed_Properties.row_outlet, Wshed_Properties.col_outlet, ...
                            CA_States.depth_tolerance, outflow_prev, idx_nan, ...
                            flags.flag_critical, Wshed_Properties.Inbank_Manning, Wshed_Properties.Overbank_Manning, ...
                            Wshed_Properties.River_Width, Wshed_Properties.River_Depth, ...
                            Qc, Qf, Qci, Qfi, C_a, ...
                            Subgrid_Properties, flags.flag_inflow, SubgridTables);
                
                    else
                        % Original baseline local inertial solver
                        [flow_rate.qout_left_t,flow_rate.qout_right_t,flow_rate.qout_up_t,flow_rate.qout_down_t,outlet_states.outlet_flow,depths.d_t,CA_States.I_tot_end_cell,outflow_bates,Hf,Qc,Qf,Qci,Qfi,C_a] = ...
                            Local_Inertial_Model_D4( ...
                            flags.flag_numerical_scheme, ...
                            Reservoir_Data.x_index,Reservoir_Data.y_index, ...
                            Reservoir_Data.k1,Reservoir_Data.h1,Reservoir_Data.k2,Reservoir_Data.k3, ...
                            Reservoir_Data.h2,Reservoir_Data.k4, ...
                            Reservoir_Data.y_ds1_index,Reservoir_Data.x_ds1_index, ...
                            Reservoir_Data.y_ds2_index,Reservoir_Data.x_ds2_index, ...
                            flags.flag_reservoir, Elevation_Properties.elevation_cell, ...
                            depths.d_tot, depths.d_p, ...
                            LULC_Properties.roughness, roughness_squared, ...
                            Wshed_Properties.cell_area, time_step, Wshed_Properties.Resolution, ...
                            outlet_index, outlet_type, slope_outlet, ...
                            Wshed_Properties.row_outlet, Wshed_Properties.col_outlet, ...
                            CA_States.depth_tolerance, outflow_prev, idx_nan, ...
                            flags.flag_critical, flags.flag_subgrid, ...
                            Wshed_Properties.Inbank_Manning, Wshed_Properties.Overbank_Manning, ...
                            Wshed_Properties.River_Width, Wshed_Properties.River_Depth, ...
                            Qc, Qf, Qci, Qfi, C_a, ...
                            Subgrid_Properties, flags.flag_overbanks, flags.flag_inflow, SubgridTables);
                    end
                
                end
            end
        end
        
%         zzz_time_tesk(k,1) = toc;
% 
%         if k == 1000
%             ttt = 1;
%             profile off;
%         end

        %% Outflows become inflows directly into qin_t (avoid extra full-size inflow arrays)
        flow_rate.qin_t(:) = 0;

        % Cardinal directions
        flow_rate.qin_t(:,2:nx)   = flow_rate.qin_t(:,2:nx)   + flow_rate.qout_right_t(:,1:nx-1);
        flow_rate.qin_t(:,1:nx-1) = flow_rate.qin_t(:,1:nx-1) + flow_rate.qout_left_t(:,2:nx);
        flow_rate.qin_t(2:ny,:)   = flow_rate.qin_t(2:ny,:)   + flow_rate.qout_down_t(1:ny-1,:);
        flow_rate.qin_t(1:ny-1,:) = flow_rate.qin_t(1:ny-1,:) + flow_rate.qout_up_t(2:ny,:);

        % Diagonal directions
        if flags.flag_D8 == 1
            flow_rate.qin_t(2:ny,1:nx-1)   = flow_rate.qin_t(2:ny,1:nx-1)   + flow_rate.qout_sw_t(1:ny-1,2:nx);
            flow_rate.qin_t(1:ny-1,1:nx-1) = flow_rate.qin_t(1:ny-1,1:nx-1) + flow_rate.qout_nw_t(2:ny,2:nx);
            flow_rate.qin_t(1:ny-1,2:nx)   = flow_rate.qin_t(1:ny-1,2:nx)   + flow_rate.qout_ne_t(2:ny,1:nx-1);
            flow_rate.qin_t(2:ny,2:nx)     = flow_rate.qin_t(2:ny,2:nx)     + flow_rate.qout_se_t(1:ny-1,1:nx-1);
        end

        % Keep only valid values
        %         flow_rate.qin_t(~isfinite(flow_rate.qin_t)) = 0;
        %         flow_rate.qin_t(idx3) = 0;

        % Including Source Terms (inflow hydrograph or water taken from the
        % domain)
        % flow_rate.qin_t = flow_rate.qin_t + BC_States.inflow/(time_step/60); % mm/h
        %         flow_rate.qin_t(~isfinite(flow_rate.qin_t)) = 0; % No flow at invalid cells
        %         flow_rate.qin_t(idx3) = 0; % No flow at these cells

        % Mass Balance Equation (outflow already taken)
        if flags.flag_CA == 1
            depths.d_t = depths.d_t + flow_rate.qin_t * time_step / 60;
        else
            % For inertial / diffusive / kinematic routing, mass balance is already
            % handled internally by the hydraulic solver.
        end

        % Water Balance Error
        water_balance_error_volume = abs(sum(sum(depths.d_t(depths.d_t<0))))*Wshed_Properties.Resolution^2*0.001; % m3
        water_balance_error_mm = water_balance_error_volume/Wshed_Properties.drainage_area*1000; % mm

%         if water_balance_error_volume > running_control.volume_error % We need to define better this parameter
%             factor_time = 1;
%             catch_index = catch_index + 5;
%             error('Mass balance error too high.')
%         else
%             catch_index = 1;
%             factor_time = factor_time + 1;
%             running_control.max_time_step = min(running_control.max_time_step*(1+0.05*factor_time),max_dt);
%         end
        % depths.d_t = max(depths.d_t,0);  % Taking away negative masses


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
        running_control.pos_save = ceil((t*60)/running_control.time_step_change);
        running_control.time_save = (running_control.pos_save - 1)*running_control.time_step_change/60;
        running_control.delta_time_save = running_control.time_save - running_control.time_save_previous;
        running_control.time_save_previous = running_control.time_save;
        running_control.actual_record_timestep = ceil((t*60)/running_control.time_step_change);

        % Refreshing time-step script
        refreshing_timestep;

        % -------------------------------------------------------------------------
        % Print diagnostics only every chosen % of simulation progress
        % This does NOT change saving logic; it only changes console printing
        % -------------------------------------------------------------------------
        do_report_now = false;
        if progress_percent >= next_report_percent || progress_percent >= 100
            do_report_now = true;
            next_report_percent = next_report_percent + running_control.report_every_percent;
        end

        %% Human Instability Calculations
        Human_Instability_Module

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

        % -------------------------------------------------------------------------
        % Simulation progress [%]
        % -------------------------------------------------------------------------
        progress_percent = (t / running_control.routing_time) * 100;

        %% Saving Output Maps

        save_output_maps;

        % Clearing stored values in-place (reuse memory)
        flow_rate.qout_left_t(:)  = 0;
        flow_rate.qout_right_t(:) = 0;
        flow_rate.qout_up_t(:)    = 0;
        flow_rate.qout_down_t(:)  = 0;
        flow_rate.qin_t(:)        = 0;

        if flags.flag_D8 == 1
            flow_rate.qout_ne_t(:) = 0;
            flow_rate.qout_se_t(:) = 0;
            flow_rate.qout_sw_t(:) = 0;
            flow_rate.qout_nw_t(:) = 0;
        end

        % Previous Depths and Moisture
        depths.d_p = depths.d_t;
        Soil_Properties.I_p = Soil_Properties.I_t;

        if flags.flag_inertial == 1
            % Saving previous outflows
            outflow_prev = outflow_bates; % Corrected previous outflow
        end

        % Saving Results in time_observations - Only Valid for Calibration
        save_automatic_cabralition_outputs;

        % Refreshing Time-step
        t = running_control.time_calculation_routing/60 + t;
        % time_step_save(k,2) = running_control.time_calculation_routing;
        % time_step_save(k,1) = t;
        k = k + 1;

        % -------------------------------------------------------------------------
        % Simulation progress [%] based on updated time
        % -------------------------------------------------------------------------
        progress_percent = (t / running_control.routing_time) * 100;

        % -------------------------------------------------------------------------
        % Print diagnostics only every chosen % of simulation progress
        % This does NOT change saving logic; it only changes console printing
        % -------------------------------------------------------------------------
        do_report_now = false;
        if progress_percent >= next_report_percent || progress_percent >= 100
            do_report_now = true;
            next_report_percent = next_report_percent + running_control.report_every_percent;
        end

        if do_report_now


            %% Mass Balance Check
            mass_balance_check

            % ------------------------------------------------------------------------
            % Convert module mass-balance errors to catchment-equivalent diagnostics
            % -------------------------------------------------------------------------
            % errors is assumed to be:
            % [Interception, Snow, ET/ETR, Infiltration, Groundwater] in m3 for current step

            dt_h = max(time_step/60, eps);   % current step duration in hours

            % Current time-step error converted to equivalent depth over basin [mm]
            errors_mm_step = (errors ./ Wshed_Properties.drainage_area) * 1000;

            % Current time-step error expressed as rate over basin [mm/h]
            errors_mm_h = errors_mm_step ./ dt_h;

            % Update cumulative errors [m3]
            cum_errors_m3 = cum_errors_m3 + errors;

            % Cumulative equivalent depth over basin [mm]
            cum_errors_mm = (cum_errors_m3 ./ Wshed_Properties.drainage_area) * 1000;

            % Totals
            total_error_m3_current    = sum(errors);
            total_error_mm_h_current  = sum(errors_mm_h);
            total_error_mm_cumulative = sum(cum_errors_mm);
            total_error_m3_cumulative = sum(cum_errors_m3);

            % Save history for later diagnostics
            mass_balance_history.count = mass_balance_history.count + 1;
            ih = mass_balance_history.count;

            mass_balance_history.time_h(ih,1)         = t;
            mass_balance_history.errors_m3(ih,:)      = errors;
            mass_balance_history.errors_mm_h(ih,:)    = errors_mm_h;
            mass_balance_history.cum_errors_m3(ih,:)  = cum_errors_m3;
            mass_balance_history.cum_errors_mm(ih,:)  = cum_errors_mm;

            %             % Plotting Data
            %             if flags.flag_dashboard == 0
            %                 HydroPol2D_running_dashboard_plot( ...
            %                     Maps, velocities.velocity_raster, DEM_raster, gauges, BC_States, ...
            %                     time_step, Wshed_Properties.Resolution, C_a, k);
            %             end

            % Determine output based on water quality flag
            if flags.flag_waterquality == 1
                perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3 = [
                    (t) / running_control.routing_time * 100, ... % Percent completed
                    (toc / ((t) / running_control.routing_time) - toc) / 3600, ... % Time remaining (hours)
                    time_step * 60, ... % Time step in seconds
                    max(max(depths.d_t(~isinf(depths.d_t)))), ... % Max depth [mm in your internal storage]
                    max(max(Hydro_States.f)), ... % Max Infiltration Rate
                    max(max(WQ_States.P_conc)), ... % Max Water Quality Pollutant concentration
                    volume_error];  % Volume error [m3]

                % Print formatted output for water quality
                fprintf('==== Water Quality Stats ====\n');
                fprintf('Percentage Complete: %.2f%%\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(1));
                fprintf('Time Remaining: %.2f hours\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(2));
                fprintf('Time Step Duration: %.0f seconds\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(3));
                fprintf('Max Depth: %.2f m\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(4));
                fprintf('Max Inf. Rate: %.2f mm/h\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(5));
                fprintf('Max Water Quality Concentration: %.2e Cmg/L\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(6));
                fprintf('Volume Error: %.3f m³\n', perc_duremain_tsec_dtmm_infmmhr_CmgL_dtmWQ_VolErrorm3(7));
            else
                perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3 = [
                    (t) / running_control.routing_time * 100, ... % Percentage complete
                    t / 60 / 24, ... % Time in days
                    (toc / ((t) / running_control.routing_time) - toc) / 3600, ... % Time remaining (hours)
                    time_step * 60, ... % Time step in seconds
                    1/1000 * max(max(depths.d_t(~isinf(depths.d_t)))), ... % Max depth [m]
                    max(max(Hydro_States.f)), ... % Max Inf Rate
                    velocities.max_velocity, ... % Max Velocity
                    volume_error];  % Volume error [m3]

                % Print formatted output for general model stats
                fprintf('---- General Model Stats ----\n');
                fprintf('Percentage Complete: %.2f%%\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(1));
                fprintf('Time Elapsed: %.2f days\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(2));
                fprintf('Time Remaining: %.2f hours\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(3));
                fprintf('Time Step Duration: %.0f seconds\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(4));
                fprintf('Max Depth: %.2f m\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(5));
                fprintf('Max Inf Rate: %.2f mm/h\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(6));
                fprintf('Max Velocity: %.2f m/s\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(7));
                fprintf('Volume Error: %.3f m³\n', perc_t_duremain_tsec_dtmm_infmmhr_CmgL_vel_VolErrorm3(8));
            end

            fprintf('\n==== Mass Balance Diagnostics Over Catchment ====\n');
            fprintf('Current time-step duration: %.4f h\n', dt_h);
            fprintf('Drainage area: %.4f m²\n', Wshed_Properties.drainage_area);
            fprintf('Negative-depth correction this step: %.4f mm\n\n', water_balance_error_mm);

            for ii = 1:numel(module_names)
                fprintf('%s\n', module_names{ii});
                fprintf('   Current error rate: %+12.4f mm/h\n', errors_mm_h(ii));
                fprintf('   Cumulative error:   %+12.4f mm\n',   cum_errors_mm(ii));
                fprintf('   Step error:         %+12.4f m³\n',   errors(ii));
                fprintf('   Cumulative error:   %+12.4f m³\n\n', cum_errors_m3(ii));
            end

            fprintf('TOTAL ACROSS MODULES\n');
            fprintf('   Current error rate: %+12.4f mm/h\n', total_error_mm_h_current);
            fprintf('   Cumulative error:   %+12.4f mm\n',   total_error_mm_cumulative);
            fprintf('   Step error:         %+12.4f m³\n',   total_error_m3_current);
            fprintf('   Cumulative error:   %+12.4f m³\n',   total_error_m3_cumulative);

            fprintf('\n');
            fprintf('===========================================\n');
            fprintf('End of Report at Time-Step %d\n', k);
            fprintf('===========================================\n\n');

        end

        % Water Quality Instability
        if tmin_wq < 0 || isnan(tmin_wq) || isinf(tmin_wq)
            error('Instability. in the Water Quality Model.')
        end

        

    catch ME % In case an error occurs in the model
        disp(ME.message)
        ME.stack.file
        ME.stack.line
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

%% Trim preallocated mass-balance history
if isfield(mass_balance_history,'count')
    nkeep = mass_balance_history.count;
    mass_balance_history.time_h        = mass_balance_history.time_h(1:nkeep,:);
    mass_balance_history.errors_m3     = mass_balance_history.errors_m3(1:nkeep,:);
    mass_balance_history.errors_mm_h   = mass_balance_history.errors_mm_h(1:nkeep,:);
    mass_balance_history.cum_errors_m3 = mass_balance_history.cum_errors_m3(1:nkeep,:);
    mass_balance_history.cum_errors_mm = mass_balance_history.cum_errors_mm(1:nkeep,:);
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
tempDir = fullfile(Paths.Temp);
save(fullfile(tempDir, ['save_map_hydro_' num2str(store) '.mat']), 'Maps', '-v7.3');


%% Returning Variables to CPU
if flags.flag_GPU == 1

    % ---- Gather structs (deep) ----
    BC_States            = gather_deep(BC_States);
    CA_States            = gather_deep(CA_States);
    Courant_Parameters   = gather_deep(Courant_Parameters);
    depths               = gather_deep(depths);
    Elevation_Properties = gather_deep(Elevation_Properties);
    GIS_data             = gather_deep(GIS_data);
    Hydro_States         = gather_deep(Hydro_States);
    Inflow_Parameters    = gather_deep(Inflow_Parameters);
    LULC_Properties      = gather_deep(LULC_Properties);
    Rainfall_Parameters  = gather_deep(Rainfall_Parameters);
    recording_parameters = gather_deep(recording_parameters);
    running_control      = gather_deep(running_control);
    Soil_Properties      = gather_deep(Soil_Properties);
    Wshed_Properties     = gather_deep(Wshed_Properties);
    outlet_states        = gather_deep(outlet_states);
    velocities           = gather_deep(velocities);
    gauges               = gather_deep(gauges);

    % IMPORTANT: gather flags LAST (or keep it as CPU-only always)
    flags                = gather_deep(flags);

    %     if flags.flag_obs_gauges == 1
    %         gauges = gather_deep(gauges);
    %         gauges.labels_observed_string = extra_parameters.gauges.labels_observed_string; % keep label text
    %     end

    if flags.flag_human_instability > 0
        Human_Instability = gather_deep(Human_Instability);
    end

    if flags.flag_waterquality == 1
        WQ_States = gather_deep(WQ_States);
    end

    if flags.flag_spatial_rainfall == 1
        Spatial_Rainfall_Parameters = gather_deep(Spatial_Rainfall_Parameters);
    end

    if flags.flag_ETP == 1
        ETP_Parameters = gather_deep(ETP_Parameters);
    end

    % ---- Gather "plain" arrays (also safe via gather_deep) ----
    elevation                = gather_deep(elevation);
    idx_nan                  = gather_deep(idx_nan);
    idx_outlet               = gather_deep(idx_outlet);
    k                        = gather_deep(k);
    nx                       = gather_deep(nx);
    outlet_index             = gather_deep(outlet_index);
    outlet_runoff_volume     = gather_deep(outlet_runoff_volume);
    outlet_type              = gather_deep(outlet_type);
    slope_outlet             = gather_deep(slope_outlet);
    t                        = gather_deep(t);
    t_previous               = gather_deep(t_previous);
    time_calculation_routing  = gather_deep(time_calculation_routing);
    time_step                = gather_deep(time_step);
    time_step_model          = gather_deep(time_step_model);
    tmin_wq                  = gather_deep(tmin_wq);
    C                        = gather_deep(C);
    min_soil_moisture        = gather_deep(min_soil_moisture);
    cumulative_infiltration  = gather_deep(cumulative_infiltration);
    max_GW_depth             = gather(max_GW_depth);
    cumulative_recharge      = gather(cumulative_recharge);

    if flags.flag_waterquality == 1
        idx_nan_5 = gather_deep(idx_nan_5);
        Out_Conc  = gather_deep(Out_Conc);
    end
end


leftovers = find_gpu_leftovers_in_workspace();
% disp(leftovers);

compare_flows;
post_processing
%% Save Workspace for post-processing
try
    save('modeled_results.mat', '-v7.3')
end

function x = gather_deep(x)
%GATHER_DEEP  Recursively gather gpuArray content inside structs/cells/arrays.
% Safe for strings/chars/logicals/datetimes/function handles/etc.

% 1) If it's a gpuArray, gather it immediately
if isa(x, 'gpuArray')
    x = gather(x);
    return
end

% 2) Struct: recurse through all fields (handles nested structs too)
if isstruct(x)
    fn = fieldnames(x);
    for ii = 1:numel(x)           % handle struct arrays
        for k = 1:numel(fn)
            f = fn{k};
            x(ii).(f) = gather_deep(x(ii).(f));
        end
    end
    return
end

% 3) Cell: recurse each cell element
if iscell(x)
    for i = 1:numel(x)
        x{i} = gather_deep(x{i});
    end
    return
end

% 4) Tables/timetables: recurse variable-wise (optional but safe)
if istable(x)
    vnames = x.Properties.VariableNames;
    for k = 1:numel(vnames)
        x.(vnames{k}) = gather_deep(x.(vnames{k}));
    end
    return
end

% 5) Everything else: leave as-is (char/string/logical/datetime/function_handle/etc.)
end

function leftovers = find_gpu_leftovers_in_workspace()
w = evalin('caller','whos');
leftovers = {};
for i = 1:numel(w)
    name = w(i).name;
    val  = evalin('caller', name);
    if contains_gpu(val)
        leftovers{end+1} = name; %#ok<AGROW>
    end
end
end

function tf = contains_gpu(x)
if isa(x,'gpuArray')
    tf = true; return
elseif isstruct(x)
    tf = false;
    fn = fieldnames(x);
    for ii = 1:numel(x)
        for k = 1:numel(fn)
            if contains_gpu(x(ii).(fn{k})), tf = true; return; end
        end
    end
elseif iscell(x)
    tf = false;
    for i = 1:numel(x)
        if contains_gpu(x{i}), tf = true; return; end
    end
elseif istable(x)
    tf = false;
    vn = x.Properties.VariableNames;
    for k = 1:numel(vn)
        if contains_gpu(x.(vn{k})), tf = true; return; end
    end
else
    tf = false;
end
end
