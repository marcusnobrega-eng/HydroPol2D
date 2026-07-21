function [B_t,P_conc,Out_Conc,tmin_wq,tot_W_out,mass_lost,Tot_Washed] = build_up_wash_off(C_3,C_4,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,B_t,time_step,nx_max,ny_max,cell_area,outlet_index,idx_nan_5,flag_wq_model,mass_lost,Tot_Washed,Bmin,Bmax,min_Bt)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 1 July, 2022                      %
%          Goal: Solve the Build-up and Wash-off Model            %

%% Pollutant Rate Flux
% --------------- Vector of Flows  % ---------------
% 3 dim array with pages 1, 2, 3, 4, and 5 being left, right, up, down, and
% outlet flow rates in mm/hr
q_out_t = cat(3,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow); % mm/h
large_timestep = time_step*60; % seconds
idx_bt_1 = ((B_t*1000/cell_area < min_Bt) + isnan(B_t) + isinf(B_t)) > 0; % Values below Min B_t
idx_bt = repmat(idx_bt_1,1,1,5);
B_t_extra = B_t;
B_begin = B_t;
Br_kg = Bmin/1000*cell_area; % kg
Bm_kg = Bmax/1000*cell_area; % kg
outlet_mask = logical(outlet_index) & (outlet_flow > 0);
if any(outlet_mask(:))
    outlet_water_volume_L = max(nansum(outlet_flow(outlet_mask)) * cell_area * time_step/60, 0);
else
    outlet_water_volume_L = 0;
end

% --------------- Choosing Which Wash-off Equation % ---------------
if flag_wq_model == 1 % Rating Curve
    B_t_extra(B_t*1000/cell_area <= Bmin) = 0; % 100 g/m2
    B_t_extra(B_t*1000/cell_area >= Bmax) = Bm_kg; % 100 g/m2
    f_bt = (1 + (max(B_t_extra - Br_kg,0)));
    W_out_t = C_3.*(q_out_t/1000/3600*cell_area).^(C_4).*f_bt; % kg/hr (rating curve)
else
    %     W_out_t = C_3.*q_out_t.^(C_4).*B_t; % kg/hr (mass based curve)
    W_out_t = C_3.*(q_out_t/1000/3600*cell_area).^(C_4).*B_t; % kg/hr (mass based curve in terms of velocity)
end
W_out_t(idx_bt) = 0; % No flux if mass is too low (below 1 g/m2)
% --------------- Total Outflows per Cell  % ---------------
tot_W_out = sum(W_out_t,3); % kg/hr
tot_q_out = sum(q_out_t,3); % mm/hr
small_number = 1e-16;

% --------------- Replacing Nans to 0  % ---------------%
tot_W_out(idx_nan_5(:,:,1)) = 0;
tot_q_out(idx_nan_5(:,:,1)) = 0;
W_out_t(idx_nan_5) = 0;
% --------------- Assigning Washoff Rates for Each Direction   % ---------------%
W_out_left_t = W_out_t(:,:,1);
W_out_right_t = W_out_t(:,:,2);
W_out_up_t = W_out_t(:,:,3);
W_out_down_t = W_out_t(:,:,4);
W_out_outlet_t = W_out_t(:,:,5);
% ---------------% Creating Pollutant Inflow Matrices  % ---------------%
W_in_left_t = [zeros(ny_max,1),W_out_right_t(:,1:(nx_max-1))]; % left
W_in_right_t = [W_out_left_t(:,(2:(nx_max))) zeros(ny_max,1)]; % right
W_in_up_t = [zeros(1,nx_max) ; W_out_down_t(1:(ny_max-1),:)]; % up
W_in_down_t = [W_out_up_t(2:(ny_max),:) ; zeros(1,nx_max)]; % down
W_in_t = cat(3,W_in_left_t,W_in_right_t,W_in_up_t,W_in_down_t); % Total Inflow
% --------------- Net Pollutant Rate % ---------------%
% At the same time, pollutants are entering and leaving the cells. The
% matrix dW measures the difference between outflows and inflows
dW = (tot_W_out - sum(W_in_t,3)); % Outflow_pol - Inflow_pol (kg/hr)
%% Minimum stable time step for the explicit pollutant update
% dW is outflow minus inflow. Only a positive dW can deplete a cell; a
% receiving cell with zero mass must not force a zero internal time step.
tmin_wq = pollutant_stable_timestep(B_t, dW, min_Bt, cell_area, large_timestep);
%% Overall Mass Balance at Cells
% B_t = B_t + dW*time_step/60; % Refreshing Pollutant Mass
B_t_mid = B_t - dW*time_step/60; % Refreshing Pollutant Mass (1/2 of the time-step)
if min(min(B_t_mid)) < 0 || tmin_wq < time_step*60 % Break time-step internally
    n_steps_decimal = time_step*60/(tmin_wq); % decimal time-step
    n_steps_floor = floor(n_steps_decimal); % integer lowest value
    steps  = (n_steps_floor+1);
    outflow_mass = 0; % Starting to measure outflow pollutant mass (kg)
    outlet_mass = 0; % kg captured by the outlet sink during the time-step
    for i = 1:(steps)
        % Break time step internally. Flow rates won't change but B_t will
        if i == (steps)
            dt = (n_steps_decimal - n_steps_floor)*tmin_wq; % sec
        else
            dt = time_step*60/(steps); % sec
        end
        % New Pollutant Flux Rates
        if flag_wq_model == 1 % Rating Curve
            B_t_extra(B_t*1000/cell_area <= Bmin) = 0; %
            B_t_extra(B_t*1000/cell_area >= Bmax) = Bm_kg; %
            f_bt = (1 + (max(B_t_extra - Br_kg,0)));
            W_out_t = C_3.*(q_out_t/1000/3600*cell_area).^(C_4).*f_bt; % kg/hr (rating curve)
        else
            W_out_t = C_3.*(q_out_t/1000/3600*cell_area).^(C_4).*B_t; % kg/hr (mass based curve in terms of velocity)
        end
        idx_bt_1 = ((B_t*1000/cell_area < min_Bt) + isnan(B_t) + isinf(B_t)) > 0; % Values below 1 g / m2 or inf or nan
        idx_bt = repmat(idx_bt_1,1,1,5);
        % --------------- Total Outflows per Cell Valid for This time-step  % ---------------
        tot_W_out = sum(W_out_t,3); % kg/hr
        tot_W_out(idx_bt_1) = 0; % No flux if mass is too low (below 1 g/m2)
        W_out_t(idx_bt) = 0; % No flux if mass is too low (below 1 g/m2)
        % --------------- Assigning Washoff Rates for Each Direction   % ---------------%
        W_out_left_t = W_out_t(:,:,1);
        W_out_right_t = W_out_t(:,:,2);
        W_out_up_t = W_out_t(:,:,3);
        W_out_down_t = W_out_t(:,:,4);
        W_out_outlet_t = W_out_t(:,:,5);
        % ---------------% Creating Pollutant Inflow Matrices  % ---------------%
        W_in_left_t = [zeros(ny_max,1),W_out_right_t(:,1:(nx_max-1))]; % left
        W_in_right_t = [W_out_left_t(:,(2:(nx_max))) zeros(ny_max,1)]; % right
        W_in_up_t = [zeros(1,nx_max) ; W_out_down_t(1:(ny_max-1),:)]; % up
        W_in_down_t = [W_out_up_t(2:(ny_max),:) ; zeros(1,nx_max)]; % down
        W_in_t = cat(3,W_in_left_t,W_in_right_t,W_in_up_t,W_in_down_t); % outlet
        % --------------- Net Pollutant Rate % ---------------%
        % At the same time, pollutants are entering and leaving the cells. The
        % matrix dW measures the difference between inflows and outflows
        dW = (tot_W_out - sum(W_in_t,3)); % Inflow pol - Outflow pol (kg/hr)
        outflow_mass = outflow_mass + tot_W_out*dt/3600; % kg of pollutant that left the cell
        % dW = max((sum(W_in_t,3) - tot_W_out),-B_t/(time_step/60) + small_number); % Inflow pol - Outflow pol (kg/hr)
        B_t = B_t - dW*dt/3600; % Mass Balance for the incremental time-step
        if any(outlet_mask(:))
            % Outlet washoff is already included in dW. Record that face
            % flux only; clearing B_t here would remove the remaining cell
            % store a second time.
            outlet_capture = W_out_outlet_t(outlet_mask) * dt / 3600;
            outlet_capture = outlet_capture(isfinite(outlet_capture) & outlet_capture > 0);
            outlet_mass = outlet_mass + sum(outlet_capture);
        end
    end
    % Average Pollutant Flux in the Time-Step
    tot_W_out = outflow_mass/(time_step/60); % kg/hr
    % Average Pollutant Flux a the outlet in the Time-step
    W_out_outlet_t = outlet_mass/(time_step/60); % kg/hr
    % Rounding Negative Values
    mass_lost = sum(sum(B_t(B_t<0))) + mass_lost;
    B_t(B_t<0) = 0;
else
    B_t = B_t_mid; % We can use the same time-step from the hydrodynamic model
    if any(outlet_mask(:))
        % The outlet-face flux is already included in dW. Do not clear the
        % residual cell store after the conservative update.
        outlet_capture = W_out_outlet_t(outlet_mask) * time_step / 60;
        outlet_capture = outlet_capture(isfinite(outlet_capture) & outlet_capture > 0);
        outlet_mass = sum(outlet_capture);
    else
        outlet_mass = 0;
    end
    % Rounding Negative Values
    mass_lost = sum(sum(B_t(B_t<0))) + mass_lost;
    B_t(B_t<0) = 0;
    % No need to calculate any average since we used the whole time-step
end
%% Minimum stable time step after the update
tmin_wq = pollutant_stable_timestep(B_t, dW, min_Bt, cell_area, large_timestep);
% ---------------% Imposing Constraint at B_t  % ---------------% %
%%% This constratint rounds any inf or nan at B_t to 0
% Since we are changing time-steps each 60 seconds (usually), we could
% experience instability within this period. To this end, we force any
% weird number to 0
% Another important consideration would be round very tiny pollutant
% massess to zero. So let's assume we have a minimum value in kg/m2
% If we have a minimum of 0.01 kg/ha, we would have 10^-6 kg/m2 as a
% minimum bound. Therefore, if we round B_t to the 6th decimal place, we
% sole this problem.

B_t(abs(B_t) < 1e-12) = 0;
% B_t(B_t*1000/cell_area < min_Bt) = 0; % smaller than 0.0001 g/m2
% Adding a constraint in q_out to avoid huge numbers
% ---------------% Imposing Constraint at Flow Rate  % ---------------% %
tot_q_out = max(tot_q_out,10); % 10 mm/hr
% ---------------% Pollutant Concentration % ---------------% %
P_conc = max(10^6*(tot_W_out)./(tot_q_out*cell_area),0); % mg/L
% ---------------% Outlet Concentration % ---------------% %
if outlet_water_volume_L > 0
    Out_Conc = max(1e6*outlet_mass/outlet_water_volume_L,0); % mg/L
else
    Out_Conc = 0;
end
% Tot_Washed = (max(B_begin - B_t,0)) + Tot_Washed; % Actually, total received
Tot_Washed = tot_W_out*time_step/60 + Tot_Washed; % (kg) Actually, total received

end

function tmin_wq = pollutant_stable_timestep(B_t, dW, min_Bt, cell_area, fallback_s)
% Return the depletion-limited explicit step in seconds. Cells with no
% available pollutant or net inflow cannot become negative in this update.
mass_floor_kg = max(min_Bt / 1000 .* cell_area, 1e-15);
eligible = isfinite(B_t) & isfinite(dW) & B_t > mass_floor_kg & dW > 0;
if ~any(eligible, 'all')
    tmin_wq = fallback_s;
    return
end

tmin_wq = 3600 * min(B_t(eligible) ./ dW(eligible));
if ~isfinite(tmin_wq) || tmin_wq <= 0
    tmin_wq = fallback_s;
end
end

