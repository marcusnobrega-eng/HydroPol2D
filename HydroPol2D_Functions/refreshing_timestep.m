% Refreshing the Time-step
% Goal: Calculate the new time-step of the model considering water quantity
% and water quality
% Developer: Marcus Nobrega, Ph.D.

% In case we are not in the first time-step and we are in a time that we
% need to change the time-step

if running_control.delta_time_save > 0 || k == 1 % First time-step
    if flags.flag_CA ~= 1 % If it is not inertial (CA model)
        flow_depth = depths.d_t; % Depth in which velocities will be calculated (mm)
        flow_depth(depths.d_t < CA_States.depth_tolerance) = 1e12; % Smaller depths won't alter velocities
        velocities.vel_left = (flow_rate.qout_left_t/1000/3600)*Wshed_Properties.Resolution^2./(Wshed_Properties.Resolution*flow_depth/1000); % m/s
        velocities.vel_right = (flow_rate.qout_right_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
        velocities.vel_up = (flow_rate.qout_up_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
        velocities.vel_down = (flow_rate.qout_down_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s

        if flags.flag_D8 == 1
            velocities.vel_ne = (flow_rate.qout_ne_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
            velocities.vel_se = (flow_rate.qout_se_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
            velocities.vel_sw = (flow_rate.qout_sw_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
            velocities.vel_nw = (flow_rate.qout_nw_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth/1000); % m/s
        end

    else % Inertial or diffusive or kinematic
        flow_depth = Hf*1000; % Depth in which velocities will be calculated (mm)
        flow_depth(flow_depth < CA_States.depth_tolerance) = 1e12; % Smaller depths won't alter velocities
        velocities.vel_right = (flow_rate.qout_right_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth(:,:,1)/1000); % m/s
        velocities.vel_left = -[zeros(ny,1), velocities.vel_right(:,2:end)];
        velocities.vel_up = (flow_rate.qout_up_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth(:,:,2)/1000); % m/s
        velocities.vel_down = -[velocities.vel_up(2:end,:); zeros(1,nx)];

        if flags.flag_D8 == 1
            velocities.vel_ne = (flow_rate.qout_ne_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth(:,:,6)/1000); % m/s
            velocities.vel_se = (flow_rate.qout_se_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth(:,:,7)/1000); % m/s
            velocities.vel_sw = (flow_rate.qout_sw_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth(:,:,8)/1000); % m/s
            velocities.vel_nw = (flow_rate.qout_nw_t/1000/3600)*Wshed_Properties.Resolution./(flow_depth(:,:,9)/1000); % m/s
        end

    end
    % Maximum Velocity
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

    % Velocity components
    if flags.flag_inertial == 1 || flags.flag_diffusive == 1 || flags.flag_kinematic == 1
        velocities.right_component = 1/2*(velocities.vel_right - velocities.vel_left);
        velocities.up_component = 1/2*(velocities.vel_up - velocities.vel_down);
    else
        velocities.right_component = velocities.vel_right - velocities.vel_left;
        velocities.up_component = velocities.vel_up - velocities.vel_down;
    end

    % D8 model
    if flags.flag_D8 == 1
        if flags.flag_CA == 1
            velocities.right_component = velocities.right_component + ...
                sqrt(2)*(1/2*(velocities.vel_ne + velocities.vel_sw)) + ...
                sqrt(2)*(1/2*(velocities.vel_se + velocities.vel_nw));
            velocities.up_component = velocities.up_component + ...
                sqrt(2)*(1/2*(velocities.vel_ne + velocities.vel_sw)) + ...
                -sqrt(2)*(1/2*(velocities.vel_se + velocities.vel_nw));
        else
            velocities.right_component = velocities.right_component + ...
                sqrt(2)*(velocities.vel_ne - velocities.vel_sw) + ...
                sqrt(2)*(velocities.vel_se - velocities.vel_nw);
            velocities.up_component = velocities.up_component + ...
                sqrt(2)*(velocities.vel_ne - velocities.vel_sw) + ...
                -sqrt(2)*(velocities.vel_se - velocities.vel_nw);
        end

    end
    velocities.total_velocity = sqrt((velocities.right_component).^2 + (velocities.up_component).^2);
    velocities.max_velocity = max(max(velocities.total_velocity));


    % Old max velocity
    if k > 1
        old_velocity = velocities.max_velocity;
    end

    % Courant Stability
    if flags.flag_D8 == 1
        Courant_Parameters.factor_grid = sqrt(1/2);
    else
        Courant_Parameters.factor_grid = 1;
    end
    if flags.flag_GPU == 1
        Courant_Parameters.factor_grid = gpuArray(Courant_Parameters.factor_grid);
    end

    if velocities.max_velocity > 0
        % Old Solution considering flow velocity
        new_timestep = (Courant_Parameters.factor_grid*Wshed_Properties.Resolution/velocities.max_velocity); % seconds
        Courant_Parameters.time_step_factor = max(Courant_Parameters.alfa_max - Courant_Parameters.slope_alfa*(max(velocities.max_velocity - Courant_Parameters.v_threshold,0)),Courant_Parameters.alfa_min);

        new_timestep = new_timestep*Courant_Parameters.time_step_factor;
        new_timestep = double(min(new_timestep,running_control.max_time_step));

        % Bates time-step
        wave_celerity = max(sqrt(9.81*max(depths.d_tot,depths.d_t)/1000),1e-10); % Using d_t, which is the depth at the end of the time-step
        wave_celerity(isinf(wave_celerity)) = nan;
        if flags.flag_adaptive_timestepping == 1
            new_timestep = (min(min((Courant_Parameters.alfa_min*Wshed_Properties.Resolution./(wave_celerity))))); % alpha of 0.4
        else
            new_timestep = (min(min((Courant_Parameters.alfa_min*Wshed_Properties.Resolution./(max_vel+wave_celerity))))); % alpha of 0.4
        end
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

    if flags.flag_waterquality == 1
        % Adding Water Quality Minimum Time-step
        Courant_Parameters.alfa_wq = 1; % Reduction factor of water quality
        new_timestep = min(Courant_Parameters.alfa_wq*tmin_wq,new_timestep);
    end

    % Rounding time-step to min_timestep or max_timestep with the required
    % precision. This is not very interesting, maybe we should delete
    % it
    running_control.time_calculation_routing = new_timestep;
    running_control.time_calculation_routing = max(running_control.time_step_increments*floor(running_control.time_calculation_routing/running_control.time_step_increments),running_control.min_time_step);
    running_control.time_calculation_routing = min(running_control.time_calculation_routing,running_control.max_time_step);
    time_step = running_control.time_calculation_routing/60; % new time-step for the next run
    if running_control.time_calculation_routing == running_control.min_time_step % Check if we reached the minimum time-step
        unstable = 1; % If this is 1, it means we possibly had an unstable period at least
    end
elseif  k == 1
    running_control.time_calculation_routing = time_step*60; % sec
else
    running_control.time_calculation_routing = time_step*60; % sec
end

velocities.max_velocity_previous = velocities.max_velocity; % Assigning right velocities