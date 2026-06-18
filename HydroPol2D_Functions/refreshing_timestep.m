% Refreshing the Time-step
% Goal: Calculate the new time-step of the model considering water quantity
% and water quality
% Developer: Marcus Nobrega, Ph.D.
%
% Updated revision:
% - Fixes depth_timestep inconsistency
% - Uses velocity magnitudes for CFL control
% - Fixes velocity_raster logic
% - Fixes undefined max_vel
% - Separates dt_adv and dt_wave, then takes min(...)
% - Corrects D8 projection factor to 1/sqrt(2)
% - Removes dead negative-max-velocity logic

% In case we are not in the first time-step and we are in a time that we
% need to change the time-step
if running_control.delta_time_save > 0 || k == 1

    %% --------------------------------------------------------------------
    % Store previous maximum velocity BEFORE recomputing current one
    % ---------------------------------------------------------------------
    if isfield(velocities, 'max_velocity')
        old_velocity = velocities.max_velocity;
    else
        old_velocity = 0;
    end

    %% --------------------------------------------------------------------
    % Compute face velocities and depth used in timestep calculation
    % ---------------------------------------------------------------------
    if flags.flag_CA == 1
        % Cellular Automata model
        flow_depth = depths.d_t;                  % mm
        depth_timestep = flow_depth / 1000;       % m

        % Prevent unrealistically high velocities in nearly dry cells
        safe_depth = flow_depth;
        safe_depth(safe_depth < CA_States.depth_tolerance) = 1e12; % mm

        % Face velocities (signed according to your original convention)
        velocities.vel_left  = (flow_rate.qout_left_t  / 1000 / 3600) .* ...
            Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
        velocities.vel_right = (flow_rate.qout_right_t / 1000 / 3600) .* ...
            Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
        velocities.vel_up    = (flow_rate.qout_up_t    / 1000 / 3600) .* ...
            Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
        velocities.vel_down  = (flow_rate.qout_down_t  / 1000 / 3600) .* ...
            Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s

        if flags.flag_D8 == 1
            velocities.vel_ne = (flow_rate.qout_ne_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
            velocities.vel_se = (flow_rate.qout_se_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
            velocities.vel_sw = (flow_rate.qout_sw_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
            velocities.vel_nw = (flow_rate.qout_nw_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth / 1000); % m/s
        end

    else
        % Inertial / Diffusive / Kinematic models
        flow_depth = Hf * 1000;                        % mm
        depth_timestep = max(flow_depth / 1000, [], 3); % m

        % Prevent unrealistically high velocities in nearly dry faces
        safe_depth = flow_depth;
        safe_depth(safe_depth < CA_States.depth_tolerance) = 1e12; % mm

        % Face velocities
        velocities.vel_right = (flow_rate.qout_right_t / 1000 / 3600) .* ...
            Wshed_Properties.Resolution ./ (safe_depth(:,:,1) / 1000); % m/s
        velocities.vel_left  = -[zeros(ny,1), velocities.vel_right(:,2:end)];

        velocities.vel_up    = (flow_rate.qout_up_t / 1000 / 3600) .* ...
            Wshed_Properties.Resolution ./ (safe_depth(:,:,2) / 1000); % m/s
        velocities.vel_down  = -[velocities.vel_up(2:end,:); zeros(1,nx)];

        if flags.flag_D8 == 1
            velocities.vel_ne = (flow_rate.qout_ne_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth(:,:,6) / 1000); % m/s
            velocities.vel_se = (flow_rate.qout_se_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth(:,:,7) / 1000); % m/s
            velocities.vel_sw = (flow_rate.qout_sw_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth(:,:,8) / 1000); % m/s
            velocities.vel_nw = (flow_rate.qout_nw_t / 1000 / 3600) .* ...
                Wshed_Properties.Resolution ./ (safe_depth(:,:,9) / 1000); % m/s
        end
    end

    %% --------------------------------------------------------------------
    % Maximum speed per direction (use absolute value, not signed max)
    % ---------------------------------------------------------------------
    velocities.max_velocity_left  = max(abs(velocities.vel_left(:)));
    velocities.max_velocity_right = max(abs(velocities.vel_right(:)));
    velocities.max_velocity_up    = max(abs(velocities.vel_up(:)));
    velocities.max_velocity_down  = max(abs(velocities.vel_down(:)));

    if flags.flag_D8 == 1
        velocities.max_velocity_ne = max(abs(velocities.vel_ne(:)));
        velocities.max_velocity_se = max(abs(velocities.vel_se(:)));
        velocities.max_velocity_sw = max(abs(velocities.vel_sw(:)));
        velocities.max_velocity_nw = max(abs(velocities.vel_nw(:)));
    end

    %% --------------------------------------------------------------------
    % Velocity raster = maximum local speed magnitude among all directions
    % ---------------------------------------------------------------------
    velocities.velocity_raster = max(abs(velocities.vel_left), abs(velocities.vel_right));
    velocities.velocity_raster = max(velocities.velocity_raster, abs(velocities.vel_up));
    velocities.velocity_raster = max(velocities.velocity_raster, abs(velocities.vel_down));

    if flags.flag_D8 == 1
        velocities.velocity_raster = max(velocities.velocity_raster, abs(velocities.vel_ne));
        velocities.velocity_raster = max(velocities.velocity_raster, abs(velocities.vel_se));
        velocities.velocity_raster = max(velocities.velocity_raster, abs(velocities.vel_sw));
        velocities.velocity_raster = max(velocities.velocity_raster, abs(velocities.vel_nw));

        velocities.velocity_vector = [ ...
            velocities.max_velocity_left, ...
            velocities.max_velocity_right, ...
            velocities.max_velocity_up, ...
            velocities.max_velocity_down, ...
            velocities.max_velocity_ne, ...
            velocities.max_velocity_se, ...
            velocities.max_velocity_sw, ...
            velocities.max_velocity_nw];
    else
        velocities.velocity_vector = [ ...
            velocities.max_velocity_left, ...
            velocities.max_velocity_right, ...
            velocities.max_velocity_up, ...
            velocities.max_velocity_down];
    end

    %% --------------------------------------------------------------------
    % Cell-centered velocity components
    % ---------------------------------------------------------------------
    if flags.flag_inertial == 1 || flags.flag_diffusive == 1 || ...
       flags.flag_kinematic == 1 || flags.flag_full_momentum == 1
    
        velocities.right_component = 0.5 * (velocities.vel_right - velocities.vel_left);
        velocities.up_component    = 0.5 * (velocities.vel_up    - velocities.vel_down);
    
    else
    
        velocities.right_component = velocities.vel_right - velocities.vel_left;
        velocities.up_component    = velocities.vel_up    - velocities.vel_down;
    
    end

    %% --------------------------------------------------------------------
    % D8 contributions projected into x/y components
    % Correct projection of diagonal velocity onto Cartesian axes = 1/sqrt(2)
    % ---------------------------------------------------------------------
    if flags.flag_D8 == 1
        c45 = 1 / sqrt(2);

        if flags.flag_CA == 1
            velocities.right_component = velocities.right_component + ...
                c45 * (0.5 * (velocities.vel_ne + velocities.vel_sw)) + ...
                c45 * (0.5 * (velocities.vel_se + velocities.vel_nw));

            velocities.up_component = velocities.up_component + ...
                c45 * (0.5 * (velocities.vel_ne + velocities.vel_sw)) - ...
                c45 * (0.5 * (velocities.vel_se + velocities.vel_nw));
        else
            velocities.right_component = velocities.right_component + ...
                c45 * (velocities.vel_ne - velocities.vel_sw) + ...
                c45 * (velocities.vel_se - velocities.vel_nw);

            velocities.up_component = velocities.up_component + ...
                c45 * (velocities.vel_ne - velocities.vel_sw) - ...
                c45 * (velocities.vel_se - velocities.vel_nw);
        end
    end

    %% --------------------------------------------------------------------
    % Total cell-centered velocity magnitude
    % ---------------------------------------------------------------------
    velocities.total_velocity = hypot(velocities.right_component, velocities.up_component);

    % Identifying stagnant areas
    idx_stagnant = velocities.total_velocity < 0.01; % Areas with very low total velocities

    % Maximum total velocity magnitude
    velocities.max_velocity = max(velocities.total_velocity(:));

    %% --------------------------------------------------------------------
    % Courant factor depending on connectivity
    % ---------------------------------------------------------------------
    if flags.flag_D8 == 1
        Courant_Parameters.factor_grid = sqrt(1/2);
    else
        Courant_Parameters.factor_grid = 1;
    end

    if flags.flag_GPU == 1
        Courant_Parameters.factor_grid = gpuArray(Courant_Parameters.factor_grid);
    end

    %% --------------------------------------------------------------------
    % Candidate 1: advective timestep based on maximum velocity
    % ---------------------------------------------------------------------
    if velocities.max_velocity > 0
        Courant_Parameters.slope_alfa = 2;
        Courant_Parameters.v_threshold = 5;

        dt_adv = Courant_Parameters.factor_grid * Wshed_Properties.Resolution / ...
            max(velocities.max_velocity, 1e-10); % seconds

        Courant_Parameters.time_step_factor = max( ...
            Courant_Parameters.alfa_max - ...
            Courant_Parameters.slope_alfa * max(velocities.max_velocity - Courant_Parameters.v_threshold, 0), ...
            Courant_Parameters.alfa_min);

        dt_adv = dt_adv * Courant_Parameters.time_step_factor;
        dt_adv = double(min(dt_adv, running_control.max_time_step));
    else
        dt_adv = running_control.max_time_step;
    end

    %% --------------------------------------------------------------------
    % Candidate 2: Bates-like shallow-water timestep
    % NOTE: keeping stagnant masking as requested by user
    % ---------------------------------------------------------------------
    wave_celerity = sqrt(9.81 * max(depth_timestep, 0)); % m/s
    wave_celerity = max(wave_celerity, 1e-10);
    wave_celerity(isinf(wave_celerity)) = nan;
    wave_celerity(idx_stagnant) = nan; % kept intentionally as requested

    if flags.flag_adaptive_timestepping == 1
        dt_wave = nanmin(Courant_Parameters.alfa_min * Wshed_Properties.Resolution ./ wave_celerity, [], 'all');
    else
        max_vel = velocities.velocity_raster; % fixed bug: use local speed raster
        dt_wave = nanmin(Courant_Parameters.alfa_min * Wshed_Properties.Resolution ./ ...
            (max_vel + wave_celerity), [], 'all');
    end

    if isnan(dt_wave)
        dt_wave = running_control.max_time_step;
    end

    %% --------------------------------------------------------------------
    % Candidate 3: Full-momentum SWE CFL timestep
    % Active only when flags.flag_full_momentum == 1
    %
    % SERGHEI-like CFL:
    %   dt <= CFL * dx / max( |u| + sqrt(g*h), |v| + sqrt(g*h) )
    %
    % Important:
    % - Do NOT remove stagnant cells here. Even if u = v = 0, the gravity wave
    %   speed sqrt(g*h) still restricts the timestep.
    % - Use cell water depth, not face Hf, for the SWE wave celerity.
    % ---------------------------------------------------------------------
    if flags.flag_full_momentum == 1
    
        g_swe = 9.81;
    
        % Cell-centered water depth [m]
        h_swe = max(depths.d_t / 1000, 0);
    
        % Dry tolerance [m]
        if isfield(CA_States, 'depth_tolerance')
            h_dry_swe = max(CA_States.depth_tolerance / 1000, 1e-8);
        else
            h_dry_swe = 1e-8;
        end
    
        % Cell-centered velocity components [m/s]
        u_swe = abs(velocities.right_component);
        v_swe = abs(velocities.up_component);
    
        % Remove invalid values
        u_swe(~isfinite(u_swe)) = 0;
        v_swe(~isfinite(v_swe)) = 0;
    
        % Gravity wave speed [m/s]
        c_swe = sqrt(g_swe * h_swe);
    
        % Full SWE signal speed per cell [m/s]
        speed_x = u_swe + c_swe;
        speed_y = v_swe + c_swe;
        speed_swe = max(speed_x, speed_y);
    
        % Ignore dry cells only. Do NOT ignore stagnant wet cells.
        speed_swe(h_swe <= h_dry_swe) = NaN;
        speed_swe(~isfinite(speed_swe)) = NaN;
    
        % CFL number for explicit full momentum.
        % SERGHEI uses CFL <= 0.5. Use the model's alfa_min but cap at 0.5.
        CFL_full_momentum = min(Courant_Parameters.alfa_min, 0.5);
    
        dt_full_momentum = nanmin( ...
            CFL_full_momentum * Wshed_Properties.Resolution ./ speed_swe, ...
            [], 'all');
    
        if isnan(dt_full_momentum) || isinf(dt_full_momentum)
            dt_full_momentum = running_control.max_time_step;
        end
    
        % Store diagnostics
        Courant_Parameters.dt_full_momentum = dt_full_momentum;
        Courant_Parameters.CFL_full_momentum = CFL_full_momentum;
        Courant_Parameters.max_speed_full_momentum = max(speed_swe(:), [], 'omitnan');
    
    else
    
        dt_full_momentum = running_control.max_time_step;
    
    end
    %% --------------------------------------------------------------------
    % Final timestep = minimum of all active stability constraints
    % ---------------------------------------------------------------------
    if flags.flag_full_momentum == 1
        % new_timestep = min([dt_adv, dt_wave, dt_full_momentum]);
        new_timestep = dt_full_momentum;
    else
        if flags.flag_timestep == 1 % 1 we use all terms
            new_timestep = min(dt_adv, dt_wave);
        else
            new_timestep = dt_wave;
        end
    end

    %% --------------------------------------------------------------------
    % Water quality restriction
    % ---------------------------------------------------------------------
    if flags.flag_waterquality == 1
        Courant_Parameters.alfa_wq = 1; % Reduction factor of water quality
        new_timestep = min(Courant_Parameters.alfa_wq * tmin_wq, new_timestep);
    end

    %% --------------------------------------------------------------------
    % Safety fallback
    % ---------------------------------------------------------------------
    if isnan(new_timestep) || isinf(new_timestep)
        new_timestep = running_control.max_time_step;
    end

    %% --------------------------------------------------------------------
    % Snap timestep to allowed increment range
    % ---------------------------------------------------------------------
    running_control.time_calculation_routing = new_timestep;

    running_control.time_calculation_routing = max( ...
        running_control.time_step_increments * ...
        floor(running_control.time_calculation_routing / running_control.time_step_increments), ...
        running_control.min_time_step);

    running_control.time_calculation_routing = min( ...
        running_control.time_calculation_routing, ...
        running_control.max_time_step);

    time_step = running_control.time_calculation_routing / 60; % new time-step for the next run

    if running_control.time_calculation_routing == running_control.min_time_step
        unstable = 1; % Reached minimum timestep; may indicate unstable period
    end

elseif k == 1
    running_control.time_calculation_routing = time_step * 60; % sec
else
    running_control.time_calculation_routing = time_step * 60; % sec
end

% Store current max velocity for next call
velocities.max_velocity_previous = velocities.max_velocity;