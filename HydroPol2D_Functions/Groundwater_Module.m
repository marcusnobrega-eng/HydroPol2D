%% ========================================================================
%  GROUNDWATER RECHARGE AND FLOW MODULE
%  - Initializes recharge
%  - Simulates unsaturated zone storage
%  - Runs linear reservoir recharge
%  - Simulates groundwater flow using Boussinesq 2D
%  - Tracks groundwater states and model errors
% =========================================================================

%% ------------------------------------------------------------------------
% 1. INITIALIZE CURRENT RECHARGE
%    (Only for the first time step)
% -------------------------------------------------------------------------
if k == 1
    current_recharge = zeros(size(DEM_raster.Z));
    current_recharge(isnan(DEM_raster.Z)) = nan;
end

%% ------------------------------------------------------------------------
% 2. CALCULATE GROUNDWATER DEPTH AND UNSATURATED ZONE STORAGE
% -------------------------------------------------------------------------
GW_Depth = (BC_States.h_t - (elevation - Soil_Properties.Soil_Depth)); % Groundwater depth [m]
UZ_max_storage = (Soil_Properties.Soil_Depth - GW_Depth) .* ...
                 (Soil_Properties.teta_sat - Soil_Properties.teta_i);  % Max storage in unsaturated zone [m water]

%% ------------------------------------------------------------------------
% 3. RUN LINEAR RESERVOIR MODEL FOR RECHARGE
% -------------------------------------------------------------------------
[recharge_rate, Soil_Moisture, cumulative_recharge] = simulate_groundwater_recharge( ...
    0 * Hydro_States.f / 1000 / 3600, ...                % No initial infiltration [m/s]
    Soil_Properties.I_t / 1000, ...                      % Initial soil moisture [m]
    Soil_Properties.k, ...                               % Hydraulic conductivity [m/s]
    time_step * 60, ...                                  % Time step [s]
    min_soil_moisture / 1000, ...                        % Minimum soil moisture [m]
    UZ_max_storage, ...                                  % Unsaturated zone maximum storage [m]
    LULC_Properties.idx_imp, ...                         % Impervious surface mask
    current_recharge ...                                 % Recharge memory [mm]
);

% Update soil moisture for the next time step [mm]
Soil_Properties.I_t = Soil_Moisture * 1000;

%% ------------------------------------------------------------------------
% 4. GROUNDWATER FLOW SIMULATION (BOUSSINESQ 2D EXPLICIT MODEL)
% -------------------------------------------------------------------------
if flags.flag_baseflow == 1
    
    % Run Boussinesq 2D Explicit Model
    [BC_States.h_t, ~, ~, q_exf, q_river, error] = Boussinesq_2D_explicit( ...
        time_step * 60, ...                                % Time step [s]
        Wshed_Properties.Resolution, ...                   % X cell size [m]
        Wshed_Properties.Resolution, ...                   % Y cell size [m]
        BC_States.h_t, ...                                 % Current groundwater head [m]
        (Elevation_Properties.elevation_cell - Soil_Properties.Soil_Depth), ... % Bedrock elevation [m]
        Soil_Properties.Sy, ...                            % Specific yield [-]
        recharge_rate, ...                                 % Recharge rate [m/s]
        Soil_Properties.ksat_gw / 1000 / 3600, ...          % GW saturated hydraulic conductivity [m/s]
        idx_rivers, ...                                    % River cells indices
        LULC_Properties.River_K_coeff * Soil_Properties.ksat / 1000 / 3600, ... % River conductance [m/s]
        (Elevation_Properties.elevation_cell + depths.d_t / 1000), ...          % River stage [m]
        Elevation_Properties.elevation_cell, ...           % River bed elevation [m]
        0.5, ...                                           % Porosity [-]
        Soil_Properties.Soil_Depth, ...                    % Soil depth [m]
        Wshed_Properties.domain, ...                       % Watershed domain mask
        [], [], Wshed_Properties.perimeter ...             % Optional boundary conditions
    );
    
    % Update previous groundwater head state
    BC_States.h_0 = BC_States.h_t;
    
    % Update surface water storage (excess saturation + river exchanges)
    depths.d_t = depths.d_t + time_step * 60 * (q_river * 1000 + q_exf * 1000); % [mm]
end

%% ------------------------------------------------------------------------
% 5. TRACK MAXIMUM GROUNDWATER DEPTH
% -------------------------------------------------------------------------
if k == 1
    max_GW_depth = zeros(size(DEM_raster.Z));
    max_GW_depth(isnan(max_GW_depth)) = nan;
end

% Update maximum groundwater depth observed during simulation [m]
max_GW_depth = max(max_GW_depth, BC_States.h_t - (elevation - Soil_Properties.Soil_Depth));

%% ------------------------------------------------------------------------
% 6. UPDATE EFFECTIVE RECHARGE MEMORY
% -------------------------------------------------------------------------
current_recharge = cumulative_recharge - (1/1000) * q_exf * (time_step * 60); % [mm]

%% ------------------------------------------------------------------------
% 7. STORE MODEL ERROR
% -------------------------------------------------------------------------
errors(5) = error; % Store Boussinesq model error metric

