function [SWE_t, H_snow_t, M_snow, P_snow, P_rain, rho_snow, E_s, mass_balance_error_t] = Snow_Model_Function(SWE_prev, H_snow_prev, T_air, T_min, P, wind, lat, DOY, H_snow_t0, alpha, epsilon, C_e, DDF, T_thesh, rho_snow_init, rho_max, k_t, k_swe, k_D)
%% ═══════════════════════════════════════════════════════════════════════
%  Function: Snow_Model_Function
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 04/02/2025
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Perform a single time-step update of the snowpack properties using a
%      physically based approach. This function partitions precipitation into 
%      snow and rain, computes solar and longwave radiation, updates snow water 
%      equivalent (SWE) and snow depth, calculates snowmelt via degree-day and 
%      radiation methods, estimates sublimation, and checks the mass balance.
%
%  ➤ Inputs:
%      • SWE_prev      - Previous Snow Water Equivalent (mm) [Nx-by-Ny]
%      • H_snow_prev   - Previous snow depth (mm) [Nx-by-Ny]
%      • T_air         - Air temperature (°C) [Nx-by-Ny]
%      • T_min         - Minimum air temperature (°C) [Nx-by-Ny]
%      • P             - Precipitation (mm) [Nx-by-Ny]
%      • wind          - Wind speed (m/s) [Nx-by-Ny]
%      • lat           - Latitude (degrees) [Nx-by-Ny]
%      • DOY           - Day of year (scalar)
%      • H_snow_t0     - Initial snow depth reference (mm) [Nx-by-Ny]
%      • alpha         - Snow albedo parameter
%      • epsilon       - Snow emissivity
%      • C_e           - Sublimation coefficient
%      • DDF           - Degree-day factor (mm/°C/day)
%      • T_thesh       - Temperature threshold for snow/rain partitioning (°C)
%      • rho_snow_init - Initial snow density (kg/m³)
%      • rho_max       - Maximum snow density (kg/m³)
%      • k_t           - Snow compaction rate due to temperature
%      • k_swe         - Snow compaction rate due to SWE
%      • k_D           - Compaction rate due to snowpack thickness parameter D
%
%  ➤ Outputs:
%      • SWE_t               - Updated Snow Water Equivalent (mm) [Nx-by-Ny]
%      • H_snow_t            - Updated snow depth (mm) [Nx-by-Ny]
%      • M_snow              - Computed snowmelt (mm) [Nx-by-Ny]
%      • P_snow              - Precipitation partitioned as snow (mm) [Nx-by-Ny]
%      • P_rain              - Precipitation partitioned as rain (mm) [Nx-by-Ny]
%      • rho_snow            - Updated snow density (kg/m³) [Nx-by-Ny]
%      • E_s                 - Sublimation from the snowpack (mm) [Nx-by-Ny]
%      • mass_balance_error_t- Mass balance error for this time-step (mm)
%
%  ➤ Local (Auxiliary) Functions:
%      • compute_radiation                  - Calculates incoming solar and longwave radiation.
%      • compute_net_radiation              - Computes net radiation at the snow surface.
%      • refined_precipitation_partitioning - Partitions precipitation into snow and rain.
%      • update_snow_density                - Updates snow density and computes new snow depth.
%      • compute_snowmelt                   - Estimates snowmelt from degree-day and radiation.
%      • compute_sublimation                - Estimates sublimation based on wind and humidity.
%
%  ➤ Notes:
%      • This function uses a combination of energy balance and degree-day methods
%        to simulate snowmelt.
%      • Precipitation is partitioned based on air temperature using a logistic/linear
%        function to determine the snow fraction.
%      • The mass balance error is computed by comparing the initial SWE plus added
%        snowfall against the sum of snowmelt, sublimation, and final SWE.
%      • Auxiliary functions are included below for radiation, snow density, snowmelt,
%        and sublimation computations.
% ═══════════════════════════════════════════════════════════════════════

    % Compute solar and longwave radiation
    [S_down, L_down] = compute_radiation(T_air, lat, DOY);

    % Partition precipitation into snow or rain
    P_snow = refined_precipitation_partitioning(P, T_air);

    % Compute Rainfall
    P_rain = P - P_snow;

    % Update SWE
    SWE_t = SWE_prev + P_snow;

    % Compute net radiation
    Q_net = compute_net_radiation(S_down, L_down, T_air, alpha, epsilon);

    % Update snow density and depth
    [rho_snow, H_snow_t] = update_snow_density(SWE_prev, rho_snow_init, T_air, H_snow_t0, k_t, k_swe, k_D, rho_max);

    % Compute snowmelt
    M_snow = compute_snowmelt(T_air, Q_net, DDF, alpha, SWE_t);

    % Update SWE_t
    SWE_t = SWE_t - M_snow;

    % Compute sublimation
    E_s = compute_sublimation(wind, T_air, T_min, SWE_prev, C_e, SWE_t);

    % Update SWE
    SWE_t = SWE_t - E_s;

    % Mass balance check
    total_initial_SWE = nansum(SWE_prev(:));
    total_precip = nansum(P_snow(:));
    total_melt = nansum(M_snow(:));
    total_sublimation = nansum(E_s(:));
    total_final_SWE = nansum(SWE_t(:));

    mass_balance_error_t = total_initial_SWE + total_precip - total_melt - total_sublimation - total_final_SWE;
end

%% AUXILIARY FUNCTIONS

% Compute solar and longwave radiation without humidity
function [S_down, L_down] = compute_radiation(T_air, lat, DOY)
    sigma = 5.67e-8; % Stefan-Boltzmann constant (W/m²K⁴)
    T_kelvin = T_air + 273.15; % Convert to Kelvin

    % Compute solar zenith angle
    decl = 23.45 * sind(360 * (DOY - 81) / 365); % Solar declination
    lat_rad = deg2rad(lat);
    decl_rad = deg2rad(decl);
    theta_z = acos(sin(lat_rad) .* sin(decl_rad) + cos(lat_rad) .* cos(decl_rad)); % Zenith angle
    
    % Shortwave radiation using clear-sky transmissivity
    S0 = 1367; % Solar constant (W/m²)
    tau = 0.75; % Fixed clear-sky transmissivity
    S_down = tau .* S0 .* cos(theta_z); % Shortwave radiation
    S_down(S_down < 0) = 0;

    % Longwave radiation using Brutsaert’s formula (based on temperature)
    epsilon_a = 0.642 + 0.035 * sqrt(max(T_kelvin - 273.15,0)); % Approximate emissivity
    L_down = squeeze(epsilon_a) .* sigma .* squeeze(T_kelvin).^4; % Downward longwave radiation
end

% Compute net radiation
function Q_net = compute_net_radiation(S_down, L_down, T_air, alpha, epsilon)
    % Albedo for snow
    % alpha = params.alpha; % Snow albedo
    % epsilon = params.epsilon; % Emissivity of snow
    sigma = 5.67e-8; % Stefan-Boltzmann constant (W/m²K⁴)

    % Assume snow temperature is at freezing point (0°C) for simplicity
    T_snow = 0; % Snow surface temperature (°C)
    T_snow_kelvin = T_snow + 273.15; % Convert to Kelvin

    % Outgoing longwave radiation from snow surface
    L_up = epsilon * sigma * T_snow_kelvin^4; % Longwave emission from snow surface

    % Compute net radiation
    Q_net = S_down * (1 - alpha) - L_up + L_down - epsilon * sigma * (T_snow_kelvin^4 - (T_air + 273.15).^4);
end


% Refined precipitation partitioning (using logistic function for snow fraction)
function P_snow = refined_precipitation_partitioning(P, T_air)
    % Partition precipitation into snow and rain using a logistic function
    % Inputs:
    %   P        - Precipitation (mm)
    %   T_air    - Air temperature (°C)
    %   a        - Parameter for the logistic function that determines the steepness
    % Outputs:
    %   P_snow   - Precipitation that falls as snow (mm)
    
    % Logistic Function
    % Calculate snow fraction using logistic function
    % snow_fraction = 1 ./ (1 + exp(-a * (T_air - 0))); % Threshold at 0°C
    % P_snow = P .* snow_fraction; % Snow precipitation

    % Linear Function
    % Compute snow fraction using a linear function
    Tthreshold_min = 4; % Celcius 
    Tthreshold_max = 7; % Celcius
    snow_fraction = (Tthreshold_max - T_air) / (Tthreshold_max - Tthreshold_min);
    snow_fraction(snow_fraction > 1) = 1; snow_fraction(snow_fraction < 0) = 0;

    % Compute snowfall
    P_snow = P .* snow_fraction;

end

% Update snow density and depth (compaction included)
function [rho_snow, H_snow] = update_snow_density(SWE, rho_snow, T_air, D, k_t, k_swe, k_D, rho_max)
    rho_snow = rho_snow + k_t .* T_air + k_swe .* SWE + k_D .* D;
    rho_snow = min(rho_snow, rho_max);
    
    H_snow = SWE ./ rho_snow;
    H_snow(SWE == 0) = 0;
end

% Compute snowmelt using Degree-Day or Energy Balance
function M_snow = compute_snowmelt(T_air, radiation, DDF, alpha, SWE_t)
    M_snow_temp = max(DDF .* T_air, 0); % Degree-day melt
    M_snow_rad = (1 - alpha) .* radiation / 334; % Radiation-driven melt (mm/day)
    M_snow = max(M_snow_temp + M_snow_rad, 0); % Ensure no negative melt
    M_snow = min(M_snow, SWE_t); % Ensure it does not melt more than the available
    M_snow(isnan(T_air)) = nan;
end

function E_s = compute_sublimation(wind, T_air, T_min, SWE, C_e, SWE_t)
    % Compute sublimation rate based on physically-based specific humidity

    % Estimate saturation vapor pressure (hPa)
    e_s = 6.112 .* exp((17.67 .* T_air) ./ (T_air + 243.5));

    e_a = 10*(0.61 * exp((17.27 * T_min) ./ (T_min + 237.3))); % Actual vapor pressure [hPa] (alternative if RH is unavailable)

    % Compute specific humidity values
    q_s = 0.622 .* e_s ./ (1013.25 - e_s); % Snow surface specific humidity
    q_a = 0.622 .* e_a ./ (1013.25 - e_a); % Air specific humidity

    % Compute sublimation flux
    E_s = C_e .* wind .* (q_s - q_a) .* SWE;
    E_s = max(E_s, 0); % Ensure no negative sublimation

    % Ensure E_s not larger than SWE_t
    E_s = min(E_s, SWE_t);
end
    