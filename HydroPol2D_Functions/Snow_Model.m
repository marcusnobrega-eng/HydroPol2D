clc; clear; close all;

%% Create Synthetic Catchment
Nx = 50; Ny = 50; T = 10; % Grid size and time steps
[x, y] = meshgrid(1:Nx, 1:Ny);

% Generate synthetic inputs
T_air = repmat(-5 + 1 * ones(Nx, Ny, T), [1, 1, 1]); % Air temperature (-5 to 5°C)
P = repmat(5 * ones(Nx, Ny, T), [1, 1, 1]); % Precipitation (0 to 5 mm)
wind = repmat(2 + 3 * ones(Nx, Ny, T), [1, 1, 1]); % Wind speed (2 to 5 m/s)
lat = 30 + 10 * (y / Ny); % Latitude varying from 30° to 40°
DOY = linspace(1, 365, T); % Days of year
H_snow_t = 0.1 + 0.05 * ones(Nx, Ny); % Snowpack thickness parameter [m]
Albedo = ones(Nx,Ny)*0.85; % Albedo value [-]

% Model parameters
params.alpha = Albedo;
params.epsilon = 0.98; % Emissivity of snow
params.C_e = 0.001; % Sublimation coefficient
params.DDF = 4; % mm/°C/day, degree-day factor
params.T_thresh = 0; % Threshold temperature for snow/rain partitioning
params.rho_snow_init = 100; % Initial snow density (kg/m³)
params.rho_max = 400; % Max snow density (kg/m³)
params.k_t = 0.1; % Snow compaction rate due to temperature
params.k_swe = 0.001; % Snow compaction rate due to SWE
params.k_D = 0.02; % Compaction rate due to D
params.snow_fraction_a = 0.2; % Snow fraction parameter for logistic function

%% Run Snow Model
[SWE, H_snow, mass_balance_error] = snow_model_2D(T_air, P, wind, lat, DOY, H_snow_t, params);

%% Visualization
for t = [1, round(T/2), T] % Show results for initial, mid, and final time steps
    figure;
    subplot(1,2,1);
    imagesc(SWE(:,:,t));
    colorbar; title(['SWE at t = ' num2str(t)]);
    xlabel('X'); ylabel('Y'); axis equal tight;
    
    subplot(1,2,2);
    imagesc(P(:,:,t));
    colorbar; title(['Precipitation at t = ' num2str(t)]);
    xlabel('X'); ylabel('Y'); axis equal tight;
end


function [SWE, H_snow, mass_balance_error] = snow_model_2D(T_air, P, wind, lat, DOY, H_snow_t, params)
    % Spatially distributed snow model with automatic radiation and mass balance check
    % Inputs:
    %   T_air  - Air temperature (time, Nx, Ny) in °C
    %   P      - Precipitation (time, Nx, Ny) in mm
    %   wind   - Wind speed (time, Nx, Ny) in m/s
    %   lat    - Latitude (Nx, Ny) in degrees
    %   DOY    - Day of year (time vector)
    %   H_snow_t      - Snowpack thickness parameter (Nx, Ny)
    %   params - Struct with model parameters
    % Outputs:
    %   SWE    - Snow water equivalent (mm)
    %   H_snow - Snow depth (m)
    %   mass_balance_error - Mass balance check over time

    [Nx, Ny, T] = size(T_air);
    
    % Initialize variables
    SWE = zeros(Nx, Ny, T);
    H_snow = zeros(Nx, Ny, T);
    rho_snow = ones(Nx, Ny) * params.rho_snow_init;
    mass_balance_error = zeros(T, 1);
    initial_SWE = zeros(Nx, Ny);

    for t = 1:T
        if t == 1
            initial_SWE = SWE(:,:,1);
        else
            initial_SWE = SWE(:,:,t-1);
        end

        % Compute solar and longwave radiation without RH
        [S_down, L_down] = compute_radiation(T_air(:,:,t), lat, DOY(t), params);

        % Partition precipitation into snow or rain using the refined method
        P_snow = refined_precipitation_partitioning(P(:,:,t), T_air(:,:,t), params.snow_fraction_a);

        % Compute net radiation
        Q_net = compute_net_radiation(S_down, L_down, T_air(:,:,t), params);

        % Update snow density and depth
        [rho_snow, H_snow_t] = update_snow_density(H_snow(:,:,max(t-1,1)), SWE(:,:,max(t-1,1)), rho_snow, T_air(:,:,t), H_snow_t, params);

        % Compute snowmelt
        M_snow = compute_snowmelt(T_air(:,:,t), Q_net, params, SWE(:,:,max(t-1,1)));

        % Compute sublimation without RH (assume empirical rate)
        E_s = compute_sublimation(wind(:,:,t), T_air(:,:,t), SWE(:,:,max(t-1,1)), params);

        % Update SWE
        SWE_t = SWE(:,:,max(t-1,1)) + P_snow - M_snow - E_s;
        SWE_t = max(SWE_t, 0);

        % Store results
        SWE(:,:,t) = SWE_t;
        H_snow(:,:,t) = H_snow_t;

        % Mass balance check
        total_initial_SWE = sum(initial_SWE(:));
        total_precip = sum(P_snow(:));
        total_melt = sum(M_snow(:));
        total_sublimation = sum(E_s(:));
        total_final_SWE = sum(SWE_t(:));

        mass_balance_error(t) = total_initial_SWE + total_precip - total_melt - total_sublimation - total_final_SWE;
    end
end

%% AUXILIARY FUNCTIONS

% Compute solar and longwave radiation without humidity
function [S_down, L_down] = compute_radiation(T_air, lat, DOY, params)
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
function Q_net = compute_net_radiation(S_down, L_down, T_air, params)
    % Albedo for snow
    alpha = params.alpha; % Snow albedo
    epsilon = params.epsilon; % Emissivity of snow
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
function P_snow = refined_precipitation_partitioning(P, T_air, a)
    % Partition precipitation into snow and rain using a logistic function
    % Inputs:
    %   P        - Precipitation (mm)
    %   T_air    - Air temperature (°C)
    %   a        - Parameter for the logistic function that determines the steepness
    % Outputs:
    %   P_snow   - Precipitation that falls as snow (mm)
    
    % Calculate snow fraction using logistic function
    snow_fraction = 1 ./ (1 + exp(-a * (T_air - 0))); % Threshold at 0°C
    P_snow = P .* snow_fraction; % Snow precipitation
end

% Update snow density and depth (compaction included)
function [rho_snow, H_snow] = update_snow_density(H_snow, SWE, rho_snow, T_air, D, params)
    k_t = params.k_t; 
    k_swe = params.k_swe; 
    rho_snow = rho_snow + k_t .* T_air + k_swe .* SWE + params.k_D .* D;
    rho_snow = min(rho_snow, params.rho_max);
    
    H_snow = SWE ./ rho_snow;
    H_snow(SWE == 0) = 0;
end

% Compute snowmelt using Degree-Day or Energy Balance
function M_snow = compute_snowmelt(T_air, radiation, params, SWE_t)
    M_snow_temp = max(params.DDF .* T_air, 0); % Degree-day melt
    M_snow_rad = (1 - params.alpha) .* radiation / 334; % Radiation-driven melt (mm/day)
    M_snow = max(M_snow_temp + M_snow_rad, 0); % Ensure no negative melt
    M_snow = min(M_snow, SWE_t); % Ensure it does not melt more than the available
end

function E_s = compute_sublimation(wind, T_air, SWE, params)
    % Estimate saturation vapor pressure (hPa)
    e_s = 6.112 .* exp((17.67 .* T_air) ./ (T_air + 243.5));

    % Assume constant specific humidity for mid-latitude snow-covered areas (g/kg)
    q_a = 2 / 1000; % Convert g/kg to kg/kg

    % Estimate actual vapor pressure (hPa)
    e_a = q_a .* (1013.25 ./ (0.622 + q_a));

    % Estimate relative humidity (RH) in %
    RH = (e_a ./ e_s) * 100;
    RH = max(0, min(RH, 100)); % Ensure RH is within 0-100%

    % Compute specific humidity values
    q_s = 0.622 .* e_s ./ (1013.25 - e_s); % Snow surface specific humidity
    q_a = 0.622 .* e_a ./ (1013.25 - e_a); % Air specific humidity

    % Compute sublimation flux
    E_s = params.C_e .* wind .* (q_s - q_a) .* SWE;
    E_s = max(E_s, 0); % Ensure no negative sublimation
end
    