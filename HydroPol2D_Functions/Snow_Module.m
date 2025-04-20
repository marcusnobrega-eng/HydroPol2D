%% ═══════════════════════════════════════════════════════════════════════
%  Module: Snow_Module
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 04/02/2025
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Estimate snow water equivalent, snowpack, snowmelt, and snow density.
%      This module partitions precipitation between snow and rain based on 
%      temperature, computes snowmelt using a degree-day approach, updates 
%      snowpack properties, and adjusts the water balance by updating the 
%      effective water depth.
%
%  ➤ Inputs:
%      • flags: Structure with flags controlling processes, including:
%              - flag_snow_modeling (activates snow modeling)
%      • BC_States: Structure containing:
%              - Eff_Rainfall: Effective rainfall (mm)
%              - Average_Daily_Temperature: Air temperature (°C)
%              - wind: Wind speed (m/s)
%              - min_temp: Minimum temperature (°C)
%      • Wshed_Properties:
%              - pixel_latitude: Latitude values for the watershed pixels
%      • day_of_year: Current day of the year (DOY)
%      • Snow_Properties.SWE_t, Snow_Properties.H_snow_t: Existing snow water 
%              equivalent (SWE) and snowpack depth (mm)
%      • DEM_raster: Digital Elevation Model structure (provides cellsize)
%      • depths:
%              - d_p: Ponded water depth (mm)
%              - d_t: Total water depth (mm), updated based on snowmelt and rain
%      • k: Time-step counter (for initialization of maximum snow depth)
%
%  ➤ Outputs:
%      • Snow_Properties updated:
%              - SWE_t: Snow Water Equivalent (mm)
%              - H_snow_t: Snowpack depth (mm)
%              - M_snow: Snowmelt (mm)
%              - P_snow: Snow precipitation (mm)
%              - P_rain: Rain precipitation (mm)
%              - rho_snow: Snow density (kg/m³)
%              - E_s: Canopy evaporation/sublimation from snow (mm)
%      • errors(2): Mass balance error in snow (m³), adjusted by DEM cell area
%      • depths.d_t: Updated total water depth on the grid (mm)
%      • max_Hsnow: Maximum snow depth recorded (mm)
%
%  ➤ Notes:
%      • Snow modeling is activated when flag_snow_modeling is set to 1.
%      • The Snow_Model_Function is called to update snow properties based on
%        current conditions (temperature, wind, precipitation, etc.).
%      • A mass balance error is computed; a warning is issued if this error 
%        exceeds 100.
%      • The effective water depth (depths.d_t) is computed as the sum of ponded
%        water depth, snowmelt, and rain, thereby updating the mass balance.
%      • If snow modeling is inactive, the effective rainfall is simply added to 
%        the ponded water depth.
% ═══════════════════════════════════════════════════════════════════════

flags.flag_snow_modeling = 1;


if flags.flag_snow_modeling == 1

    % Snow Parameters
    Snow_Properties.alpha = 0.8;
    Snow_Properties.epsilon = 0.98; % Emissivity of snow
    Snow_Properties.C_e = 0.001; % Sublimation coefficient
    Snow_Properties.DDF = 2; % mm/°C/day, degree-day factor
    Snow_Properties.T_thresh = 0; % Threshold temperature for snow/rain partitioning
    Snow_Properties.rho_snow_init = 100; % Initial snow density (kg/m³)
    Snow_Properties.rho_max = 400; % Max snow density (kg/m³)
    Snow_Properties.k_t = 0.1; % Snow compaction rate due to temperature
    Snow_Properties.k_swe = 0.001; % Snow compaction rate due to SWE
    Snow_Properties.k_D = 0.02; % Compaction rate due to D
    Snow_Properties.snow_fraction_a = 0.2; % Snow fraction parameter for logistic function

    % Snow Model
    P = BC_States.Eff_Rainfall; % mm
    lat = Wshed_Properties.pixel_latitude;
    T_air = BC_States.Average_Daily_Temperature;
    wind = BC_States.wind;
    DOY = day_of_year;
    T_min = BC_States.min_temp;
    [Snow_Properties.SWE_t, Snow_Properties.H_snow_t, Snow_Properties.M_snow, Snow_Properties.P_snow, Snow_Properties.P_rain, Snow_Properties.rho_snow, Snow_Properties.E_s, mass_balance_error_snow] = Snow_Model_Function(Snow_Properties.SWE_t, Snow_Properties.H_snow_t, T_air, T_min, P, wind, lat, DOY, Snow_Properties.H_snow_t, Snow_Properties.alpha, Snow_Properties.epsilon, Snow_Properties.C_e, Snow_Properties.DDF, Snow_Properties.T_thresh, Snow_Properties.rho_snow_init, Snow_Properties.rho_max, Snow_Properties.k_t, Snow_Properties.k_swe, Snow_Properties.k_D);

    if mass_balance_error_snow > 100
        warning('Mass balance error in snow too large.')
    end
    
    errors(2) = mass_balance_error_snow*DEM_raster.cellsize^2; % m3 

    % Total Snow Water Equivalent
    if k == 1
        max_Hsnow = zeros(size(DEM_raster.Z)); max_Hsnow(isnan(max_Hsnow)) = nan;
    end
    max_Hsnow = max(max_Hsnow,Snow_Properties.H_snow_t); % Max snow depth [mm]

    % Depth mass balance
    depths.d_t = depths.d_p + Snow_Properties.M_snow + Snow_Properties.P_rain; % mm (considering what actually snow melted)

else
    % Depth mass balance
    depths.d_t = depths.d_p + BC_States.Eff_Rainfall; % mm
end