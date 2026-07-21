% Snow accumulation, melt, and sublimation.
% This script is called from Hydrological_Model and shares its workspace.

if flags.flag_snow_modeling == 1
    required_fields = {'SWE_t','H_snow_t','rho_snow','alpha','epsilon','C_e', ...
        'DDF','T_snow_all','T_rain_all','rho_snow_init','rho_max','k_t','k_swe','k_D'};
    missing_fields = required_fields(~isfield(Snow_Properties, required_fields));
    if ~isempty(missing_fields)
        error('Snow_Properties is incomplete. Missing: %s', strjoin(missing_fields, ', '));
    end

    dt_s = time_step * 60;
    forcing_fields = {'Eff_Rainfall','Average_Daily_Temperature','min_temp','wind'};
    missing_forcing = forcing_fields(~isfield(BC_States, forcing_fields));
    if ~isempty(missing_forcing)
        error(['Snow modeling requires internal meteorological forcing. Missing BC_States fields: %s. ', ...
            'Set flag_ETP = 1 and flag_input_ETP_map = 0, and provide ETP_input_data.xlsx.'], ...
            strjoin(missing_forcing, ', '));
    end
    active_cells = isfinite(Snow_Properties.SWE_t);
    has_invalid_climate = any(~isfinite(BC_States.Average_Daily_Temperature(active_cells))) || ...
        any(~isfinite(BC_States.min_temp(active_cells))) || ...
        any(~isfinite(BC_States.wind(active_cells)));
    if has_invalid_climate
        error(['Snow modeling received incomplete internal meteorological forcing. ', ...
            'Check the ETP_input_data.xlsx time coverage and station values.']);
    end
    P = BC_States.Eff_Rainfall;
    lat = Wshed_Properties.pixel_latitude;
    T_air = BC_States.Average_Daily_Temperature;
    wind = BC_States.wind;
    T_min = BC_States.min_temp;

    [Snow_Properties.SWE_t, Snow_Properties.H_snow_t, ...
        Snow_Properties.M_snow, Snow_Properties.P_snow, ...
        Snow_Properties.P_rain, Snow_Properties.rho_snow, ...
        Snow_Properties.E_s, mass_balance_error_snow] = ...
        Snow_Model_Function(Snow_Properties.SWE_t, Snow_Properties.H_snow_t, ...
        Snow_Properties.rho_snow, T_air, T_min, P, wind, lat, day_of_year, ...
        Snow_Properties.alpha, Snow_Properties.epsilon, Snow_Properties.C_e, ...
        Snow_Properties.DDF, Snow_Properties.T_snow_all, Snow_Properties.T_rain_all, ...
        Snow_Properties.rho_snow_init, Snow_Properties.rho_max, ...
        Snow_Properties.k_t, Snow_Properties.k_swe, Snow_Properties.k_D, dt_s);

    errors(2) = mass_balance_error_snow * DEM_raster.cellsize^2 / 1000;

    if k == 1
        max_Hsnow = zeros(size(DEM_raster.Z), 'like', Snow_Properties.H_snow_t);
        max_Hsnow(isnan(DEM_raster.Z)) = NaN;
    end
    max_Hsnow = max(max_Hsnow, Snow_Properties.H_snow_t);

    % Only rainfall and melt reach the surface-water store in this step.
    depths.d_t = depths.d_p + Snow_Properties.M_snow + Snow_Properties.P_rain;
else
    depths.d_t = depths.d_t + BC_States.Eff_Rainfall;
    if k == 1
        Snow_Properties.E_s = 0 * depths.d_p;
        Snow_Properties.SWE_t = 0 * depths.d_p;
        Snow_Properties.H_snow_t = 0 * depths.d_p;
    end
end
