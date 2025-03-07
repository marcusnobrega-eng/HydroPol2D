% Evaporation / Evapotranspiration Module

% Current Reservoir UZ + Surface Storage
S_UZ_ETR = nansum(nansum(C_a.*Soil_Properties.I_t/1000));
S_p_ETR = nansum(nansum(C_a.*depths.d_t/1000));
S_ETR_0 = S_UZ_ETR + S_p_ETR;

% Evapotranspiration Moduel
BC_States.delta_E = zeros(size(elevation,1),size(elevation,2));
Hydro_States.ETR = zeros(size(elevation,1),size(elevation,2));
% Limiting evapotranspiration in soils with low moisture
if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1
    BC_States.delta_E = zeros(size(Hydro_States.ETR));
    % Real evapotranspiration data
    Hydro_States.ETR  = input_evaporation + input_transpiration; % % Real evapotranspiration [mm/day] using the input data
    % Areas with little soil moisture
    Hydro_States.idx_ETR = and(Soil_Properties.I_t <= min_soil_moisture,depths.d_t == 0); % Adopting the minimum depth
    Hydro_States.ETR(isnan(Hydro_States.ETR)) = 0; % Fixing nan values to 0
    Hydro_States.ETR((idx_nan)) = nan; % Cells outside of the domain remain nan
    % Areas with surface water
    idx_open_water =  depths.d_t > 0;
    Hydro_States.ETR(idx_open_water) =  min((depths.d_t(idx_open_water))/(time_step/60/24),Hydro_States.ETP(idx_open_water)); % Limiting to the available surface water (using only evaporation)
    % From the remainder, evaporation occurs
    BC_States.delta_E(idx_open_water) = (time_step/60/24)*min((max(depths.d_t(idx_open_water) - (time_step/60/24)*Hydro_States.ETR(idx_open_water),0))/(time_step/60/24),Hydro_States.Ep(idx_open_water));
    idx_soil_available = and(depths.d_t == 0, Soil_Properties.I_t > min_soil_moisture); % Cells with no surface depth and with available soil moisture
    Hydro_States.ETR(idx_soil_available) =  min((Soil_Properties.I_t(idx_soil_available) - min_soil_moisture(idx_soil_available))/(time_step/60/24),Hydro_States.ETR(idx_soil_available)); % Limiting to the available soil moisture
    Hydro_States.ETR(LULC_Properties.idx_imp) = 0; % No evapotranspiration in impervious areas

    % Soil Matrix Mass Balance
    Soil_Properties.I_t(idx_soil_available) = Soil_Properties.Itp(idx_soil_available) - Hydro_States.ETR(idx_soil_available)*(time_step/60/24); % mm in the soil (taking away ETR prior)

    Hydro_States.ETP = zeros(size(Soil_Properties.I_t)); % Setting ETP to 0 to restart again
    BC_States.delta_E(idx_soil_available) = 0; % mm (already taken from the soil)
elseif flags.flag_ETP == 1
    % Areas with surface water
    idx_open_water =  depths.d_t > 0 & ~idx_nan;
    Hydro_States.ETR  = Hydro_States.ETP; % Real evapotranspiration [mm/24h]
    Hydro_States.ETR(LULC_Properties.idx_imp) = 0; % No evapotranspiration in impervious cells

    % Open waters ETR
    Hydro_States.ETR(idx_open_water) =  min((depths.d_t(idx_open_water))/(time_step/60/24),Hydro_States.ETR(idx_open_water)); % Limiting to the available surface water (using only evaporation)

    % ETR influence on surface depth
    depths.d_t(idx_open_water) = depths.d_t(idx_open_water) - (time_step/60/24)*(Hydro_States.ETR(idx_open_water)); % Considering ETR

    % % From the remainder, evaporation occurs
    BC_States.delta_E(idx_open_water) = (time_step/60/24)*min((depths.d_t(idx_open_water)/(time_step/60/24)),Hydro_States.Ep(idx_open_water)); % mm in a time-step

    % Influence of E on surface depth
    depths.d_t(idx_open_water) = depths.d_t(idx_open_water) - BC_States.delta_E(idx_open_water); % Considering E

    % No water and no soil moisture
    Hydro_States.idx_ETR = Soil_Properties.I_t <= min_soil_moisture; % Adopting the minimum depth
    Hydro_States.ETR(Hydro_States.idx_ETR & ~idx_open_water & ~idx_nan) = 0*min_soil_moisture(Hydro_States.idx_ETR & ~idx_open_water); % No ETR in these cases

    % Areas with soil moisture and no surface depth
    idx_real_ETR = ~Hydro_States.idx_ETR & ~idx_open_water & ~idx_nan;
    Hydro_States.ETR(idx_real_ETR) = min(Hydro_States.ETR(idx_real_ETR), (max(Soil_Properties.I_t(idx_real_ETR) - min_soil_moisture(idx_real_ETR),0))/(time_step/60/24)); % mm/day
    Hydro_States.ETR(LULC_Properties.idx_imp) = 0; % No evapotranspiration in impervious areas
    % Hydro_States.ETR(isnan(Hydro_States.ETP)) = 0;
    Hydro_States.ETR(idx_nan) = nan;

    % Soil Matrix Mass Balance
    Soil_Properties.I_t(~idx_open_water) = Soil_Properties.I_t(~idx_open_water) - Hydro_States.ETR(~idx_open_water)*(time_step/60/24); % mm in the soil (cells that are not open water)
else
    BC_States.delta_E = zeros(size(elevation,1),size(elevation,2)); BC_States.delta_E(idx_nan) = nan;
end

% Final Reservoir UZ + Surface Storage
S_UZ_ETR = nansum(nansum(C_a.*Soil_Properties.I_t/1000));
S_p_ETR = nansum(nansum(C_a.*depths.d_t/1000));
S_ETR_t = S_UZ_ETR + S_p_ETR;

% ETR and E fluxes
flux_E_ETR = (-1)*(nansum(nansum(C_a.*BC_States.delta_E/1000)) +  nansum(nansum(C_a.*Hydro_States.ETR/1000*(time_step/60/24))));

% Mass balance error
dS_ETR = (S_ETR_t - S_ETR_0);
error = dS_ETR - flux_E_ETR;



