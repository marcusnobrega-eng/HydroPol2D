%% ═══════════════════════════════════════════════════════════════════════
%  Module: evapotranspirationModule
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 02/24/2025  (adjust as needed)
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Estimate evaporation and evapotranspiration (ETR) by updating the soil 
%      and surface water storage. This module calculates the water storage in 
%      the unsaturated zone (UZ) and surface, applies ETR corrections based on 
%      soil moisture availability and open water conditions, updates the soil 
%      moisture mass balance, and computes the mass balance error.
%
%  ➤ Inputs:
%      • elevation         - Elevation matrix (used for domain size and NaN mask)
%      • C_a               - Cell area (m²) for converting mm to m³
%      • Soil_Properties   - Structure with soil moisture properties:
%             - I_t: Soil moisture storage (mm)
%             - Itp: Soil moisture storage prior to ETR update (mm)
%             - (min_soil_moisture is used as a lower bound)
%      • depths            - Structure containing:
%             - d_t: Surface water depth (mm)
%             - d_p: Ponded water depth (mm)
%      • flags             - Structure with control flags:
%             - flag_ETP: Activates evapotranspiration processes
%             - flag_input_ETP_map: Use input-based evapotranspiration data
%      • Hydro_States      - Structure with hydrological states:
%             - ETR: Evapotranspiration (mm/day) [to be updated]
%             - ETP: Potential evapotranspiration (mm/day)
%             - Ep: Evaporation component (mm/day)
%             - idx_ETR: Logical index for cells with minimal soil moisture
%      • BC_States         - Structure with boundary conditions:
%             - delta_E: Evaporation flux (mm/day) [to be updated]
%      • input_evaporation & input_transpiration:
%             - Real evapotranspiration data (mm/day)
%      • LULC_Properties   - Structure with land use/land cover masks:
%             - idx_imp: Logical mask for impervious areas
%      • time_step         - Time step (minutes)
%      • DEM_raster        - Digital Elevation Model structure (provides cellsize)
%      • idx_nan           - Logical mask for cells with invalid/missing data
%
%  ➤ Outputs:
%      • Updated Soil_Properties.I_t (soil moisture in the soil matrix)
%      • Updated depths.d_t (surface water depth, mm)
%      • Hydro_States.ETR  - Computed evapotranspiration (mm/day)
%      • BC_States.delta_E - Evaporation flux (mm/day)
%      • S_ETR_t           - Final total water storage in UZ + surface (m³)
%      • flux_E_ETR        - Net evaporation flux (m³)
%      • errors(3)         - Mass balance error for ETR (m³), adjusted by DEM cell area
%
%  ➤ Process Summary:
%      1. Compute initial storage S_ETR_0 from soil and surface water.
%      2. Based on flag_ETP and flag_input_ETP_map, process:
%           - Real evapotranspiration using input data,
%           - Open water limitations,
%           - Soil moisture restrictions.
%      3. Update soil moisture mass balance accordingly.
%      4. Compute final storage S_ETR_t, evaporation flux, and mass balance error.
%
%  ➤ Notes:
%      • In impervious areas (LULC_Properties.idx_imp), ETR is set to zero.
%      • Cells with insufficient soil moisture or open water have special treatment.
%      • The module calculates a mass balance error to verify conservation of water,
%        with the error stored in errors(3).
%      • Adjustments are applied to ensure that soil moisture does not fall below 
%        a small threshold (set to 1e-16 when needed).
% ═══════════════════════════════════════════════════════════════════════


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

% Check on I_t values
% if nanmin(nanmin(Soil_Properties.I_t)) <= 0
%     Soil_Properties.I_t(Soil_Properties.I_t <= 0) = 1e-16;
% end

% Final Reservoir UZ + Surface Storage
S_UZ_ETR = nansum(nansum(C_a.*Soil_Properties.I_t/1000));
S_p_ETR = nansum(nansum(C_a.*depths.d_t/1000));
S_ETR_t = S_UZ_ETR + S_p_ETR;

% ETR and E fluxes
flux_E_ETR = (-1)*(nansum(nansum(C_a.*BC_States.delta_E/1000)) +  nansum(nansum(C_a.*Hydro_States.ETR/1000*(time_step/60/24))));

% Mass balance error
dS_ETR = (S_ETR_t - S_ETR_0);
error = dS_ETR - flux_E_ETR;

errors(3) = error;



