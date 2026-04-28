%% ═══════════════════════════════════════════════════════════════════════
%  Module: evapotranspirationModule
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Updated: 03/09/2026
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Compute actual evapotranspiration (ETR) and open-water evaporation
%      while updating soil moisture and surface water storage consistently.
%
%  ➤ Main improvement in this revision:
%      The previous version only computed ETR when
%          flags.flag_ETP == 1  AND  flags.flag_input_ETP_map == 1
%      which meant that internally computed ETP/Ep cases had no actual ET sink.
%
%      This revised version supports BOTH cases:
%
%      1) Input ET maps available:
%         - Actual ET demand is prescribed by:
%               input_evaporation + input_transpiration
%
%      2) Internal meteorological ETP mode:
%         - Actual soil ET demand is derived from:
%               Hydro_States.ETP
%         - Open-water evaporation is driven by:
%               Hydro_States.Ep
%
%  ➤ Important mass-balance fix in this revision:
%      The ET masks and ET limitations are now computed from FROZEN pre-ET
%      storage states:
%          • I_t_before_ET
%          • d_t_before_ET
%
%      This prevents mismatch between:
%          • the storage basis used to compute the ET fluxes, and
%          • the storage basis used to update the states.
%
%      In other words:
%      - surface evaporation is computed from d_t_before_ET and removed from
%        depths.d_t using exactly the same quantity BC_States.delta_E
%      - soil ETR is computed from I_t_before_ET and removed from
%        Soil_Properties.I_t using exactly the same quantity
%        Hydro_States.ETR * dt_days
%
%  ➤ Units:
%      • Hydro_States.ETP, Hydro_States.Ep, Hydro_States.ETR : mm/day
%      • BC_States.delta_E                                   : mm over current time step
%      • Soil_Properties.I_t                                 : mm
%      • depths.d_t                                          : mm
% ═══════════════════════════════════════════════════════════════════════


%% ---------------------------------------------------------------------
%  1) Initial storage before ET update
% ----------------------------------------------------------------------
if flags.flag_ETP == 1 || k == 1
S_UZ_ETR = nansum(nansum(C_a .* Soil_Properties.I_t ./ 1000));
S_p_ETR  = nansum(nansum(C_a .* depths.d_t          ./ 1000));
S_ETR_0  = S_UZ_ETR + S_p_ETR;
S_ETR_t = S_ETR_0;


%---------------------------------------------------------------------
%  2) Initialize output fields
% ----------------------------------------------------------------------
BC_States.delta_E = zeros(size(elevation,1), size(elevation,2));   % mm over current step
Hydro_States.ETR  = zeros(size(elevation,1), size(elevation,2));   % mm/day

BC_States.delta_E(idx_nan) = nan;
Hydro_States.ETR(idx_nan)  = nan;

% Time step in days
dt_days = time_step / 60 / 24;

end
%% ---------------------------------------------------------------------
%  3) Only process ET if ET is enabled
% ----------------------------------------------------------------------
if flags.flag_ETP == 1

    % ------------------------------------------------------------------
    % Freeze pre-ET storages for internal consistency
    % ------------------------------------------------------------------
    I_t_before_ET = Soil_Properties.I_t;
    d_t_before_ET = depths.d_t;


    % ------------------------------------------------------------------
    % Define masks FROM PRE-ET STORAGE STATES
    % ------------------------------------------------------------------
    idx_open_water = d_t_before_ET > 0 & ~idx_nan;

    idx_soil_available = d_t_before_ET == 0 & ...
        I_t_before_ET > min_soil_moisture & ...
        ~idx_nan;

    idx_dry = d_t_before_ET == 0 & ...
        I_t_before_ET <= min_soil_moisture & ...
        ~idx_nan;

    % ------------------------------------------------------------------
    % Build actual ET demand field for soil cells
    % ------------------------------------------------------------------
    if flags.flag_input_ETP_map == 1
        % CASE A: Actual ET prescribed by input maps
        ETR_demand = input_evaporation + input_transpiration;   % mm/day
    else
        % CASE B: Actual ET derived from internally computed ETP
        ETR_demand = Hydro_States.ETP;                          % mm/day
    end

    % Clean demand field
    ETR_demand(isnan(ETR_demand)) = 0;
    ETR_demand(idx_nan) = nan;

    % Impervious cells do not extract ET from soil
    ETR_demand(LULC_Properties.idx_imp) = 0;

    % Initialize valid cells explicitly
    Hydro_States.ETR(~idx_nan)  = 0;
    BC_States.delta_E(~idx_nan) = 0;

    % ------------------------------------------------------------------
    % 3.1) OPEN WATER CELLS
    %      Evaporation removed from surface storage only
    % ------------------------------------------------------------------
    if any(idx_open_water(:))

        % Available equivalent evaporation rate from ponded storage [mm/day]
        max_open_water_evap_rate = d_t_before_ET(idx_open_water) ./ max(dt_days, eps);

        % Actual open-water evaporation depth over this time step [mm]
        BC_States.delta_E(idx_open_water) = dt_days .* min( ...
            max_open_water_evap_rate, ...
            Hydro_States.Ep(idx_open_water));

        % This flux is from surface water, not from soil ETR
        Hydro_States.ETR(idx_open_water) = 0;

        % Update surface storage using the exact same removed depth
        depths.d_t(idx_open_water) = d_t_before_ET(idx_open_water) ...
            - BC_States.delta_E(idx_open_water);

        depths.d_t(idx_open_water) = max(depths.d_t(idx_open_water), 0);
    end

    % ------------------------------------------------------------------
    % 3.2) SOIL CELLS WITHOUT PONDED WATER
    %      ETR removed from soil moisture only
    % ------------------------------------------------------------------
    if any(idx_soil_available(:))

        % Maximum extractable ET rate based on PRE-ET soil storage [mm/day]
        max_soil_ET_rate = (I_t_before_ET(idx_soil_available) ...
            - min_soil_moisture(idx_soil_available)) ./ max(dt_days, eps);

        max_soil_ET_rate = max(max_soil_ET_rate, 0);

        % Actual ET rate [mm/day]
        Hydro_States.ETR(idx_soil_available) = min( ...
            ETR_demand(idx_soil_available), ...
            max_soil_ET_rate);

        % Update soil storage using the exact same pre-ET storage basis
        Soil_Properties.I_t(idx_soil_available) = ...
            I_t_before_ET(idx_soil_available) ...
            - Hydro_States.ETR(idx_soil_available) .* dt_days;

        % Enforce lower bound
        Soil_Properties.I_t(idx_soil_available) = max( ...
            Soil_Properties.I_t(idx_soil_available), ...
            min_soil_moisture(idx_soil_available));

        % No surface evaporation term for these cells
        BC_States.delta_E(idx_soil_available) = 0;
    end

    % ------------------------------------------------------------------
    % 3.3) DRY CELLS
    % ------------------------------------------------------------------
    Hydro_States.idx_ETR = idx_dry;

    if any(idx_dry(:))
        Hydro_States.ETR(idx_dry)  = 0;
        BC_States.delta_E(idx_dry) = 0;
    end

    % ------------------------------------------------------------------
    % 3.4) Final cleanup in invalid cells
    % ------------------------------------------------------------------
    Hydro_States.ETR(idx_nan)  = nan;
    BC_States.delta_E(idx_nan) = nan;

    % ------------------------------------------------------------------
    % 3.5) Optional convention from original code
    % ------------------------------------------------------------------
    if flags.flag_input_ETP_map == 1
        Hydro_States.ETP = zeros(size(Soil_Properties.I_t));
        Hydro_States.ETP(idx_nan) = nan;
    end
end


%% ---------------------------------------------------------------------
%  4) Safety check on soil moisture
% ----------------------------------------------------------------------
% if nanmin(nanmin(Soil_Properties.I_t)) <= 0
%     Soil_Properties.I_t(Soil_Properties.I_t <= 0) = 1e-16;
% end


%% ---------------------------------------------------------------------
%  5) Final storage after ET update
% ----------------------------------------------------------------------

if flags.flag_ETP == 1
    S_UZ_ETR = nansum(nansum(C_a .* Soil_Properties.I_t ./ 1000));
    S_p_ETR  = nansum(nansum(C_a .* depths.d_t          ./ 1000));
    S_ETR_t  = S_UZ_ETR + S_p_ETR;
end

%% ---------------------------------------------------------------------
%  6) Flux accounting
% ----------------------------------------------------------------------
% Surface evaporation:
%   BC_States.delta_E is already the removed depth over the current step [mm]
%
% Soil evapotranspiration:
%   Hydro_States.ETR is a rate [mm/day], so convert using dt_days

if flags.flag_ETP == 1 || k == 1
removed_surface_m3 = nansum(nansum(C_a .* BC_States.delta_E ./ 1000));
removed_soil_m3    = nansum(nansum(C_a .* Hydro_States.ETR .* dt_days ./ 1000));

% Keep original sign convention for compatibility
flux_E_ETR = (-1) * (removed_surface_m3 + removed_soil_m3);


%% ---------------------------------------------------------------------
%  7) Mass-balance error
% ----------------------------------------------------------------------
dS_ETR = S_ETR_t - S_ETR_0;
error  = dS_ETR - flux_E_ETR;


end

if k == 1 && flags.flag_ETP ~= 1
    errors(3) = 0;
else
    errors(3) = error;
end


% Optional debug version, easier to inspect:
% errors(3) = dS_ETR + removed_surface_m3 + removed_soil_m3;