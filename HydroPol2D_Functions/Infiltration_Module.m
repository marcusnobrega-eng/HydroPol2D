%% =========================================================================
%  INFILTRATION_MODULE
%  =========================================================================
%  HydroPol2D – Hydrological Infiltration Component
%  Developed by Marcus N. Gomes Jr., PhD
%
%  Updated formulation:
%  -------------------------------------------------------------------------
%  Computes surface infiltration using a Darcy-based formulation coupled
%  with the van Genuchten–Mualem soil hydraulic model.
%
%  Important conceptual update:
%
%     • theta_bucket is used to compute suction pressure head h(theta_bucket)
%     • theta_top is used to compute the conductivity at the bottom of the
%       near-surface layer, K(theta_top)
%     • the soil entry conductivity is an arithmetic interface average
%       between saturated surface conductivity Ks and K(theta_top)
%
%  The near-surface representative conductivity is:
%
%     Ksoil = wK * Ks + (1 - wK) * K(theta_top)
%
%  where:
%
%     wK = Soil_Properties.Ksurf_weight
%
%  Default:
%
%     wK = 0.5
%
%  This represents:
%
%     z = 0      : wet/saturated surface, K = Ks
%     z = Ltop   : bottom of near-surface layer, K = K(theta_top)
%
%  The infiltration capacity is computed as:
%
%     C = Ksoil * [ min(h0 - h(theta_bucket), dh_max) / Ltop + 1 ]
%
%  where theta_top is a within-step near-surface wetting predictor:
%
%     theta_top = theta_bucket + w_wet / Ltop
%
%  and:
%
%     w_wet = min(surface water availability,
%                 remaining unsaturated-zone storage,
%                 top-layer storage deficit)
%
%  This predictor does NOT add water to the soil bucket. The soil bucket is
%  updated only after actual infiltration f is computed.
%
%  -------------------------------------------------------------------------
%  Inputs through struct fields:
%
%  Soil_Properties:
%     I_t            - Current soil storage [mm]
%     Soil_Depth     - Total soil depth [m]
%     theta_i        - Initial water content [-]
%     theta_sat      - Saturated water content [-]
%     theta_r        - Residual water content [-]
%     ksat           - Saturated hydraulic conductivity [mm/h]
%     alpha_vg       - van Genuchten alpha parameter [1/m]
%     n_vg           - van Genuchten n parameter [-]
%     l_vg           - Mualem pore connectivity parameter [-]
%     Ltop           - Effective surface conductance / wetting length [m]
%     dh_max         - Maximum head difference cap [m]
%     Ksurf_weight   - Weight on saturated surface conductivity [-]
%
%  Hydro_States:
%     i_a            - Available infiltration supply [mm/h]
%     f              - Infiltration rate [mm/h]
%
%  depths:
%     d_t            - Surface water depth [mm]
%
%  BC_States:
%     h_t            - Groundwater head [m]
%
%  LULC_Properties:
%     idx_imp        - Impervious mask logical
%
%  flags:
%     flag_infiltration
%     flag_subgrid
%     flag_overbanks
%
%  -------------------------------------------------------------------------
%  Outputs:
%
%     Hydro_States.f        - Infiltration flux [mm/h]
%     Soil_Properties.I_t   - Updated soil storage [mm]
%     depths.d_t            - Updated surface water depth [mm]
%     cumulative_infiltration
%     errors(4)             - Mass balance diagnostic [m3]
%
%  =========================================================================


% -----------------------------
% Area / initial storage accounting
% -----------------------------

if flags.flag_infiltration == 1 || k == 1

    Coarse_Area = Wshed_Properties.Resolution^2;

    % Current unsaturated-zone storage
    S_UZ_inf_0 = nansum(nansum(Coarse_Area .* Soil_Properties.I_t / 1000));

    if flags.flag_subgrid == 1 && flags.flag_overbanks == 1

        Vol = ((Wshed_Properties.Resolution - Wshed_Properties.River_Width) .* ...
                Wshed_Properties.Resolution .* ...
                max((depths.d_t / 1000 - Wshed_Properties.River_Depth), 0) + ...
               (Wshed_Properties.Resolution .* ...
                Wshed_Properties.River_Width .* ...
                depths.d_t / 1000));  % [m3 per cell]

        S_p_inf_0 = nansum(nansum(Vol));

        % Equivalent coarse-cell surface-water depth [mm]
        eff_depth = 1000 * (Vol / Coarse_Area);

    else

        eff_depth = depths.d_t;

        S_p_inf_0 = nansum(nansum(Coarse_Area .* depths.d_t / 1000));

    end

    S_inf_0 = S_UZ_inf_0 + S_p_inf_0;

end


% -----------------------------
% Infiltration
% -----------------------------

if flags.flag_infiltration == 1

    % ---------------------------------------------------------------------
    % Default hydraulic parameters if not provided
    % ---------------------------------------------------------------------

    if ~isfield(Soil_Properties, 'l_vg') || isempty(Soil_Properties.l_vg)
        Soil_Properties.l_vg = 0.5;        % Mualem pore connectivity [-]
    end

    if ~isfield(Soil_Properties, 'Ltop') || isempty(Soil_Properties.Ltop)
        Soil_Properties.Ltop = 0.05;       % Effective surface conductance length [m]
    end

    if ~isfield(Soil_Properties, 'dh_max') || isempty(Soil_Properties.dh_max)
        Soil_Properties.dh_max = 1.0;      % Maximum driving head cap [m]
    end

    if ~isfield(Soil_Properties, 'Ksurf_weight') || isempty(Soil_Properties.Ksurf_weight)
        Soil_Properties.Ksurf_weight = 0.5; % Arithmetic surface-interface weight [-]
    end

    has_layered_soil = isfield(Soil_Properties, 'Layers') && ...
        isstruct(Soil_Properties.Layers) && ...
        isfield(Soil_Properties.Layers, 'near_surface_storage_mm');


    % ---------------------------------------------------------------------
    % Time
    % ---------------------------------------------------------------------

    dt_h = time_step / 60;     % [h]
    dt_s = time_step * 60;     % [s], kept for compatibility if needed elsewhere


    % ---------------------------------------------------------------------
    % Available surface-water supply
    % ---------------------------------------------------------------------

    % Supply rate from effective surface depth [mm/h]
    Hydro_States.i_a = eff_depth ./ max(dt_h, 1e-12);

    % Ponding head at surface [m]
    h0 = max(eff_depth, 0) / 1000;


    % ---------------------------------------------------------------------
    % Unsaturated-zone capacity as function of groundwater depth
    % ---------------------------------------------------------------------

    % Saturated thickness above soil bottom [m]
    GW_Depth = BC_States.h_t - (elevation - Soil_Properties.Soil_Depth);

    % Depth to water table from land surface [m]
    zwt = Soil_Properties.Soil_Depth - GW_Depth;
    zwt = max(zwt, 0);

    if has_layered_soil
        Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
        UZ_max_storage = Soil_Properties.Layers.total_vadose_capacity_mm ./ 1000;
    else
        % Maximum unsaturated-zone storage [m water]
        UZ_max_storage = zwt .* ...
            (Soil_Properties.theta_sat - Soil_Properties.theta_r);
    end

    % Current storage [m water]
    S = max(Soil_Properties.I_t, 0) / 1000;

    % Remaining unsaturated-zone storage capacity [m water]
    if has_layered_soil
        S_rem = layered_soil_remaining_capacity_mm(Soil_Properties, idx_nan) ./ 1000;
    else
        S_rem = max(UZ_max_storage - S, 0);
    end

    % If the water table is at or above the surface, no unsaturated capacity exists
    idx_noUZ = (zwt <= 0) | (UZ_max_storage <= 0);


    % ---------------------------------------------------------------------
    % van Genuchten parameters
    % ---------------------------------------------------------------------

    if has_layered_soil
        theta_r = Soil_Properties.Layers.theta_r_near_surface;
        theta_s = Soil_Properties.Layers.theta_sat_near_surface;
        n = Soil_Properties.Layers.n_vg_near_surface;
        a = Soil_Properties.Layers.alpha_vg_near_surface;   % [1/m]
        Ks = Soil_Properties.Layers.ksat_near_surface / 1000 / 3600;   % [m/s]
        near_thick = max(Soil_Properties.Layers.near_surface_thickness_m, 1e-6);
    else
        theta_r = Soil_Properties.theta_r;
        theta_s = Soil_Properties.theta_sat;
        n = Soil_Properties.n_vg;
        a = Soil_Properties.alpha_vg;   % [1/m]
        Ks = Soil_Properties.ksat / 1000 / 3600;   % [m/s]
        near_thick = [];
    end

    m = 1 - 1 ./ n;

    ell = Soil_Properties.l_vg;     % [-]

    Ltop   = max(Soil_Properties.Ltop, 1e-6);  % [m]
    if has_layered_soil
        Ltop = min(Ltop, near_thick);
        Ltop = max(Ltop, 1e-6);
    end
    dh_max = Soil_Properties.dh_max;           % [m]

    % In this simplified implementation:
    % Lmix is assumed equal to Ltop.
    Lmix = Ltop;


    % ---------------------------------------------------------------------
    % Bucket-average moisture state for pressure head h(theta_bucket)
    % ---------------------------------------------------------------------

    if has_layered_soil
        S_near = max(Soil_Properties.Layers.near_surface_storage_mm, 0) ./ 1000;
        theta_bucket = theta_r + (S_near ./ near_thick);
    else
        zwt_safe = max(zwt, 1e-6);
        theta_bucket = theta_r + (S ./ zwt_safe);
    end

    theta_bucket = max(theta_bucket, theta_r);
    theta_bucket = min(theta_bucket, theta_s);

    Se_bucket = (theta_bucket - theta_r) ./ max(theta_s - theta_r, 1e-12);
    Se_bucket = min(max(Se_bucket, 1e-6), 1);

    % Pressure head h(theta_bucket) [m]
    % Negative under unsaturated conditions.
    h_bucket = -(1 ./ a) .* ((Se_bucket .^ (-1 ./ m) - 1) .^ (1 ./ n));

    % At saturation, pressure head is set to zero.
    h_bucket(Se_bucket >= 0.999999) = 0;


    % ---------------------------------------------------------------------
    % Near-surface moisture predictor for K(theta_top)
    % ---------------------------------------------------------------------
    %
    % This predictor estimates how wet the bottom of the near-surface layer
    % becomes during the current time step if surface water is available.
    %
    % Important:
    %   w_wet is NOT added to Soil_Properties.I_t here.
    %   It is only used to compute theta_top for hydraulic conductivity.
    %
    % Actual mass is added later using:
    %   Soil_Properties.I_t = Soil_Properties.I_t + Hydro_States.f * dt_h
    %
    % ---------------------------------------------------------------------

    % Surface water available before infiltration during this step [m]
    w_avail = max(eff_depth, 0) / 1000;

    if has_layered_soil
        % Water required to fill the represented near-surface storage [m water]
        w_top_deficit = max( ...
            Soil_Properties.Layers.near_surface_capacity_mm - ...
            Soil_Properties.Layers.near_surface_storage_mm, 0) ./ 1000;
    else
        % Water required to wet the top layer from theta_bucket to saturation [m water]
        w_top_deficit = Lmix .* max(theta_s - theta_bucket, 0);
    end

    % Wetting amount used only for near-surface conductivity prediction [m water]
    w_wet = min(w_avail, S_rem);
    w_wet = min(w_wet, w_top_deficit);
    w_wet = max(w_wet, 0);

    % Predicted near-surface water content [-]
    theta_top = theta_bucket + w_wet ./ Lmix;

    theta_top = max(theta_top, theta_r);
    theta_top = min(theta_top, theta_s);

    Se_top = (theta_top - theta_r) ./ max(theta_s - theta_r, 1e-12);
    Se_top = min(max(Se_top, 1e-6), 1);


    % ---------------------------------------------------------------------
    % Mualem–van Genuchten conductivity at bottom of near-surface layer
    % ---------------------------------------------------------------------
    %
    % z = 0:
    %   surface is assumed wet/saturated when water is available,
    %   so K_surface = Ks.
    %
    % z = Ltop:
    %   bottom of near-surface layer has conductivity Ktop = K(theta_top).
    %
    % Effective bucket-scale near-surface conductivity:
    %
    %   Ksoil = wK * Ks + (1 - wK) * Ktop
    %
    % This arithmetic interface conductivity is intentionally NOT a strict
    % harmonic resistance average. It is a bucket-scale closure representing
    % rapid wetting of the surface boundary while the bottom of the layer may
    % remain less conductive.
    %
    % ---------------------------------------------------------------------

    term_top = 1 - Se_top .^ (1 ./ m);
    term_top = min(max(term_top, 0), 1);

    Kr_top = Se_top .^ ell .* (1 - term_top .^ m) .^ 2;
    Kr_top = min(max(Kr_top, 0), 1);

    % Conductivity at bottom of top layer [m/s]
    Ktop = Ks .* Kr_top;

    % Saturated conductivity at the wet surface [m/s]
    Ksurf = Ks + zeros(size(Ktop), 'like', Ktop);

    % Surface-interface arithmetic weight [-]
    wK = Soil_Properties.Ksurf_weight;

    if isscalar(wK)
        wK = wK + zeros(size(Ktop), 'like', Ktop);
    end

    wK = min(max(wK, 0), 1);

    % Use arithmetic surface-interface conductivity only when water is
    % available at the boundary. Otherwise, infiltration is supply-limited
    % to zero anyway, and Ktop is retained for consistency.
    idx_wet_boundary = (Hydro_States.i_a > 0) | (h0 > 0);

    Ksoil = Ktop;

    Ksoil(idx_wet_boundary) = ...
        wK(idx_wet_boundary) .* Ksurf(idx_wet_boundary) + ...
        (1 - wK(idx_wet_boundary)) .* Ktop(idx_wet_boundary);

    % Safety bounds
    Ksoil = min(max(Ksoil, 0), Ksurf);


    % ---------------------------------------------------------------------
    % Darcy infiltration capacity
    % ---------------------------------------------------------------------
    %
    % Conductivity:
    %     Ksoil = arithmetic interface average between Ks and K(theta_top)
    %
    % Pressure head:
    %     h_bucket = h(theta_bucket)
    %
    % Capacity:
    %     C = Ksoil * [min(h0 - h(theta_bucket), dh_max) / Ltop + 1]
    %
    % ---------------------------------------------------------------------

    dh = h0 - h_bucket;     % [m]
    dh = max(dh, 0);
    dh = min(dh, dh_max);

    grad = (dh ./ Ltop) + 1;   % dimensionless

    C = Ksoil .* grad;         % [m/s]

    % Convert capacity to [mm/h]
    C = C * 1000 * 3600;


    % ---------------------------------------------------------------------
    % Final infiltration rate
    % ---------------------------------------------------------------------
    %
    % Limits:
    %   1. available surface-water supply
    %   2. hydraulic infiltration capacity
    %   3. remaining unsaturated-zone storage capacity
    %   4. impervious mask
    %   5. no-unsaturated-zone mask
    %
    % ---------------------------------------------------------------------

    f_store_lim = (S_rem * 1000) ./ max(dt_h, 1e-12);  % [mm/h]

    Hydro_States.f = min(Hydro_States.i_a, C);
    Hydro_States.f = min(Hydro_States.f, f_store_lim);

    % Impervious areas: no infiltration
    Hydro_States.f(LULC_Properties.idx_imp) = 0;

    % No unsaturated-zone capacity: no infiltration
    Hydro_States.f(idx_noUZ) = 0;

    % Safety
    Hydro_States.f(isnan(Hydro_States.f)) = 0;
    Hydro_States.f = max(Hydro_States.f, 0);


    % ---------------------------------------------------------------------
    % Update soil storage and compute infiltrated volume
    % ---------------------------------------------------------------------

    Soil_Properties.I_p = Soil_Properties.I_t;

    % Infiltrated depth per coarse cell [mm]
    inf_volume_cell = Hydro_States.f * dt_h;

    if has_layered_soil
        [Soil_Properties, inf_volume_cell] = add_layered_infiltration( ...
            Soil_Properties, inf_volume_cell, idx_nan);
        Hydro_States.f = inf_volume_cell ./ max(dt_h, 1e-12);
    else
        % Update soil storage [mm]
        Soil_Properties.I_t = Soil_Properties.I_t + inf_volume_cell;

        % Enforce soil storage bound [mm]
        Soil_Properties.I_t = min(Soil_Properties.I_t, UZ_max_storage * 1000);
        Soil_Properties.I_t = max(Soil_Properties.I_t, 0);
    end

    % Total infiltrated volume over domain [m3], for diagnostics
    inf_volume = 1 / 1000 * ...
        nansum(nansum(inf_volume_cell .* Coarse_Area));

    % Remove infiltrated depth from surface water depth
    depths.d_t = depths.d_t - inf_volume_cell .* Coarse_Area ./ C_a;
    depths.d_t(depths.d_t < 0) = 0;

end


% -----------------------------
% Subgrid inbank/overbank transitions
% -----------------------------

if flags.flag_subgrid == 1 && flags.flag_overbanks

    % Inbank -> Overbank
    idx = depths.d_t / 1000 > Wshed_Properties.River_Depth & ...
          depths.d_p / 1000 <= Wshed_Properties.River_Depth & ...
          (Wshed_Properties.River_Width > 0);

    if sum(sum(idx)) > 0

        depths.d_t(idx) = inbank_to_overbank( ...
            Wshed_Properties.River_Width(idx), ...
            Wshed_Properties.River_Depth(idx), ...
            depths.d_t(idx) / 1000, ...
            Wshed_Properties.Resolution);

        C_a(idx) = Wshed_Properties.Resolution^2;

    end


    % Overbank -> Inbank
    idx = depths.d_t / 1000 <= Wshed_Properties.River_Depth & ...
          depths.d_p / 1000 >= Wshed_Properties.River_Depth & ...
          (Wshed_Properties.River_Width > 0);

    if sum(sum(idx)) > 0

        depths.d_t(idx) = overbank_to_inbank( ...
            Wshed_Properties.River_Width(idx), ...
            Wshed_Properties.River_Depth(idx), ...
            depths.d_t(idx) / 1000, ...
            Wshed_Properties.Resolution);

        C_a(idx) = Wshed_Properties.River_Width(idx) * ...
                   Wshed_Properties.Resolution;

    end

end


% -----------------------------
% Final storage accounting + error
% -----------------------------

if flags.flag_infiltration == 1 || k == 1

    S_UZ_inf_t = nansum(nansum(Coarse_Area .* Soil_Properties.I_t / 1000));

    if flags.flag_subgrid == 1 && flags.flag_overbanks == 1

        S_p_inf_t = ...
            nansum(nansum( ...
                (Wshed_Properties.Resolution - Wshed_Properties.River_Width) .* ...
                 Wshed_Properties.Resolution .* ...
                 max((depths.d_t / 1000 - Wshed_Properties.River_Depth), 0))) + ...
            nansum(nansum( ...
                 Wshed_Properties.Resolution .* ...
                 Wshed_Properties.River_Width .* ...
                 depths.d_t / 1000));

    else

        S_p_inf_t = nansum(nansum(Coarse_Area .* depths.d_t / 1000));

    end

    S_inf_t = S_UZ_inf_t + S_p_inf_t;

    dS_inf = S_inf_t - S_inf_0;

    errors(4) = dS_inf;  % [m3]

end


% -----------------------------
% Cumulative infiltration
% -----------------------------

if k == 1

    cumulative_infiltration = ones(size(DEM_raster.Z));
    cumulative_infiltration(isnan(cumulative_infiltration)) = nan;

end

if flags.flag_infiltration == 1

    % Hydro_States.f is [mm/h]
    % time_step/60 is [h]
    % Therefore:
    %   Hydro_States.f * (time_step/60) gives [mm]
    %   division by 1000 gives [m]
    cumulative_infiltration = cumulative_infiltration + ...
        Hydro_States.f * (time_step / 60) / 1000;

else

    Hydro_States.f = 0 * elevation;
    inf_volume = 0 * elevation;

end


%% Local Functions

function [h] = inbank_to_overbank(B, H, h_p, R)

    h = 1000 * (B .* H + (R - B) .* h_p) ./ R;

end


function [h] = overbank_to_inbank(B, H, h_p, R)

    h = 1000 * (B .* h_p + (R - B) .* (H - h_p)) ./ B;

end
