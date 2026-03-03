%% =========================================================================
%  INFILTRATION_MODULE
%  =========================================================================
%  HydroPol2D – Hydrological Infiltration Component
%  Developed by Marcus N. Gomes Jr., PhD

%  Description:
%  -------------------------------------------------------------------------
%  Computes surface infiltration using a Darcy-based formulation coupled
%  with the van Genuchten–Mualem soil hydraulic model. The infiltration
%  capacity is controlled by:
%
%     • Representative soil moisture state (bucket storage)
%     • van Genuchten pressure head h(θ)
%     • Mualem hydraulic conductivity K(θ)
%     • Surface ponding head (h0)
%     • Effective top conductance length (Ltop)
%     • Maximum driving head cap (dh_max) for numerical stability
%     • Unsaturated storage capacity limited by groundwater depth
%
%  The infiltration flux is computed as:
%
%     C = K(θ) * [ (min(h0 - h(θ), dh_max) / Ltop) + 1 ]
%
%  Final infiltration is limited by:
%     (i) available surface water supply,
%    (ii) soil storage capacity remaining,
%   (iii) impervious area mask.
%
%  This formulation replaces the previous Green–Ampt implementation
%  with a reduced Richards-consistent flux approximation while preserving:
%
%     • Mass conservation
%     • HydroPol2D output structure
%     • Subgrid inbank/overbank corrections
%     • Cumulative infiltration accounting
%
%  -------------------------------------------------------------------------
%  Inputs (through struct fields):
%
%  Soil_Properties:
%     I_t           - Current soil storage [mm]
%     Soil_Depth    - Total soil depth [m]
%     theta_i        - Initial water content [-]
%     theta_sat      - Saturated water content [-]
%     theta_r        - Residual water content [-]
%     ksat          - Saturated hydraulic conductivity [mm/h]
%     alpha_vg      - van Genuchten α parameter [1/m]
%     n_vg          - van Genuchten n parameter [-]
%     l_vg          - Mualem pore connectivity parameter [-]
%     Ltop          - Effective surface conductance length [m]
%     dh_max        - Maximum head difference cap [m]
%
%  Hydro_States:
%     i_a           - Available infiltration supply [mm/h]
%     f             - Infiltration rate [mm/h]
%
%  depths:
%     d_t           - Surface water depth [mm]
%
%  BC_States:
%     h_t           - Groundwater head [m]
%
%  LULC_Properties:
%     idx_imp       - Impervious mask (logical)
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
%     errors(4)             - Mass balance diagnostic [m³]
%
%  -------------------------------------------------------------------------
%  Notes:
%
%  • When groundwater reaches the surface (zwt → 0), infiltration is set
%    to zero and excess water is routed as surface flow.
%
%  • The parameter Ltop represents an effective hydraulic resistance
%    thickness of the near-surface layer (typically 0.05–0.10 m).
%
%  • The parameter dh_max prevents unrealistically large suction-driven
%    gradients under very dry conditions.
%
%  • This module maintains compatibility with existing HydroPol2D
%    groundwater and routing components.
%
%  -------------------------------------------------------------------------
%  Version:
%     v2.0 – Darcy–van Genuchten infiltration formulation
%
%
%  =========================================================================
% -----------------------------
% Area / initial storage accounting 
% -----------------------------
Coarse_Area = Wshed_Properties.Resolution^2;

% Current UZ + depth Storage
S_UZ_inf_0 = nansum(nansum(Coarse_Area .* Soil_Properties.I_t/1000));

if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    Vol = ((Wshed_Properties.Resolution - Wshed_Properties.River_Width).*Wshed_Properties.Resolution.*max((depths.d_t/1000 - Wshed_Properties.River_Depth),0) + ...
           (Wshed_Properties.Resolution.*Wshed_Properties.River_Width.*depths.d_t/1000)); % m3 per cell

    S_p_inf_0 = nansum(nansum(Vol));
    eff_depth = 1000*(Vol / Coarse_Area); % Equivalent coarse depth (mm)
else
    eff_depth = depths.d_t;
    S_p_inf_0 = nansum(nansum(Coarse_Area .* depths.d_t/1000));
end

S_inf_0 = S_UZ_inf_0 + S_p_inf_0;

% -----------------------------
% Infiltration
% -----------------------------
if flags.flag_infiltration == 1
    
    % ---------------------------------------------------------------------
    % Default hydraulic parameters (only if not provided in struct)
    % ---------------------------------------------------------------------
    if ~isfield(Soil_Properties,'l_vg') || isempty(Soil_Properties.l_vg)
        Soil_Properties.l_vg = 0.5;        % Mualem pore connectivity [-]
    end

    if ~isfield(Soil_Properties,'Ltop') || isempty(Soil_Properties.Ltop)
        Soil_Properties.Ltop = 0.05;       % Effective surface conductance length [m]
    end

    if ~isfield(Soil_Properties,'dh_max') || isempty(Soil_Properties.dh_max)
        Soil_Properties.dh_max = 1.0;      % Max driving head cap [m]
    end

    % Time
    dt_h = (time_step/60);     % hours
    dt_s = time_step * 60;     % seconds

    % Inflow rate from ponded depth (mm/h) – consistent with your existing definition
    Hydro_States.i_a = eff_depth ./ dt_h;

    % Ponding head at surface (m)
    h0 = max(eff_depth, 0) / 1000;

    % ---------------------------------------------------------
    % Unsaturated-zone capacity as function of groundwater depth
    % (same idea you already use when baseflow is on, but always applied)
    % ---------------------------------------------------------
    GW_Depth = (BC_States.h_t - (elevation - Soil_Properties.Soil_Depth));  % [m] saturated thickness
    zwt = Soil_Properties.Soil_Depth - GW_Depth;                           % [m] depth to WT from surface
    zwt = max(zwt, 0);

    % Max UZ storage [m water]
    UZ_max_storage = zwt .* (Soil_Properties.theta_sat - Soil_Properties.theta_r);

    % Current storage [m water]
    S = max(Soil_Properties.I_t, 0) / 1000;

    % Remaining storage capacity [m water]
    S_rem = max(UZ_max_storage - S, 0);

    % If fully saturated / no unsat thickness -> no infiltration
    idx_noUZ = (zwt <= 0) | (UZ_max_storage <= 0);

    % ---------------------------------------------------------
    % Map bucket storage -> representative theta
    % theta = theta_r + S/zwt (bounded)
    % ---------------------------------------------------------
    zwt_safe = max(zwt, 1e-6);
    theta = Soil_Properties.theta_r + (S ./ zwt_safe);

    theta = max(theta, Soil_Properties.theta_r);
    theta = min(theta, Soil_Properties.theta_sat);

    % ---------------------------------------------------------
    % van Genuchten: h(theta)
    % Required fields:
    %   Soil_Properties.theta_r
    %   Soil_Properties.alpha_vg [1/m]
    %   Soil_Properties.n_vg
    % m = 1 - 1/n
    % ---------------------------------------------------------
    theta_r = Soil_Properties.theta_r;
    theta_s = Soil_Properties.theta_sat;

    Se = (theta - theta_r) ./ max(theta_s - theta_r, 1e-12);
    Se = min(max(Se, 1e-6), 1);

    n = Soil_Properties.n_vg;
    m = 1 - 1./n;
    a = Soil_Properties.alpha_vg;   % [1/m]

    % Pressure head h(theta) [m], negative unsaturated
    h_soil = -(1./a) .* ((Se.^(-1./m) - 1).^(1./n));
    h_soil(Se >= 0.999999) = 0;

    % ---------------------------------------------------------
    % Mualem–van Genuchten: K(theta)
    % Required fields:
    %   Soil_Properties.ksat [mm/h]
    %   Soil_Properties.l_vg [-]  (typical 0.5)
    % ---------------------------------------------------------
    Ks = Soil_Properties.ksat / 1000 / 3600;   % [m/s]
    ell = Soil_Properties.l_vg;                % [-]

    term = (1 - Se.^(1./m));
    Kr = Se.^ell .* (1 - term.^m).^2;
    Kr = min(max(Kr, 0), 1);

    Ksoil = Ks .* Kr;                          % [m/s]

    % ---------------------------------------------------------
    % Darcy infiltration capacity with Ltop and dh_max
    % fcap = Ksoil * ( (min(h0 - h_soil, dh_max))/Ltop + 1 )
    % ---------------------------------------------------------
    Ltop   = Soil_Properties.Ltop;     % [m] e.g., 0.05
    dh_max = Soil_Properties.dh_max;   % [m] e.g., 1.0

    dh = h0 - h_soil;                  % [m]
    dh = max(dh, 0);
    dh = min(dh, dh_max);

    grad = (dh ./ max(Ltop, 1e-6)) + 1;   % dimensionless
    C = Ksoil .* grad;            % [m/s]

    % Convert to mm/h
    C = C * 1000 * 3600;      % [mm/h]

    % ---------------------------------------------------------
    % Final infiltration: supply limit + capacity limit + storage limit
    % Supply limit is Hydro_States.i_a by definition.
    % Storage limit: cannot exceed remaining capacity this step.
    % ---------------------------------------------------------
    f_store_lim = (S_rem * 1000) ./ max(dt_h, 1e-12);  % [mm/h]

    Hydro_States.f = min(Hydro_States.i_a, C);
    Hydro_States.f = min(Hydro_States.f, f_store_lim);

    % Impervious areas: no infiltration
    Hydro_States.f(LULC_Properties.idx_imp) = 0;

    % No UZ capacity (WT at/above surface): no infiltration
    Hydro_States.f(idx_noUZ) = 0;

    % Safety
    Hydro_States.f(isnan(Hydro_States.f)) = 0;
    Hydro_States.f = max(Hydro_States.f, 0);

    % ---------------------------------------------------------
    % Update soil storage and compute infiltrated volume
    % ---------------------------------------------------------
    Soil_Properties.I_p = Soil_Properties.I_t;  % keep previous (used elsewhere sometimes)

    % Infiltrated depth per coarse cell [mm]
    inf_volume_cell = Hydro_States.f * dt_h;

    % Update soil storage [mm]
    Soil_Properties.I_t = Soil_Properties.I_t + inf_volume_cell;

    % Enforce soil storage bound [mm]
    Soil_Properties.I_t = min(Soil_Properties.I_t, UZ_max_storage * 1000);
    Soil_Properties.I_t = max(Soil_Properties.I_t, 0);

    % Total infiltrated volume (domain) [m3] (for diagnostics only)
    inf_volume = 1/1000 * nansum(nansum((inf_volume_cell .* Coarse_Area)));

    % Remove infiltrated depth from surface water depth (your original area correction)
    depths.d_t = depths.d_t - inf_volume_cell .* Coarse_Area ./ C_a;
    depths.d_t(depths.d_t < 0) = 0;

end

% -----------------------------
% Subgrid inbank/overbank transitions (unchanged)
% -----------------------------
if flags.flag_subgrid == 1 && flags.flag_overbanks
    % Inbank -> Overbank
    idx = depths.d_t/1000 > Wshed_Properties.River_Depth & depths.d_p/1000 <= Wshed_Properties.River_Depth & (Wshed_Properties.River_Width > 0);
    if sum(sum(idx)) > 0
        depths.d_t(idx) = inbank_to_overbank(Wshed_Properties.River_Width(idx), Wshed_Properties.River_Depth(idx), depths.d_t(idx)/1000, Wshed_Properties.Resolution);
        C_a(idx) = Wshed_Properties.Resolution^2;
    end

    % Overbank -> Inbank
    idx = depths.d_t/1000 <= Wshed_Properties.River_Depth & depths.d_p/1000 >= Wshed_Properties.River_Depth & (Wshed_Properties.River_Width > 0);
    if sum(sum(idx)) > 0
        depths.d_t(idx) = overbank_to_inbank(Wshed_Properties.River_Width(idx), Wshed_Properties.River_Depth(idx), depths.d_t(idx)/1000, Wshed_Properties.Resolution);
        C_a(idx) = Wshed_Properties.River_Width(idx) * Wshed_Properties.Resolution;
    end
end

% -----------------------------
% Final storage accounting + error (unchanged)
% -----------------------------
S_UZ_inf_t = nansum(nansum(Coarse_Area .* Soil_Properties.I_t/1000));

if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    S_p_inf_t = nansum(nansum((Wshed_Properties.Resolution - Wshed_Properties.River_Width).*Wshed_Properties.Resolution.*max((depths.d_t/1000 - Wshed_Properties.River_Depth),0))) + ...
               nansum(nansum(Wshed_Properties.Resolution .* Wshed_Properties.River_Width .* depths.d_t/1000));
else
    S_p_inf_t = nansum(nansum(Coarse_Area .* depths.d_t/1000));
end

S_inf_t = S_UZ_inf_t + S_p_inf_t;

dS_inf = (S_inf_t - S_inf_0);
errors(4) = dS_inf; % m3

% -----------------------------
% Cumulative infiltration (unchanged interface)
% -----------------------------
if k == 1
    cumulative_infiltration = ones(size(DEM_raster.Z));
    cumulative_infiltration(isnan(cumulative_infiltration)) = nan;
end

if flags.flag_infiltration == 1
    cumulative_infiltration = cumulative_infiltration + Hydro_States.f/1000/3600 * (time_step/60);
else
    Hydro_States.f = 0 * elevation;
    inf_volume = 0 * elevation;
end

%%%% Local Functions %%%%

function [h] = inbank_to_overbank(B,H,h_p,R)
    h = 1000*(B.*H + (R - B).*h_p) ./ R;
end

function [h] = overbank_to_inbank(B,H,h_p,R)
    h = 1000*(B.*h_p + (R - B).*(H - h_p)) ./ B;
end