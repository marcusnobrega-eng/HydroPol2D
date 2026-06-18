%% =========================================================================
%  GROUNDWATER MODULE (UPDATED PHYSICS)
%  =========================================================================
%  HydroPol2D – Recharge + Groundwater Flow Coupling
%
%  Description:
%  -------------------------------------------------------------------------
%  This module couples the surface–vadose–groundwater system by:
%    1) Computing unsaturated-zone (UZ) maximum storage as a function of
%       dynamic groundwater depth (water table position)
%    2) Computing recharge with a reduced Richards-consistent Darcy closure
%       (van Genuchten pressure head h(θ) + Mualem K(θ))
%    3) Routing groundwater using the explicit 2D Boussinesq solver
%    4) Updating surface water depth with exfiltration and river exchanges
%    5) Ensuring mass consistency when WT rises and UZ capacity shrinks
%
%  Key Interfaces (units):
%    Hydro_States.f      [mm/h]  infiltration flux from surface to vadose
%    Soil_Properties.I_t [mm]    vadose storage (bucket water depth)
%    recharge_rate       [m/s]   flux to groundwater (can be negative = upflux)
%    depths.d_t          [mm]    surface water depth
%    BC_States.h_t       [m]     groundwater head
%
%  Notes:
%  -------------------------------------------------------------------------
%  • This replaces the older "wetting front" consistency update.
%  • If the water table rises and UZ capacity decreases, any vadose storage
%    exceeding the new capacity is converted to surface water (saturation-
%    excess), preserving water mass without double counting groundwater.
%  • `simulate_groundwater_recharge` is assumed updated to the Darcy–vG
%    recharge closure, while preserving its input/output signature.
%
%  Updated:
%    02/21/2026 – vG/Mualem Darcy recharge + consistent WT–UZ coupling
% =========================================================================

%% ------------------------------------------------------------------------
% 1) INITIALIZE CURRENT RECHARGE MEMORY (mm)
% -------------------------------------------------------------------------
if k == 1 || ~exist('cumulative_recharge', 'var') || isempty(cumulative_recharge)
    cumulative_recharge = zeros(size(DEM_raster.Z));
    cumulative_recharge(isnan(DEM_raster.Z)) = nan;
end

%% ------------------------------------------------------------------------
% INITIALIZE SHARED ARRAYS
% -------------------------------------------------------------------------
if k == 1 || ~exist('zero_matrix', 'var') || isempty(zero_matrix)
    z_terrain = Elevation_Properties.elevation_cell;
    zero_matrix = zeros(size(Elevation_Properties.elevation_cell));
    zero_matrix(isnan(Elevation_Properties.elevation_cell)) = nan;
end

if k == 1 || ~exist('max_GW_depth', 'var') || isempty(max_GW_depth)
    max_GW_depth = zeros(size(DEM_raster.Z));
    max_GW_depth(isnan(DEM_raster.Z)) = nan;
end

dt_s = time_step * 60;                    % [s]

recharge_rate = zero_matrix;              % [m/s]
q_exf   = zero_matrix;                    % [m/s] exfiltration to surface
q_river = zero_matrix;                    % [m/s] exchange with river
error   = 0;

has_layered_soil = isfield(Soil_Properties, 'Layers') && ...
    isstruct(Soil_Properties.Layers) && ...
    isfield(Soil_Properties.Layers, 'near_surface_storage_mm');

%% ------------------------------------------------------------------------
% EXIT EARLY IF GROUNDWATER / RECHARGE MODELING IS DISABLED
% -------------------------------------------------------------------------
if flags.flag_groundwater_modeling == 0
    errors(5) = 0;
    return
end

if ~exist('GW_States', 'var') || ~isstruct(GW_States)
    GW_States = struct();
end
[flags, GW_States] = initialize_groundwater_scheduler( ...
    flags, GW_States, zero_matrix, idx_nan);
GW_States.cell_area_m2 = ones(size(zero_matrix), 'like', zero_matrix) .* Wshed_Properties.cell_area;
GW_States.cell_area_m2(idx_nan) = nan;

%% ------------------------------------------------------------------------
% 2) COMPUTE WATER TABLE POSITION AND UNSATURATED ZONE STORAGE CAPACITY
% -------------------------------------------------------------------------
% Saturated thickness above bedrock (as used in your original code) [m]
GW_Depth = BC_States.h_t - ...
    (z_terrain - Soil_Properties.Soil_Depth);  % [m]

% Depth to water table from surface [m]
zwt = Soil_Properties.Soil_Depth - GW_Depth;  % [m]
zwt = max(zwt, 0);
zwt = min(zwt, Soil_Properties.Soil_Depth);

if has_layered_soil
    Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
    UZ_max_storage = Soil_Properties.Layers.total_vadose_capacity_mm ./ 1000;
else
    % Maximum UZ storage [m water]
    % Soil_Properties.I_t is above-residual storage, so capacity must be
    % theta_sat - theta_r, not theta_sat - theta_i.
    UZ_max_storage = zwt .* ...
        (Soil_Properties.theta_sat - Soil_Properties.theta_r);
end

%% ------------------------------------------------------------------------
% 3) COMPUTE RECHARGE
% -------------------------------------------------------------------------
% IMPORTANT:
% Infiltration_Module may or may not be active.
%
% If infiltration is active, Infiltration_Module already added
% Hydro_States.f * dt_h to Soil_Properties.I_t and removed it from depths.d_t.
%
% If infiltration is inactive, recharge can still occur from existing
% vadose-zone storage Soil_Properties.I_t.
%
% Therefore, do NOT pass Hydro_States.f again into the recharge function.
% The recharge function should drain the current vadose storage only.
% -------------------------------------------------------------------------

if has_layered_soil
    [recharge_down_rate, Soil_Properties, cumulative_recharge] = ...
        simulate_layered_groundwater_recharge( ...
        Soil_Properties, dt_s, idx_nan, cumulative_recharge);
else
    inf_ms = zero_matrix;                     % [m/s] prevents double-counting f

    S0_m   = Soil_Properties.I_t / 1000;      % [m water], current UZ storage

    % For the current above-residual storage convention, the minimum admissible
    % UZ storage is zero unless you intentionally define a non-drainable store.
    minS_m = zeros(size(S0_m), 'like', S0_m); % [m water]

    [recharge_down_rate, Soil_Moisture, cumulative_recharge] = simulate_groundwater_recharge( ...
        inf_ms, ...
        S0_m, ...
        Soil_Properties.alpha_vg, ...
        dt_s, ...
        minS_m, ...
        UZ_max_storage, ...
        LULC_Properties.idx_imp, ...
        cumulative_recharge ...
        );

    % Update vadose storage state for next modules [mm]
    Soil_Properties.I_t = Soil_Moisture * 1000;
end

recharge_down_rate(isnan(recharge_down_rate)) = nan;

%% ------------------------------------------------------------------------
% 3b) COMPUTE SIMPLE LAYERED CAPILLARY RISE AND ACCUMULATE NET EXCHANGE
% -------------------------------------------------------------------------

% Bedrock elevation [m]
z_bed = z_terrain - Soil_Properties.Soil_Depth;

capillary_rate = zero_matrix;              % [m/s], positive upward to vadose
if has_layered_soil && flags.flag_capillary_rise == 1
    Sy_for_cap = local_groundwater_specific_yield(Soil_Properties);
    saturated_thickness_m = max(BC_States.h_t - z_bed, 0);
    pending_net_m = GW_States.pending_net_exchange_m;
    pending_net_m(~isfinite(pending_net_m)) = 0;
    available_gw_depth_m = max(saturated_thickness_m .* Sy_for_cap + pending_net_m, 0);

    [capillary_rate, Soil_Properties] = simulate_layered_capillary_rise( ...
        Soil_Properties, BC_States, z_terrain, dt_s, idx_nan, available_gw_depth_m);
end

recharge_rate = recharge_down_rate - capillary_rate; % net vadose-GW exchange [m/s]
recharge_rate(idx_nan) = nan;

net_exchange_depth_m = recharge_rate .* dt_s;
downward_depth_m = max(recharge_down_rate, 0) .* dt_s;
capillary_depth_m = max(capillary_rate, 0) .* dt_s;

GW_States.pending_net_exchange_m = GW_States.pending_net_exchange_m + net_exchange_depth_m;
GW_States.pending_recharge_m = GW_States.pending_recharge_m + downward_depth_m;
GW_States.pending_capillary_m = GW_States.pending_capillary_m + capillary_depth_m;
GW_States.elapsed_s = GW_States.elapsed_s + dt_s;
GW_States.total_recharge_m3 = GW_States.total_recharge_m3 + ...
    nansum(nansum(downward_depth_m .* GW_States.cell_area_m2));
GW_States.total_capillary_m3 = GW_States.total_capillary_m3 + ...
    nansum(nansum(capillary_depth_m .* GW_States.cell_area_m2));
GW_States.last_surface_recharge_rate_m_s = recharge_down_rate;
GW_States.last_capillary_rate_m_s = capillary_rate;
GW_States.last_net_exchange_rate_m_s = recharge_rate;

%% ------------------------------------------------------------------------
% 4) GROUNDWATER UPDATE
% -------------------------------------------------------------------------
% IMPORTANT:
%
% flags.flag_baseflow controls only lateral groundwater propagation.
%
% flags.flag_groundwater_modeling controls recharge and local groundwater
% table update.
%
% Therefore:
%
%   flag_groundwater_modeling = 1, flag_baseflow = 0
%       -> recharge + local recharge-to-water-table update only
%
%   flag_groundwater_modeling = 1, flag_baseflow = 1
%       -> recharge + Boussinesq lateral groundwater propagation
%
% -------------------------------------------------------------------------

dt_stable_s = estimate_groundwater_stable_timestep( ...
    Soil_Properties, BC_States, z_bed, Wshed_Properties.Resolution, ...
    Wshed_Properties.Resolution, flags, idx_nan);

target_dt_s = min(GW_States.target_dt_s, dt_stable_s);
target_dt_s = max(target_dt_s, GW_States.min_dt_s);

Sy_head = local_groundwater_specific_yield(Soil_Properties);
head_change_m = abs(GW_States.pending_net_exchange_m) ./ max(Sy_head, 1e-9);
head_change_due = max(head_change_m(:), [], 'omitnan') >= GW_States.max_head_change_m;

is_final_step = false;
if exist('t', 'var') && exist('running_control', 'var') && ...
        isstruct(running_control) && isfield(running_control, 'routing_time')
    is_final_step = t >= (running_control.routing_time - 1e-9);
end

if flags.flag_groundwater_async == 0
    groundwater_update_due = true;
else
    groundwater_update_due = GW_States.elapsed_s >= target_dt_s || ...
        head_change_due || is_final_step;
end

groundwater_update_ran = false;
dt_gw_s = max(GW_States.elapsed_s, eps);
recharge_for_gw = zero_matrix;
if groundwater_update_due
    recharge_for_gw = GW_States.pending_net_exchange_m ./ dt_gw_s;
    recharge_for_gw(idx_nan) = nan;
end

if groundwater_update_due && flags.flag_baseflow == 1

    %% --------------------------------------------------------------------
    % 4a) FULL GROUNDWATER PROPAGATION: BOUSSINESQ SOLVER
    %% --------------------------------------------------------------------

    [BC_States.h_t, ~, ~, q_exf, q_river, error] = Boussinesq_2D_explicit( ...
        dt_gw_s, ...
        Wshed_Properties.Resolution, ...
        Wshed_Properties.Resolution, ...
        BC_States.h_t, ...
        z_bed, ...
        Soil_Properties.Sy, ...
        recharge_for_gw, ...
        Soil_Properties.ksat_gw / 1000 / 3600, ...
        idx_rivers, ...
        LULC_Properties.River_K_coeff * Soil_Properties.ksat / 1000 / 3600, ...
        z_terrain + depths.d_t / 1000, ...
        z_terrain, ...
        flags.groundwater_courant, ...
        Soil_Properties.Soil_Depth, ...
        Wshed_Properties.domain, ...
        [], [], Wshed_Properties.perimeter ...
        );

    % Update previous groundwater head
    BC_States.h_0 = BC_States.h_t;

    % Add groundwater exfiltration to surface water [mm]
    depths.d_t = depths.d_t + dt_gw_s * (q_exf * 1000);
    depths.d_t(depths.d_t < 0) = 0;
    groundwater_update_ran = true;

elseif groundwater_update_due

    %% --------------------------------------------------------------------
    % 4b) LOCAL RECHARGE-TO-WATER-TABLE UPDATE
    %% --------------------------------------------------------------------
    % No lateral groundwater propagation.
    %
    % Recharge still raises the local water table even if
    % flag_infiltration == 0, because recharge drains existing vadose-zone
    % storage.
    %
    % Use groundwater specific yield for water-table movement. Vadose
    % storage capacity still uses theta_sat - theta_r, but saturated
    % groundwater storage changes are governed by Sy.
    %% --------------------------------------------------------------------

    Sy_local = local_groundwater_specific_yield(Soil_Properties);

    % Recharge depth over the current time step [m water]
    recharge_depth_m = recharge_for_gw .* dt_gw_s;

    % Local water-table rise [m]
    dh_wt = recharge_depth_m ./ Sy_local;

    % Candidate groundwater table elevation [m]
    h_candidate = BC_States.h_t + dh_wt;

    % Do not allow water table below bedrock
    h_candidate = max(h_candidate, z_bed);

    % If the groundwater table rises above the land surface, only the excess
    % groundwater storage above the surface becomes surface water.
    excess_head_m  = max(h_candidate - z_terrain, 0);   % [m head]
    excess_water_m = excess_head_m .* Sy_local;         % [m water]

    % Clamp groundwater table at the land surface
    BC_States.h_t = min(h_candidate, z_terrain);

    % Update previous groundwater head
    BC_States.h_0 = BC_States.h_t;

    % Add only true above-surface groundwater excess to surface water [mm]
    depths.d_t = depths.d_t + excess_water_m * 1000;
    depths.d_t(depths.d_t < 0) = 0;
    groundwater_update_ran = true;

end

if groundwater_update_ran
    GW_States.last_groundwater_dt_s = dt_gw_s;
    GW_States.last_groundwater_recharge_rate_m_s = recharge_for_gw;
    GW_States.last_q_exf_m_s = q_exf;
    GW_States.last_q_river_m_s = q_river;
    GW_States.total_exfiltration_m3 = GW_States.total_exfiltration_m3 + ...
        nansum(nansum(max(q_exf, 0) .* dt_gw_s .* GW_States.cell_area_m2));
    GW_States.n_updates = GW_States.n_updates + 1;
    GW_States.pending_net_exchange_m(:) = 0;
    GW_States.pending_recharge_m(:) = 0;
    GW_States.pending_capillary_m(:) = 0;
    GW_States.elapsed_s = 0;
else
    GW_States.last_groundwater_dt_s = 0;
    GW_States.last_groundwater_recharge_rate_m_s = zero_matrix;
    GW_States.last_q_exf_m_s = zero_matrix;
    GW_States.last_q_river_m_s = zero_matrix;
end

%% ------------------------------------------------------------------------
% 4c) CONSISTENT UZ CAPACITY CLAMP AFTER WATER-TABLE UPDATE
% -------------------------------------------------------------------------
% If the water table rises, the unsaturated thickness decreases.
% Therefore, the maximum vadose storage decreases.
%
% Any vadose water above the new storage capacity becomes saturation-excess
% surface water.
%
% IMPORTANT:
% Soil_Properties.I_t is above-residual storage.
% Therefore the capacity must use theta_sat - theta_r, not theta_sat - theta_i.
% -------------------------------------------------------------------------

if groundwater_update_ran && has_layered_soil
    layer_options = struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005);
    [Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
        Soil_Properties, LULC_Properties, BC_States, z_terrain, idx_nan, layer_options);
    [Soil_Properties, excess_mm] = clamp_layered_soil_storage(Soil_Properties, idx_nan);

    if any(excess_mm(:) > 0)
        depths.d_t = depths.d_t + excess_mm;
        depths.d_t(depths.d_t < 0) = 0;
    end
elseif groundwater_update_ran
    zwt_new = z_terrain - BC_States.h_t;                  % [m]
    zwt_new = max(zwt_new, 0);
    zwt_new = min(zwt_new, Soil_Properties.Soil_Depth);

    UZ_max_storage_new = zwt_new .* ...
        (Soil_Properties.theta_sat - Soil_Properties.theta_r);   % [m water]

    S_m = Soil_Properties.I_t / 1000;                            % [m water]

    excess_m = max(S_m - UZ_max_storage_new, 0);                 % [m water]

    if any(excess_m(:) > 0)

        % Remove excess from vadose storage
        Soil_Properties.I_t = (S_m - excess_m) * 1000;
        Soil_Properties.I_t = max(Soil_Properties.I_t, 0);

        % Add excess vadose water to the surface
        depths.d_t = depths.d_t + excess_m * 1000;
        depths.d_t(depths.d_t < 0) = 0;

    end
end

%% ------------------------------------------------------------------------
% 5) TRACK MAXIMUM GROUNDWATER DEPTH (saturated thickness) [m]
% -------------------------------------------------------------------------

max_GW_depth = max(max_GW_depth, BC_States.h_t - ...
    (elevation - Soil_Properties.Soil_Depth));

%% ------------------------------------------------------------------------
% 7) STORE MODEL ERROR (Boussinesq solver diagnostic)
% -------------------------------------------------------------------------
errors(5) = error;

function Sy_local = local_groundwater_specific_yield(Soil_Properties)
if isfield(Soil_Properties, 'Sy') && ~isempty(Soil_Properties.Sy)
    Sy_local = Soil_Properties.Sy;
else
    Sy_local = Soil_Properties.theta_sat - Soil_Properties.theta_r;
end
Sy_local = max(Sy_local, 1e-6);
end
