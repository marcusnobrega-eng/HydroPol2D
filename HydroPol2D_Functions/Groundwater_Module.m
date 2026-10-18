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
seepage_discharge_m3_s = 0;

has_layered_soil = isfield(Soil_Properties, 'Layers') && ...
    isstruct(Soil_Properties.Layers) && ...
    isfield(Soil_Properties.Layers, 'near_surface_storage_mm');
has_prescribed_recharge = isfield(flags, 'flag_prescribed_recharge') && ...
    flags.flag_prescribed_recharge == 1 && ...
    isfield(BC_States, 'gw_prescribed_recharge') && ...
    isstruct(BC_States.gw_prescribed_recharge) && ...
    isfield(BC_States.gw_prescribed_recharge, 'time_min') && ...
    ~isempty(BC_States.gw_prescribed_recharge.time_min);
surface_coupling_enabled = ~(isfield(flags, 'flag_groundwater_surface_coupling') && ...
    flags.flag_groundwater_surface_coupling == 0);

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

if has_prescribed_recharge
    Soil_Properties.I_t(~idx_nan) = 0;
end

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

perv_frac = ones(size(zwt), 'like', zwt);
if isfield(LULC_Properties, 'frac_perv') && ~isempty(LULC_Properties.frac_perv)
    perv_frac = min(max(LULC_Properties.frac_perv, 0), 1);
elseif isfield(LULC_Properties, 'idx_imp') && ~isempty(LULC_Properties.idx_imp)
    perv_frac = double(~LULC_Properties.idx_imp);
end
perv_frac(idx_nan) = nan;

if has_layered_soil
    Soil_Properties = sync_layered_soil_storage(Soil_Properties, idx_nan);
    UZ_max_storage = Soil_Properties.Layers.total_vadose_capacity_mm ./ 1000;
else
    % Maximum UZ storage [m water]
    % Soil_Properties.I_t is above-residual storage, so capacity must be
    % theta_sat - theta_r, not theta_sat - theta_i.
    UZ_max_storage = perv_frac .* zwt .* ...
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

% Bedrock elevation [m]
z_bed = z_terrain - Soil_Properties.Soil_Depth;

if has_prescribed_recharge
    recharge_down_rate = zero_matrix;
    capillary_rate = zero_matrix;
    recharge_rate = zero_matrix;
    downward_depth_m = zero_matrix;
    capillary_depth_m = zero_matrix;
    net_exchange_depth_m = zero_matrix;
    GW_States.elapsed_s = GW_States.elapsed_s + dt_s;
    GW_States.last_surface_recharge_rate_m_s = zero_matrix;
    GW_States.last_capillary_rate_m_s = zero_matrix;
    GW_States.last_net_exchange_rate_m_s = zero_matrix;
else
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

if isfield(flags,'flag_deep_aquifer') && flags.flag_deep_aquifer == 1
    assert(isfield(BC_States,'deep_aquifer_subbasin_id') && ...
        isfield(BC_States,'deep_aquifer_channel_weight'), ...
        'Deep aquifer requires aligned subbasin and channel-weight rasters.');
    if ~isfield(GW_States,'deep_aquifer_storage_m3') && ...
            isfield(BC_States,'deep_aquifer_initial_storage_m3')
        GW_States.deep_aquifer_storage_m3 = ...
            BC_States.deep_aquifer_initial_storage_m3;
        GW_States.deep_aquifer_discharge_m3_s = ...
            BC_States.deep_aquifer_initial_discharge_m3_s;
        GW_States.cumulative_deep_aquifer_recharge_m3 = ...
            zeros(size(GW_States.deep_aquifer_storage_m3));
        GW_States.cumulative_deep_aquifer_discharge_m3 = ...
            zeros(size(GW_States.deep_aquifer_storage_m3));
    end
    [GW_States,recharge_down_rate,deep_surface_addition_mm] = ...
        hp2d_deep_aquifer_step(GW_States,recharge_down_rate, ...
        BC_States.deep_aquifer_subbasin_id,BC_States.deep_aquifer_channel_weight, ...
        Soil_Properties.deep_recharge_fraction, ...
        Soil_Properties.deep_aquifer_response_time_days, ...
        Wshed_Properties.cell_area,dt_s);
    depths.d_t = depths.d_t + deep_surface_addition_mm;
end

%% ------------------------------------------------------------------------
% 3b) COMPUTE SIMPLE LAYERED CAPILLARY RISE AND ACCUMULATE NET EXCHANGE
% -------------------------------------------------------------------------

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
end

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

Sy_head = local_groundwater_specific_yield(Soil_Properties);
prescribed_input_mm_h = 0;
if has_prescribed_recharge
    [prescribed_rate,GW_States,prescribed_input_mm_h] = compute_prescribed_recharge_rate( ...
        BC_States.gw_prescribed_recharge,flags,GW_States,dt_s,zero_matrix,idx_nan,t);
    prescribed_depth_m=max(prescribed_rate,0).*dt_s;
    GW_States.pending_net_exchange_m=GW_States.pending_net_exchange_m+prescribed_depth_m;
    GW_States.pending_recharge_m=GW_States.pending_recharge_m+prescribed_depth_m;
    GW_States.total_recharge_m3=GW_States.total_recharge_m3+ ...
        nansum(prescribed_depth_m.*GW_States.cell_area_m2,'all');
    GW_States.last_surface_recharge_rate_m_s=prescribed_rate;
    GW_States.last_net_exchange_rate_m_s=prescribed_rate;
end

% River exchange is a local two-store process. Apply its exact solution at
% every surface step, then include it in the local groundwater storage update.
if isfield(flags,'flag_river_groundwater_exchange') && ...
        flags.flag_river_groundwater_exchange==1
    effective_head=BC_States.h_t+ ...
        GW_States.pending_net_exchange_m./max(Sy_head,1e-9);
    [effective_after,depths.d_t,q_river,river_exchange_m3]= ...
        structured_river_exchange(effective_head,depths.d_t,z_bed,z_terrain, ...
        Soil_Properties,idx_rivers,dt_s,GW_States.cell_area_m2,idx_nan);
    GW_States.pending_net_exchange_m=GW_States.pending_net_exchange_m+ ...
        (effective_after-effective_head).*Sy_head;
    if ~isfield(GW_States,'total_river_exchange_m3')
        GW_States.total_river_exchange_m3=0;
    end
    GW_States.total_river_exchange_m3= ...
        GW_States.total_river_exchange_m3+river_exchange_m3;
end

% Recharge, capillary rise, and river exchange alter local storage now. They
% do not wait for the slower lateral Boussinesq redistribution.
h_before_local=BC_States.h_t;
h_candidate=BC_States.h_t+GW_States.pending_net_exchange_m./max(Sy_head,1e-9);
h_candidate=max(h_candidate,z_bed);
local_excess_m=max(h_candidate-z_terrain,0).*Sy_head;
BC_States.h_t=min(h_candidate,z_terrain);
BC_States.h_0=BC_States.h_t;
if surface_coupling_enabled
    depths.d_t=depths.d_t+local_excess_m.*1000;
    depths.d_t=max(depths.d_t,0);
end
GW_States.pending_net_exchange_m(:)=0;
q_exf=local_excess_m./dt_s;
q_exf(idx_nan)=nan;
step_exfiltration_volume_m3=nansum(local_excess_m.*GW_States.cell_area_m2,'all');

if ~isfield(GW_States,'lateral_reference_head_m') || ...
        ~isequal(size(GW_States.lateral_reference_head_m),size(BC_States.h_t))
    GW_States.lateral_reference_head_m=BC_States.h_t;
end
scale=flags.groundwater_head_tolerance_m+ ...
    flags.groundwater_relative_tolerance.*max(BC_States.h_t-z_bed,0);
normalized_head_change=abs(BC_States.h_t-GW_States.lateral_reference_head_m)./max(scale,realmin);
normalized_head_change(idx_nan)=0;
head_change_due=max(normalized_head_change(:),[],'omitnan')>=1;

if flags.flag_baseflow==1
    dt_stable_s=estimate_groundwater_stable_timestep( ...
        Soil_Properties,BC_States,z_bed,Wshed_Properties.Resolution, ...
        Wshed_Properties.Resolution,idx_nan);
    next_lateral_dt_s=min(dt_stable_s,GW_States.suggested_dt_s);
else
    next_lateral_dt_s=inf;
end

is_final_step = false;
if exist('t', 'var') && exist('running_control', 'var') && ...
        isstruct(running_control) && isfield(running_control, 'routing_time')
    is_final_step = t >= (running_control.routing_time - 1e-9);
end

groundwater_update_due=flags.flag_baseflow==1 && ...
    (GW_States.elapsed_s>=next_lateral_dt_s || head_change_due || is_final_step);
groundwater_update_ran=true;
lateral_update_ran=false;
dt_gw_s = max(GW_States.elapsed_s, eps);
recharge_for_gw = zero_matrix;

dirichlet_mask = [];
h_dirichlet = [];
if isfield(BC_States, 'gw_dirichlet_mask') && isfield(BC_States, 'gw_dirichlet_head') && ...
        ~isempty(BC_States.gw_dirichlet_mask) && ~isempty(BC_States.gw_dirichlet_head)
    dirichlet_mask = logical(BC_States.gw_dirichlet_mask);
    h_dirichlet = BC_States.gw_dirichlet_head;
end

seepage_mask = [];
h_seepage = [];
if isfield(flags, 'flag_groundwater_seepage_face') && flags.flag_groundwater_seepage_face == 1 && ...
        ~isempty(dirichlet_mask) && ~isempty(h_dirichlet)
    seepage_mask = dirichlet_mask;
    h_seepage = h_dirichlet;
    dirichlet_mask = [];
    h_dirichlet = [];
end

gw_ksat_m_s = Soil_Properties.ksat_gw / 1000 / 3600;
if ~isempty(dirichlet_mask) && any(dirichlet_mask(:))
    dirichlet_ksat_mult = max(double(flags.gw_dirichlet_ksat_multiplier), 1);
    if dirichlet_ksat_mult ~= 1
        gw_ksat_m_s(dirichlet_mask) = gw_ksat_m_s(dirichlet_mask) .* dirichlet_ksat_mult;
    end
end

if groundwater_update_due

    %% --------------------------------------------------------------------
    % 4a) FULL GROUNDWATER PROPAGATION: BOUSSINESQ SOLVER
    %% --------------------------------------------------------------------

    initial_step_s=min(GW_States.suggested_dt_s,dt_gw_s);
    minimum_step_s=max(min(dt_s,1),eps);
    [BC_States.h_t,~,~,lateral_q_exf,error,seepage_discharge_m3_s,adaptive] = ...
        hp2d_boussinesq_adaptive(dt_gw_s,initial_step_s,minimum_step_s, ...
        flags.groundwater_relative_tolerance,flags.groundwater_head_tolerance_m, ...
        Wshed_Properties.Resolution,Wshed_Properties.Resolution,BC_States.h_t, ...
        z_bed,Soil_Properties.Sy,recharge_for_gw,gw_ksat_m_s,z_terrain, ...
        Soil_Properties.Soil_Depth,Wshed_Properties.domain,dirichlet_mask, ...
        h_dirichlet,seepage_mask,h_seepage,Wshed_Properties.perimeter,0.25);

    % Update previous groundwater head
    BC_States.h_0 = BC_States.h_t;

    % Add groundwater exfiltration to surface water [mm]
    if surface_coupling_enabled
        depths.d_t = depths.d_t + dt_gw_s * (lateral_q_exf * 1000);
        depths.d_t(depths.d_t < 0) = 0;
    end
    q_exf=q_exf+lateral_q_exf;
    step_exfiltration_volume_m3=step_exfiltration_volume_m3+ ...
        nansum(max(lateral_q_exf,0).*dt_gw_s.*GW_States.cell_area_m2,'all');
    GW_States.suggested_dt_s=adaptive.suggested_dt_s;
    GW_States.adaptive_accepted_steps=GW_States.adaptive_accepted_steps+adaptive.accepted_steps;
    GW_States.adaptive_rejected_steps=GW_States.adaptive_rejected_steps+adaptive.rejected_steps;
    GW_States.lateral_reference_head_m=BC_States.h_t;
    GW_States.n_updates=GW_States.n_updates+1;
    GW_States.elapsed_s=0;
    GW_States.pending_recharge_m(:)=0;
    GW_States.pending_capillary_m(:)=0;
    lateral_update_ran=true;
elseif flags.flag_baseflow~=1
    GW_States.elapsed_s=0;
    GW_States.pending_recharge_m(:)=0;
    GW_States.pending_capillary_m(:)=0;
end

if groundwater_update_ran
    recharge_volume_m3=nansum(max(GW_States.last_surface_recharge_rate_m_s,0).* ...
        dt_s.*GW_States.cell_area_m2,'all');
    GW_States.last_groundwater_dt_s = double(lateral_update_ran).*dt_gw_s;
    GW_States.last_groundwater_recharge_rate_m_s = GW_States.last_net_exchange_rate_m_s;
    GW_States.last_q_exf_m_s = q_exf;
    GW_States.last_q_river_m_s = q_river;
    GW_States.last_groundwater_recharge_volume_m3 = recharge_volume_m3;
    GW_States.total_exfiltration_m3 = GW_States.total_exfiltration_m3 + ...
        step_exfiltration_volume_m3;
    storage_change_rate = (BC_States.h_t - h_before_local) .* ...
        local_groundwater_specific_yield(Soil_Properties) ./ max(dt_s, eps);
    GW_States.last_storage_change_rate_m_s = mean(storage_change_rate(~idx_nan), 'omitnan');
    GW_States.last_boundary_seepage_discharge_m3_s = max(seepage_discharge_m3_s, 0);
    if isfield(Wshed_Properties, 'drainage_area') && ...
            isfinite(Wshed_Properties.drainage_area) && Wshed_Properties.drainage_area > 0
        drainage_area_m2 = Wshed_Properties.drainage_area;
    else
        drainage_area_m2 = nansum(GW_States.cell_area_m2(~idx_nan), 'all');
    end
    GW_States.last_boundary_seepage_discharge_mm_h = ...
        GW_States.last_boundary_seepage_discharge_m3_s / max(drainage_area_m2, eps) * 1000 * 3600;
    [mass_balance_q_m3_s, mass_balance_q_mm_h, active_area_m2] = ...
        compute_mass_balance_seepage_discharge( ...
        GW_States.last_net_exchange_rate_m_s,BC_States.h_t,h_before_local,z_bed,z_terrain, ...
        Soil_Properties,GW_States.cell_area_m2,idx_nan,dt_s,flags);
    GW_States.last_mass_balance_seepage_discharge_m3_s = max(mass_balance_q_m3_s, 0);
    GW_States.last_mass_balance_seepage_discharge_mm_h = max(mass_balance_q_mm_h, 0);
    GW_States.last_mass_balance_active_area_m2 = active_area_m2;
    if isfield(flags, 'flag_groundwater_mass_balance_discharge') && ...
            flags.flag_groundwater_mass_balance_discharge == 1
        GW_States.last_seepage_discharge_m3_s = GW_States.last_mass_balance_seepage_discharge_m3_s;
        GW_States.last_seepage_discharge_m_s = GW_States.last_mass_balance_seepage_discharge_m3_s;
        GW_States.last_seepage_discharge_mm_h = GW_States.last_mass_balance_seepage_discharge_mm_h;
    else
        GW_States.last_seepage_discharge_m3_s = GW_States.last_boundary_seepage_discharge_m3_s;
        GW_States.last_seepage_discharge_m_s = GW_States.last_boundary_seepage_discharge_m3_s;
        GW_States.last_seepage_discharge_mm_h = GW_States.last_boundary_seepage_discharge_mm_h;
    end
    GW_States.last_prescribed_input_mm_h = prescribed_input_mm_h;
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

if groundwater_update_ran && has_layered_soil && ~has_prescribed_recharge
    layer_options = struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005);
    [Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
        Soil_Properties, LULC_Properties, BC_States, z_terrain, idx_nan, layer_options);
    [Soil_Properties, excess_mm] = clamp_layered_soil_storage(Soil_Properties, idx_nan);

    if any(excess_mm(:) > 0)
        depths.d_t = depths.d_t + excess_mm;
        depths.d_t(depths.d_t < 0) = 0;
    end
elseif groundwater_update_ran && ~has_prescribed_recharge
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

function [head,surface_depth_mm,rate_m_s,net_volume_m3] = structured_river_exchange( ...
        head,surface_depth_mm,bedrock,surface,soil,river_mask,dt_s,cell_area,idx_nan)
assert(isfield(soil,'riverbed_hydraulic_conductivity_m_s') && ...
    isfield(soil,'cell_riverbed_area_m2') && isfield(soil,'riverbed_thickness_m'), ...
    'HydroPol2D:MissingRiverbedParameters', ...
    'River-groundwater exchange requires riverbed conductivity, area, and thickness.');
kbed=double(soil.riverbed_hydraulic_conductivity_m_s);
area=max(double(soil.cell_riverbed_area_m2),0);
thickness=double(soil.riverbed_thickness_m);
assert(isfinite(thickness) && thickness>0 && all(kbed(~idx_nan)>=0), ...
    'HydroPol2D:InvalidRiverbedParameters', ...
    'Riverbed conductivity must be nonnegative and thickness positive.');
active=logical(river_mask) & area>0 & ~idx_nan;
stage=surface+max(surface_depth_mm,0)./1000;
sy=local_groundwater_specific_yield(soil);
groundwater_available=max(head-bedrock,0).*sy.*cell_area;
surface_available=max(surface_depth_mm,0)./1000.*cell_area;
[volume,~]=hp2d_river_groundwater_exchange_step(kbed,thickness,area,head,stage, ...
    groundwater_available,surface_available,sy.*cell_area,cell_area,dt_s,active);
head=head-volume./max(sy.*cell_area,eps);
surface_depth_mm=surface_depth_mm+volume./max(cell_area,eps).*1000;
surface_depth_mm=max(surface_depth_mm,0);
rate_m_s=volume./max(cell_area,eps)./dt_s;
rate_m_s(idx_nan)=nan;
net_volume_m3=nansum(volume,'all');
end

function [q_m3_s, q_mm_h, active_area_m2] = compute_mass_balance_seepage_discharge( ...
        recharge_rate_m_s, h_now, h_prev, z_bed, z_surface, Soil_Properties, ...
        cell_area_m2, idx_nan, dt_s, flags)
active_depth_threshold_m = max(double(flags.groundwater_active_depth_threshold_m), 0);
active_mask = ~idx_nan & isfinite(h_now) & isfinite(z_bed) & ...
    (max(h_now - z_bed, 0) > active_depth_threshold_m);
active_area_m2 = nansum(nansum(double(active_mask) .* cell_area_m2));
if ~any(active_mask(:)) || ~isfinite(active_area_m2) || active_area_m2 <= 0
    q_m3_s = 0;
    q_mm_h = 0;
    return
end

h_half = 0.5 * (h_now + h_prev);
Sy_mb = storage_dependent_specific_yield(Soil_Properties, z_surface, h_half);
dh_dt = (h_now - h_prev) ./ max(dt_s, eps);
cell_flux_m_s = recharge_rate_m_s - Sy_mb .* dh_dt;
q_flux_m_s = mean(cell_flux_m_s(active_mask), 'omitnan');
if ~isfinite(q_flux_m_s)
    q_flux_m_s = 0;
end
q_flux_m_s = max(q_flux_m_s, 0);
q_m3_s = q_flux_m_s * active_area_m2;
q_mm_h = q_flux_m_s * 1000 * 3600;
end

function Sy_mb = storage_dependent_specific_yield(Soil_Properties, z_surface, h_total)
has_vg = isfield(Soil_Properties, 'theta_sat') && isfield(Soil_Properties, 'theta_r') && ...
    isfield(Soil_Properties, 'alpha_vg') && isfield(Soil_Properties, 'n_vg') && ...
    ~isempty(Soil_Properties.theta_sat) && ~isempty(Soil_Properties.theta_r) && ...
    ~isempty(Soil_Properties.alpha_vg) && ~isempty(Soil_Properties.n_vg);
if ~has_vg
    Sy_mb = local_groundwater_specific_yield(Soil_Properties);
    return
end

effective_porosity = max(Soil_Properties.theta_sat - Soil_Properties.theta_r, 0);
alpha_vg = max(abs(Soil_Properties.alpha_vg), eps);
n_vg = max(Soil_Properties.n_vg, 1 + 1e-6);
suction_head = max(z_surface - h_total, 0);
Sy_mb = effective_porosity .* ...
    (1 - (1 + (alpha_vg .* suction_head) .^ n_vg) .^ (-(n_vg + 1) ./ n_vg));
Sy_mb = max(Sy_mb, 1e-6);
end

function [recharge_rate_m_s, GW_States, input_mm_h] = compute_prescribed_recharge_rate( ...
        forcing, flags, GW_States, dt_s, template, idx_nan, t_min)
% Compute delayed groundwater recharge from a prescribed net-input series.
% The BSM workbook provides net rainfall-like input to a linear UZ reservoir;
% HydroPol uses positive downward recharge, so the equivalent equations are:
%   dS/dt = P - k S
%   R     = k S

recharge_rate_m_s = zeros(size(template), 'like', template);
recharge_rate_m_s(idx_nan) = nan;
input_mm_h = 0;

if nargin < 7 || isempty(t_min)
    t_min = 0;
end

if isempty(forcing.time_min) || isempty(forcing.input_m_s)
    return
end

input_rate_m_s = interp1( ...
    double(forcing.time_min(:)), ...
    double(forcing.input_m_s(:)), ...
    double(t_min), ...
    'previous', 0);
input_mm_h = input_rate_m_s * 1000 * 3600;

storage_mm = GW_States.prescribed_uz_storage_mm;
storage_mm(~isfinite(storage_mm)) = 0;

damping_factor = double(flags.prescribed_recharge_damping_factor);
k_per_s = max(double(flags.prescribed_recharge_k_per_s), 0);
input_map_m_s = damping_factor .* input_rate_m_s .* ones(size(template), 'like', template);
input_map_m_s(idx_nan) = nan;

recharge_rate_m_s = k_per_s .* max(storage_mm ./ 1000, 0);
recharge_rate_m_s(idx_nan) = nan;

storage_next_m = max(storage_mm ./ 1000 + dt_s .* (input_map_m_s - recharge_rate_m_s), 0);
GW_States.prescribed_uz_storage_mm = storage_next_m .* 1000;
GW_States.prescribed_uz_storage_mm(idx_nan) = nan;
end
