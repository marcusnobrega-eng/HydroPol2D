%% =========================================================================
%  LAYERED_SUBSURFACE_MODULE (flag_layered_subsurface = 1)
%  =========================================================================
%  Replaces Evaporation_Evapotranspiration_Module (soil part), Infiltration_Module
%  and Groundwater_Module with the layered soil-over-regolith column: slot geometry
%  per cell down to min(DTB, max depth), one water table splitting the slot it is in,
%  Noah-MP MMF rise/fall rules, gravity drainage, a Darcy flux at the table, layered
%  transmissivity, and Feddes root-water stress. Same order as HydroPol2D-Python
%  core.layered_vertical.advance_layered_subsurface:
%    1) infiltration into the soil slots   2) open-water evaporation on wet cells
%    3) soil ET (Feddes, or the rest of the demand from the saturated root zone)
%    4) drain cascade                      5) Darcy flux at the table
%    6) deep-aquifer split of the recharge 7) water-table update (exfiltration)
%    8) river-groundwater exchange         9) lateral flow on the scheduler
%   10) exfiltration to the surface
%  State: Subsurface.water_mm [cells x slots] (water above theta_r in the unsaturated
%  part of each slot) and BC_States.h_t. Soil_Properties.Layers and I_t mirror slots
%  1-3 for outputs. Docs: HydroPol2D-Python docs/layered-subsurface.md.
% =========================================================================

dt_s = time_step*60;
dt_h = time_step/60;
dt_days = time_step/60/24;
[lay_ny, lay_nx] = size(elevation);
lay_valid = Subsurface.valid;
lay_area = Wshed_Properties.Resolution^2;
lay_cfg = Subsurface.config;
% Infiltration_Module defaults (that module does not run in layered mode).
if ~isfield(Soil_Properties, 'l_vg') || isempty(Soil_Properties.l_vg), Soil_Properties.l_vg = 0.5; end
if ~isfield(Soil_Properties, 'Ltop') || isempty(Soil_Properties.Ltop), Soil_Properties.Ltop = 0.05; end
if ~isfield(Soil_Properties, 'dh_max') || isempty(Soil_Properties.dh_max), Soil_Properties.dh_max = 1.0; end
if ~isfield(Soil_Properties, 'Ksurf_weight') || isempty(Soil_Properties.Ksurf_weight)
    Soil_Properties.Ksurf_weight = 0.5;
end

if k == 1 || ~exist('cumulative_recharge', 'var') || isempty(cumulative_recharge)
    cumulative_recharge = zeros(size(DEM_raster.Z));
    cumulative_recharge(isnan(DEM_raster.Z)) = nan;
end
if k == 1 || ~exist('cumulative_infiltration', 'var') || isempty(cumulative_infiltration)
    cumulative_infiltration = zeros(size(DEM_raster.Z));
    cumulative_infiltration(isnan(DEM_raster.Z)) = nan;
end
if k == 1 || ~exist('max_GW_depth', 'var') || isempty(max_GW_depth)
    max_GW_depth = zeros(size(DEM_raster.Z));
    max_GW_depth(isnan(DEM_raster.Z)) = nan;
end
lay_zero = zeros(size(elevation), 'like', elevation);
lay_zero(idx_nan) = nan;
if ~exist('GW_States', 'var') || ~isstruct(GW_States)
    GW_States = struct();
end
[flags, GW_States] = initialize_groundwater_scheduler(flags, GW_States, lay_zero, idx_nan);
GW_States.cell_area_m2 = lay_area + lay_zero;

% Pervious fraction (as Infiltration_Module).
lay_pervious = ones(size(elevation), 'like', elevation);
if isfield(LULC_Properties, 'frac_perv') && ~isempty(LULC_Properties.frac_perv)
    lay_pervious = min(max(LULC_Properties.frac_perv, 0), 1);
elseif isfield(LULC_Properties, 'idx_imp') && ~isempty(LULC_Properties.idx_imp)
    lay_pervious = double(~LULC_Properties.idx_imp);
end

% ET demands over this step (mm): open water from Ep; soil ET as in the ET module,
% before the pervious fraction (applied inside the column step).
lay_open_water_mm = zeros(size(elevation), 'like', elevation);
lay_soil_demand_mm = zeros(size(elevation), 'like', elevation);
if flags.flag_ETP == 1
    lay_open_water_mm = max(Hydro_States.Ep, 0).*dt_days;
    if flags.flag_input_ETP_map == 1
        lay_rate = input_evaporation + input_transpiration;
    else
        if isfield(LULC_Properties, 'Kc') && ~isempty(LULC_Properties.Kc)
            lay_rate = Hydro_States.ETP.*LULC_Properties.Kc;
        else
            lay_rate = Hydro_States.ETP;
        end
    end
    lay_rate(~isfinite(lay_rate)) = 0;
    lay_soil_demand_mm = max(lay_rate, 0).*dt_days;
    lay_open_water_mm(~isfinite(lay_open_water_mm)) = 0;
end

lay_storage_before_m3 = sum(hp2d_column_storage_mm(Subsurface.column, BC_States.h_t(:), ...
    Subsurface.water_mm), 'omitnan')/1000*lay_area;
lay_surface_before = depths.d_t(:);
lay_surface_before(~lay_valid) = 0;
lay_surface_before_m3 = sum(max(lay_surface_before, 0))/1000*lay_area;

lay_P = struct('valid', lay_valid, 'pervious', lay_pervious(:), ...
    'root_depth_m', LULC_Properties.root_depth_m(:), ...
    'ksat_near_surface_mm_h', Soil_Properties.Layers.ksat_near_surface(:), ...
    'l_vg', Soil_Properties.l_vg, 'ltop_m', Soil_Properties.Ltop, ...
    'dh_max_m', Soil_Properties.dh_max, 'ksurf_weight', Soil_Properties.Ksurf_weight, ...
    'feddes', Subsurface.feddes, 'drainage_substeps', Subsurface.drainage_substeps, ...
    'infiltration', flags.flag_infiltration == 1, ...
    'evapotranspiration', flags.flag_ETP == 1, ...
    'groundwater', flags.flag_groundwater_modeling == 1);
for lay_name = {'l_vg', 'ltop_m', 'dh_max_m', 'ksurf_weight'}
    if ~isscalar(lay_P.(lay_name{1})), lay_P.(lay_name{1}) = lay_P.(lay_name{1})(:); end
end
lay_S = struct('water_mm', Subsurface.water_mm, 'head_m', BC_States.h_t(:), ...
    'surface_depth_mm', max(lay_surface_before, 0), ...
    'open_water_evaporation_mm', lay_open_water_mm(:), ...
    'soil_et_demand_mm', lay_soil_demand_mm(:));
lay_S.head_m(~lay_valid) = 0;
lay_S = hp2d_layered_vertical_step(Subsurface.column, lay_S, lay_P, dt_s);

lay_water = lay_S.water_mm;
lay_head = lay_S.head_m;
lay_surface = lay_S.surface_depth_mm;
lay_exfiltration = zeros(size(lay_head), 'like', lay_head);
lay_deep_added_mm = zeros(size(lay_head), 'like', lay_head);
lay_deep_recharge_mm = zeros(size(lay_head), 'like', lay_head);
lay_river_m3 = 0;
lay_boundary_out_m3 = 0;
lay_recharge = lay_S.recharge_mm;
lay_capillary = lay_S.capillary_mm;

if flags.flag_groundwater_modeling == 1
    lay_to_table = lay_recharge;
    if isfield(flags, 'flag_deep_aquifer') && flags.flag_deep_aquifer == 1
        if ~isfield(GW_States, 'deep_aquifer_storage_m3') && ...
                isfield(BC_States, 'deep_aquifer_initial_storage_m3')
            GW_States.deep_aquifer_storage_m3 = BC_States.deep_aquifer_initial_storage_m3;
            GW_States.deep_aquifer_discharge_m3_s = BC_States.deep_aquifer_initial_discharge_m3_s;
            GW_States.cumulative_deep_aquifer_recharge_m3 = zeros(size(GW_States.deep_aquifer_storage_m3));
            GW_States.cumulative_deep_aquifer_discharge_m3 = zeros(size(GW_States.deep_aquifer_storage_m3));
        end
        lay_rate_map = reshape(lay_recharge/1000/dt_s, lay_ny, lay_nx);
        lay_rate_map(idx_nan) = nan;
        [GW_States, lay_shallow_rate, lay_deep_add] = hp2d_deep_aquifer_step(GW_States, ...
            lay_rate_map, BC_States.deep_aquifer_subbasin_id, BC_States.deep_aquifer_channel_weight, ...
            Soil_Properties.deep_recharge_fraction, Soil_Properties.deep_aquifer_response_time_days, ...
            Wshed_Properties.cell_area, dt_s);
        lay_shallow = lay_shallow_rate(:)*dt_s*1000;
        lay_shallow(~lay_valid | ~isfinite(lay_shallow)) = 0;
        lay_deep_recharge_mm = lay_recharge - lay_shallow;
        lay_to_table = lay_shallow;
        lay_deep_added_mm = lay_deep_add(:);
        lay_deep_added_mm(~isfinite(lay_deep_added_mm)) = 0;
        lay_surface = lay_surface + lay_deep_added_mm;
    end
    [lay_head, lay_water, lay_exfiltration, lay_unmet] = hp2d_water_table_update( ...
        Subsurface.column, lay_head, lay_water, lay_to_table - lay_capillary);
    lay_capillary = lay_capillary - lay_unmet;

    if isfield(flags, 'flag_river_groundwater_exchange') && flags.flag_river_groundwater_exchange == 1
        lay_kbed = double(Soil_Properties.riverbed_hydraulic_conductivity_m_s(:)) + zeros(size(lay_head));
        lay_bed_area = max(double(Soil_Properties.cell_riverbed_area_m2(:)), 0) + zeros(size(lay_head));
        lay_active = logical(idx_rivers(:)) & lay_bed_area > 0 & lay_valid;
        lay_stage = elevation(:) + max(lay_surface, 0)/1000;
        lay_sy = max(hp2d_effective_specific_yield(Subsurface.column, lay_head), 1e-3);
        lay_gw_available = hp2d_drainable_storage_mm(Subsurface.column, lay_head)/1000*lay_area;
        lay_surface_available = max(lay_surface, 0)/1000*lay_area;
        lay_volume = hp2d_river_groundwater_exchange_step(lay_kbed, ...
            double(Soil_Properties.riverbed_thickness_m), lay_bed_area, lay_head, lay_stage, ...
            lay_gw_available, lay_surface_available, lay_sy*lay_area, lay_area + 0*lay_head, ...
            dt_s, lay_active);
        lay_volume(~isfinite(lay_volume)) = 0;
        lay_surface = max(lay_surface + lay_volume/lay_area*1000, 0);
        [lay_head, lay_water, lay_exf2, lay_unmet2] = hp2d_water_table_update( ...
            Subsurface.column, lay_head, lay_water, -lay_volume/lay_area*1000);
        lay_exfiltration = lay_exfiltration + lay_exf2;
        lay_surface = max(lay_surface - lay_unmet2, 0);
        lay_river_m3 = sum(lay_volume);
        if ~isfield(GW_States, 'total_river_exchange_m3'), GW_States.total_river_exchange_m3 = 0; end
        GW_States.total_river_exchange_m3 = GW_States.total_river_exchange_m3 + lay_river_m3;
        GW_States.last_q_river_m_s = reshape(lay_volume/lay_area/dt_s, lay_ny, lay_nx);
    end

    GW_States.elapsed_s = GW_States.elapsed_s + dt_s;
    lay_final = exist('t', 'var') && exist('running_control', 'var') && ...
        isstruct(running_control) && isfield(running_control, 'routing_time') && ...
        t >= running_control.routing_time - 1e-9;
    lay_target_dt = max(min(lay_cfg.groundwater_target_dt_s, GW_States.suggested_dt_s), ...
        lay_cfg.groundwater_min_dt_s);
    lay_due = GW_States.elapsed_s >= lay_target_dt || lay_final;
    if flags.flag_baseflow == 1 && lay_due
        lay_dir_mask = []; lay_dir_head = [];
        if isfield(BC_States, 'gw_dirichlet_mask') && isfield(BC_States, 'gw_dirichlet_head') && ...
                ~isempty(BC_States.gw_dirichlet_mask) && ~isempty(BC_States.gw_dirichlet_head)
            lay_dir_mask = logical(BC_States.gw_dirichlet_mask);
            lay_dir_head = double(BC_States.gw_dirichlet_head);
        end
        lay_R = hp2d_layered_lateral_advance(Subsurface.column, reshape(lay_head, lay_ny, lay_nx), ...
            lay_water, GW_States.elapsed_s, Wshed_Properties.Resolution, ...
            Wshed_Properties.Resolution, lay_dir_mask, lay_dir_head, lay_cfg.groundwater_courant);
        lay_head = lay_R.head_m(:);
        lay_water = lay_R.water_mm;
        lay_exfiltration = lay_exfiltration + lay_R.exfiltration_mm(:);
        lay_boundary_out_m3 = -lay_R.boundary_inflow_m3;
        GW_States.n_updates = GW_States.n_updates + 1;
        GW_States.adaptive_accepted_steps = GW_States.adaptive_accepted_steps + lay_R.substeps;
        GW_States.last_groundwater_dt_s = GW_States.elapsed_s;
        GW_States.elapsed_s = 0;
    elseif flags.flag_baseflow ~= 1 && lay_due
        GW_States.elapsed_s = 0;
    end
    lay_exfiltration(~lay_valid) = 0;
    if ~(isfield(flags, 'flag_groundwater_surface_coupling') && flags.flag_groundwater_surface_coupling == 0)
        lay_surface = lay_surface + lay_exfiltration;
    else
        lay_boundary_out_m3 = lay_boundary_out_m3 + sum(lay_exfiltration)/1000*lay_area;
    end
end

% Write back the state and the legacy views.
Subsurface.water_mm = lay_water;
lay_head_map = reshape(lay_head, lay_ny, lay_nx);
lay_head_map(idx_nan) = nan;
BC_States.h_t = lay_head_map;
BC_States.h_0 = lay_head_map;
lay_surface_map = reshape(max(lay_surface, 0), lay_ny, lay_nx);
lay_surface_map(idx_nan) = depths.d_t(idx_nan);
depths.d_t = lay_surface_map;
Soil_Properties = hp2d_sync_layered_soil(Soil_Properties, Subsurface, lay_head, idx_nan);

Hydro_States.f = reshape(lay_S.infiltration_mm/max(dt_h, 1e-12), size(elevation));
Hydro_States.f(idx_nan) = nan;
Hydro_States.ETR = reshape(lay_S.soil_et_mm/max(dt_days, eps), size(elevation));
Hydro_States.ETR(idx_nan) = nan;
BC_States.delta_E = reshape(lay_S.surface_evaporation_mm, size(elevation));
BC_States.delta_E(idx_nan) = nan;
cumulative_infiltration = cumulative_infiltration + reshape(lay_S.infiltration_mm, size(elevation))/1000;
cumulative_recharge = cumulative_recharge + reshape(max(lay_recharge, 0), size(elevation));
inf_volume = sum(lay_S.infiltration_mm)/1000*lay_area;
GW_States.total_recharge_m3 = GW_States.total_recharge_m3 + sum(lay_recharge)/1000*lay_area;
GW_States.total_capillary_m3 = GW_States.total_capillary_m3 + sum(lay_capillary)/1000*lay_area;
GW_States.total_exfiltration_m3 = GW_States.total_exfiltration_m3 + sum(lay_exfiltration)/1000*lay_area;
GW_States.last_surface_recharge_rate_m_s = reshape(lay_recharge/1000/dt_s, size(elevation));
GW_States.last_capillary_rate_m_s = reshape(lay_capillary/1000/dt_s, size(elevation));
GW_States.last_net_exchange_rate_m_s = reshape((lay_recharge - lay_capillary)/1000/dt_s, size(elevation));
GW_States.last_q_exf_m_s = reshape(lay_exfiltration/1000/dt_s, size(elevation));
GW_States.last_boundary_outflow_m3 = lay_boundary_out_m3;
if ~isfield(GW_States, 'total_boundary_outflow_m3'), GW_States.total_boundary_outflow_m3 = 0; end
GW_States.total_boundary_outflow_m3 = GW_States.total_boundary_outflow_m3 + lay_boundary_out_m3;
max_GW_depth = max(max_GW_depth, BC_States.h_t - BC_States.z0);

% Module closure (m3): column + surface change against ET, evaporation, deep-aquifer
% recharge and return, river exchange stays inside, lateral boundary outflow.
lay_storage_after_m3 = sum(hp2d_column_storage_mm(Subsurface.column, lay_head, lay_water), ...
    'omitnan')/1000*lay_area;
lay_surface_after_m3 = sum(max(lay_surface, 0).*lay_valid)/1000*lay_area;
errors(3) = 0;
errors(4) = 0;
errors(5) = (lay_storage_after_m3 + lay_surface_after_m3) - ...
    (lay_storage_before_m3 + lay_surface_before_m3) + ...
    (sum(lay_S.soil_et_mm) + sum(lay_S.surface_evaporation_mm) + sum(lay_deep_recharge_mm) ...
    - sum(lay_deep_added_mm))/1000*lay_area + lay_boundary_out_m3;
