function [ledger, step_residual_m3] = hp2d_update_system_mass_ledger( ...
    ledger, flags, BC_States, Hydro_States, Snow_Properties, Soil_Properties, GW_States, ...
    depths, outlet_states, Wshed_Properties, Elevation_Properties, time_step, ...
    C_a, SubgridTables)
%HP2D_UPDATE_SYSTEM_MASS_LEDGER Non-intrusive event water-balance ledger.
%
% The legacy mass_balance_check script is a progress-report diagnostic. This
% helper is deliberately separate: it evaluates the same water stores and
% fluxes at every accepted model step, without changing any model state.

coarse_cell_area = ones(size(Elevation_Properties.elevation_cell), 'like', depths.d_t) ...
    .* Wshed_Properties.cell_area;
coarse_cell_area(~isfinite(Elevation_Properties.elevation_cell)) = NaN;

use_exact_subgrid = flags.flag_subgrid == 1 && flags.flag_overbanks ~= 1 && ...
    ~isempty(SubgridTables) && isfield(SubgridTables, 'sfincs_exact') && ...
    SubgridTables.sfincs_exact;

% Fluxes for the step that has just been accepted.
if ~isfield(ledger, 'initial_storage_m3')
    % Backward-compatible initialization for callers that constructed the
    % optional diagnostic ledger before the event-scale fields were added.
    ledger.initial_storage_m3 = ledger.previous_storage_m3;
end
rain_m3 = nansum(nansum(BC_States.delta_p_agg ./ 1000 .* coarse_cell_area));
inflow_m3 = nansum(nansum(BC_States.inflow ./ 1000 .* coarse_cell_area));
if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    inflow_m3 = nansum(nansum(BC_States.inflow .* coarse_cell_area ./ 1000));
end

stage_m3 = 0;
if flags.flag_stage_hydrograph == 1 && isfield(BC_States, 'stage_cells') && ...
        isfield(BC_States, 'stage_depth') && isfield(BC_States, 'stage_depth_previous')
    stage_m3 = nansum(nansum(BC_States.stage_cells .* coarse_cell_area .* ...
        (BC_States.stage_depth - BC_States.stage_depth_previous)));
end

prescribed_recharge_m3 = 0;
if isfield(flags, 'flag_prescribed_recharge') && flags.flag_prescribed_recharge == 1 && ...
        isfield(BC_States, 'last_groundwater_recharge_volume_m3')
    prescribed_recharge_m3 = max(double(BC_States.last_groundwater_recharge_volume_m3), 0);
end

canopy_evap_m3 = nansum(nansum(coarse_cell_area .* Hydro_States.E_int ./ 1000));
etr_m3 = nansum(nansum(coarse_cell_area .* Hydro_States.ETR ./ 1000 ...
    .* (time_step / 60 / 24)));
open_water_evap_m3 = nansum(nansum(coarse_cell_area .* BC_States.delta_E ./ 1000));
snow_sublimation_m3 = nansum(nansum(coarse_cell_area .* Snow_Properties.E_s ./ 1000));
outlet_m3 = nansum(nansum(outlet_states.outlet_flow .* coarse_cell_area)) ...
    ./ 1000 ./ 3600 .* time_step .* 60;

% Current total storage.
canopy_storage_m3 = nansum(nansum(coarse_cell_area .* Hydro_States.S ./ 1000));
if use_exact_subgrid
    eta_storage = SubgridTables.z_zmin + max(depths.d_t ./ 1000, 0);
    surface_storage_m3 = nansum(nansum(hp2d_sfincs_cell_volume_from_zs( ...
        SubgridTables, eta_storage)));
elseif flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    surface_storage_m3 = nansum(nansum(hp2d_neal_cell_volume( ...
        max(depths.d_t ./ 1000, 0), Wshed_Properties.River_Width, ...
        Wshed_Properties.River_Depth, Wshed_Properties.Resolution)));
else
    surface_storage_m3 = nansum(nansum(coarse_cell_area .* depths.d_t ./ 1000));
end
soil_storage_m3 = nansum(nansum(coarse_cell_area .* Soil_Properties.I_t ./ 1000));
groundwater_storage_m3 = nansum(nansum(coarse_cell_area .* Soil_Properties.Sy .* ...
    (BC_States.h_t - (Elevation_Properties.elevation_cell - Soil_Properties.Soil_Depth))));
if flags.flag_groundwater_modeling == 1 && isstruct(GW_States) && ...
        isfield(GW_States, 'pending_net_exchange_m')
    groundwater_storage_m3 = groundwater_storage_m3 + nansum(nansum( ...
        coarse_cell_area .* GW_States.pending_net_exchange_m));
end
snow_storage_m3 = nansum(nansum(coarse_cell_area .* Snow_Properties.SWE_t ./ 1000));
storage_m3 = canopy_storage_m3 + surface_storage_m3 + soil_storage_m3 + ...
    groundwater_storage_m3 + snow_storage_m3;

net_flux_m3 = rain_m3 + inflow_m3 + stage_m3 + prescribed_recharge_m3 ...
    - canopy_evap_m3 - etr_m3 - open_water_evap_m3 - outlet_m3 - snow_sublimation_m3;
step_residual_m3 = (storage_m3 - ledger.previous_storage_m3) - net_flux_m3;

ledger.previous_storage_m3 = storage_m3;
ledger.final_storage_m3 = storage_m3;
ledger.cumulative_precipitation_m3 = ledger.cumulative_precipitation_m3 + rain_m3;
ledger.cumulative_boundary_inflow_m3 = ledger.cumulative_boundary_inflow_m3 + inflow_m3 + stage_m3;
ledger.cumulative_prescribed_recharge_m3 = ledger.cumulative_prescribed_recharge_m3 + prescribed_recharge_m3;
ledger.cumulative_canopy_evaporation_m3 = ledger.cumulative_canopy_evaporation_m3 + canopy_evap_m3;
ledger.cumulative_etr_m3 = ledger.cumulative_etr_m3 + etr_m3;
ledger.cumulative_open_water_evaporation_m3 = ledger.cumulative_open_water_evaporation_m3 + open_water_evap_m3;
ledger.cumulative_snow_sublimation_m3 = ledger.cumulative_snow_sublimation_m3 + snow_sublimation_m3;
ledger.cumulative_outlet_m3 = ledger.cumulative_outlet_m3 + outlet_m3;
ledger.cumulative_residual_m3 = ledger.cumulative_residual_m3 + step_residual_m3;
ledger.step_count = ledger.step_count + 1;
end
