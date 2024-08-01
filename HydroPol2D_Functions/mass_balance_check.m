%% Checking Mass Balance Calculations
% Developer: Marcus Nobrega
% Goal: Check inflows, storage, and outflows in the domain
% Not currently working for the stage_hydrograph case

% Runoff Coefficient Calculation
BC_States.outflow_volume  = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600*time_step*60 + BC_States.outflow_volume ;
if flags.flag_spatial_rainfall == 1
    if flags.flag_inflow == 1
        inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area)) + ...
            nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area;
    else
        inflow_vol = nansum(nansum(nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area));
    end
    BC_States.inflow_volume = inflow_vol + BC_States.inflow_volume; % m3
elseif flags.flag_spatial_rainfall ~= 1 && flags.flag_inflow == 1
    inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area)) + BC_States.delta_p_agg/1000*Wshed_Properties.drainage_area;
    BC_States.inflow_volume = inflow_vol +  BC_States.inflow_volume; % check future
else
    inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area)) + BC_States.delta_p_agg/1000*Wshed_Properties.drainage_area ;
    BC_States.inflow_volume = inflow_vol + BC_States.inflow_volume; % check future
end

% Storage Calculation
previous_storage = current_storage;
current_storage = nansum(nansum(Wshed_Properties.cell_area*depths.d_t/1000)); % m3
if flags.flag_infiltration == 0
    Hydro_States.f = zeros(size(outlet_states.outlet_flow,1),size(outlet_states.outlet_flow,2));
end
% dS/dt = Qin - Qout = Rain + Inflow - Outflow - infiltration
delta_storage = current_storage - previous_storage;
outflow_flux = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600 + ...
    nansum(nansum(Hydro_States.f/1000/3600*Wshed_Properties.cell_area)); % m3 per sec
flux_volumes = inflow_vol - outflow_flux*time_step*60; % dt(Qin - Qout) m3
volume_error = flux_volumes - delta_storage;

% Introducing the volume error to the inflow cells
if flags.flag_inflow == 1 && flags.flag_waterbalance == 1
    mask = Wshed_Properties.inflow_mask;
    if any(any(mask))
        depths.d_t(mask) = depths.d_t(mask) - 1000*(volume_error)/sum(sum(mask))/Wshed_Properties.cell_area;
    end
end
