%% Checking Mass Balance Calculations
% Developer: Marcus Nobrega
% Goal: Check inflows, storage, and outflows in the domain
% Not currently working for the stage_hydrograph case

% Storage Calculation
if flags.flag_subgrid ~= 1
    C_a = ones(ny,nx)*Wshed_Properties.cell_area; % m2
    C_a(isnan(Elevation_Properties.elevation_cell)) = nan;
end

% Runoff Coefficient Calculation
BC_States.outflow_volume  = nansum(nansum(outlet_states.outlet_flow.*C_a))/1000/3600*time_step*60 + BC_States.outflow_volume ;
if flags.flag_stage_hydrograph == 1
    % inflow_stage = sum(sum(stage_cells*Wshed_Properties.cell_area.*(stage_depth - stage_depth_previous))); % m3
    inflow_stage = sum(sum(stage_cells*Wshed_Properties.cell_area.*(stage_depth - stage_depth_previous))) + ...
                    sum(sum(CA_States.I_tot_end_cell(stage_cells)));
else
    inflow_stage = 0;
end
if flags.flag_spatial_rainfall == 1
    if flags.flag_inflow == 1
        if flags.flag_subgrid == 1
            inflow_vol = nansum(nansum(BC_States.inflow.*(C_a)./Wshed_Properties.cell_area/1000*Wshed_Properties.cell_area)) + ...
            nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area + inflow_stage;
        else
        inflow_vol = nansum(nansum(BC_States.inflow/1000*Wshed_Properties.cell_area)) + ...
            nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area + inflow_stage;
        end
    else
        inflow_vol = nansum(nansum(nansum(nansum(BC_States.delta_p_agg))/1000*Wshed_Properties.cell_area)) + inflow_stage;
    end
    BC_States.inflow_volume = inflow_vol + BC_States.inflow_volume; % m3
elseif flags.flag_spatial_rainfall ~= 1 && flags.flag_inflow == 1
    inflow_vol = nansum(nansum(BC_States.inflow.*(C_a)./Wshed_Properties.cell_area/1000*Wshed_Properties.cell_area)) + nansum(nansum(BC_States.delta_p_agg/1000.*C_a)) + inflow_stage;
    BC_States.inflow_volume = inflow_vol +  BC_States.inflow_volume; % check future
else
    inflow_vol = nansum(nansum(BC_States.inflow/1000.*C_a)) + nansum(nansum(BC_States.delta_p_agg/1000.*C_a)) + inflow_stage ;
    BC_States.inflow_volume = inflow_vol  + BC_States.inflow_volume; % check future
end

if flags.flag_subgrid == 1
    current_storage = nansum(nansum((Wshed_Properties.Resolution - Wshed_Properties.River_Width).*Wshed_Properties.Resolution.*max((depths.d_t/1000 - Wshed_Properties.River_Depth),0))) + ...
                      nansum(nansum(Wshed_Properties.Resolution.*Wshed_Properties.River_Width.*depths.d_t/1000)); % m3
else
    current_storage = nansum(nansum(C_a.*depths.d_t/1000)); % m3
end

if flags.flag_infiltration == 0
    inf_volume = 0;
end
% dS/dt = Qin - Qout = Rain + Inflow - Outflow - infiltration
delta_storage = current_storage - previous_storage;
outflow_vol = nansum(nansum(outlet_states.outlet_flow.*C_a))/1000/3600*time_step*60 + ...
    inf_volume; % m3 
flux_volumes = inflow_vol - outflow_vol; % dt(Qin - Qout) m3
volume_error = flux_volumes - delta_storage;

% if inflow_vol ~= 0
%     if abs(volume_error)/inflow_vol*100 > 100
%         factor_time = 1;
%         catch_index = catch_index  + 1;
%         error('Mass balance error larger than 10%')
%     end
% end

% Previous storage
previous_storage = current_storage;

% Introducing the volume error to the inflow cells
if flags.flag_inflow == 1 && flags.flag_waterbalance == 1
    mask = Wshed_Properties.inflow_mask;
    if any(any(mask))
        depths.d_t(mask) = depths.d_t(mask) - 1000*(volume_error)/sum(sum(mask))/Wshed_Properties.cell_area;
    end
end
