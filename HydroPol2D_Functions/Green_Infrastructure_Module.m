%% =========================================================================
%  GREEN_INFRASTRUCTURE_MODULE (flag_green_infrastructure = 1)
%  =========================================================================
%  Subgrid green-infrastructure compartments (bioretention, infiltration
%  trenches, permeable pavement), as in HydroPol2D-Python
%  core.engine._advance_green_infrastructure. Runs in Hydrological_Model after
%  the host cell's snow, canopy, infiltration, ET and groundwater, and before
%  the inflow hydrograph and the routing. hp2d_green_infrastructure_step moves
%  water between the host surface (depths.d_t) and the GI stores, then:
%    - the native exfiltration request enters the soil below the cell: the
%      unsaturated soil layers top-down (add_layered_infiltration), or the
%      three soil slots of the layered column (flag_layered_subsurface); what
%      the soil cannot take goes back to the gravel;
%    - underdrain water goes to its target surface cell, or leaves the model;
%      HydroPol2D_Main_While books that volume as outlet flow at the host
%      cells after the routing (hp2d_reservoir_book_external).
%  The host share of the rain, ET and LAI is set by
%  hp2d_green_infrastructure_scale_forcing before the vertical step. GI
%  storage and ET enter the mass balance (mass_balance_check,
%  hp2d_update_system_mass_ledger).
%  State: GI_States (surface_m3, media_m3, gravel_m3, edge_q_m3_s per
%  compartment or edge, cumulative volumes and the fluxes of the last step).
% =========================================================================

gi_dt_s = time_step*60;
gi_area = Wshed_Properties.cell_area;
gi_n_cells = numel(depths.d_t);

% Host surface water and gross forcing at the GI cells.
gi_host_mm = double(gather(depths.d_t(Green_Infrastructure.host_cells)));
gi_host_mm(~isfinite(gi_host_mm)) = 0;
gi_host_m3 = max(gi_host_mm, 0)/1000*gi_area;
% Lumped forcing is a scalar, spatial forcing a map.
gi_at_cells = @(x) double(gather(x(min(Green_Infrastructure.cell_index, numel(x)))));
gi_rain_rate = gi_at_cells(BC_States.delta_p_agg)*3600/gi_dt_s;
gi_rain_rate(~isfinite(gi_rain_rate)) = 0;
gi_et_rate = zeros(Green_Infrastructure.compartment_count, 1);
gi_apply_kc = false;
if flags.flag_ETP == 1
    if flags.flag_input_ETP_map == 1
        % Prescribed ET maps: an actual demand, no crop coefficient (Python's
        % prescribed soil ET).
        if isfield(BC_States, 'current_evaporation_map') && ~isempty(BC_States.current_evaporation_map)
            gi_et_map = BC_States.current_evaporation_map + BC_States.current_transpiration_map;
        else
            gi_et_map = input_evaporation + input_transpiration;
        end
    else
        % Reference ET from the station records: the design crop coefficient applies.
        gi_et_map = Hydro_States.ETP;
        gi_apply_kc = true;
    end
    gi_et_rate = gi_at_cells(gi_et_map);
    gi_et_rate(~isfinite(gi_et_rate)) = 0;
end
gi_layered = isfield(flags, 'flag_layered_subsurface') && flags.flag_layered_subsurface == 1 && ...
    exist('Subsurface', 'var') && isstruct(Subsurface);
gi_recording = isfield(Green_Infrastructure, 'record_steps') && Green_Infrastructure.record_steps;
gi_subsurface = [];
if gi_layered, gi_subsurface = Subsurface; end
if gi_recording
    % Diagnostic run: record the GI inputs (for a replay) and the water the
    % module moves between the host surface, the soil and the GI stores.
    gi_record_in = struct('host_m3', gi_host_m3, 'rain_rate_mm_h', gi_rain_rate, ...
        'et_rate_mm_day', gi_et_rate, 'dt_s', gi_dt_s);
    gi_balance0 = [nansum(double(gather(depths.d_t(:)))), ...
        hp2d_gi_soil_storage_mm(Soil_Properties, gi_subsurface, BC_States, gi_layered)]/1000*gi_area;
    gi_store0 = sum(GI_States.surface_m3) + sum(GI_States.media_m3) + sum(GI_States.gravel_m3);
end

[GI_States, gi_F, gi_host_after_m3] = hp2d_green_infrastructure_step(Green_Infrastructure, ...
    GI_States, gi_host_m3, gi_rain_rate, gi_et_rate, gi_apply_kc, gi_dt_s);

% Native exfiltration into the soil below the cell (Python commit 2d9431f).
gi_request_cell_m3 = accumarray(Green_Infrastructure.cell_index, gi_F.native_request, [gi_n_cells 1]);
gi_request_mm = reshape(gi_request_cell_m3/gi_area*1000, size(depths.d_t));
if any(gi_request_cell_m3 > 0)
    if gi_layered
        % Layered column: the three soil slots above the water table, top-down.
        gi_valid = Subsurface.valid;
        gi_head = BC_States.h_t(:);
        gi_head(~gi_valid) = 0;
        gi_view = hp2d_soil_view(Subsurface.column, gi_head, Subsurface.water_mm, ...
            LULC_Properties.root_depth_m(:));
        [gi_water3, gi_accepted_mm] = hp2d_add_slot_infiltration(gi_view, gi_request_mm(:), gi_valid);
        Subsurface.water_mm(gi_valid, 1:3) = gi_water3(gi_valid, :);
        Soil_Properties = hp2d_sync_layered_soil(Soil_Properties, Subsurface, gi_head, idx_nan);
    elseif isfield(Soil_Properties, 'Layers') && isstruct(Soil_Properties.Layers) && ...
            isfield(Soil_Properties.Layers, 'near_surface_storage_mm')
        % Legacy soil: near-surface, root-zone and transmission layers, top-down.
        [Soil_Properties, gi_accepted_mm] = add_layered_infiltration(Soil_Properties, ...
            gi_request_mm, idx_nan);
    else
        % Single soil bucket: room left above the water table (Infiltration_Module).
        gi_perv = ones(size(depths.d_t));
        if isfield(LULC_Properties, 'frac_perv') && ~isempty(LULC_Properties.frac_perv)
            gi_perv = min(max(double(LULC_Properties.frac_perv), 0), 1);
        elseif isfield(LULC_Properties, 'idx_imp') && ~isempty(LULC_Properties.idx_imp)
            gi_perv = double(~LULC_Properties.idx_imp);
        end
        gi_zwt = Soil_Properties.Soil_Depth - (BC_States.h_t - (elevation - Soil_Properties.Soil_Depth));
        gi_zwt = min(max(gi_zwt, 0), Soil_Properties.Soil_Depth);
        gi_room_mm = max(gi_perv.*gi_zwt.*(Soil_Properties.theta_sat - Soil_Properties.theta_r)*1000 - ...
            Soil_Properties.I_t, 0);
        gi_accepted_mm = min(gi_request_mm, gi_room_mm);
        gi_accepted_mm(idx_nan) = 0;
        Soil_Properties.I_t = Soil_Properties.I_t + cast(gi_accepted_mm, 'like', Soil_Properties.I_t);
    end
    gi_accepted_cell_m3 = double(gather(gi_accepted_mm(:)))/1000*gi_area;
    gi_accepted_cell_m3(~isfinite(gi_accepted_cell_m3)) = 0;
else
    gi_accepted_cell_m3 = zeros(gi_n_cells, 1);
end
gi_fraction = zeros(gi_n_cells, 1);
gi_has = gi_request_cell_m3 > 0;
gi_fraction(gi_has) = min(1, gi_accepted_cell_m3(gi_has)./gi_request_cell_m3(gi_has));
gi_native = gi_F.native_request.*gi_fraction(Green_Infrastructure.cell_index);
GI_States.gravel_m3 = GI_States.gravel_m3 + max(gi_F.native_request - gi_native, 0);

% Host surface: capture out, overflow and surface-target underdrain in.
gi_delta_mm = (gi_host_after_m3 - gi_host_m3)/gi_area*1000;
depths.d_t(Green_Infrastructure.host_cells) = depths.d_t(Green_Infrastructure.host_cells) + ...
    cast(gi_delta_mm, 'like', depths.d_t);

% Underdrain water that leaves the model, at the compartments' host cells.
gi_external = gi_F.underdrain.*(Green_Infrastructure.target_kind(Green_Infrastructure.facility_index) == 3);
GI_States.step_external_m3_host = accumarray(Green_Infrastructure.compartment_host, gi_external, ...
    [numel(Green_Infrastructure.host_cells) 1]);
GI_States.step_total_external_m3 = gi_F.external_m3;
GI_States.step_total_et_m3 = sum(gi_F.et);
GI_States.step_dt_s = gi_dt_s;

gi_names = {'rain', 'capture', 'infiltration', 'overflow', 'et', 'underdrain'};
for gi_i = 1:numel(gi_names)
    GI_States.(['step_' gi_names{gi_i} '_m3']) = gi_F.(gi_names{gi_i});
    GI_States.(['cumulative_' gi_names{gi_i} '_m3']) = ...
        GI_States.(['cumulative_' gi_names{gi_i} '_m3']) + gi_F.(gi_names{gi_i});
end
GI_States.step_native_exfiltration_m3 = gi_native;
GI_States.cumulative_native_exfiltration_m3 = GI_States.cumulative_native_exfiltration_m3 + gi_native;
GI_States.cumulative_external_m3 = GI_States.cumulative_external_m3 + gi_F.external_m3;
GI_States.gravel_edge_q_m3_s = gi_F.gravel_edge_q_m3_s;

if gi_recording
    if gi_layered, gi_subsurface = Subsurface; end
    gi_balance1 = [nansum(double(gather(depths.d_t(:)))), ...
        hp2d_gi_soil_storage_mm(Soil_Properties, gi_subsurface, BC_States, gi_layered)]/1000*gi_area;
    gi_store1 = sum(GI_States.surface_m3) + sum(GI_States.media_m3) + sum(GI_States.gravel_m3);
    gi_record_in.coupling_error_m3 = sum(gi_balance1 - gi_balance0) + (gi_store1 - gi_store0) + ...
        gi_F.external_m3 + sum(gi_F.et) - sum(gi_F.rain);
    % The same balance from the transfers the module books (double), and
    % from the model states (above; their class is that of the inputs).
    gi_record_in.bookkeeping_error_m3 = sum(gi_delta_mm)/1000*gi_area + sum(gi_native) + ...
        (gi_store1 - gi_store0) + gi_F.external_m3 + sum(gi_F.et) - sum(gi_F.rain);
    gi_record_in.surface_gain_m3 = gi_balance1(1) - gi_balance0(1);
    gi_record_in.soil_gain_m3 = gi_balance1(2) - gi_balance0(2);
    gi_record_in.store_gain_m3 = gi_store1 - gi_store0;
    gi_record_in.native_accepted_fraction = gi_fraction(Green_Infrastructure.cell_index);
    gi_record_in.surface_m3 = GI_States.surface_m3;
    gi_record_in.media_m3 = GI_States.media_m3;
    gi_record_in.gravel_m3 = GI_States.gravel_m3;
    gi_record_in.edge_q_m3_s = GI_States.edge_q_m3_s;
    gi_record_in.host_after_m3 = gi_host_after_m3;
    for gi_i = 1:numel(gi_names)
        gi_record_in.(gi_names{gi_i}) = gi_F.(gi_names{gi_i});
    end
    gi_record_in.native_request = gi_F.native_request;
    gi_record_in.native = gi_native;
    gi_record_in.external_m3 = gi_F.external_m3;
    gi_record_in.t_min = t;
    GI_States.step_log{end + 1, 1} = gi_record_in;
end

clear gi_host_mm gi_host_m3 gi_rain_rate gi_et_rate gi_et_map gi_F gi_host_after_m3 ...
    gi_request_cell_m3 gi_request_mm gi_accepted_mm gi_accepted_cell_m3 gi_fraction gi_has ...
    gi_native gi_delta_mm gi_external gi_valid gi_head gi_view gi_water3 gi_perv gi_zwt ...
    gi_room_mm gi_record_in gi_i gi_names gi_apply_kc gi_at_cells gi_layered gi_recording ...
    gi_balance0 gi_balance1 gi_store0 gi_store1 gi_subsurface

function storage_mm = hp2d_gi_soil_storage_mm(Soil_Properties, Subsurface, BC_States, layered)
%HP2D_GI_SOIL_STORAGE_MM Soil water summed over the grid [mm x cells] (diagnostic).
if layered
    storage_mm = sum(double(gather(hp2d_column_storage_mm(Subsurface.column, BC_States.h_t(:), ...
        Subsurface.water_mm))), 'omitnan');
else
    storage_mm = nansum(double(gather(Soil_Properties.I_t(:))));
end
end
