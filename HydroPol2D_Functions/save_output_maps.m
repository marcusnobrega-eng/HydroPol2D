%% Save Output Maps
% Goal: Save model output maps
% Developer: Marcus Nobrega, Ph.D.
%
% Updated logic in this version
% -----------------------------
% 1) Maps are saved ONLY here, when a true save event happens.
% 2) update_spatial no longer writes directly into Maps.Hydro.* every step.
% 3) Current rainfall / evaporation / transpiration maps are taken from:
%       - BC_States.current_rainfall_mean_map
%       - BC_States.current_evaporation_map
%       - BC_States.current_transpiration_map
% 4) ETP / ETR maps saved here come from the current Hydro_States fields.
% 5) saver_count only advances when recording_parameters.delta_record > 0.
%
% Assumptions
% -----------
% - update_spatial now populates the current forcing/state fields.
% - This routine remains the single place where map stacks are written.

% Maps with generally coarser resolution
recording_parameters.actual_record_state = find(running_control.time_records < t_save,1,'last');
recording_parameters.delta_record = recording_parameters.actual_record_state - recording_parameters.last_record_maps;
recording_parameters.last_record_maps = recording_parameters.actual_record_state;

if recording_parameters.delta_record > 0
    ttt = 1; %#ok<NASGU>
end

set(0, 'DefaultFigureWindowStyle', 'normal')

if flags.flag_automatic_calibration ~= 1

    % =====================================================================
    % FIRST TIME-STEP
    % =====================================================================
    if k == 1

        % -----------------------------------------------------------------
        % Core hydro maps
        % -----------------------------------------------------------------
        Maps.Hydro.d(:,:,1) = depths.d_t;

        if flags.flag_human_instability == 1
            Maps.Hydro.risk(:,:,1) = Human_Instability.risk_t;
        elseif flags.flag_human_instability == 2
        elseif flags.flag_human_instability == 3
            Maps.Hydro.risk_cm(:,:,1) = Human_Instability.risk_t_cm;
            Maps.Hydro.risk_tm(:,:,1) = Human_Instability.risk_t_tm;
            Maps.Hydro.risk_am(:,:,1) = Human_Instability.risk_t_am;
            Maps.Hydro.risk_om(:,:,1) = Human_Instability.risk_t_om;
            Maps.Hydro.risk_cf(:,:,1) = Human_Instability.risk_t_cf;
            Maps.Hydro.risk_tf(:,:,1) = Human_Instability.risk_t_tf;
            Maps.Hydro.risk_af(:,:,1) = Human_Instability.risk_t_af;
            Maps.Hydro.risk_of(:,:,1) = Human_Instability.risk_t_of;
        end

        if flags.flag_infiltration == 1
            Maps.Hydro.I_t(:,:,1) = Soil_Properties.I_t;
            Maps.Hydro.C(:,:,1)   = C;
            Maps.Hydro.f(:,:,1)   = Hydro_States.f;
        end

        if flags.flag_ETP == 1
            Maps.Hydro.ETP_save(:,:,1) = Hydro_States.ETP;
            Maps.Hydro.ETR_save(:,:,1) = Hydro_States.ETR;
        end

        % -----------------------------------------------------------------
        % Current rainfall / evaporation / transpiration maps
        % These are now written here, not in update_spatial
        % -----------------------------------------------------------------
        if flags.flag_rainfall > 0 && ...
           (flags.flag_spatial_rainfall == 1 || flags.flag_input_rainfall_map == 1 || flags.flag_satellite_rainfall == 1)

            if isfield(BC_States,'current_rainfall_mean_map') && ~isempty(BC_States.current_rainfall_mean_map)
                Maps.Hydro.spatial_rainfall_maps(:,:,1) = BC_States.current_rainfall_mean_map;
            else
                Maps.Hydro.spatial_rainfall_maps(:,:,1) = zeros(size(depths.d_t));
            end

            if isfield(BC_States,'current_avg_spatial_rainfall') && ~isempty(BC_States.current_avg_spatial_rainfall)
                BC_States.average_spatial_rainfall(1,1) = BC_States.current_avg_spatial_rainfall;
            else
                BC_States.average_spatial_rainfall(1,1) = 0;
            end
        end

        if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1
            if isfield(BC_States,'current_evaporation_map') && ~isempty(BC_States.current_evaporation_map)
                Maps.Hydro.spatial_evaporation_maps(:,:,1) = BC_States.current_evaporation_map;
            else
                Maps.Hydro.spatial_evaporation_maps(:,:,1) = zeros(size(depths.d_t));
            end

            if isfield(BC_States,'current_transpiration_map') && ~isempty(BC_States.current_transpiration_map)
                Maps.Hydro.spatial_transpiration_maps(:,:,1) = BC_States.current_transpiration_map;
            else
                Maps.Hydro.spatial_transpiration_maps(:,:,1) = zeros(size(depths.d_t));
            end

            if isfield(BC_States,'current_avg_spatial_evaporation') && ~isempty(BC_States.current_avg_spatial_evaporation)
                BC_States.average_spatial_evaporation(1,1) = BC_States.current_avg_spatial_evaporation;
            else
                BC_States.average_spatial_evaporation(1,1) = 0;
            end

            if isfield(BC_States,'current_avg_spatial_transpiration') && ~isempty(BC_States.current_avg_spatial_transpiration)
                BC_States.average_spatial_transpiration(1,1) = BC_States.current_avg_spatial_transpiration;
            else
                BC_States.average_spatial_transpiration(1,1) = 0;
            end
        end

        % -----------------------------------------------------------------
        % Interception
        % -----------------------------------------------------------------
        if flags.flag_abstraction == 1
            Maps.Hydro.Abstraction(:,:,1) = Hydro_States.S;
        end

        % -----------------------------------------------------------------
        % Groundwater depth
        % -----------------------------------------------------------------
        if flags.flag_baseflow == 1
            Maps.Hydro.GWdepth_save(:,:,1) = BC_States.h_t - elevation + Soil_Properties.Soil_Depth; % GW depth
        end

        % -----------------------------------------------------------------
        % Snowpack
        % -----------------------------------------------------------------
        if flags.flag_snow_modeling == 1
            Maps.Hydro.Snowpack(:,:,1) = Snow_Properties.H_snow_t; % mm
        end

        % -----------------------------------------------------------------
        % Water quality
        % -----------------------------------------------------------------
        if flags.flag_waterquality == 1
            Maps.WQ_States.Pol_Conc_Map(:,:,1) = zeros(ny,nx);
            Maps.WQ_States.Pol_Mass_Map(:,:,1) = zeros(ny,nx);
            Maps.WQ_States.Pol_Load_Map(:,:,1) = zeros(ny,nx);
        end

    % =====================================================================
    % SAVING MAPS ONLY WHEN A TRUE SAVE EVENT HAPPENS
    % =====================================================================
    elseif recording_parameters.delta_record > 0

        t_store = recording_parameters.actual_record_state; %#ok<NASGU>

        % -----------------------------------------------------------------
        % Core hydro maps
        % -----------------------------------------------------------------
        Maps.Hydro.d(:,:,saver_count) = depths.d_t;

        if flags.flag_human_instability == 1
            Maps.Hydro.risk(:,:,saver_count) = Human_Instability.risk_t;
        elseif flags.flag_human_instability == 2
        elseif flags.flag_human_instability == 3
            Maps.Hydro.risk_cm(:,:,saver_count) = Human_Instability.risk_t_cm;
            Maps.Hydro.risk_tm(:,:,saver_count) = Human_Instability.risk_t_tm;
            Maps.Hydro.risk_am(:,:,saver_count) = Human_Instability.risk_t_am;
            Maps.Hydro.risk_om(:,:,saver_count) = Human_Instability.risk_t_om;
            Maps.Hydro.risk_cf(:,:,saver_count) = Human_Instability.risk_t_cf;
            Maps.Hydro.risk_tf(:,:,saver_count) = Human_Instability.risk_t_tf;
            Maps.Hydro.risk_af(:,:,saver_count) = Human_Instability.risk_t_af;
            Maps.Hydro.risk_of(:,:,saver_count) = Human_Instability.risk_t_of;
        end

        if flags.flag_infiltration == 1
            Maps.Hydro.I_t(:,:,saver_count) = Soil_Properties.I_t;
            Maps.Hydro.C(:,:,saver_count)   = C;
            Maps.Hydro.f(:,:,saver_count)   = Hydro_States.f;
        end

        if flags.flag_ETP == 1
            Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP;
            Maps.Hydro.ETR_save(:,:,saver_count) = Hydro_States.ETR;
        end

        % -----------------------------------------------------------------
        % Rainfall maps: save current map only when save event happens
        % -----------------------------------------------------------------
        if flags.flag_rainfall > 0 && ...
           (flags.flag_spatial_rainfall == 1 || flags.flag_input_rainfall_map == 1 || flags.flag_satellite_rainfall == 1)

            if isfield(BC_States,'current_rainfall_mean_map') && ~isempty(BC_States.current_rainfall_mean_map)
                Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = BC_States.current_rainfall_mean_map;
            end

            if isfield(BC_States,'current_avg_spatial_rainfall') && ~isempty(BC_States.current_avg_spatial_rainfall)
                BC_States.average_spatial_rainfall(t_store,1) = BC_States.current_avg_spatial_rainfall;
            end
        end

        % -----------------------------------------------------------------
        % Map-based evaporation / transpiration
        % -----------------------------------------------------------------
        if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1

            if isfield(BC_States,'current_evaporation_map') && ~isempty(BC_States.current_evaporation_map)
                Maps.Hydro.spatial_evaporation_maps(:,:,saver_count) = BC_States.current_evaporation_map;
            end

            if isfield(BC_States,'current_transpiration_map') && ~isempty(BC_States.current_transpiration_map)
                Maps.Hydro.spatial_transpiration_maps(:,:,saver_count) = BC_States.current_transpiration_map;
            end

            if isfield(BC_States,'current_avg_spatial_evaporation') && ~isempty(BC_States.current_avg_spatial_evaporation)
                BC_States.average_spatial_evaporation(t_store,1) = BC_States.current_avg_spatial_evaporation;
            end

            if isfield(BC_States,'current_avg_spatial_transpiration') && ~isempty(BC_States.current_avg_spatial_transpiration)
                BC_States.average_spatial_transpiration(t_store,1) = BC_States.current_avg_spatial_transpiration;
            end
        end

        % -----------------------------------------------------------------
        % Groundwater
        % -----------------------------------------------------------------
        if flags.flag_baseflow == 1
            Maps.Hydro.GWdepth_save(:,:,saver_count) = BC_States.h_t - elevation + Soil_Properties.Soil_Depth; % GW depth
        end

        % -----------------------------------------------------------------
        % Snowpack
        % -----------------------------------------------------------------
        if flags.flag_snow_modeling == 1
            Maps.Hydro.Snowpack(:,:,saver_count) = Snow_Properties.H_snow_t;
        end

        % -----------------------------------------------------------------
        % Interception
        % -----------------------------------------------------------------
        if flags.flag_abstraction == 1
            Maps.Hydro.Abstraction(:,:,saver_count) = Hydro_States.S;
        end

        % -----------------------------------------------------------------
        % Water quality
        % -----------------------------------------------------------------
        if flags.flag_waterquality == 1
            Maps.WQ_States.Pol_Conc_Map(:,:,saver_count) = WQ_States.P_conc;
            Maps.WQ_States.Pol_Mass_Map(:,:,saver_count) = WQ_States.B_t; % kg
            Maps.WQ_States.Pol_Load_Map(:,:,saver_count) = ...
                1/1000 * WQ_States.P_conc .* CA_States.I_tot_end_cell / (time_step*60) / Wshed_Properties.Resolution^2; % kg/s/m2
        end
        % Plotting Data
%         if flags.flag_dashboard == 0
%             HydroPol2D_running_dashboard_plot( ...
%                 Maps, velocities.velocity_raster, DEM_raster, gauges, BC_States, ...
%                 time_step, Wshed_Properties.Resolution, C_a, saver_count);
%         end

    end
end

% =========================================================================
% Hydrographs and Pollutographs with same time resolution
% =========================================================================
recording_parameters.actual_record_hydrograph = find(running_control.time_record_hydrograph <= t_save,1,'last');
recording_parameters.delta_record_hydrograph = ...
    recording_parameters.actual_record_hydrograph - recording_parameters.last_record_hydrograph;
recording_parameters.last_record_hydrograph = recording_parameters.actual_record_hydrograph;

if flags.flag_automatic_calibration ~= 1

    % ---------------------------------------------------------------------
    % First time-step
    % ---------------------------------------------------------------------
    if k == 1
        outlet_states.outlet_hydrograph(1,1) = nansum(nansum(outlet_states.outlet_flow .* C_a)) / 1000 / 3600; % m3/s
        running_control.time_hydrograph(1,1) = running_control.time_calculation_routing(k,1) / 60;
        outlet_states.depth_outlet(1,1) = mean(depths.d_t(idx_outlet));

        % Maximum flooded area
        Flooded_Area = max(sum(sum((depths.d_t > 150))) * Wshed_Properties.Resolution^2, Flooded_Area);

        % Maximum risk area
        if flags.flag_human_instability == 1
            Risk_Area = max(sum(sum((Human_Instability.risk_t == 1))) * Wshed_Properties.Resolution^2, Risk_Area);
        end

        if flags.flag_waterquality(1,1) == 1
            WQ_States.WQ_States.outet_pollutograph(1,1) = Out_Conc; %#ok<STRNU>
            WQ_States.EMC_outlet(1,1) = WQ_States.mass_outlet / WQ_States.vol_outlet; % mg/L
            WQ_States.mass_outlet_save(1,1) = WQ_States.mass_outlet;
            WQ_States.vol_outlet_save(1,1)  = WQ_States.vol_outlet;
        end

        % Saving data of input gauges
        if flags.flag_obs_gauges == 1
            for i = 1:gauges.num_obs_gauges
                gauges.x_cell = gauges.easting_obs_gauges(i);
                gauges.y_cell = gauges.northing_obs_gauges(i);
                gauges.wse_cell(1,i)       = depths.d_t(gauges.y_cell,gauges.x_cell)/1000 + Elevation_Properties.elevation_cell(gauges.y_cell,gauges.x_cell);
                gauges.depth_cell(1,i)     = depths.d_t(gauges.y_cell,gauges.x_cell)/1000;
                gauges.hydrograph_cell(1,i)= CA_States.I_tot_end_cell(gauges.y_cell,gauges.x_cell)/(running_control.time_calculation_routing);
            end

            if flags.flag_dashboard == 1
                % --- Gauge labels are CPU-only (stored in extra_parameters) ---
                if exist('extra_parameters','var') && isfield(extra_parameters,'gauges') && ...
                        isfield(extra_parameters.gauges,'labels_observed_string') && ...
                        ~isempty(extra_parameters.gauges.labels_observed_string)
                    gaugeVarNames = string(extra_parameters.gauges.labels_observed_string);
                else
                    gaugeVarNames = "G" + string(1:size(gauges.depth_cell,2));
                end

                % --- Time (always CPU double) ---
                t_now_days = double(running_control.time_hydrograph(end)/60/24);
                if flags.flag_GPU == 1
                    t_now_days = gather(t_now_days);
                end
                ax.app.gauges_time = array2table(t_now_days + date_begin, 'VariableNames', {'Time'});

                % --- Depth row (gather if needed) ---
                depthRow = gauges.depth_cell(end,:);
                if isa(depthRow,'gpuArray'); depthRow = gather(depthRow); end
                ax.app.gauges_data = array2table(depthRow, 'VariableNames', gaugeVarNames);

                ax.app.gauges.x_coord_gauges = gauges.easting_obs_gauges_absolute;
                ax.app.gauges.y_coord_gauges = gauges.northing_obs_gauges_absolute;

                mapshow(ax.ax_d, mappoint(ax.app.gauges.x_coord_gauges', ax.app.gauges.y_coord_gauges'), ...
                    'DisplayType', 'point', ...
                    'Marker', 'o', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerSize', 10, 'LineWidth', 1.1);

                if flags.flag_reservoir == 1
                    ax.app.reservoirs = Reservoir_Data;
                    boundary.index = {};

                    for i = 1:size(ax.app.reservoirs.index,1)
                        try
                            boundary.index = [boundary.index; num2str(gather(ax.app.reservoirs.index(i)))];
                            boundary.wse_cell_1(1,i) = depths.d_t(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i))/1000 + Elevation_Properties.elevation_cell(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i));
                            boundary.depth_cell_1(1,i) = depths.d_t(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i))/1000;
                            boundary.Qin_1(1,i) = flow_rate.qin_t(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i))*(Wshed_Properties.Resolution^2/1000)/3600;

                            boundary.wse_cell_2(1,i) = depths.d_t(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i))/1000 + Elevation_Properties.elevation_cell(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i));
                            boundary.depth_cell_2(1,i) = depths.d_t(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i))/1000;
                            boundary.Qin_2(1,i) = flow_rate.qin_t(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i))*(Wshed_Properties.Resolution^2/1000)/3600;

                            boundary.wse_cell_3(1,i) = depths.d_t(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds2_index(i))/1000 + Elevation_Properties.elevation_cell(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds2_index(i));
                            boundary.depth_cell_3(1,i) = depths.d_t(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds2_index(i))/1000;
                            boundary.Qin_3(1,i) = flow_rate.qin_t(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds_2index(i))*(Wshed_Properties.Resolution^2/1000)/3600;
                        catch
                            continue
                        end
                    end

                    ax.ax_list_bc.Items = boundary.index;

                    if flags.flag_GPU == 1
                        ax.app.boundary_time = array2table(gather(double(running_control.time_hydrograph(end)/60/24)) + date_begin, 'VariableNames', {'Time'});
                    else
                        ax.app.boundary_time = array2table(double(running_control.time_hydrograph(end)/60/24) + date_begin, 'VariableNames', {'Time'});
                    end

                    try
                        ax.app.boundary_data.wse_1 = array2table(subsref({boundary.wse_cell_1(end,:), gather(boundary.wse_cell_1(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                        ax.app.boundary_data.d_1   = array2table(subsref({gather(boundary.depth_cell_1(end,:)), gather(boundary.wse_cell_1(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                        ax.app.boundary_data.Qin_1 = array2table(subsref({gather(boundary.Qin_1(end,:)), gather(boundary.Qin_1(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);

                        ax.app.boundary_data.wse_2 = array2table(subsref({boundary.wse_cell_2(end,:), gather(boundary.wse_cell_2(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                        ax.app.boundary_data.d_2   = array2table(subsref({gather(boundary.depth_cell_2(end,:)), gather(boundary.wse_cell_2(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                        ax.app.boundary_data.Qin_2 = array2table(subsref({gather(boundary.Qin_2(end,:)), gather(boundary.Qin_2(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);

                        ax.app.boundary_data.wse_3 = array2table(subsref({boundary.wse_cell_3(end,:), gather(boundary.wse_cell_3(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                        ax.app.boundary_data.d_3   = array2table(subsref({gather(boundary.depth_cell_3(end,:)), gather(boundary.wse_cell_3(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                        ax.app.boundary_data.Qin_3 = array2table(subsref({gather(boundary.Qin_3(end,:)), gather(boundary.Qin_3(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index);
                    catch
                    end

                    try
                        mapshow(ax.ax_bc, mappoint(ax.app.reservoirs.x_us', ax.app.reservoirs.y_us'), ...
                            'DisplayType', 'point', 'Marker', 'x', 'MarkerEdgeColor', 'yellow', 'MarkerSize', 10, 'LineWidth', 1.1);
                        mapshow(ax.ax_bc, mappoint(ax.app.reservoirs.x_ds1', ax.app.reservoirs.y_ds1'), ...
                            'DisplayType', 'point', 'Marker', 'x', 'MarkerEdgeColor', 'green', 'MarkerSize', 10, 'LineWidth', 1.1);
                        mapshow(ax.ax_bc, mappoint(ax.app.reservoirs.x_ds2', ax.app.reservoirs.y_ds2'), ...
                            'DisplayType', 'point', 'Marker', 'x', 'MarkerEdgeColor', 'blue', 'MarkerSize', 10, 'LineWidth', 1.1);
                    catch
                    end
                end
            end
        end

    % ---------------------------------------------------------------------
    % Save hydrograph/pollutograph records only when hydrograph record advances
    % ---------------------------------------------------------------------
    elseif recording_parameters.delta_record_hydrograph > 0

        t_store = recording_parameters.actual_record_hydrograph;

        outlet_states.outlet_hydrograph(t_store,1) = nansum(nansum(outlet_states.outlet_flow .* C_a))/1000/3600; % m3/s
        running_control.time_hydrograph(t_store,1) = t;
        outlet_states.depth_outlet(t_store,1) = mean(depths.d_t(idx_outlet));

        % Maximum flooded area
        Flooded_Area = max(sum(sum((depths.d_t > 150))) * Wshed_Properties.Resolution^2, Flooded_Area);

        % Maximum risk areas
        if flags.flag_human_instability == 1
            Risk_Area = max(sum(sum((Human_Instability.risk_t == 1))) * Wshed_Properties.Resolution^2, Risk_Area);
        end

        if flags.flag_waterquality == 1
            WQ_States.EMC_outlet(t_store,1)       = WQ_States.mass_outlet / WQ_States.vol_outlet; % mg/L
            WQ_States.mass_outlet_save(t_store,1) = WQ_States.mass_outlet;
            WQ_States.vol_outlet_save(t_store,1)  = WQ_States.vol_outlet;
        end

        if flags.flag_obs_gauges == 1
            for i = 1:gauges.num_obs_gauges
                gauges.x_cell = gauges.easting_obs_gauges(i);
                gauges.y_cell = gauges.northing_obs_gauges(i);
                gauges.wse_cell(t_store,i)        = depths.d_t(gauges.y_cell,gauges.x_cell)/1000 + Elevation_Properties.elevation_cell(gauges.y_cell,gauges.x_cell);
                gauges.depth_cell(t_store,i)      = depths.d_t(gauges.y_cell,gauges.x_cell)/1000;
                gauges.hydrograph_cell(t_store,i) = CA_States.I_tot_end_cell(gauges.y_cell,gauges.x_cell)/(running_control.time_calculation_routing);
            end

            if flags.flag_dashboard == 1
                % --- Gauge labels (CPU-only) ---
                if exist('extra_parameters','var') && isfield(extra_parameters,'gauges') && ...
                        isfield(extra_parameters.gauges,'labels_observed_string') && ...
                        ~isempty(extra_parameters.gauges.labels_observed_string)
                    gaugeVarNames = string(extra_parameters.gauges.labels_observed_string);
                else
                    gaugeVarNames = "G" + string(1:size(gauges.depth_cell,2));
                end

                % --- Append last row (gather if needed) ---
                depthRow = gauges.depth_cell(end,:);
                if isa(depthRow,'gpuArray'); depthRow = gather(depthRow); end
                ax.app.gauges_data = [ax.app.gauges_data; array2table(depthRow, 'VariableNames', gaugeVarNames)];

                if flags.flag_reservoir == 1
                    Reservoir_Data = ax.app.reservoirs;
                    for i = 1:size(ax.app.reservoirs.index,1)
                        try
                            boundary.wse_cell_1(t_store,i) = depths.d_t(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i))/1000 + Elevation_Properties.elevation_cell(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i));
                            boundary.depth_cell_1(t_store,i) = depths.d_t(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i))/1000;
                            boundary.Qin_1(t_store,i) = flow_rate.qin_t(Reservoir_Data.y_index(i),Reservoir_Data.x_index(i))*(Wshed_Properties.Resolution^2/1000)/3600;

                            boundary.wse_cell_2(t_store,i) = depths.d_t(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i))/1000 + Elevation_Properties.elevation_cell(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i));
                            boundary.depth_cell_2(t_store,i) = depths.d_t(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i))/1000;
                            boundary.Qin_2(t_store,i) = flow_rate.qin_t(Reservoir_Data.y_ds1_index(i),Reservoir_Data.x_ds1_index(i))*(Wshed_Properties.Resolution^2/1000)/3600;

                            boundary.wse_cell_3(t_store,i) = depths.d_t(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds2_index(i))/1000 + Elevation_Properties.elevation_cell(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds2_index(i));
                            boundary.depth_cell_3(t_store,i) = depths.d_t(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds2_index(i))/1000;
                            boundary.Qin_3(t_store,i) = flow_rate.qin_t(Reservoir_Data.y_ds2_index(i),Reservoir_Data.x_ds_2index(i))*(Wshed_Properties.Resolution^2/1000)/3600;
                        catch
                            continue
                        end
                    end

                    try
                        ax.app.boundary_data.wse_1 = [ax.app.boundary_data.wse_1; array2table(subsref({boundary.wse_cell_1(end,:), gather(boundary.wse_cell_1(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                        ax.app.boundary_data.d_1   = [ax.app.boundary_data.d_1;   array2table(subsref({gather(boundary.depth_cell_1(end,:)), gather(boundary.wse_cell_1(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                        ax.app.boundary_data.Qin_1 = [ax.app.boundary_data.Qin_1; array2table(subsref({gather(boundary.Qin_1(end,:)), gather(boundary.Qin_1(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];

                        ax.app.boundary_data.wse_2 = [ax.app.boundary_data.wse_2; array2table(subsref({boundary.wse_cell_2(end,:), gather(boundary.wse_cell_2(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                        ax.app.boundary_data.d_2   = [ax.app.boundary_data.d_2;   array2table(subsref({gather(boundary.depth_cell_2(end,:)), gather(boundary.wse_cell_2(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                        ax.app.boundary_data.Qin_2 = [ax.app.boundary_data.Qin_2; array2table(subsref({gather(boundary.Qin_2(end,:)), gather(boundary.Qin_2(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];

                        ax.app.boundary_data.wse_3 = [ax.app.boundary_data.wse_3; array2table(subsref({boundary.wse_cell_3(end,:), gather(boundary.wse_cell_3(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                        ax.app.boundary_data.d_3   = [ax.app.boundary_data.d_3;   array2table(subsref({gather(boundary.depth_cell_3(end,:)), gather(boundary.wse_cell_3(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                        ax.app.boundary_data.Qin_3 = [ax.app.boundary_data.Qin_3; array2table(subsref({gather(boundary.Qin_3(end,:)), gather(boundary.Qin_3(end,:))}, struct('type', '{}', 'subs', {{flags.flag_GPU + 1}})), 'VariableNames', boundary.index)];
                    catch
                    end

                    if flags.flag_GPU == 1
                        ax.app.gauges_time   = [ax.app.gauges_time; array2table(gather(double(running_control.time_hydrograph(end)/60/24)) + date_begin, 'VariableNames', {'Time'})];
                        ax.app.boundary_time = ax.app.gauges_time;
                    else
                        ax.app.gauges_time   = [ax.app.gauges_time; array2table(double(running_control.time_hydrograph(end)/60/24) + date_begin, 'VariableNames', {'Time'})];
                        ax.app.boundary_time = ax.app.gauges_time;
                    end
                else
                    if flags.flag_GPU == 1
                        ax.app.gauges_time = [ax.app.gauges_time; array2table(gather(double(running_control.time_hydrograph(end)/60/24)) + date_begin, 'VariableNames', {'Time'})];
                    else
                        ax.app.gauges_time = [ax.app.gauges_time; array2table(double(running_control.time_hydrograph(end)/60/24) + date_begin, 'VariableNames', {'Time'})];
                    end
                end
            end
        end

        if flags.flag_waterquality == 1
            if flags.flag_obs_gauges == 1
                for i = 1:gauges.num_obs_gauges
                    gauges.x_cell = gauges.northing_obs_gauges(i);
                    gauges.y_cell = gauges.northing_obs_gauges(i);
                    WQ_States.pollutograph_cell(t_store,i) = WQ_States.P_conc(gauges.y_cell,gauges.x_cell); % mg/L
                end
            end
        end
    end

    % ---------------------------------------------------------------------
    % Only after a true save event do we advance saver_count and flush Maps
    % ---------------------------------------------------------------------
    if recording_parameters.delta_record > 0
        if flags.flag_dashboard == 1
            ax.timer = minutes(double(gather(running_control.time_records(recording_parameters.actual_record_state)))) + date_begin;
            ax.percentage = gather((t)/running_control.routing_time*100);
            ax = HydroPol2D_running_dashboard(ax, Maps, v0, DEM_raster, gauges, BC_States, time_step, Resolution, 1, 1, C_a);
        end

        saver_count = saver_count + 1;

        if saver_count > 12
            saver_count = 1;
            tempDir = Paths.Temp;
            save(fullfile(tempDir, ['save_map_hydro_' num2str(store) '.mat']), 'Maps', '-v7.3');
            store = store + 1;
        end
    end
end

% =========================================================================
% Saving Maximum Depths and Outlet Flows
% =========================================================================
if flags.flag_automatic_calibration ~= 1
    if k == 1
        depths.dmax_final      = depths.d_t;
        velocities.vmax_final  = velocities.velocity_raster;

        if flags.flag_waterquality == 1
            WQ_States.max_Conc = zeros(size(Elevation_Properties.elevation_cell));
        end
    else
        depths.dmax_final     = max(depths.d_t, depths.dmax_final);
        velocities.vmax_final = max(velocities.velocity_raster, velocities.vmax_final);

        if flags.flag_waterquality == 1
            WQ_States.max_Conc = max(WQ_States.max_Conc, WQ_States.P_conc);
        end
    end
end

outlet_runoff_volume = ...
    nansum(nansum(outlet_states.outlet_flow .* C_a)) * time_step / 60 / Wshed_Properties.drainage_area + outlet_runoff_volume; % mm