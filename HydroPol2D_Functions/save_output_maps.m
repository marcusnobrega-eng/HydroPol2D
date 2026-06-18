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
% 6) Lightweight time-series outputs are also snapshotted to
%       Paths.Temp/current_time_series_outputs.mat
%    whenever a hydrograph record is written.
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
        [velocity_snapshot, hazard_snapshot] = snisb_velocity_hazard_snapshot( ...
            velocities.velocity_raster, depths.d_t, idx_nan);
        Maps.Hydro.velocity(:,:,1) = velocity_snapshot;
        Maps.Hydro.hazard_dv(:,:,1) = hazard_snapshot;

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
        [velocity_snapshot, hazard_snapshot] = snisb_velocity_hazard_snapshot( ...
            velocities.velocity_raster, depths.d_t, idx_nan);
        Maps.Hydro.velocity(:,:,saver_count) = velocity_snapshot;
        Maps.Hydro.hazard_dv(:,:,saver_count) = hazard_snapshot;

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

        % -----------------------------------------------------------------
        % Preallocate lightweight time-series outputs at the first time-step
        % -----------------------------------------------------------------
        % No new user input is required. These arrays follow the existing
        % hydrograph recording time base and are later written to a temporary
        % MAT-file in Paths.Temp/current_time_series_outputs.mat.
        n_hydro_records = max(1, numel(running_control.time_record_hydrograph));

        if ~isfield(running_control,'time_hydrograph') || isempty(running_control.time_hydrograph)
            running_control.time_hydrograph = zeros(n_hydro_records,1);
        elseif size(running_control.time_hydrograph,1) < n_hydro_records
            running_control.time_hydrograph(n_hydro_records,1) = 0;
        end

        if ~isfield(outlet_states,'outlet_hydrograph') || isempty(outlet_states.outlet_hydrograph)
            if isfield(outlet_states,'outlet_flow') && isnumeric(outlet_states.outlet_flow)
                outlet_states.outlet_hydrograph = zeros(n_hydro_records,1,'like',outlet_states.outlet_flow);
            else
                outlet_states.outlet_hydrograph = zeros(n_hydro_records,1);
            end
        elseif size(outlet_states.outlet_hydrograph,1) < n_hydro_records
            outlet_states.outlet_hydrograph(n_hydro_records,1) = 0;
        end

        if ~isfield(outlet_states,'depth_outlet') || isempty(outlet_states.depth_outlet)
            outlet_states.depth_outlet = zeros(n_hydro_records,1,'like',depths.d_t);
        elseif size(outlet_states.depth_outlet,1) < n_hydro_records
            outlet_states.depth_outlet(n_hydro_records,1) = 0;
        end

        if ~isfield(BC_States,'average_spatial_rainfall_hydrograph') || isempty(BC_States.average_spatial_rainfall_hydrograph)
            if isfield(BC_States,'current_avg_spatial_rainfall') && isnumeric(BC_States.current_avg_spatial_rainfall)
                BC_States.average_spatial_rainfall_hydrograph = zeros(n_hydro_records,1,'like',BC_States.current_avg_spatial_rainfall);
            else
                BC_States.average_spatial_rainfall_hydrograph = zeros(n_hydro_records,1);
            end
        elseif size(BC_States.average_spatial_rainfall_hydrograph,1) < n_hydro_records
            BC_States.average_spatial_rainfall_hydrograph(n_hydro_records,1) = 0;
        end

        if ~isfield(BC_States,'average_spatial_evaporation_hydrograph') || isempty(BC_States.average_spatial_evaporation_hydrograph)
            if isfield(BC_States,'current_avg_spatial_evaporation') && isnumeric(BC_States.current_avg_spatial_evaporation)
                BC_States.average_spatial_evaporation_hydrograph = zeros(n_hydro_records,1,'like',BC_States.current_avg_spatial_evaporation);
            else
                BC_States.average_spatial_evaporation_hydrograph = zeros(n_hydro_records,1);
            end
        elseif size(BC_States.average_spatial_evaporation_hydrograph,1) < n_hydro_records
            BC_States.average_spatial_evaporation_hydrograph(n_hydro_records,1) = 0;
        end

        if ~isfield(BC_States,'average_spatial_transpiration_hydrograph') || isempty(BC_States.average_spatial_transpiration_hydrograph)
            if isfield(BC_States,'current_avg_spatial_transpiration') && isnumeric(BC_States.current_avg_spatial_transpiration)
                BC_States.average_spatial_transpiration_hydrograph = zeros(n_hydro_records,1,'like',BC_States.current_avg_spatial_transpiration);
            else
                BC_States.average_spatial_transpiration_hydrograph = zeros(n_hydro_records,1);
            end
        elseif size(BC_States.average_spatial_transpiration_hydrograph,1) < n_hydro_records
            BC_States.average_spatial_transpiration_hydrograph(n_hydro_records,1) = 0;
        end

        if flags.flag_waterquality == 1
            if ~isfield(WQ_States,'EMC_outlet') || isempty(WQ_States.EMC_outlet)
                if isfield(WQ_States,'mass_outlet') && isnumeric(WQ_States.mass_outlet)
                    WQ_States.EMC_outlet = zeros(n_hydro_records,1,'like',WQ_States.mass_outlet);
                else
                    WQ_States.EMC_outlet = zeros(n_hydro_records,1);
                end
            elseif size(WQ_States.EMC_outlet,1) < n_hydro_records
                WQ_States.EMC_outlet(n_hydro_records,1) = 0;
            end

            if ~isfield(WQ_States,'mass_outlet_save') || isempty(WQ_States.mass_outlet_save)
                if isfield(WQ_States,'mass_outlet') && isnumeric(WQ_States.mass_outlet)
                    WQ_States.mass_outlet_save = zeros(n_hydro_records,1,'like',WQ_States.mass_outlet);
                else
                    WQ_States.mass_outlet_save = zeros(n_hydro_records,1);
                end
            elseif size(WQ_States.mass_outlet_save,1) < n_hydro_records
                WQ_States.mass_outlet_save(n_hydro_records,1) = 0;
            end

            if ~isfield(WQ_States,'vol_outlet_save') || isempty(WQ_States.vol_outlet_save)
                if isfield(WQ_States,'vol_outlet') && isnumeric(WQ_States.vol_outlet)
                    WQ_States.vol_outlet_save = zeros(n_hydro_records,1,'like',WQ_States.vol_outlet);
                else
                    WQ_States.vol_outlet_save = zeros(n_hydro_records,1);
                end
            elseif size(WQ_States.vol_outlet_save,1) < n_hydro_records
                WQ_States.vol_outlet_save(n_hydro_records,1) = 0;
            end
        end

        if flags.flag_obs_gauges == 1
            n_obs_gauges_ts = gauges.num_obs_gauges;

            if ~isfield(gauges,'wse_cell') || isempty(gauges.wse_cell)
                gauges.wse_cell = zeros(n_hydro_records,n_obs_gauges_ts,'like',depths.d_t);
            elseif size(gauges.wse_cell,1) < n_hydro_records
                gauges.wse_cell(n_hydro_records,n_obs_gauges_ts) = 0;
            end

            if ~isfield(gauges,'depth_cell') || isempty(gauges.depth_cell)
                gauges.depth_cell = zeros(n_hydro_records,n_obs_gauges_ts,'like',depths.d_t);
            elseif size(gauges.depth_cell,1) < n_hydro_records
                gauges.depth_cell(n_hydro_records,n_obs_gauges_ts) = 0;
            end

            if ~isfield(gauges,'hydrograph_cell') || isempty(gauges.hydrograph_cell)
                gauges.hydrograph_cell = zeros(n_hydro_records,n_obs_gauges_ts,'like',depths.d_t);
            elseif size(gauges.hydrograph_cell,1) < n_hydro_records
                gauges.hydrograph_cell(n_hydro_records,n_obs_gauges_ts) = 0;
            end

            if flags.flag_waterquality == 1
                if ~isfield(WQ_States,'pollutograph_cell') || isempty(WQ_States.pollutograph_cell)
                    if isfield(WQ_States,'P_conc') && isnumeric(WQ_States.P_conc)
                        WQ_States.pollutograph_cell = zeros(n_hydro_records,n_obs_gauges_ts,'like',WQ_States.P_conc);
                    else
                        WQ_States.pollutograph_cell = zeros(n_hydro_records,n_obs_gauges_ts,'like',depths.d_t);
                    end
                elseif size(WQ_States.pollutograph_cell,1) < n_hydro_records
                    WQ_States.pollutograph_cell(n_hydro_records,n_obs_gauges_ts) = 0;
                end
            end
        end

        outlet_area_reporting = C_a;
        if flags.flag_subgrid == 1 && flags.flag_overbanks ~= 1
            outlet_area_reporting = ones(size(C_a), 'like', C_a) .* Wshed_Properties.cell_area;
            outlet_area_reporting(isnan(Elevation_Properties.elevation_cell)) = NaN;
        end
        outlet_states.outlet_hydrograph(1,1) = nansum(nansum(outlet_states.outlet_flow .* outlet_area_reporting)) / 1000 / 3600; % m3/s
        running_control.time_hydrograph(1,1) = running_control.time_calculation_routing(k,1) / 60;
        outlet_states.depth_outlet(1,1) = mean(depths.d_t(idx_outlet));

        % Current forcing diagnostics saved at hydrograph resolution
        current_avg_rainfall_ts = 0;
        if flags.flag_rainfall > 0
            if isfield(BC_States,'current_avg_spatial_rainfall') && ~isempty(BC_States.current_avg_spatial_rainfall)
                current_avg_rainfall_ts = BC_States.current_avg_spatial_rainfall;
            end

            flag_real_time_satellite_ts = isfield(flags,'flag_real_time_satellite_rainfall') && flags.flag_real_time_satellite_rainfall == 1;
            if (isempty(current_avg_rainfall_ts) || ~isfinite(current_avg_rainfall_ts)) && ...
                    flags.flag_spatial_rainfall ~= 1 && ...
                    flags.flag_input_rainfall_map ~= 1 && ...
                    flags.flag_satellite_rainfall ~= 1 && ...
                    flag_real_time_satellite_ts ~= 1 && ...
                    isfield(BC_States,'delta_p_agg') && ~isempty(BC_States.delta_p_agg)

                dt_for_ts = t - t_previous;
                if isfinite(dt_for_ts) && dt_for_ts > 0
                    current_avg_rainfall_ts = mean(BC_States.delta_p_agg,'all','omitnan') * 60 / dt_for_ts;
                else
                    current_avg_rainfall_ts = 0;
                end
            end
        end
        BC_States.average_spatial_rainfall_hydrograph(1,1) = current_avg_rainfall_ts;

        current_avg_evaporation_ts = 0;
        if flags.flag_ETP == 1 && isfield(BC_States,'current_avg_spatial_evaporation') && ~isempty(BC_States.current_avg_spatial_evaporation)
            current_avg_evaporation_ts = BC_States.current_avg_spatial_evaporation;
        end
        BC_States.average_spatial_evaporation_hydrograph(1,1) = current_avg_evaporation_ts;

        current_avg_transpiration_ts = 0;
        if flags.flag_ETP == 1 && isfield(BC_States,'current_avg_spatial_transpiration') && ~isempty(BC_States.current_avg_spatial_transpiration)
            current_avg_transpiration_ts = BC_States.current_avg_spatial_transpiration;
        end
        BC_States.average_spatial_transpiration_hydrograph(1,1) = current_avg_transpiration_ts;

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
        end

    % ---------------------------------------------------------------------
    % Save hydrograph/pollutograph records only when hydrograph record advances
    % ---------------------------------------------------------------------
    elseif recording_parameters.delta_record_hydrograph > 0

        t_store = recording_parameters.actual_record_hydrograph;

        outlet_area_reporting = C_a;
        if flags.flag_subgrid == 1 && flags.flag_overbanks ~= 1
            outlet_area_reporting = ones(size(C_a), 'like', C_a) .* Wshed_Properties.cell_area;
            outlet_area_reporting(isnan(Elevation_Properties.elevation_cell)) = NaN;
        end
        outlet_states.outlet_hydrograph(t_store,1) = nansum(nansum(outlet_states.outlet_flow .* outlet_area_reporting))/1000/3600; % m3/s
        running_control.time_hydrograph(t_store,1) = t;
        outlet_states.depth_outlet(t_store,1) = mean(depths.d_t(idx_outlet));

        % Current forcing diagnostics saved at hydrograph resolution
        current_avg_rainfall_ts = 0;
        if flags.flag_rainfall > 0
            if isfield(BC_States,'current_avg_spatial_rainfall') && ~isempty(BC_States.current_avg_spatial_rainfall)
                current_avg_rainfall_ts = BC_States.current_avg_spatial_rainfall;
            end

            flag_real_time_satellite_ts = isfield(flags,'flag_real_time_satellite_rainfall') && flags.flag_real_time_satellite_rainfall == 1;
            if (isempty(current_avg_rainfall_ts) || ~isfinite(current_avg_rainfall_ts)) && ...
                    flags.flag_spatial_rainfall ~= 1 && ...
                    flags.flag_input_rainfall_map ~= 1 && ...
                    flags.flag_satellite_rainfall ~= 1 && ...
                    flag_real_time_satellite_ts ~= 1 && ...
                    isfield(BC_States,'delta_p_agg') && ~isempty(BC_States.delta_p_agg)

                dt_for_ts = t - t_previous;
                if isfinite(dt_for_ts) && dt_for_ts > 0
                    current_avg_rainfall_ts = mean(BC_States.delta_p_agg,'all','omitnan') * 60 / dt_for_ts;
                else
                    current_avg_rainfall_ts = 0;
                end
            end
        end
        BC_States.average_spatial_rainfall_hydrograph(t_store,1) = current_avg_rainfall_ts;

        current_avg_evaporation_ts = 0;
        if flags.flag_ETP == 1 && isfield(BC_States,'current_avg_spatial_evaporation') && ~isempty(BC_States.current_avg_spatial_evaporation)
            current_avg_evaporation_ts = BC_States.current_avg_spatial_evaporation;
        end
        BC_States.average_spatial_evaporation_hydrograph(t_store,1) = current_avg_evaporation_ts;

        current_avg_transpiration_ts = 0;
        if flags.flag_ETP == 1 && isfield(BC_States,'current_avg_spatial_transpiration') && ~isempty(BC_States.current_avg_spatial_transpiration)
            current_avg_transpiration_ts = BC_States.current_avg_spatial_transpiration;
        end
        BC_States.average_spatial_transpiration_hydrograph(t_store,1) = current_avg_transpiration_ts;

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
        [~, velocities.hazardmax_final] = snisb_velocity_hazard_snapshot( ...
            velocities.velocity_raster, depths.d_t, idx_nan);

        if flags.flag_waterquality == 1
            WQ_States.max_Conc = zeros(size(Elevation_Properties.elevation_cell));
        end
    else
        depths.dmax_final     = max(depths.d_t, depths.dmax_final);
        velocities.vmax_final = max(velocities.velocity_raster, velocities.vmax_final);
        [~, hazard_now] = snisb_velocity_hazard_snapshot( ...
            velocities.velocity_raster, depths.d_t, idx_nan);
        velocities.hazardmax_final = max(hazard_now, velocities.hazardmax_final);

        if flags.flag_waterquality == 1
            WQ_States.max_Conc = max(WQ_States.max_Conc, WQ_States.P_conc);
        end
    end
end

outlet_area_reporting = C_a;
if flags.flag_subgrid == 1 && flags.flag_overbanks ~= 1
    outlet_area_reporting = ones(size(C_a), 'like', C_a) .* Wshed_Properties.cell_area;
    outlet_area_reporting(isnan(Elevation_Properties.elevation_cell)) = NaN;
end
outlet_runoff_volume = ...
    nansum(nansum(outlet_states.outlet_flow .* outlet_area_reporting)) * time_step / 60 / Wshed_Properties.drainage_area + outlet_runoff_volume; % mm

% =========================================================================
% Temporary lightweight time-series output snapshot
% =========================================================================
% This file is overwritten during the run. It is intentionally limited to
% time-series and scalar diagnostics so map chunks remain handled only by the
% existing save_map_hydro_*.mat workflow above.
if flags.flag_automatic_calibration ~= 1 && ...
        (k == 1 || recording_parameters.delta_record_hydrograph > 0)

    try
        if isfield(recording_parameters,'actual_record_hydrograph') && ...
                ~isempty(recording_parameters.actual_record_hydrograph)
            ts_end = recording_parameters.actual_record_hydrograph;
        else
            ts_end = 1;
        end
        ts_end = max(1, ts_end);

        TempTS = struct();
        TempTS.meta.created_datetime = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
        TempTS.meta.k = k;
        TempTS.meta.t_current_min = t;
        TempTS.meta.t_save_min = t_save;
        TempTS.meta.actual_record_hydrograph = ts_end;
        TempTS.meta.record_time_hydrographs_min = running_control.record_time_hydrographs;
        TempTS.meta.description = ['Temporary HydroPol2D time-series snapshot. ', ...
            'Raster map stacks are not stored here.'];

        % Time vector
        if isfield(running_control,'time_hydrograph') && ~isempty(running_control.time_hydrograph)
            n_tmp = min(ts_end, size(running_control.time_hydrograph,1));
            tmp_data = running_control.time_hydrograph(1:n_tmp,:);
            if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
            TempTS.time.time_hydrograph_min = tmp_data;
        end

        % Outlet time series
        if exist('outlet_states','var') && isstruct(outlet_states)
            if isfield(outlet_states,'outlet_hydrograph') && ~isempty(outlet_states.outlet_hydrograph)
                n_tmp = min(ts_end, size(outlet_states.outlet_hydrograph,1));
                tmp_data = outlet_states.outlet_hydrograph(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.outlet.hydrograph_m3_s = tmp_data;
            end

            if isfield(outlet_states,'depth_outlet') && ~isempty(outlet_states.depth_outlet)
                n_tmp = min(ts_end, size(outlet_states.depth_outlet,1));
                tmp_data = outlet_states.depth_outlet(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.outlet.depth_outlet_mm = tmp_data;
            end
        end

        % Catchment-average forcing diagnostics saved at hydrograph resolution
        if exist('BC_States','var') && isstruct(BC_States)
            if isfield(BC_States,'average_spatial_rainfall_hydrograph') && ~isempty(BC_States.average_spatial_rainfall_hydrograph)
                n_tmp = min(ts_end, size(BC_States.average_spatial_rainfall_hydrograph,1));
                tmp_data = BC_States.average_spatial_rainfall_hydrograph(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.forcing.catchment_average_rainfall_mm_h = tmp_data;
            end

            if isfield(BC_States,'average_spatial_evaporation_hydrograph') && ~isempty(BC_States.average_spatial_evaporation_hydrograph)
                n_tmp = min(ts_end, size(BC_States.average_spatial_evaporation_hydrograph,1));
                tmp_data = BC_States.average_spatial_evaporation_hydrograph(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.forcing.catchment_average_evaporation_mm_day = tmp_data;
            end

            if isfield(BC_States,'average_spatial_transpiration_hydrograph') && ~isempty(BC_States.average_spatial_transpiration_hydrograph)
                n_tmp = min(ts_end, size(BC_States.average_spatial_transpiration_hydrograph,1));
                tmp_data = BC_States.average_spatial_transpiration_hydrograph(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.forcing.catchment_average_transpiration_mm_day = tmp_data;
            end
        end

        % Water-quality time series
        if flags.flag_waterquality == 1 && exist('WQ_States','var') && isstruct(WQ_States)
            if isfield(WQ_States,'EMC_outlet') && ~isempty(WQ_States.EMC_outlet)
                n_tmp = min(ts_end, size(WQ_States.EMC_outlet,1));
                tmp_data = WQ_States.EMC_outlet(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.waterquality.EMC_outlet_mg_L = tmp_data;
            end

            if isfield(WQ_States,'mass_outlet_save') && ~isempty(WQ_States.mass_outlet_save)
                n_tmp = min(ts_end, size(WQ_States.mass_outlet_save,1));
                tmp_data = WQ_States.mass_outlet_save(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.waterquality.mass_outlet = tmp_data;
            end

            if isfield(WQ_States,'vol_outlet_save') && ~isempty(WQ_States.vol_outlet_save)
                n_tmp = min(ts_end, size(WQ_States.vol_outlet_save,1));
                tmp_data = WQ_States.vol_outlet_save(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.waterquality.volume_outlet = tmp_data;
            end

            if isfield(WQ_States,'pollutograph_cell') && ~isempty(WQ_States.pollutograph_cell)
                n_tmp = min(ts_end, size(WQ_States.pollutograph_cell,1));
                tmp_data = WQ_States.pollutograph_cell(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.waterquality.pollutograph_cell_mg_L = tmp_data;
            end
        end

        % Observation gauge time series
        if flags.flag_obs_gauges == 1 && exist('gauges','var') && isstruct(gauges)
            if isfield(gauges,'hydrograph_cell') && ~isempty(gauges.hydrograph_cell)
                n_tmp = min(ts_end, size(gauges.hydrograph_cell,1));
                tmp_data = gauges.hydrograph_cell(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.gauges.hydrograph_cell = tmp_data;
            end

            if isfield(gauges,'depth_cell') && ~isempty(gauges.depth_cell)
                n_tmp = min(ts_end, size(gauges.depth_cell,1));
                tmp_data = gauges.depth_cell(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.gauges.depth_cell_m = tmp_data;
            end

            if isfield(gauges,'wse_cell') && ~isempty(gauges.wse_cell)
                n_tmp = min(ts_end, size(gauges.wse_cell,1));
                tmp_data = gauges.wse_cell(1:n_tmp,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.gauges.wse_cell_m = tmp_data;
            end

            if isfield(gauges,'labels_observed_string') && ~isempty(gauges.labels_observed_string)
                TempTS.gauges.labels = gauges.labels_observed_string;
            end
        end

        % Scalar diagnostics accumulated during the run
        tmp_data = Flooded_Area;
        if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
        TempTS.diagnostics.Flooded_Area_m2 = tmp_data;

        if exist('Risk_Area','var')
            tmp_data = Risk_Area;
            if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
            TempTS.diagnostics.Risk_Area_m2 = tmp_data;
        end

        if exist('outlet_runoff_volume','var')
            tmp_data = outlet_runoff_volume;
            if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
            TempTS.diagnostics.outlet_runoff_volume_mm = tmp_data;
        end

        % Mass-balance report history, if available
        if exist('mass_balance_history','var') && isstruct(mass_balance_history) && ...
                isfield(mass_balance_history,'count') && mass_balance_history.count > 0
            ih = mass_balance_history.count;

            if isfield(mass_balance_history,'time_h')
                tmp_data = mass_balance_history.time_h(1:ih,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.mass_balance.time_h = tmp_data;
            end
            if isfield(mass_balance_history,'errors_m3')
                tmp_data = mass_balance_history.errors_m3(1:ih,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.mass_balance.errors_m3 = tmp_data;
            end
            if isfield(mass_balance_history,'errors_mm_h')
                tmp_data = mass_balance_history.errors_mm_h(1:ih,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.mass_balance.errors_mm_h = tmp_data;
            end
            if isfield(mass_balance_history,'cum_errors_m3')
                tmp_data = mass_balance_history.cum_errors_m3(1:ih,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.mass_balance.cum_errors_m3 = tmp_data;
            end
            if isfield(mass_balance_history,'cum_errors_mm')
                tmp_data = mass_balance_history.cum_errors_mm(1:ih,:);
                if isa(tmp_data,'gpuArray'); tmp_data = gather(tmp_data); end
                TempTS.mass_balance.cum_errors_mm = tmp_data;
            end
        end

        % Resolve temporary output folder without requiring any new input.
        if exist('Paths','var') && isstruct(Paths) && isfield(Paths,'Temp') && ~isempty(Paths.Temp)
            tempDir = Paths.Temp;
        elseif exist('folderName_2','var') && ~isempty(folderName_2)
            tempDir = folderName_2;
        elseif exist('temp_dir','var') && ~isempty(temp_dir)
            tempDir = temp_dir;
        else
            tempDir = pwd;
        end

        if ~exist(tempDir,'dir')
            mkdir(tempDir);
        end

        tempFile  = fullfile(tempDir, 'current_time_series_outputs.tmp.mat');
        finalFile = fullfile(tempDir, 'current_time_series_outputs.mat');

        save(tempFile, 'TempTS', '-v7.3');
        movefile(tempFile, finalFile, 'f');

    catch temp_timeseries_save_error
        warning('HydroPol2D:TempTimeSeriesSaveFailed', ...
            'Could not save temporary time-series outputs: %s', ...
            temp_timeseries_save_error.message);
    end
end

function [velocity_map, hazard_map] = snisb_velocity_hazard_snapshot(velocity_raster, depth_mm, idx_nan)
%SNISB_VELOCITY_HAZARD_SNAPSHOT Build velocity and depth*velocity maps.
% Depth is stored internally in millimetres; exported hazard uses metres.
% Cells below 0.01 m are assigned zero hazard to keep exposure semantics
% consistent across the SNISB dam-break workflow.
    depth_m = max(double(depth_mm) / 1000, 0);
    velocity_map = double(velocity_raster);
    if isscalar(velocity_map)
        velocity_map = zeros(size(depth_m)) + velocity_map;
    end
    velocity_map(~isfinite(velocity_map)) = 0;
    hazard_map = depth_m .* velocity_map;
    hazard_map(depth_m < 0.01) = 0;
    hazard_map(~isfinite(hazard_map)) = 0;
    velocity_map(idx_nan) = NaN;
    hazard_map(idx_nan) = NaN;
end
