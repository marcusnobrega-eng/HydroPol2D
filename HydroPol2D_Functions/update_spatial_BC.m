
%% ═══════════════════════════════════════════════════════════════════════
%  Function: update_spatial_BC
%
%  ---------------------------------------------------------------------
%  This version assumes:
%     - stage hydrograph is piecewise constant
%     - inflow hydrograph is piecewise constant
%     - rainfall forcings are piecewise constant
%
%  Interpretation of forcing times
%  ---------------------------------------------------------------------
%  For all time series used here, time(i) is interpreted as the END of the
%  forcing interval i. Therefore, the value at index i is assumed valid
%  from the previous forcing time up to time(i).
%
%  Optimization strategy used here
%  ---------------------------------------------------------------------
%  Instead of recomputing overlap weights every model step, we use a moving
%  cursor that only advances forward in time. This is much cheaper because:
%
%    - model time t increases monotonically
%    - forcing timestamps are static
%    - only the currently active forcing interval is used
%
%  Important unit convention for forcing variables
%  ---------------------------------------------------------------------
%  Lumped rainfall:
%     Rainfall_Parameters.intensity_rainfall is a native piecewise-constant
%     rate series [mm/h] defined on Rainfall_Parameters.time_rainfall.
%
%  Inflow hydrograph:
%     Inflow_Parameters.inflow_hydrograph_rate is a native piecewise-
%     constant rate series [mm/h] defined on Inflow_Parameters.time_inflow.
%
%  Stage hydrograph:
%     Stage_Parameters.stage is a native piecewise-constant depth series
%     [m] defined on Stage_Parameters.time_stage.
%
%  Spatial / map / satellite rainfall:
%     These are piecewise-constant raster rate fields [mm/h].
%
%  For rate forcings, update_spatial_BC computes the exact integrated depth
%  over [t_previous, t] using overlap durations.
%
%  For stage, update_spatial_BC computes the exact time-averaged stage over
%  [t_previous, t].
%% ═══════════════════════════════════════════════════════════════════════

% -------------------------------------------------------------------------
% Current model-step duration [minutes]
% -------------------------------------------------------------------------
dt_model_current    = max(t - t_previous, 0);
dt_model_current_hr = dt_model_current / 60; %#ok<NASGU>

% -------------------------------------------------------------------------
% One-time initialization
% -------------------------------------------------------------------------
if k == 1
    % ---------------------------------------------------------------------
    % Canonical grid size
    % ---------------------------------------------------------------------
    gridSize = size(Elevation_Properties.elevation_cell);
    BC_States.gridSize = gridSize;

    % ---------------------------------------------------------------------
    % Recover static grid size from memory
    % ---------------------------------------------------------------------
    gridSize = BC_States.gridSize;

    % ---------------------------------------------------------------------
    % Initialize rainfall holder
    % For spatial/map rainfall this is a raster; for lumped rainfall it can
    % later become scalar depending on the branch used.
    % ---------------------------------------------------------------------
    % BC_States.delta_p_agg = zeros(gridSize);
    % BC_States.delta_p_agg(idx_nan) = nan;

    % ---------------------------------------------------------------------
    % Ensure Hydro_States fields exist
    % ---------------------------------------------------------------------
    if ~isfield(Hydro_States,'ETP') || isempty(Hydro_States.ETP)
        Hydro_States.ETP = zeros(gridSize);
        Hydro_States.ETP(idx_nan) = nan;
    end
    if ~isfield(Hydro_States,'ETR') || isempty(Hydro_States.ETR)
        Hydro_States.ETR = zeros(gridSize);
        Hydro_States.ETR(idx_nan) = nan;
    end
    if ~isfield(Hydro_States,'Ep') || isempty(Hydro_States.Ep)
        Hydro_States.Ep = zeros(gridSize);
        Hydro_States.Ep(idx_nan) = nan;
    end

    % ---------------------------------------------------------------------
    % Keep most recent valid internal ETP/ETR in memory
    % ---------------------------------------------------------------------
    if ~isfield(Hydro_States,'last_valid_ETP') || isempty(Hydro_States.last_valid_ETP)
        Hydro_States.last_valid_ETP = Hydro_States.ETP;
    end
    if ~isfield(Hydro_States,'last_valid_ETR') || isempty(Hydro_States.last_valid_ETR)
        Hydro_States.last_valid_ETR = Hydro_States.ETR;
    end

    % ---------------------------------------------------------------------
    % Initialize piecewise-constant cursors
    % ---------------------------------------------------------------------

    % Lumped rainfall cursor
    if ~isfield(BC_States,'cursor_lumped_rain') || isempty(BC_States.cursor_lumped_rain)
        BC_States.cursor_lumped_rain = 1;
    end

    % Inflow cursor
    if ~isfield(BC_States,'cursor_inflow') || isempty(BC_States.cursor_inflow)
        BC_States.cursor_inflow = 1;
    end

    % Stage cursor
    if ~isfield(BC_States,'cursor_stage') || isempty(BC_States.cursor_stage)
        BC_States.cursor_stage = 1;
    end

    % Spatial rainfall cursor
    if flags.flag_spatial_rainfall == 1 && ...
            flags.flag_satellite_rainfall == 0 && ...
            flags.flag_input_rainfall_map == 0

        if ~isfield(Spatial_Rainfall_Parameters,'cursor_rain') || isempty(Spatial_Rainfall_Parameters.cursor_rain)
            Spatial_Rainfall_Parameters.cursor_rain = 1;
        end
    end

    % Input rainfall map cursor
    if isfield(flags,'flag_input_rainfall_map') && flags.flag_input_rainfall_map == 1
        if ~isfield(Input_Rainfall,'cursor_rain') || isempty(Input_Rainfall.cursor_rain)
            Input_Rainfall.cursor_rain = 1;
        end
    end

    % Satellite rainfall cursor
    if flags.flag_satellite_rainfall == 1 || flags.flag_real_time_satellite_rainfall == 1
        if ~isfield(Input_Rainfall,'cursor_rain') || isempty(Input_Rainfall.cursor_rain)
            Input_Rainfall.cursor_rain = 1;
        end
    end

    % ---------------------------------------------------------------------
    % Initialize cache fields for rainfall products
    % ---------------------------------------------------------------------
    if flags.flag_spatial_rainfall == 1
        if ~isfield(Spatial_Rainfall_Parameters,'cache_index_rain')
            Spatial_Rainfall_Parameters.cache_index_rain = [];
        end
        if ~isfield(Spatial_Rainfall_Parameters,'cache_map_rain')
            Spatial_Rainfall_Parameters.cache_map_rain = [];
        end
    end

    if isfield(flags,'flag_input_rainfall_map') && flags.flag_input_rainfall_map == 1
        if ~isfield(Input_Rainfall,'cache_index_rain')
            Input_Rainfall.cache_index_rain = [];
        end
        if ~isfield(Input_Rainfall,'cache_map_rain')
            Input_Rainfall.cache_map_rain = [];
        end
    end

    % ---------------------------------------------------------------------
    % Ensure map-based ETP forcing fields exist
    % ---------------------------------------------------------------------
    if ~isfield(BC_States,'delta_e_agg') || isempty(BC_States.delta_e_agg)
        BC_States.delta_e_agg = zeros(gridSize);
        BC_States.delta_e_agg(idx_nan) = nan;
    end
    if ~isfield(BC_States,'delta_tr_agg') || isempty(BC_States.delta_tr_agg)
        BC_States.delta_tr_agg = zeros(gridSize);
        BC_States.delta_tr_agg(idx_nan) = nan;
    end

    if flags.flag_input_ETP_map == 1
        if ~isfield(Input_Evaporation,'cache_index_et')
            Input_Evaporation.cache_index_et = [];
        end
        if ~isfield(Input_Evaporation,'cache_map_et')
            Input_Evaporation.cache_map_et = [];
        end
        if ~isfield(Input_Transpiration,'cache_index_et')
            Input_Transpiration.cache_index_et = [];
        end
        if ~isfield(Input_Transpiration,'cache_map_et')
            Input_Transpiration.cache_map_et = [];
        end

        if ~isfield(Input_Evaporation,'cursor_et') || isempty(Input_Evaporation.cursor_et)
            Input_Evaporation.cursor_et = 1;
        end
        if ~isfield(Input_Transpiration,'cursor_et') || isempty(Input_Transpiration.cursor_et)
            Input_Transpiration.cursor_et = 1;
        end
    end

    % ---------------------------------------------------------------------
    % Ensure rainfall / inflow fields exist
    % ---------------------------------------------------------------------
    if ~isfield(BC_States,'delta_p_agg') || isempty(BC_States.delta_p_agg)
        BC_States.delta_p_agg = 0;
    end
    if ~isfield(BC_States,'inflow') || isempty(BC_States.inflow)
        BC_States.inflow = zeros(gridSize);
        BC_States.inflow(idx_nan) = nan;
    end
    if ~isfield(BC_States,'delta_inflow_agg') || isempty(BC_States.delta_inflow_agg)
        BC_States.delta_inflow_agg = zeros(Inflow_Parameters.n_stream_gauges,1);
    end

    % ---------------------------------------------------------------------
    % Current map holders
    % ---------------------------------------------------------------------
    if ~isfield(BC_States,'current_rainfall_mean_map') || isempty(BC_States.current_rainfall_mean_map)
        BC_States.current_rainfall_mean_map = zeros(gridSize);
        BC_States.current_rainfall_mean_map(idx_nan) = nan;
    end
    if ~isfield(BC_States,'current_avg_spatial_rainfall') || isempty(BC_States.current_avg_spatial_rainfall)
        BC_States.current_avg_spatial_rainfall = 0;
    end

    if ~isfield(BC_States,'current_evaporation_map') || isempty(BC_States.current_evaporation_map)
        BC_States.current_evaporation_map = zeros(gridSize);
        BC_States.current_evaporation_map(idx_nan) = nan;
    end
    if ~isfield(BC_States,'current_transpiration_map') || isempty(BC_States.current_transpiration_map)
        BC_States.current_transpiration_map = zeros(gridSize);
        BC_States.current_transpiration_map(idx_nan) = nan;
    end
    if ~isfield(BC_States,'current_avg_spatial_evaporation') || isempty(BC_States.current_avg_spatial_evaporation)
        BC_States.current_avg_spatial_evaporation = 0;
    end
    if ~isfield(BC_States,'current_avg_spatial_transpiration') || isempty(BC_States.current_avg_spatial_transpiration)
        BC_States.current_avg_spatial_transpiration = 0;
    end

    % ---------------------------------------------------------------------
    % Build rainfall interpolation support only once
    % ---------------------------------------------------------------------
    if flags.flag_spatial_rainfall == 1 && ...
            flags.flag_satellite_rainfall == 0 && ...
            flags.flag_input_rainfall_map == 0

        if ~isfield(Spatial_Rainfall_Parameters,'x_coordinate') || isempty(Spatial_Rainfall_Parameters.x_coordinate)
            Spatial_Rainfall_Parameters.x_coordinate = ...
                Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1);
        end

        if ~isfield(Spatial_Rainfall_Parameters,'y_coordinate') || isempty(Spatial_Rainfall_Parameters.y_coordinate)
            Spatial_Rainfall_Parameters.y_coordinate = ...
                Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2);
        end

        if ~isfield(Spatial_Rainfall_Parameters,'x_grid') || isempty(Spatial_Rainfall_Parameters.x_grid)
            Spatial_Rainfall_Parameters.x_grid = ...
                GIS_data.xulcorner + Wshed_Properties.Resolution * (1:size(DEM_raster.Z,2))';
        end

        if ~isfield(Spatial_Rainfall_Parameters,'y_grid') || isempty(Spatial_Rainfall_Parameters.y_grid)
            Spatial_Rainfall_Parameters.y_grid = ...
                GIS_data.yulcorner - Wshed_Properties.Resolution * (1:size(DEM_raster.Z,1))';
        end
    end
end

% -------------------------------------------------------------------------
% Recover static grid size after initialization
% -------------------------------------------------------------------------
gridSize = BC_States.gridSize;

%% ═══════════════════════════════════════════════════════════════════════
%  Stage-Hydrograph Boundary Condition (exact time-average over [t_previous,t])
% ═══════════════════════════════════════════════════════════════════════
% if flags.flag_stage_hydrograph == 1
%     if k == 1
%         stage_depth = zeros(Stage_Parameters.n_stage_gauges,1);
%     else
%         stage_depth_previous = stage_depth;
%     end
%     dt_step = max(t - t_previous, 0);
% 
%     if dt_step > 0
%         for z = 1:Stage_Parameters.n_stage_gauges
%             [stage_depth(z,1), BC_States.cursor_stage] = localAverageStateOverStep( ...
%                 Stage_Parameters.time_stage, ...
%                 Stage_Parameters.stage(:,z), ...
%                 t_previous, ...
%                 t, ...
%                 BC_States.cursor_stage);
%         end
%     else
%         % fallback for zero-duration step
%         BC_States.cursor_stage = localAdvanceEndStampedCursor( ...
%             Stage_Parameters.time_stage, ...
%             t, ...
%             BC_States.cursor_stage);
% 
%         iz_stage = localClampIndex(BC_States.cursor_stage, 1, size(Stage_Parameters.stage,1));
% 
%         for z = 1:Stage_Parameters.n_stage_gauges
%             stage_depth(z,1) = Stage_Parameters.stage(iz_stage,z);
%         end
%     end
% end
% 
% if k == 1
%    stage_depth_previous = stage_depth;
% end

% Stage-Hydrograph Boundary Condition - benchmark mode
if flags.flag_stage_hydrograph == 1

    stage_depth = zeros(Stage_Parameters.n_stage_gauges,1);

    for z = 1:Stage_Parameters.n_stage_gauges

        z2 = find(Stage_Parameters.time_stage <= t,1,'last');

        if isempty(z2)
            z2 = 1;
        end

        stage_depth(z,1) = Stage_Parameters.stage(z2,z);

    end

    % Apply Stage Boundary Condition
    for i = 1:Stage_Parameters.n_stage_gauges

        stage_cells = Wshed_Properties.stage_mask(:,:,i);

        % Stage is assumed to be water depth in meters.
        % Convert to mm for model state.
        depths.d_t(logical(stage_cells)) = 1000 * stage_depth(i,1);

    end

    stage_depth_previous = stage_depth;

end
%% Apply Stage Boundary Condition
if flags.flag_stage_hydrograph == 1
    for i = 1:Stage_Parameters.n_stage_gauges
        stage_cells = Wshed_Properties.stage_mask(:,:,i);

        % Stage is assumed to be water depth in meters.
        % Convert to mm for model state.
        depths.d_t(logical(stage_cells)) = ...
            1000 * stage_depth(i,1) * stage_cells(stage_cells == 1);
    end
end

%% ═══════════════════════════════════════════════════════════════════════
%  Inflows for current model step (exact integral over [t_previous,t])
% ═══════════════════════════════════════════════════════════════════════
if flags.flag_inflow > 0

    BC_States.delta_inflow_agg = zeros(Inflow_Parameters.n_stream_gauges,1);
    dt_step = max(t - t_previous, 0);

    if dt_step > 0
        for z = 1:Inflow_Parameters.n_stream_gauges
            [avg_inflow_rate_mm_h, BC_States.cursor_inflow] = localAverageRateOverStep( ...
                Inflow_Parameters.time_inflow, ...
                Inflow_Parameters.inflow_hydrograph_rate(:,z), ...
                t_previous, ...
                t, ...
                BC_States.cursor_inflow);

            BC_States.delta_inflow_agg(z,1) = avg_inflow_rate_mm_h * dt_step / 60; % [mm]
        end
    end
end

%% Inflows Boundary Condition
if flags.flag_inflow == 1
    BC_States.inflow = zeros(gridSize);

    for i = 1:Inflow_Parameters.n_stream_gauges
        BC_States.inflow = BC_States.inflow + ...
            BC_States.delta_inflow_agg(i) * Wshed_Properties.inflow_cells(:,:,i); % mm
    end

    if flags.flag_subgrid == 1
        BC_States.inflow = BC_States.inflow .* Wshed_Properties.cell_area ./ C_a;
    end

    BC_States.inflow(idx_nan) = nan;
end

%% ═══════════════════════════════════════════════════════════════════════
%  Aggregating Precipitation to the New Time-step
%  All rainfall is treated as piecewise constant
% ═══════════════════════════════════════════════════════════════════════
if flags.flag_rainfall > 0

    % Reset rainfall forcing for the current step
    % Only reset scalar/lumped rainfall. Map branches overwrite the whole array.
    if flags.flag_spatial_rainfall ~= 1 && ...
            flags.flag_input_rainfall_map ~= 1 && ...
            flags.flag_satellite_rainfall ~= 1 && ...
            flags.flag_real_time_satellite_rainfall ~= 1
    
        BC_States.delta_p_agg = 0;
    end

    % ---------------------------------------------------------------------
    % Lumped rainfall (piecewise constant depth already on model-step grid)
    % BC_States.delta_p is produced in preprocessing with
    % accumulated_incremental(...), so it is already a depth [mm] for the
    % active model-step interval.
    % ---------------------------------------------------------------------
if flags.flag_spatial_rainfall ~= 1 && ...
        flags.flag_input_rainfall_map ~= 1 && ...
        flags.flag_satellite_rainfall ~= 1 && ...
        flags.flag_real_time_satellite_rainfall ~= 1

    dt_step = max(t - t_previous, 0);

    % Default value for storage/diagnostics.
    % This prevents undefined avg_rain_rate_mm_h when dt_step == 0.
    avg_rain_rate_mm_h = 0;

    if dt_step > 0

        [avg_rain_rate_mm_h, BC_States.cursor_lumped_rain] = localAverageRateOverStep( ...
            Rainfall_Parameters.time_rainfall, ...
            Rainfall_Parameters.intensity_rainfall(:), ...
            t_previous, ...
            t, ...
            BC_States.cursor_lumped_rain);

        if isnan(avg_rain_rate_mm_h)
            warning('No rainfall data for this time. Assuming 0 mm/h.')
            avg_rain_rate_mm_h = 0;
        end

        BC_States.delta_p_agg = avg_rain_rate_mm_h * dt_step / 60; % [mm]

    else

        BC_States.delta_p_agg = 0;

    end

    % -----------------------------------------------------------------
    % Current rainfall diagnostic for temporary time-series output.
    % For lumped rainfall, the catchment-average rainfall rate is simply
    % the lumped rainfall rate [mm/h].
    % -----------------------------------------------------------------
    BC_States.current_avg_spatial_rainfall = avg_rain_rate_mm_h;

    if isfield(BC_States,'current_rainfall_mean_map') && ...
            ~isempty(BC_States.current_rainfall_mean_map)

        BC_States.current_rainfall_mean_map(:) = avg_rain_rate_mm_h;
        BC_States.current_rainfall_mean_map(idx_nan) = nan;

    end
        % ---------------------------------------------------------------------
        % Spatial rainfall from gauges (piecewise constant rate field in mm/h)
        % ---------------------------------------------------------------------
    elseif flags.flag_spatial_rainfall == 1 && ...
            flags.flag_input_rainfall_map ~= 1 && ...
            flags.flag_satellite_rainfall ~= 1 && ...
            flags.flag_real_time_satellite_rainfall ~= 1

        [avg_rain_map, Spatial_Rainfall_Parameters.cursor_rain, Spatial_Rainfall_Parameters] = ...
            localAverageSpatialGaugeRainOverStep( ...
            Spatial_Rainfall_Parameters.rainfall_spatial_duration, ...
            Spatial_Rainfall_Parameters.rainfall_raingauges, ...
            t_previous, ...
            t, ...
            Spatial_Rainfall_Parameters.cursor_rain, ...
            Spatial_Rainfall_Parameters, ...
            gridSize, ...
            idx_nan, ...
            @Rainfall_Interpolator);

        BC_States.delta_p_agg = avg_rain_map * (dt_model_current / 60);

        BC_States.current_rainfall_mean_map = localGather(avg_rain_map);
        BC_States.current_avg_spatial_rainfall = mean(avg_rain_map,'all','omitnan');

%         rainfall_spatial_aggregation(:,:,1) = BC_States.current_rainfall_mean_map; %#ok<NASGU>
%         Rainfall_Parameters.index_aggregation = 1;

        % ---------------------------------------------------------------------
        % Input rainfall maps (piecewise constant rate field in mm/h)
        % ---------------------------------------------------------------------
    elseif flags.flag_input_rainfall_map == 1 && ...
            flags.flag_satellite_rainfall ~= 1 && ...
            flags.flag_real_time_satellite_rainfall ~= 1

        [avg_rain_map, Input_Rainfall.cursor_rain, Input_Rainfall] = ...
            localAverageInputRainMapOverStep( ...
            Input_Rainfall.time, ...
            t_previous, ...
            t, ...
            Input_Rainfall.cursor_rain, ...
            Input_Rainfall, ...
            DEM_raster, ...
            GIS_data, ...
            flags, ...
            idx_nan, ...
            k);

        BC_States.delta_p_agg = avg_rain_map * (dt_model_current / 60);

        BC_States.current_rainfall_mean_map = localGather(avg_rain_map);
        BC_States.current_avg_spatial_rainfall = mean(avg_rain_map,'all','omitnan');

%         rainfall_spatial_aggregation(:,:,1) = BC_States.current_rainfall_mean_map; %#ok<NASGU>
%         Rainfall_Parameters.index_aggregation = 1;

        % ---------------------------------------------------------------------
        % Satellite rainfall (piecewise constant rate field in mm/h)
        % ---------------------------------------------------------------------
    elseif flags.flag_input_rainfall_map == 0 && ...
            flags.flag_satellite_rainfall == 1 && ...
            flags.flag_real_time_satellite_rainfall ~= 1

        recording_parameters.actual_record_state_rainfall = ...
            find(running_control.time_records < t,1,'last');

        recording_parameters.delta_record_rainfall = ...
            recording_parameters.actual_record_state_rainfall - ...
            recording_parameters.last_record_maps_rainfall;

        recording_parameters.last_record_maps_rainfall = ...
            recording_parameters.actual_record_state_rainfall;

        Input_Rainfall.cursor_rain = localAdvanceEndStampedCursor( ...
            Input_Rainfall.time, ...
            t, ...
            Input_Rainfall.cursor_rain);

        % Keep satellite call structure unchanged except for removing overlap logic
        product = 'PDIRNow1hourly';

        [rainfall_raster, register, ~] = Satellite_rainfall_processing( ...
            [], [], register, product, date_begin, date_end, ...
            flags.flag_satellite_rainfall, ...
            flags.flag_real_time_satellite_rainfall, ...
            DEM_raster);

        input_rainfall = rainfall_raster.Z;
        input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;
        input_rainfall(idx_nan) = nan;

        % Convert rate field [mm/h] to rainfall depth over current step [mm]
        BC_States.delta_p_agg = input_rainfall * (dt_model_current / 60);

        % Store current rate map for diagnostics
        BC_States.current_rainfall_mean_map = localGather(input_rainfall);
        BC_States.current_avg_spatial_rainfall = mean(input_rainfall,'all','omitnan');

%         rainfall_spatial_aggregation(:,:,1) = BC_States.current_rainfall_mean_map; %#ok<NASGU>
%         Rainfall_Parameters.index_aggregation = 1;

        % ---------------------------------------------------------------------
        % Real-time satellite rainfall
        % Kept as placeholder exactly as in your current structure
        % ---------------------------------------------------------------------
    elseif flags.flag_input_rainfall_map == 0 && ...
            flags.flag_satellite_rainfall ~= 1 && ...
            flags.flag_real_time_satellite_rainfall == 1

        if flags.flag_real_time_satellite_rainfall == 1
            % Placeholder kept exactly as in your original code
        end
    end
end

% -------------------------------------------------------------------------
% Optional sub-grid rainfall correction kept commented exactly as before
% -------------------------------------------------------------------------
% BC_States.delta_p_agg = BC_States.delta_p_agg .* Wshed_Properties.cell_area ./ C_a;
% BC_States.delta_p_agg(idx_nan) = nan;

%% ═══════════════════════════════════════════════════════════════════════
%  Aggregating ETP for next time-step (internal meteorological mode)
% ═══════════════════════════════════════════════════════════════════════
if flags.flag_ETP == 1 && flags.flag_input_ETP_map ~= 1

    z1 = find(ETP_Parameters.climatologic_spatial_duration <= t_previous,1,'last');
    z2 = find(ETP_Parameters.climatologic_spatial_duration <= t,1,'last');

    if isempty(z1) && isempty(z2)

        Hydro_States.ETP = zeros(gridSize); Hydro_States.ETP(idx_nan) = nan;
        Hydro_States.ETR = zeros(gridSize); Hydro_States.ETR(idx_nan) = nan;
        Hydro_States.Ep  = zeros(gridSize); Hydro_States.Ep(idx_nan)  = nan;

    else

        if ~isempty(z1) && z1 == z2 && z2 == length(ETP_Parameters.climatologic_spatial_duration)

            Hydro_States.ETP = zeros(gridSize); Hydro_States.ETP(idx_nan) = nan;

        elseif (isempty(z1) && z2 > 0) || (z2 > z1 && z2 < length(ETP_Parameters.climatologic_spatial_duration))

            if ~isempty(ETP_Parameters.time_ETP)
                day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
            else
                day_of_year = day(extra_parameters_ETP.time_ETP(z2,1),'dayofyear');
            end

            [Hydro_States.ETP, Hydro_States.Ep, ~, ~, ~] = ...
                ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2), ...
                Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations, ...
                ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations, ...
                ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);

            Hydro_States.ETP(isnan(Hydro_States.ETP)) = 0; Hydro_States.ETP(idx_nan) = nan;
            Hydro_States.Ep(isnan(Hydro_States.Ep))   = 0; Hydro_States.Ep(idx_nan)  = nan;

            if nansum(Hydro_States.ETP,'all') == 0
                warning('No ETP data from current forcing. Using last valid internal ETP/ETR in memory.')
                Hydro_States.ETP = Hydro_States.last_valid_ETP;
                Hydro_States.ETR = Hydro_States.last_valid_ETR;

                if nansum(Hydro_States.ETP,'all') == 0
                    warning('No valid internal ETP/ETR available. Assuming 0.')
                    Hydro_States.ETP = zeros(gridSize); Hydro_States.ETP(idx_nan) = nan;
                    Hydro_States.ETR = zeros(gridSize); Hydro_States.ETR(idx_nan) = nan;
                end
            else
                Hydro_States.last_valid_ETP = Hydro_States.ETP;
                Hydro_States.last_valid_ETR = Hydro_States.ETR;
            end

        elseif z1 == 1 && z2 == 1 && k == 1

            if ~isempty(ETP_Parameters.time_ETP)
                day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
            else
                day_of_year = day(extra_parameters_ETP.time_ETP(z2,1),'dayofyear');
            end

            [Hydro_States.ETP, Hydro_States.Ep, ~, ~, ~] = ...
                ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2), ...
                Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations, ...
                ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations, ...
                ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);

            Hydro_States.ETP(isnan(Hydro_States.ETP)) = 0; Hydro_States.ETP(idx_nan) = nan;
            Hydro_States.Ep(isnan(Hydro_States.Ep))   = 0; Hydro_States.Ep(idx_nan)  = nan;

            if nansum(Hydro_States.ETP,'all') == 0
                Hydro_States.ETP = Hydro_States.last_valid_ETP;
                Hydro_States.ETR = Hydro_States.last_valid_ETR;
            else
                Hydro_States.last_valid_ETP = Hydro_States.ETP;
                Hydro_States.last_valid_ETR = Hydro_States.ETR;
            end
        end
    end
end

%% ═══════════════════════════════════════════════════════════════════════
%  Aggregating ETP for next time-step (input evaporation/transpiration maps)
%  Exact overlap-weighted average over [t_previous,t]
% ═══════════════════════════════════════════════════════════════════════
if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1

    dt_step = max(t - t_previous, 0);

    if ~isfield(Input_Evaporation,'cursor_et') || isempty(Input_Evaporation.cursor_et)
        Input_Evaporation.cursor_et = 1;
    end
    if ~isfield(Input_Transpiration,'cursor_et') || isempty(Input_Transpiration.cursor_et)
        Input_Transpiration.cursor_et = 1;
    end

    if dt_step > 0

        [avg_evaporation_map, Input_Evaporation.cursor_et, Input_Evaporation] = ...
            localAverageInputETMapOverStep( ...
            Input_Evaporation.time, ...
            Input_Evaporation.labels_Directory, ...
            t_previous, ...
            t, ...
            Input_Evaporation.cursor_et, ...
            Input_Evaporation, ...
            DEM_raster, ...
            GIS_data, ...
            flags, ...
            idx_nan, ...
            k);

        [avg_transpiration_map, Input_Transpiration.cursor_et, Input_Transpiration] = ...
            localAverageInputETMapOverStep( ...
            Input_Transpiration.time, ...
            Input_Transpiration.labels_Directory, ...
            t_previous, ...
            t, ...
            Input_Transpiration.cursor_et, ...
            Input_Transpiration, ...
            DEM_raster, ...
            GIS_data, ...
            flags, ...
            idx_nan, ...
            k);

        avg_evaporation_map(idx_nan)   = nan;
        avg_transpiration_map(idx_nan) = nan;

        avg_evaporation_map   = max(avg_evaporation_map, 0);
        avg_transpiration_map = max(avg_transpiration_map, 0);

        BC_States.current_evaporation_map   = localGather(avg_evaporation_map);
        BC_States.current_transpiration_map = localGather(avg_transpiration_map);

        BC_States.current_avg_spatial_evaporation   = mean(avg_evaporation_map,'all','omitnan');
        BC_States.current_avg_spatial_transpiration = mean(avg_transpiration_map,'all','omitnan');

        % ETP maps are assumed in mm/day
        BC_States.delta_e_agg  = avg_evaporation_map   * (dt_step / 1440);
        BC_States.delta_tr_agg = avg_transpiration_map * (dt_step / 1440);

    else
        if k == 1
            BC_States.current_evaporation_map   = zeros(gridSize);
            BC_States.current_transpiration_map = zeros(gridSize);

            BC_States.current_evaporation_map(idx_nan)   = nan;
            BC_States.current_transpiration_map(idx_nan) = nan;

            BC_States.current_avg_spatial_evaporation   = 0;
            BC_States.current_avg_spatial_transpiration = 0;
        end
    end
end

%% ===========================
%  Local helper functions
% ===========================

function A = localGather(A)
% Gather GPU array to CPU only when needed.
if isa(A,'gpuArray')
    A = gather(A);
end
end

function idx = localClampIndex(idx, lo, hi)
% Clamp scalar index to valid bounds.
idx = max(lo, min(hi, idx));
end

function dt_default = localDefaultDtFromTimes(times, fallback_dt)
% Retained for compatibility. Not required by piecewise-constant forcing,
% but kept so the function remains self-contained.
if numel(times) >= 2
    d = diff(times(:));
    d = d(d > 0 & isfinite(d));
    if isempty(d)
        dt_default = fallback_dt;
    else
        dt_default = median(d,'omitnan');
    end
else
    dt_default = fallback_dt;
end

if isempty(dt_default) || ~isfinite(dt_default) || dt_default <= 0
    dt_default = fallback_dt;
end
end

function out = localComputeOverlapWeights_EndStamped(times, t0, t1, dt_default)
% Retained for compatibility / future tests.
% Not used by the new piecewise-constant stage/inflow/rainfall logic.

out.idx = [];
out.overlap_min = [];
out.wfrac = [];

if isempty(times) || ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
    return
end

times = double(times(:));
t0 = double(t0);
t1 = double(t1);
dt_default = double(dt_default);

if ~isfinite(dt_default) || dt_default <= 0
    return
end

valid_time = isfinite(times);
if ~any(valid_time)
    return
end

orig_idx = find(valid_time);
times_valid = times(valid_time);

[times_valid, sort_idx] = sort(times_valid, 'ascend');
orig_idx = orig_idx(sort_idx);

nTimes = numel(times_valid);
if nTimes == 0
    return
end

t_end = times_valid;
t_start = zeros(nTimes,1);

if nTimes == 1
    t_start(1) = t_end(1) - dt_default;
else
    t_start(1)        = t_end(1) - dt_default;
    t_start(2:nTimes) = t_end(1:nTimes-1);
end

dt_full = t_end - t_start;
bad = dt_full <= 0 | ~isfinite(dt_full);
if any(bad)
    t_start(bad) = t_end(bad) - dt_default;
    dt_full = t_end - t_start;
end

overlap_min = zeros(nTimes,1);

for i = 1:nTimes
    a = max(t0, t_start(i));
    b = min(t1, t_end(i));
    overlap_min(i) = max(0, b - a);
end

keep = find(overlap_min > 0);
if isempty(keep)
    return
end

dt_i = dt_full(keep);
dt_i(dt_i <= 0 | ~isfinite(dt_i)) = dt_default;

out.idx = orig_idx(keep).';
out.overlap_min = overlap_min(keep).';
out.wfrac = out.overlap_min ./ dt_i.';
end

function [spatial_rainfall, Spatial_Rainfall_Parameters] = localGetSpatialGaugeRainCached( ...
    iz, Spatial_Rainfall_Parameters, gridSize, idx_nan, RainfallInterpolatorFcn)
% Cached interpolation for gauge-based rainfall.
%
% This is still very useful even with piecewise-constant forcing:
% if the rainfall cursor has not advanced, the map does not need to be
% interpolated again.

useCache = false;

if isfield(Spatial_Rainfall_Parameters,'cache_index_rain') && ...
        isfield(Spatial_Rainfall_Parameters,'cache_map_rain')   && ...
        ~isempty(Spatial_Rainfall_Parameters.cache_index_rain)  && ...
        ~isempty(Spatial_Rainfall_Parameters.cache_map_rain)

    if Spatial_Rainfall_Parameters.cache_index_rain == iz
        useCache = true;
        spatial_rainfall = Spatial_Rainfall_Parameters.cache_map_rain;
    end
end

if ~useCache
    rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges( ...
        iz, 1:Spatial_Rainfall_Parameters.n_raingauges)';

    idx_rainfall = isnan(rainfall);

    rainfall_valid = rainfall(~idx_rainfall);
    x_coordinate   = Spatial_Rainfall_Parameters.x_coordinate(~idx_rainfall);
    y_coordinate   = Spatial_Rainfall_Parameters.y_coordinate(~idx_rainfall);

    if isempty(x_coordinate)
        spatial_rainfall = zeros(gridSize);
        spatial_rainfall(idx_nan) = nan;
    else
        spatial_rainfall = RainfallInterpolatorFcn( ...
            x_coordinate, ...
            y_coordinate, ...
            rainfall_valid, ...
            Spatial_Rainfall_Parameters.x_grid, ...
            Spatial_Rainfall_Parameters.y_grid);

        spatial_rainfall(idx_nan) = nan;
    end

    Spatial_Rainfall_Parameters.cache_index_rain = iz;
    Spatial_Rainfall_Parameters.cache_map_rain   = spatial_rainfall;
end
end

function [input_rainfall_Z, Input_Rainfall] = localGetInputRainMapCached( ...
    iz, Input_Rainfall, DEM_raster, GIS_data, flags, idx_nan, k)
% Cached reading/alignment for input rainfall rasters.

useCache = false;

if isfield(Input_Rainfall,'cache_index_rain') && ...
        isfield(Input_Rainfall,'cache_map_rain')   && ...
        ~isempty(Input_Rainfall.cache_index_rain)  && ...
        ~isempty(Input_Rainfall.cache_map_rain)

    if Input_Rainfall.cache_index_rain == iz
        useCache = true;
        input_rainfall_Z = Input_Rainfall.cache_map_rain;
    end
end

if ~useCache
    input_rainfall_Z = localReadAlignRain( ...
        Input_Rainfall.labels_Directory{iz}, ...
        DEM_raster, ...
        GIS_data, ...
        flags, ...
        k);

    % ---------------------------------------------------------------------
    % Quality control kept from your current code
    % ---------------------------------------------------------------------
    if max(input_rainfall_Z(:), [], 'omitnan') > 300
        warning('Rainfall values larger than 300 mm/h. Using previous raster for this forcing interval.')
        prevIdx = localClampIndex(iz - 1, 1, numel(Input_Rainfall.labels_Directory));

        input_rainfall_Z = localReadAlignRain( ...
            Input_Rainfall.labels_Directory{prevIdx}, ...
            DEM_raster, ...
            GIS_data, ...
            flags, ...
            k);

    elseif any(isnan(input_rainfall_Z(~idx_nan)))
        warning('Rainfall has NaNs inside the domain. Using previous raster for this forcing interval.')
        prevIdx = localClampIndex(iz - 1, 1, numel(Input_Rainfall.labels_Directory));

        input_rainfall_Z = localReadAlignRain( ...
            Input_Rainfall.labels_Directory{prevIdx}, ...
            DEM_raster, ...
            GIS_data, ...
            flags, ...
            k);
    end

    input_rainfall_Z((idx_nan == 0) & isnan(input_rainfall_Z)) = 0;
    input_rainfall_Z(idx_nan) = nan;

    Input_Rainfall.cache_index_rain = iz;
    Input_Rainfall.cache_map_rain   = input_rainfall_Z;
end
end

function Z = localReadAlignRain(labelCell, DEM_raster, GIS_data, flags, k)
% Read rainfall raster and align it to the model DEM grid if needed.

fp = string(labelCell);
if strlength(fp) == 0
    error('Input rainfall map path is empty for this timestep.');
end

[A, rR] = readgeoraster(fp);

if isprop(rR,'CoordinateSystemType') && strcmpi(rR.CoordinateSystemType,'geographic')
    if ~isprop(rR,'GeographicCRS') || isempty(rR.GeographicCRS)
        rR.GeographicCRS = geocrs(4326);
    end
end

if flags.flag_resample
    targetRes = GIS_data.resolution_resample;
else
    targetRes = DEM_raster.cellsize;
end

needResample = false;

if isa(rR,'map.rasterref.GeographicCellsReference')
    if abs(rR.CellExtentInLatitude - targetRes) > 1e-12
        needResample = true;
    end
else
    if abs(rR.CellExtentInWorldX - targetRes) > 1e-9
        needResample = true;
    end
end

% Force alignment if raster size differs from DEM size
if ~isequal(size(A,1), DEM_raster.size(1)) || ~isequal(size(A,2), DEM_raster.size(2))
    needResample = true;
end

if needResample
    reset_cache = (k == 1);
    A = raster_cutter(DEM_raster, rR, A, reset_cache);
    Z = A.Z;
else
    if isstruct(A) && isfield(A,'Z')
        Z = A.Z;
    else
        Z = A;
    end
end

Z = double(Z);
Z = max(Z, 0);
end

function Z = localReadAlignET(labelCell, DEM_raster, GIS_data, flags, k)
% Read ET raster and align to DEM if needed.

fp = string(labelCell);
if strlength(fp) == 0
    error('Input ETP map path is empty for this timestep.');
end

[A, rR] = readgeoraster(fp);

if isprop(rR,'CoordinateSystemType') && strcmpi(rR.CoordinateSystemType,'geographic')
    if ~isprop(rR,'GeographicCRS') || isempty(rR.GeographicCRS)
        rR.GeographicCRS = geocrs(4326);
    end
end

if flags.flag_resample
    targetRes = GIS_data.resolution_resample;
else
    targetRes = DEM_raster.cellsize;
end

needResample = false;

if isa(rR,'map.rasterref.GeographicCellsReference')
    if abs(rR.CellExtentInLatitude - targetRes) > 1e-12
        needResample = true;
    end
else
    if abs(rR.CellExtentInWorldX - targetRes) > 1e-9
        needResample = true;
    end
end

if ~isequal(size(A,1), DEM_raster.size(1)) || ~isequal(size(A,2), DEM_raster.size(2))
    needResample = true;
end

if needResample
    reset_cache = (k == 1);
    A = raster_cutter(DEM_raster, rR, A, reset_cache);
    Z = A.Z;
else
    if isstruct(A) && isfield(A,'Z')
        Z = A.Z;
    else
        Z = A;
    end
end

Z = double(Z);
Z = max(Z, 0);
end

function A = localCastToModelPrecision(Z, flags)
% Retained for compatibility. Not used actively in this version unless you
% decide to re-enable casting/moving to GPU here.

if flags.flag_single == 1
    A = single(Z);
    if flags.flag_GPU == 1
        A = gpuArray(A);
    end
else
    A = Z;
    if flags.flag_GPU == 1
        A = gpuArray(A);
    end
end
end

function M = localAccumMap(X, w)
% Retained for compatibility. Not needed by the new piecewise-constant
% rainfall logic, but kept so the file remains self-contained.

M = X * w;
end

function cursor = localAdvanceEndStampedCursor(times, t_now, cursor)
% Move cursor forward for end-stamped piecewise-constant forcing.
%
% Interpretation:
%   times(i) is the END of interval i.
%
% Example:
%   times = [30 60 90]
%   then:
%       interval 1 ends at 30
%       interval 2 ends at 60
%       interval 3 ends at 90
%
% If current model time t_now is within interval i, this function returns i.
%
% This is cheap because model time only moves forward.

times = double(times(:));
n = numel(times);

if isempty(times) || n == 0
    cursor = [];
    return
end

if isempty(cursor) || ~isfinite(cursor) || cursor < 1
    cursor = 1;
end

cursor = min(max(1, cursor), n);

while cursor < n && t_now > times(cursor)
    cursor = cursor + 1;
end
end

function [avg_rate, cursor] = localAverageRateOverStep(times, values, t0, t1, cursor)
% Exact time-averaged rate over [t0,t1] for an end-stamped piecewise-constant series.
%
% times(i) = END of interval i
% values(i) = constant value during interval i
%
% Returns avg_rate over [t0,t1].

avg_rate = 0;

times  = double(times(:));
values = double(values(:));

n = min(numel(times), numel(values));
times  = times(1:n);
values = values(1:n);

if isempty(times) || isempty(values) || ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
    return
end

dt_total = t1 - t0;
accum = 0;

% move cursor close to t0
cursor = localAdvanceEndStampedCursor(times, t0, cursor);
if isempty(cursor)
    cursor = 1;
end

dt_force = localEstimateConstantDt(times);

for i = cursor:n
    if i == 1
        t_start_i = times(1) - dt_force;
    else
        t_start_i = times(i-1);
    end
    t_end_i = times(i);

    overlap = max(0, min(t1, t_end_i) - max(t0, t_start_i));

    if overlap > 0 && isfinite(values(i))
        accum = accum + values(i) * overlap;
    end

    if t_end_i >= t1
        cursor = i;
        break
    end
end

avg_rate = accum / dt_total;
end

function [avg_state, cursor] = localAverageStateOverStep(times, values, t0, t1, cursor)
% Exact time-averaged state over [t0,t1] for an end-stamped piecewise-constant series.
%
% times(i) = END of interval i
% values(i) = constant value during interval i
%
% Returns average state over [t0,t1].

avg_state = 0;

times  = double(times(:));
values = double(values(:));

n = min(numel(times), numel(values));
times  = times(1:n);
values = values(1:n);

if isempty(times) || isempty(values) || ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
    return
end

dt_total = t1 - t0;
accum = 0;

cursor = localAdvanceEndStampedCursor(times, t0, cursor);
if isempty(cursor)
    cursor = 1;
end

dt_force = localEstimateConstantDt(times);

for i = cursor:n
    if i == 1
        t_start_i = times(1) - dt_force;
    else
        t_start_i = times(i-1);
    end
    t_end_i = times(i);

    overlap = max(0, min(t1, t_end_i) - max(t0, t_start_i));

    if overlap > 0 && isfinite(values(i))
        accum = accum + values(i) * overlap;
    end

    if t_end_i >= t1
        cursor = i;
        break
    end
end

avg_state = accum / dt_total;
end

function dt_force = localEstimateConstantDt(times)
% Estimate constant forcing dt from end-stamped time vector.
times = double(times(:));

if numel(times) >= 2
    dt_force = times(2) - times(1);
else
    dt_force = 0;
end

if ~isfinite(dt_force) || dt_force <= 0
    dt_force = 0;
end
end

function [input_et_Z, Input_ET] = localGetInputETMapCached( ...
    iz, Input_ET, DEM_raster, GIS_data, flags, idx_nan, k)

useCache = false;

if isfield(Input_ET,'cache_index_et') && ...
        isfield(Input_ET,'cache_map_et')   && ...
        ~isempty(Input_ET.cache_index_et)  && ...
        ~isempty(Input_ET.cache_map_et)

    if Input_ET.cache_index_et == iz
        useCache = true;
        input_et_Z = Input_ET.cache_map_et;
    end
end

if ~useCache
    input_et_Z = localReadAlignET( ...
        Input_ET.labels_Directory{iz}, ...
        DEM_raster, ...
        GIS_data, ...
        flags, ...
        k);

    input_et_Z((idx_nan == 0) & isnan(input_et_Z)) = 0;
    input_et_Z(idx_nan) = nan;

    Input_ET.cache_index_et = iz;
    Input_ET.cache_map_et   = input_et_Z;
end
end

function [avg_rain_map, cursor, Input_Rainfall] = ...
    localAverageInputRainMapOverStep(times, t0, t1, cursor, Input_Rainfall, DEM_raster, GIS_data, flags, idx_nan, k)

times = double(times(:));
gridSize = size(DEM_raster.Z);

% Default output
avg_rain_map = zeros(gridSize);

if isempty(times) || t1 <= t0
    avg_rain_map(idx_nan) = nan;
    return
end

dt_total = t1 - t0;

% Move cursor to active interval
cursor = localAdvanceEndStampedCursor(times, t0, cursor);
if isempty(cursor)
    cursor = 1;
end

% Use precomputed forcing dt if available
if isfield(Input_Rainfall,'dt_force') && ~isempty(Input_Rainfall.dt_force)
    dt_force = Input_Rainfall.dt_force;
else
    dt_force = localEstimateConstantDt(times);
end

% Bounds of current forcing interval
if cursor == 1
    t_start_i = times(1) - dt_force;
else
    t_start_i = times(cursor-1);
end
t_end_i = times(cursor);

% ---- FAST PATH: entire model step lies inside one forcing interval ----
if t0 >= t_start_i && t1 <= t_end_i
    [avg_rain_map, Input_Rainfall] = localGetInputRainMapCached( ...
        cursor, Input_Rainfall, DEM_raster, GIS_data, flags, idx_nan, k);
    return
end

% ---- SLOW PATH: only if the model step spans >1 forcing interval ----
accum_map = zeros(gridSize);

for i = cursor:numel(times)
    if i == 1
        t_start_i = times(1) - dt_force;
    else
        t_start_i = times(i-1);
    end
    t_end_i = times(i);

    overlap = max(0, min(t1, t_end_i) - max(t0, t_start_i));

    if overlap > 0
        [rain_map_i, Input_Rainfall] = localGetInputRainMapCached( ...
            i, Input_Rainfall, DEM_raster, GIS_data, flags, idx_nan, k);

        accum_map = accum_map + overlap .* rain_map_i;
    end

    if t_end_i >= t1
        cursor = i;
        break
    end
end

avg_rain_map = accum_map / dt_total;
avg_rain_map(idx_nan) = nan;
end

function [avg_rain_map, cursor, Spatial_Rainfall_Parameters] = ...
    localAverageSpatialGaugeRainOverStep(times, rainfall_raingauges, t0, t1, cursor, Spatial_Rainfall_Parameters, gridSize, idx_nan, RainfallInterpolatorFcn)

avg_rain_map = zeros(gridSize);
avg_rain_map(idx_nan) = nan;

times = double(times(:));
if isempty(times) || t1 <= t0
    return
end

dt_total = t1 - t0;
cursor = localAdvanceEndStampedCursor(times, t0, cursor);
if isempty(cursor)
    cursor = 1;
end

dt_force = localEstimateConstantDt(times);

accum_map = zeros(gridSize);
accum_map(idx_nan) = 0;

for i = cursor:numel(times)
    if i == 1
        t_start_i = times(1) - dt_force;
    else
        t_start_i = times(i-1);
    end
    t_end_i = times(i);

    overlap = max(0, min(t1, t_end_i) - max(t0, t_start_i));

    if overlap > 0
        [rain_map_i, Spatial_Rainfall_Parameters] = localGetSpatialGaugeRainCached( ...
            i, Spatial_Rainfall_Parameters, gridSize, idx_nan, RainfallInterpolatorFcn);

        accum_map = accum_map + rain_map_i * overlap;
    end

    if t_end_i >= t1
        cursor = i;
        break
    end
end

avg_rain_map = accum_map / dt_total;
avg_rain_map(idx_nan) = nan;
end

function [avg_et_map, cursor, Input_ET] = ...
    localAverageInputETMapOverStep(times, labels_Directory, t0, t1, cursor, Input_ET, DEM_raster, GIS_data, flags, idx_nan, k)

gridSize = size(DEM_raster.Z);

avg_et_map = zeros(gridSize);
avg_et_map(idx_nan) = nan;

times = double(times(:));
if isempty(times) || t1 <= t0
    return
end

dt_total = t1 - t0;
cursor = localAdvanceEndStampedCursor(times, t0, cursor);
if isempty(cursor)
    cursor = 1;
end

dt_force = localEstimateConstantDt(times);

accum_map = zeros(gridSize);
accum_map(idx_nan) = 0;

for i = cursor:numel(times)
    if i == 1
        t_start_i = times(1) - dt_force;
    else
        t_start_i = times(i-1);
    end
    t_end_i = times(i);

    overlap = max(0, min(t1, t_end_i) - max(t0, t_start_i));

    if overlap > 0
        [et_map_i, Input_ET] = localGetInputETMapCached( ...
            i, Input_ET, DEM_raster, GIS_data, flags, idx_nan, k);

        accum_map = accum_map + et_map_i * overlap;
    end

    if t_end_i >= t1
        cursor = i;
        break
    end
end

avg_et_map = accum_map / dt_total;
avg_et_map(idx_nan) = nan;
end