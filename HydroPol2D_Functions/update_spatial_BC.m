% Update Spatial Boundary Conditions
% Goal: Update boundary conditions and save spatial maps of these, if
% required
% Developer: Marcus Nobrega, Ph.D.

%% Stage-Hydrograph Boundary Condition
if flags.flag_stage_hydrograph == 1
    for z = 1:Stage_Parameters.n_stage_gauges
        z1 = find(Stage_Parameters.time_stage > t_previous,1,'first'); % begin of the time-step
        z2 = find(Stage_Parameters.time_stage <= t,1,'last'); % end of the time-step
        if isempty(z1)
            z1 = 1;
        end
        if isempty(z2) || z2 < z1
            z2 = z1;
        end
        if k == 1
            stage_depth_previous = 0;
        else
            stage_depth_previous = stage_depth;
        end
        stage_depth(z,1) = Stage_Parameters.stage(z2,z);
    end
end

% Stage Boundary Condition
if flags.flag_stage_hydrograph == 1
    for i = 1:Stage_Parameters.n_stage_gauges
        stage_cells = Wshed_Properties.stage_mask(:,:,i);
        depths.d_t(logical(Wshed_Properties.stage_mask(:,:,i))) = 1000*stage_depth(i,1)*stage_cells(stage_cells == 1); % mm
    end
end


%% Inflows for next time-step
% Agregating Inflows to the New Time-step
if flags.flag_inflow > 0
    for z = 1:Inflow_Parameters.n_stream_gauges
        z1 = find(BC_States.time_deltainflow > t_previous,1,'first'); % begin of the time-step
        z2 = find(BC_States.time_deltainflow <= t,1,'last'); % end of the time-step
        if isempty(z1)
            z1 = 1;
        end
        if isempty(z2) || z2 < z1
            z2 = z1;
        end
        if time_step >= time_step_model
            BC_States.delta_inflow_agg(z,1) = mean(BC_States.delta_inflow(z,z1:z2))/(time_step_model*60)*time_step*60;
        else
            BC_States.delta_inflow_agg(z,1) = BC_States.delta_inflow(z,z1)/(time_step_model*60)*time_step*60;
        end

    end
end


%% Inflows Boundary Condition
if flags.flag_inflow == 1
    BC_States.inflow = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2)); % This is to solve spatially, don't delete
    for i = 1:Inflow_Parameters.n_stream_gauges
        BC_States.inflow = BC_States.inflow + BC_States.delta_inflow_agg(i)*Wshed_Properties.inflow_cells(:,:,i); % mm
        if flags.flag_subgrid == 1
            % In this case, we need to correct the depths to ensure correct
            % mass balance
            BC_States.inflow = BC_States.inflow.*Wshed_Properties.cell_area./(C_a); % Channels can have smaller area
        end
    end
end



%% Agregating Precipitation to the New Time-step
if flags.flag_rainfall > 0
    if flags.flag_spatial_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
        z1 = find(running_control.time_deltap > t_previous,1,'first'); % begin of the time-step
        z2 = find(running_control.time_deltap <= t,1,'last'); % end of the time-step
        if isempty(z2) | z2 < z1
            z2 = z1;
        end
        if time_step >= time_step_model
            BC_States.delta_p_agg = mean(BC_States.delta_p(1,z1:z2))/(time_step_model*60)*time_step*60; % mm
        else
            BC_States.delta_p_agg = BC_States.delta_p(1,z1)/(time_step_model*60)*time_step*60;  % mm
        end
        if isnan(BC_States.delta_p_agg)
            warning('No rainfall data for this time. Assuming as 0 mm/h.')
            BC_States.delta_p_agg = 0;
        end
    elseif flags.flag_spatial_rainfall == 1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
        % Spatial Rainfall
        % Code for Spatial-Varying Rainfall
        % Times of full rainfall dataset
        z1 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg <= t_previous,1,'last');
        z2 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg <= t,1,'last');

        % Times for aggregated rainfall dataset
        zz1 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t_previous,1,'last');
        zz2 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t,1,'last');

        if Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2)  == Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) % No difference between them, same time-step
            Rainfall_Parameters.index_aggregation = 1; % We save only 1 map
        elseif zz1 == zz2 && z2 ~= z1 % We save more than 1 map
            Rainfall_Parameters.index_aggregation = Rainfall_Parameters.index_aggregation + (z2-z1);
        elseif zz2 > zz1
            Rainfall_Parameters.index_aggregation = 1;
        end

        if z1 ~= z2 || z2 == length(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg)
            Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1); % Coordinates (easting) of each rain gauge
            Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2); % Coordinates (northing) of each rain gauge
            Spatial_Rainfall_Parameters.x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,2)]'; % Pixel easting coordinates
            Spatial_Rainfall_Parameters.y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,1)]'; % Pixel northing coordinates

            % Spatial Rainfall
            if z2 == length(Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg)
                spatial_rainfall = zeros(size(Elevation_Properties.elevation_cell,1),size(Elevation_Properties.elevation_cell,2));

            else
                rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(z2,1:Spatial_Rainfall_Parameters.n_raingauges)'; % Values of rainfall at t for each rain gauge
                % idx_rainfall = logical(isnan(rainfall) | rainfall == 0);
                idx_rainfall = logical(isnan(rainfall));
                rainfall(idx_rainfall) = []; % Taking out nans
                Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = []; % Taking out nans
                Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = []; % Taking out nans
                if isempty(Spatial_Rainfall_Parameters.x_coordinate) % No rainfall
                    spatial_rainfall = zeros(size(Elevation_Properties.elevation_cell));
                    spatial_rainfall(idx_nan) = nan;
                else
                    [spatial_rainfall] = Rainfall_Interpolator(Spatial_Rainfall_Parameters.x_coordinate,Spatial_Rainfall_Parameters.y_coordinate,rainfall,Spatial_Rainfall_Parameters.x_grid,Spatial_Rainfall_Parameters.y_grid); % Interpolated Values
                    spatial_rainfall(idx_nan) = nan;
                end

                if nansum(nansum(spatial_rainfall)) > 0
                    sprintf('Rainfall interpolated')
                end
            end
            if zz2 > zz1 % Saving Maps
                % Maps.Hydro.spatial_rainfall_maps(:,:,zz2) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                % zzz = Maps.Hydro.spatial_rainfall_maps(:,:,zz2);
                Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = mean(rainfall_spatial_aggregation,3);
                zzz = Maps.Hydro.spatial_rainfall_maps(:,:,saver_count);
                BC_States.average_spatial_rainfall(zz2,1) = mean(zzz(zzz>=0));
            end
            BC_States.delta_p_agg = spatial_rainfall/3600*time_step*60; % Matrix of delta P for each pixel
            % Saving Maps for Spatial Aggregation
            rainfall_spatial_aggregation(:,:,Rainfall_Parameters.index_aggregation) = gather(spatial_rainfall); % Saving high resolution map
        end
    elseif flags.flag_input_rainfall_map == 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1
        % Input Rainfall Maps
        % Times of full rainfall dataset
        z1_input = find(Input_Rainfall.time <= t_previous,1,'last');
        z2_input = find(Input_Rainfall.time <= t,1,'last');

        % Times for aggregated rainfall dataset
        zz1 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t_previous,1,'last');
        zz2 = find(Spatial_Rainfall_Parameters.rainfall_spatial_duration <= t,1,'last');

        if Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg(2)  == Spatial_Rainfall_Parameters.rainfall_spatial_duration(2) % No difference between them, same time-step
            Rainfall_Parameters.index_aggregation = 1; % We save only 1 map
        elseif zz1 == zz2 && z1_input ~= z2_input % We save more than 1 map
            Rainfall_Parameters.index_aggregation = Rainfall_Parameters.index_aggregation + (z2_input-z1_input);
        elseif zz2 > zz1
            Rainfall_Parameters.index_aggregation = 1;
        end

        if z2_input > z1_input
            % Read Geotiff
            % input_rainfall = GRIDobj(Input_Rainfall.labels_Directory{z2_input}{1});
            try
                [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{z2_input}{1}),'CoordinateSystemType','geographic');                % Filtering Non-physical values and extreme wrong values
                if max(max(input_rainfall)) > 300 % Very unlikely to have a value like this
                    warning('Rainfall values larger than 300 mm/h. We are neglecting this raster and collecting the previous one for this time-step.')
                    [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{z2_input - 1}{1}),'CoordinateSystemType','geographic');                
                elseif isnan(sum(input_rainfall(~idx_nan)))
                    warning('Rainfall has nan values inside of the domain. Using the previous rainfall as input for this time-step.')
                    [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{z2_input - 1}{1}),'CoordinateSystemType','geographic');                                    
                end
                % No negative rainfall
                input_rainfall = max(input_rainfall,0); % No negative values                  
                % Reproject the coordinates from EPSG:4326 to EPSG:3857
                rR.GeographicCRS=geocrs(4326);
                if rR.CellExtentInLatitude ~= GIS_data.resolution_resample
                    input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,1);
                end
            catch
                [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{z2_input}{1}));
                if rR.CellExtentInWorldX ~= GIS_data.resolution_resample
                    input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,0);
                end
            end

            if flags.flag_single == 1
                input_rainfall = single(input_rainfall.Z); % Only the values
                if flags.flag_GPU == 1
                    input_rainfall = gpuArray(input_rainfall); % Only the values
                end
            else
                input_rainfall = (input_rainfall.Z); % Only the values
                if flags.flag_GPU == 1
                    input_rainfall = gpuArray(input_rainfall); % Only the values
                end
            end
            if zz2 > zz1 % Saving Maps
                % Maps.Hydro.spatial_rainfall_maps(:,:,zz2) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                % zzz = Maps.Hydro.spatial_rainfall_maps(:,:,zz2);
                Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = mean(rainfall_spatial_aggregation,3);
                zzz = Maps.Hydro.spatial_rainfall_maps(:,:,saver_count);
                BC_States.average_spatial_rainfall(zz2,1) = mean(zzz(zzz>=0));
            end
            input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;
            % Maps.Hydro.spatial_rainfall_maps(:,:,z2_input) = input_rainfall;
            BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
            % BC_States.average_spatial_rainfall(z2_input,1) = mean(input_rainfall(input_rainfall>=0));
            rainfall_spatial_aggregation(:,:,Rainfall_Parameters.index_aggregation) = gather(input_rainfall); % Saving high resolution map
        else
            if k == 1 % Very first time-step
                try
                [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{1,:}),'CoordinateSystemType','geographic');
                % Filtering Non-physical values and extreme wrong values
                if max(max(input_rainfall)) > 300   % Very unlikely to have a value like this
                    error('Rainfall values larger than 300 mm/h. We are neglecting this raster and collecting the previous one for this time-step.')
                elseif isnan(sum(input_rainfall(~idx_nan)))
                    error('Rainfall has nan values inside of the domain. Please fix it.')
                end
                input_rainfall = max(input_rainfall,0); % No negative values
                
                % Reproject the coordinates from EPSG:4326 to EPSG:3857
                rR.GeographicCRS=geocrs(4326);
                if flags.flag_resample
                    if rR.CellExtentInLatitude ~= GIS_data.resolution_resample
                        input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,1);
                    end
                else
                    if rR.CellExtentInLatitude ~= DEM_raster.cellsize
                        input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,1);
                    end
                end
            catch
                [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{1,:}));
                % Filtering Non-physical values and extreme wrong values
                if max(max(input_rainfall)) > 300 % Very unlikely to have a value like this
                    error('Rainfall values larger than 300 mm/h. We are neglecting this raster and collecting the previous one for this time-step.')
                end
                % No negative rainfall
                input_rainfall = max(input_rainfall,0); % No negative values                
                if flags.flag_resample
                    if rR.CellExtentInWorldX ~= GIS_data.resolution_resample
                        input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,0);
                    end
                else
                    if rR.CellExtentInWorldX ~= DEM_raster.cellsize
                        input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,0);
                    end
                end
            end

               input_rainfall = input_rainfall.Z; % Only the values                
               BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
            end
        end
    elseif flags.flag_input_rainfall_map == 0 && flags.flag_satellite_rainfall == 1 && flags.flag_real_time_satellite_rainfall ~= 1
        % Satellite Rainfall
        if flags.flag_satellite_rainfall == 1
            % Maps with generally coarser resolution
            recording_parameters.actual_record_state_rainfall = find(running_control.time_records < t,1,'last');
            recording_parameters.delta_record_rainfall = recording_parameters.actual_record_state_rainfall - recording_parameters.last_record_maps_rainfall;
            recording_parameters.last_record_maps_rainfall = recording_parameters.actual_record_state_rainfall;

            % Input Rainfall Maps
            % Times of full rainfall dataset
            z1_input = find(Input_Rainfall.time <= t_previous,1,'last');
            z2_input = find(Input_Rainfall.time <= t,1,'last');
            if z2_input > z1_input
                product = 'PDIRNow1hourly';
                [rainfall_raster, register,~] = Satellite_rainfall_processing([],[],register,product,date_begin,date_end,flags.flag_satellite_rainfall,flags.flag_real_time_satellite_rainfall,DEM_raster);
                input_rainfall = rainfall_raster.Z;
                input_rainfall((idx_nan == 0 & isnan(input_rainfall))) = 0;

                if recording_parameters.delta_record_rainfall > 0
                    record_map_indice = recording_parameters.last_record_maps_rainfall;
                    % Maps.Hydro.spatial_rainfall_maps(:,:,record_map_indice) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                    % zzz = Maps.Hydro.spatial_rainfall_maps(:,:,record_map_indice);
                    Maps.Hydro.spatial_rainfall_maps(:,:,saver_count) = mean(rainfall_spatial_aggregation,3); % Average in whole duration
                    zzz = Maps.Hydro.spatial_rainfall_maps(:,:,saver_count);
                    BC_States.average_spatial_rainfall(record_map_indice,1) = mean(zzz(zzz>=0));
                    Rainfall_Parameters.index_aggregation = 0; % Returning to first rainfall
                end

                Rainfall_Parameters.index_aggregation = Rainfall_Parameters.index_aggregation + (z2_input-z1_input);
                % Map from previous timestep
                rainfall_spatial_aggregation(:,:,Rainfall_Parameters.index_aggregation) = input_rainfall;

            end
        end
    elseif flags.flag_input_rainfall_map == 0 && flags.flag_satellite_rainfall ~= 1 && flags.flag_real_time_satellite_rainfall == 1
        if flags.flag_real_time_satellite_rainfall == 1
            %%%%%%%%%%%%%%%%%
        end
    end
end

% Correcting Rainfall Volumes to the sub-grid model
BC_States.delta_p_agg = BC_States.delta_p_agg.*Wshed_Properties.cell_area./C_a; % Here we assume that all rainfall goes directly to the channel

%% Aggregating ETP for next time-step

if flags.flag_ETP == 1 && flags.flag_input_ETP_map ~= 1
    z1 = find(ETP_Parameters.climatologic_spatial_duration <= t_previous,1,'last');
    z2 = find(ETP_Parameters.climatologic_spatial_duration <= t,1,'last');
    if isempty(z1) && isempty(z2)
        % No data, No ETP.
        Hydro_States.ETP = zeros(size(Elevation_Properties.elevation_cell));
    else
        if ~isempty(z1) && z1 == z2 && z2 == length(ETP_Parameters.climatologic_spatial_duration)
            % Outside of ETP data maximum duration
            Hydro_States.ETP = zeros(size(Elevation_Properties.elevation_cell));
            % Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
            Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP; % Saving ETP Maps
            ETR_save(:,:,z2) = Hydro_States.ETR ;
        elseif  (isempty(z1) && z2 > 0) || z2 > z1 && z2 < length(ETP_Parameters.climatologic_spatial_duration)
            if flags.flag_GPU == 1 || flags.flag_single == 1
                day_of_year = day(extra_parameters.ETP.time_ETP(z2,1),'dayofyear');
            else
                day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
            end
            [Hydro_States.ETP] = ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2),Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations,ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations,ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);
            if nansum(nansum(Hydro_States.ETP)) == 0
                if saver_count==1
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,12);
                else
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,saver_count-1);
                end
                Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,z2-1);
                Hydro_States.ETR  = ETR_save(:,:,z2-1);
                if nansum(nansum(Hydro_States.ETP)) == 0
                    warning('No ETP and ETR data. Assuming it equals 0')
                    Hydro_States.ETP = zeros(size(DEM));
                    Hydro_States.ETP(idx_nan) = nan;
                    Hydro_States.ETR  = zeros(size(DEM));
                    Hydro_States.ETR (idx_nan) = nan;
                end
            end
            % Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
            Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP; % Saving ETP Maps
            ETR_save(:,:,z2) = Hydro_States.ETR ;
        elseif z1 == 1 && z2 == 1 && k == 1
            % First data. We assume a constant ETP using
            % 1st data
            if flags.flag_GPU == 1 || flags.flag_single == 1
                day_of_year = day(extra_parameters.ETP.time_ETP(z2,1),'dayofyear');
            else
                day_of_year = day(ETP_Parameters.time_ETP(z2,1),'dayofyear');
            end
            [Hydro_States.ETP] = ETP_model(z2,day_of_year,ETP_Parameters.coordinates_stations(:,1),ETP_Parameters.coordinates_stations(:,2),Spatial_Rainfall_Parameters.x_grid',Spatial_Rainfall_Parameters.y_grid',ETP_Parameters.maxtemp_stations,ETP_Parameters.mintemp_stations,ETP_Parameters.avgtemp_stations,ETP_Parameters.u2_stations,ETP_Parameters.ur_stations,ETP_Parameters.G_stations,ETP_Parameters.DEM_etp,ETP_Parameters.lat,ETP_Parameters.Krs,ETP_Parameters.alfa_albedo_input,idx_nan);
            if nansum(nansum(Hydro_States.ETP)) == 0
                if saver_count==1
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,12);
                else
                    Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,saver_count-1);
                end
                Hydro_States.ETP = Maps.Hydro.ETP_save(:,:,z2-1);
                Hydro_States.ETR  = ETR_save(:,:,z2-1);
            end
            % Maps.Hydro.ETP_save(:,:,z2) = Hydro_States.ETP; % Saving ETP Maps
            Maps.Hydro.ETP_save(:,:,saver_count) = Hydro_States.ETP; % Saving ETP Maps
            if z2 == 1
                ETR_save(:,:,z2) = Hydro_States.ETP ; % this might be incorrect
            else
                ETR_save(:,:,z2) = Hydro_States.ETR ;
            end
        end
    end
end

% Reading first ETP tiles
if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1 % Input Rainfall Maps
    z1 = find(ETP_Parameters.climatologic_spatial_duration <= t_previous,1,'last');
    z2 = find(ETP_Parameters.climatologic_spatial_duration <= t,1,'last');
    if z2>z1

    % input_rainfall = GRIDobj(string(Input_Rainfall.labels_Directory{1,:}));
    try
        [input_evaporation,rR] = readgeoraster(string(Input_Evaporation.labels_Directory{z2,:}),'CoordinateSystemType','geographic');
        [input_transpiration,rR] = readgeoraster(string(Input_Transpiration.labels_Directory{z2,:}),'CoordinateSystemType','geographic');

        % Reproject the coordinates from EPSG:4326 to EPSG:3857
        rR.GeographicCRS=geocrs(4326);
        if flags.flag_resample
            if rR.CellExtentInLatitude ~= GIS_data.resolution_resample
                input_evaporation = raster_cutter(DEM_raster,rR,input_evaporation,1);
                input_transpiration = raster_cutter(DEM_raster,rR,input_transpiration,1);
            end
        else
            if rR.CellExtentInLatitude ~= DEM_raster.cellsize
                input_evaporation = raster_cutter(DEM_raster,rR,input_evaporation,1);
                input_transpiration = raster_cutter(DEM_raster,rR,input_transpiration,1);
            end
        end
    catch
        [input_evaporation,rR] = readgeoraster(string(Input_Evaporation.labels_Directory{1,:}));
        [input_transpiration,rR] = readgeoraster(string(Input_Transpiration.labels_Directory{1,:}));
        if flags.flag_resample
            if rR.CellExtentInWorldX ~= GIS_data.resolution_resample
                input_evaporation = raster_cutter(DEM_raster,rR,input_evaporation,0);
                input_transpiration = raster_cutter(DEM_raster,rR,input_transpiration,0);
            end
        else
            if rR.CellExtentInWorldX ~= DEM_raster.cellsize
                input_evaporation = raster_cutter(DEM_raster,rR,input_evaporation,0);
                input_transpiration = raster_cutter(DEM_raster,rR,input_transpiration,0);
            end
        end
    end
    input_evaporation = input_evaporation.Z; % Only the values
    input_transpiration = input_transpiration.Z; % Only the values
    Maps.Hydro.spatial_evaporation_maps(:,:,1) = input_evaporation;
    Maps.Hydro.spatial_transpiration_maps(:,:,1) = input_transpiration;


    % Spatial_Rainfall_Parameters.rainfall_spatial_duration = Input_Rainfall.time';
    BC_States.average_spatial_evaporation(1,1) = mean(input_evaporation(input_evaporation>=0));
    BC_States.average_spatial_transpiration(1,1) = mean(input_transpiration(input_transpiration>=0));
end
end
