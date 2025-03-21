function [d_db,r_db,i_db,v_db,S_p,A,RA,dates,t_d,system_output,save_counter] = HydroPol2D_real_time_dashboard(ax_d,ax_r,ax_i,ax_v,ax_date, ax_iter, d_t, r_t,i_t,v_t,d_db,r_db,i_db,v_db,k,S_p,A,RA,DEM_raster,register_data,dates,t_d,system_output,ax_system_output,flag_message,gauges_data,ax_list,flag_unit_forecast,save_counter)
%% Real-Time dashboard

    if k == 1
        [Spectrum,depth_ramp,terrain_ramp] = coloramps(); % Run coloramp function
        %Store the inditial data into the multiarrays, depths, rainfall,
        %infiltration, and velocity, (human risk is missing)...
        if flag_unit_forecast ~=1 
            d_db(:,:,end)=d_t/1000;
            r_db(:,:,end)=r_t;
            i_db(:,:,end)=i_t;
            v_db(:,:,end)=v_t;
        end
        %% Creating the custom basemap
        basemapName = "openstreetmap";
        url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png"; 
        copyright = char(uint8(169));
        attribution = copyright + "OpenStreetMap contributors";
        addCustomBasemap(basemapName,url,"Attribution",attribution)
        % Getting lat and lon from the study area
        
        [lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
        time_zone = timezone(mean(lon));
        latlim = [lat(1) lat(2)];
        lonlim = [lon(1) lon(2)];
        % Retriving the basemap image
        [A,RA,attribA] = readBasemapImage(basemapName,latlim,lonlim);

        %% Creating a Shapefile from the DEM as reference
        
        % Creating a binary mask
        binaryMask = ~isnan(DEM_raster.Z);
        boundaries = bwboundaries(binaryMask);
        % Pre-allocate arrays to store combined X and Y coordinates
        combinedX = [];
        combinedY = [];
        
        % Combine all boundary coordinates into a single array
        for k = 1:numel(boundaries)
            boundary = boundaries{k};
            X = boundary(:, 2);
            Y = boundary(:, 1);
            combinedX = [combinedX; X; NaN]; % Add NaN to separate polygons
            combinedY = [combinedY; Y; NaN]; % Add NaN to separate polygons
        end
        % Remove the trailing NaNs at the end (optional)
        combinedX = combinedX(1:end-1);
        combinedY = combinedY(1:end-1);
        % making the geostruct to alocate the shapefile
        % cheking if CRS of the project is on WGS84
        web_mercator_crs = projcrs(3857);
        no_plot=0;
        if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= web_mercator_crs.Name
            no_plot = 1;
        else
           S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
           S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)]; % Calculate bounding box for each polygon
           S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
           S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
        end
        
        % Setting up the monitors
        xmax = size(DEM_raster.Z,1);
        xend = xmax;
        ymax = size(DEM_raster.Z,2);
        yend = ymax;
        xbegin = 1;
        ybegin = 1;
        x_grid = DEM_raster.georef.SpatialRef.YWorldLimits(2) - DEM_raster.georef.SpatialRef.CellExtentInWorldX*[xbegin:1:xend];
        y_grid = DEM_raster.georef.SpatialRef.XWorldLimits(1) + DEM_raster.georef.SpatialRef.CellExtentInWorldX*[ybegin:1:yend];
        DEM_s1=size(DEM_raster.Z,1); DEM_s2=size(DEM_raster.Z,2);
    
        bm_d = mapshow(ax_d, A, RA, "AlphaData",0.45);hold(ax_d, 'on');
        bm_r = mapshow(ax_r, A, RA, "AlphaData",0.45);hold(ax_r, 'on');
        bm_i = mapshow(ax_i, A, RA, "AlphaData",0.45);hold(ax_i, 'on');
        bm_v = mapshow(ax_v, A, RA, "AlphaData",0.45);hold(ax_v, 'on');
        F_d = d_db([1:1:DEM_s1],[1:1:DEM_s2],end);
        F_r = r_db([1:1:DEM_s1],[1:1:DEM_s2],end);
        F_i = i_db([1:1:DEM_s1],[1:1:DEM_s2],end);
        F_v = v_db([1:1:DEM_s1],[1:1:DEM_s2],end);
        monitor_d = mesh(ax_d,y_grid,x_grid,F_d,'FaceColor','interp','LineStyle',"none");
        colormap(ax_d,depth_ramp);colorbar(ax_d);
        monitor_r = mesh(ax_r,y_grid,x_grid,F_r,'FaceColor','interp');
        colormap(ax_r,linspecer);colorbar(ax_r);
        monitor_i = mesh(ax_i,y_grid,x_grid,F_i,'FaceColor','interp');
        colormap(ax_i,linspecer);colorbar(ax_i);
        monitor_v = mesh(ax_v,y_grid,x_grid,F_v,'FaceColor','interp');
        colormap(ax_v,depth_ramp);colorbar(ax_v);
        shp_d = mapshow(ax_d,S_p,'FaceColor','n');hold(ax_d, 'on');
        shp_r = mapshow(ax_r,S_p,'FaceColor','n');hold(ax_r, 'on');
        shp_i = mapshow(ax_i,S_p,'FaceColor','n');hold(ax_i, 'on');
        shp_v = mapshow(ax_v,S_p,'FaceColor','n');hold(ax_v, 'on');
    
        dates(:) = register_data; 
        set(ax_date,'Value',dates(end));
        set(ax_system_output,'Value',system_output);
        set(ax_iter,'Text',num2str(1));
    
        set(ax_d,'FontName','Garamond');
        ax_d.XAxis.Exponent = 0; ax_d.XAxis.TickLabelFormat = '%.0f';
        ax_d.YAxis.Exponent = 0; ax_d.YAxis.TickLabelFormat = '%.0f';
        set(ax_r,'FontName','Garamond');
        ax_r.XAxis.Exponent = 0; ax_r.XAxis.TickLabelFormat = '%.0f';
        ax_r.YAxis.Exponent = 0; ax_r.YAxis.TickLabelFormat = '%.0f';
        set(ax_i,'FontName','Garamond');
        ax_i.XAxis.Exponent = 0; ax_i.XAxis.TickLabelFormat = '%.0f';
        ax_i.YAxis.Exponent = 0; ax_i.YAxis.TickLabelFormat = '%.0f';
        set(ax_v,'FontName','Garamond');
        ax_v.XAxis.Exponent = 0; ax_v.XAxis.TickLabelFormat = '%.0f';
        ax_v.YAxis.Exponent = 0; ax_v.YAxis.TickLabelFormat = '%.0f';
        % UTM Coordinates
        % Convert the string arrays to character vectors
        char_vector_cell = {};
        for i = 1:numel(gauges_data.labels_observed_string)
            char_vector_cell = [char_vector_cell; char(gauges_data.labels_observed_string{i})];
        end

        % Assign the modified cell array to ax_list.Items
        ax_list.Items = char_vector_cell;
        app.gauges_data = gauges_data;

        % TimerFcn to update plot
        % Create timer
        if flag_unit_forecast ~=1
            t_period = 1;
        else
            t_period = 1;
        end
        % Initializing variables
        forecast = 1;
        cumulator = 0; dates_end = dates(end);
        cumulator_2 = 1; dates_end = 'first try';
        t_d = timer('TimerFcn', {@updatePlot, ...
            monitor_d, d_db, bm_d, shp_d ...
            monitor_r, r_db, bm_r, shp_r, ...
            monitor_i, i_db, bm_i, shp_i, ...
            monitor_v, v_db, bm_v, shp_v, ...
            A,RA,S_p, y_grid, x_grid,DEM_s1,DEM_s2,dates,ax_date,ax_iter, ...
            ax_system_output,system_output,flag_unit_forecast,forecast,time_zone, ...
            cumulator,dates_end,cumulator_2},'Period', t_period ,'ExecutionMode','fixedSpacing'); 
        idx_rt =[];
        idx_loop = [];
        start(t_d);
    else
        if flag_message==1 
            t_d.TimerFcn{29} = system_output;
            save_counter = save_counter;
        elseif isempty(r_t) && flag_message~=1 % if there is no rainfall
            for i=length(d_db(1,1,:)):-1:2 % move all values of the matrix
                d_db(:,:,end-i+1) = d_db(:,:,end-i+2); % one layer up
                i_db(:,:,end-i+1) = i_db(:,:,end-i+2);
                v_db(:,:,end-i+1) = v_db(:,:,end-i+2);
            end
            d_db(:,:,end) = d_t/1000; % update for t+1 the other variables
            i_db(:,:,end) = i_t;
            v_db(:,:,end) = v_t;
            d_db(d_db<0)=NaN;d_db(d_db==0)=NaN;d_db(d_db<0.1)=NaN; % to plot just pixels with
            i_db(i_db<0)=NaN;i_db(i_db==0)=NaN; % relevant data and be able
            v_db(v_db<0)=NaN;v_db(v_db==0)=NaN; % to see the basemap
            t_d.TimerFcn{11} = i_db; % update all variables in t_d
            t_d.TimerFcn{15} = v_db;
            t_d.TimerFcn{3} = d_db;
            if flag_unit_forecast == 1
                save_counter = save_counter + 1;
                if save_counter == length(size(d_db,3))
                    save('test_1','i_db','v_db','d_db','r_db','dates');
                    save_counter = 0;
                end
            else
                save_counter = save_counter;
            end
        else
            for i=length(r_db(1,1,:)):-1:2 % move all values of the matrix
                r_db(:,:,end-i+1) = r_db(:,:,end-i+2); % one layer up
                dates(end-i+1) = dates(end-i+2);
            end
            r_db(:,:,end) = r_t; % update the rainfall by t+1 
            r_db(r_db<0)=NaN;r_db(r_db==0)=NaN;
            dates(end) = register_data;
            t_d.TimerFcn{7} = r_db; % update all variables in t_d
            t_d.TimerFcn{25} = dates;
            save_counter=save_counter;
        end
    end

    
    function updatePlot(~, ~, ...
            monitor_d, d_db, bm_d, shp_d, ...
            monitor_r, r_db, bm_r, shp_r, ...
            monitor_i, i_db, bm_i, shp_i, ...
            monitor_v, v_db, bm_v, shp_v, ...
            A, RA, S_p, y_grid, x_grid,DEM_s1,DEM_s2,dates,ax_date,ax_iter, ...
            ax_system_output,system_output,flag_unit_forecast,forecast,time_zone, ...
            cumulator,dates_end,cumulator_2)  
        
          if isempty(idx_rt)
            idx_rt = 1; t_d.TimerFcn{33} = 0;
          end

          % Plot update code
          set(monitor_d, 'CData',d_db([1:1:DEM_s1],[1:1:DEM_s2],idx_rt),'YData', x_grid,'XData',y_grid)
          set(monitor_r, 'CData',r_db([1:1:DEM_s1],[1:1:DEM_s2],idx_rt),'YData', x_grid,'XData',y_grid)
          set(monitor_i, 'CData',i_db([1:1:DEM_s1],[1:1:DEM_s2],idx_rt),'YData', x_grid,'XData',y_grid)
          set(monitor_v, 'CData',v_db([1:1:DEM_s1],[1:1:DEM_s2],idx_rt),'YData', x_grid,'XData',y_grid)
          
          temp = dates(idx_rt);
          
          set(ax_date,'Value',temp);
          set(ax_system_output,'Value',system_output(idx_rt/idx_rt));
          set(ax_iter,'Text',num2str(idx_rt));

          idx_rt = idx_rt + 1;
          % when the iteration number is bigger than the array and data is
          % loaded, change the array through save and loading workspace
          if idx_rt > size(d_db,3) && forecast ~= 1
            idx_rt = 1;
            if flag_unit_forecast==1
                % Saving the data for forecast and divided in 120/24 pieces
                if dates(1) == dates_end
                    if t_d.TimerFcn{33} == 5
                        t_d.TimerFcn{33} = 1;
                    else
                        t_d.TimerFcn{33} = t_d.TimerFcn{33} +1;
                    end
                    save(strcat('data_f_',num2str(forecast),'_',num2str(t_d.TimerFcn{33})),'i_db','v_db','d_db','r_db','dates');
                    t_d.TimerFcn{34} = dates(end);
                end
                try
                    % Updating the variables to plot, this way reduce
                    % memory usage
                    if isempty(idx_loop)
                        idx_loop = forecast;
                    end
                    load(strcat('data_f_',num2str(idx_loop),'_',num2str(t_d.TimerFcn{35})),'i_db','v_db','d_db','r_db','dates');
                    t_d.TimerFcn{11} = i_db; % update all variables in t_d
                    t_d.TimerFcn{15} = v_db;
                    t_d.TimerFcn{3} = d_db;
                    t_d.TimerFcn{7} = r_db;
                    t_d.TimerFcn{25} = dates;
                    if idx_loop == 0 && time_zone > 0
                        t_d.TimerFcn{29} = strcat("Main forecast computing from~", num2str(24 - time_zone), " hrs update");
                    elseif idx_loop == 0 && time_zone < 0
                        t_d.TimerFcn{29} = strcat("Main forecast computing from~", num2str(0 - time_zone), " hrs update");
                    else
                        t_d.TimerFcn{29} = strcat("Main forecast computing from~", num2str(idx_loop - time_zone), " hrs update");
                    end

                    t_d.TimerFcn{35} = t_d.TimerFcn{35} + 1;
                    if t_d.TimerFcn{35} > 5
                        % for the next update
                        idx_loop = idx_loop + 6;
                        % retun for the first group of data
                        if idx_loop > 18
                            idx_loop = 0;
                        end
                        t_d.TimerFcn{35} = 1;
                    end
                catch
                    % returning the loop of cumulator_2 to show the first set of data
                    t_d.TimerFcn{35} = 1;
                    % for the next update
                    idx_loop = idx_loop + 6;
                    % retun for the first group of data
                    if idx_loop > 18
                        idx_loop = 0;
                    end
                    if idx_loop == 0 && time_zone > 0
                        t_d.TimerFcn{29} = strcat("Forecast from~", num2str(24 - time_zone), " hrs update, is not ready to display");
                    elseif idx_loop == 0 && time_zone < 0
                        t_d.TimerFcn{29} = strcat("Forecast from~", num2str(0 - time_zone), " hrs update, is not ready to display");
                    else
                        t_d.TimerFcn{29} = strcat("Forecast from~", num2str(idx_loop - time_zone), " hrs update, is not ready to display");
                    end
                end
            end
          elseif idx_rt > size(d_db,3) && forecast == 1
              idx_rt = 1;
          end
          drawnow
    end
end






