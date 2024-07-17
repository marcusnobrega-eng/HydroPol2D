function [ax] = HydroPol2D_running_dashboard(ax,Maps,v_t,DEM_raster,gauges,first_time,layer)
    % Runnning dashboard
    try 
        if first_time ==1
            disp('Starting HydroPol2D_running_dashboard');
            app = HydroPol2D_Monitor_base();
            app.UIFigure.Visible = 'on';
            disp('App created and set to visible');
            ax.ax_d = app.UIAxes; ax.ax_r = app.UIAxes2; ax.ax_date = app.DateTimeTextArea; ax.ax_iter = app.iter;
            ax.ax_v = app.UIAxes_2; ax.ax_system_output = app.SystemOutput;
            ax.ax_list = app.gauges_list;
            ax.app = app;
            app.flags = ax.flags;

            [Spectrum,depth_ramp,terrain_ramp] = coloramps(); % Run coloramp function
            %Store the inditial data into the multiarrays, depths, rainfall,
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
            ax.x_grid = DEM_raster.georef.SpatialRef.YWorldLimits(2) - DEM_raster.georef.SpatialRef.CellExtentInWorldX*[xbegin:1:xend];
            ax.y_grid = DEM_raster.georef.SpatialRef.XWorldLimits(1) + DEM_raster.georef.SpatialRef.CellExtentInWorldX*[ybegin:1:yend];
            ax.DEM_s1=size(DEM_raster.Z,1); ax.DEM_s2=size(DEM_raster.Z,2);
            
            % Plot or ignore according flags
            bm_d = mapshow(ax.ax_d, A, RA, "AlphaData",0.35);hold(ax.ax_d, 'on');
            bm_r = mapshow(ax.ax_r, A, RA, "AlphaData",0.35);hold(ax.ax_r, 'on');
            bm_v = mapshow(ax.ax_v, A, RA, "AlphaData",0.35);hold(ax.ax_v, 'on');
            F_d = Maps.Hydro.d([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
            if ax.flags.flag_spatial_rainfall == 1
                F_r = Maps.Hydro.spatial_rainfall_maps([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
            else
                mask = isnan(DEM_raster.Z);
                if ax.flags.flag_rainfall == 1
                    F_r = ones(ax.DEM_s1,ax.DEM_s2); % ADD Rainfall here
                    F_r(mask) = nan;
                else
                    F_r = zeros(ax.DEM_s1,ax.DEM_s2);
                    F_r(mask) = nan;
                end                
            end
            F_v = double(v_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1));
            
            if ax.flags.flag_infiltration == 1
                ax.ax_i = app.UIAxes_4;
                bm_i = mapshow(ax.ax_i, A, RA, "AlphaData",0.35);hold(ax.ax_i, 'on');
                F_i = Maps.Hydro.I_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                ax.monitor_i = pcolor(ax.ax_i,ax.y_grid,ax.x_grid,F_i);
                set(ax.monitor_i,'EdgeColor', 'none');
                shp_i = mapshow(ax.ax_i,S_p,'FaceColor','n');hold(ax.ax_i, 'on');
                set(ax.ax_i,'FontName','Garamond');
                ax.ax_i.XAxis.Exponent = 0; ax.ax_i.XAxis.TickLabelFormat = '%.0f';
                ax.ax_i.YAxis.Exponent = 0; ax.ax_i.YAxis.TickLabelFormat = '%.0f';
                ax.ax_i.YAxis.TickLabelRotation = 90;
                ax.ax_i.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
                colormap(ax.ax_i,depth_ramp);colorbar(ax.ax_i);


                ax.ax_C = app.UIAxes_3;
                bm_C = mapshow(ax.ax_C, A, RA, "AlphaData",0.35);hold(ax.ax_C, 'on');
                F_C = Maps.Hydro.C([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                ax.monitor_C = pcolor(ax.ax_C,ax.y_grid,ax.x_grid,F_C);
                set(ax.monitor_C,'EdgeColor', 'none');
                shp_C = mapshow(ax.ax_C,S_p,'FaceColor','n');hold(ax.ax_C, 'on');
                set(ax.ax_C,'FontName','Garamond');
                ax.ax_C.XAxis.Exponent = 0; ax.ax_C.XAxis.TickLabelFormat = '%.0f';
                ax.ax_C.YAxis.Exponent = 0; ax.ax_C.YAxis.TickLabelFormat = '%.0f';
                ax.ax_C.YAxis.TickLabelRotation = 90;
                ax.ax_C.YAxis.TickValues = ax.ax_C.YAxis.TickValues(1:2:end);
                colormap(ax.ax_C,depth_ramp);colorbar(ax.ax_C);

                ax.ax_f = app.UIAxes_5;
                bm_i = mapshow(ax.ax_f, A, RA, "AlphaData",0.35);hold(ax.ax_f, 'on');
                F_f = Maps.Hydro.f([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                ax.monitor_f = pcolor(ax.ax_f,ax.y_grid,ax.x_grid,F_f);
                set(ax.monitor_f,'EdgeColor', 'none');
                shp_f = mapshow(ax.ax_f,S_p,'FaceColor','n');hold(ax.ax_f, 'on');
                set(ax.ax_f,'FontName','Garamond');
                ax.ax_f.XAxis.Exponent = 0; ax.ax_f.XAxis.TickLabelFormat = '%.0f';
                ax.ax_f.YAxis.Exponent = 0; ax.ax_f.YAxis.TickLabelFormat = '%.0f';
                ax.ax_f.YAxis.TickLabelRotation = 90;
                ax.ax_f.YAxis.TickValues = ax.ax_f.YAxis.TickValues(1:2:end);
                colormap(ax.ax_f,depth_ramp);colorbar(ax.ax_f);
            end
            if ax.flags.flag_reservoir == 1
                ax.ax_bc = app.UIAxes_12;
                ax.ax_list_bc = app.boundary_list;
                bm_bc = mapshow(ax.ax_bc, A, RA, "AlphaData",0.35);hold(ax.ax_bc, 'on');
                ax.monitor_bc = pcolor(ax.ax_bc,ax.y_grid,ax.x_grid,F_d);
                set(ax.monitor_bc,'EdgeColor', 'none');
                colormap(ax.ax_bc,depth_ramp);colorbar(ax.ax_bc);
                shp_bc = mapshow(ax.ax_bc,S_p,'FaceColor','n');hold(ax.ax_bc, 'on');
                set(ax.ax_bc,'FontName','Garamond');
                ax.ax_bc.XAxis.Exponent = 0; ax.ax_bc.XAxis.TickLabelFormat = '%.0f';
                ax.ax_bc.YAxis.Exponent = 0; ax.ax_bc.YAxis.TickLabelFormat = '%.0f';
                ax.ax_bc.YAxis.TickLabelRotation = 90;
                ax.ax_bc.YAxis.TickValues = ax.ax_bc.YAxis.TickValues(1:2:end);
            end
            if ax.flags.flag_ETP == 1
                ax.ax_ETR = app.UIAxes_7;
                bm_etr = mapshow(ax.ax_ETR, A, RA, "AlphaData",0.35);hold(ax.ax_ETR, 'on');
                F_ETR = Maps.Hydro.ETR_save([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                ax.monitor_ETR = pcolor(ax.ax_ETR,ax.y_grid,ax.x_grid,F_ETR);
                set(ax.monitor_ETR,'EdgeColor', 'none');
                colormap(ax.ax_ETR,depth_ramp);colorbar(ax.ax_ETR);
                shp_ETR = mapshow(ax.ax_ETR,S_p,'FaceColor','n');hold(ax.ax_ETR, 'on');
                set(ax.ax_ETR,'FontName','Garamond');
                ax.ax_ETR.XAxis.Exponent = 0; ax.ax_ETR.XAxis.TickLabelFormat = '%.0f';
                ax.ax_ETR.YAxis.Exponent = 0; ax.ax_ETR.YAxis.TickLabelFormat = '%.0f';
                ax.ax_ETR.YAxis.TickLabelRotation = 90;
                ax.ax_ETR.YAxis.TickValues = ax.ax_ETR.YAxis.TickValues(1:2:end);
                
                ax.ax_ETP = app.UIAxes_6;
                bm_etp = mapshow(ax.ax_ETP, A, RA, "AlphaData",0.35);hold(ax.ax_ETP, 'on');
                F_ETP = Maps.Hydro.ETP_save([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                ax.monitor_ETP = pcolor(ax.ax_ETP,ax.y_grid,ax.x_grid,F_ETP);
                set(ax.monitor_ETP,'EdgeColor', 'none');
                colormap(ax.ax_ETP,depth_ramp);colorbar(ax.ax_ETP);
                shp_ETP = mapshow(ax.ax_ETP,S_p,'FaceColor','n');hold(ax.ax_ETP, 'on');
                set(ax.ax_ETP,'FontName','Garamond');
                ax.ax_ETP.XAxis.Exponent = 0; ax.ax_ETP.XAxis.TickLabelFormat = '%.0f';
                ax.ax_ETP.YAxis.Exponent = 0; ax.ax_ETP.YAxis.TickLabelFormat = '%.0f';
                ax.ax_ETP.YAxis.TickLabelRotation = 90;
                ax.ax_ETP.YAxis.TickValues = ax.ax_ETP.YAxis.TickValues(1:2:end);
            end

            F_r(F_r==0) = nan;
 
            ax.monitor_d = pcolor(ax.ax_d,ax.y_grid,ax.x_grid,F_d);
            set(ax.monitor_d,'EdgeColor', 'none');
            colormap(ax.ax_d,depth_ramp);colorbar(ax.ax_d);
            ax.monitor_r = pcolor(ax.ax_r,ax.y_grid,ax.x_grid,F_r);
            set(ax.monitor_r,'EdgeColor', 'none');
            colormap(ax.ax_r,depth_ramp);colorbar(ax.ax_r);
            ax.monitor_v = pcolor(ax.ax_v,ax.y_grid,ax.x_grid,F_v);
            set(ax.monitor_v,'EdgeColor', 'none');
            colormap(ax.ax_v,depth_ramp);colorbar(ax.ax_v);

            shp_d = mapshow(ax.ax_d,S_p,'FaceColor','n');hold(ax.ax_d, 'on');
            shp_r = mapshow(ax.ax_r,S_p,'FaceColor','n');hold(ax.ax_r, 'on');
            shp_v = mapshow(ax.ax_v,S_p,'FaceColor','n');hold(ax.ax_v, 'on');
    
            system_output = 'Initializing the system...';
            set(ax.ax_date,'Value',num2str(1));
            set(ax.ax_system_output,'Value',system_output);
            set(ax.ax_iter,'Text',num2str(1));
    
            set(ax.ax_d,'FontName','Garamond');
            ax.ax_d.XAxis.Exponent = 0; ax.ax_d.XAxis.TickLabelFormat = '%.0f';
            ax.ax_d.YAxis.Exponent = 0; ax.ax_d.YAxis.TickLabelFormat = '%.0f';
            ax.ax_d.YAxis.TickLabelRotation = 90;
            ax.ax_d.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
            set(ax.ax_r,'FontName','Garamond');
            ax.ax_r.XAxis.Exponent = 0; ax.ax_r.XAxis.TickLabelFormat = '%.0f';
            ax.ax_r.YAxis.Exponent = 0; ax.ax_r.YAxis.TickLabelFormat = '%.0f';
            ax.ax_r.YAxis.TickLabelRotation = 90;
            ax.ax_r.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
            set(ax.ax_v,'FontName','Garamond');
            ax.ax_v.XAxis.Exponent = 0; ax.ax_v.XAxis.TickLabelFormat = '%.0f';
            ax.ax_v.YAxis.Exponent = 0; ax.ax_v.YAxis.TickLabelFormat = '%.0f';
            ax.ax_v.YAxis.TickLabelRotation = 90;
            ax.ax_v.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
            % UTM Coordinates
            % Convert the string arrays to character vectors
            char_vector_cell = {};
            for i = 1:numel(gauges.label_observed_string)
                char_vector_cell = [char_vector_cell; char(gauges.label_observed_string{i})];
            end
            ax.gauges = char_vector_cell;
            % Assign the modified cell array to ax_list.Items
            ax.ax_list.Items = char_vector_cell;

            drawnow
        else
            % Update plots
            idx_g = Maps.Hydro.d([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
            idx_g(idx_g == 0) = NaN;  
            set(ax.monitor_d, 'CData', idx_g/1000);
            if ax.flags.flag_reservoir == 1
                set(ax.monitor_bc, 'CData', idx_g/1000);
            end
            if ax.flags.flag_ETP == 1
                set(ax.monitor_ETP, 'CData', Maps.Hydro.ETP_save([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
                set(ax.monitor_ETR, 'CData', Maps.Hydro.ETR_save([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
            end
            % set(ax.monitor_d, 'CData', idx_g, 'YData', ax.x_grid, 'XData', ax.y_grid)
            if ax.flags.flag_spatial_rainfall == 1
                idx_g = Maps.Hydro.spatial_rainfall_maps([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);                        
                idx_g(idx_g == 0) = NaN;
            else
                mask = isnan(DEM_raster.Z);
                if ax.flags.flag_rainfall == 1
                    idx_g = ones(ax.DEM_s1,ax.DEM_s2); % ADD Rainfall here
                    idx_g(mask) = nan;
                    idx_g(idx_g == 0) = nan;
                else
                    idx_g = zeros(ax.DEM_s1,ax.DEM_s2); % ADD Rainfall here
                    idx_g(mask) = nan;
                    idx_g(idx_g == 0) = nan;
                end      
            end
            set(ax.monitor_r, 'CData', idx_g);
            % set(ax.monitor_r, 'CData', idx_g, 'YData', ax.x_grid, 'XData', ax.y_grid)
            if ax.flags.flag_infiltration == 1
                set(ax.monitor_i, 'CData', Maps.Hydro.I_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
                set(ax.monitor_C, 'CData', Maps.Hydro.C([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
                idx_g = Maps.Hydro.f([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
                idx_g(idx_g == 0) = NaN;
                set(ax.monitor_f, 'CData', idx_g);
                idx_g = Maps.Hydro.I_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
                idx_g(idx_g == 0) = NaN;  
                set(ax.monitor_i, 'CData', idx_g);
            end

            
            % set(ax.monitor_i, 'CData', idx_g, 'YData', ax.x_grid, 'XData', ax.y_grid)
            v_t(v_t == 0) = NaN;
            set(ax.monitor_v, 'CData', v_t);
            % set(ax.monitor_v, 'CData', v_t, 'YData', ax.x_grid, 'XData', ax.y_grid)


            set(ax.ax_date,'Value',datestr(ax.timer));
            set(ax.ax_system_output,'Value',strcat('Main model execution at', 32,num2str(ax.percentage),'%'));
            set(ax.ax_iter, 'Text', num2str(layer));

            drawnow
            % pause(0.5)
        end
    catch e
        disp('Error occurred');
        disp(e.message);
        disp(e.stack(1));
    end
end


