function [ax] = HydroPol2D_running_dashboard(ax,Maps,v_t,DEM_raster,gauges,BC_States,time_step,Resolution,first_time,layer)
    % Runnning dashboard
    mask = isnan(DEM_raster.Z);
    
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

            [Spectrum,depth_ramp,terrain_ramp,blue_ramp,Depths_RAS,pallete,Depth_RAS,Terrain_RAS,Velocity_RAS,WSE_RAS] = coloramps();        
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
            mask = isnan(DEM_raster.Z);            
            F_d(mask) = nan;
            if ax.flags.flag_spatial_rainfall == 1
                F_r = Maps.Hydro.spatial_rainfall_maps([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
            else
                if ax.flags.flag_rainfall == 1
                    F_r = BC_States.delta_p_agg/(time_step/60)*ones(ax.DEM_s1,ax.DEM_s2); % ADD Rainfall here
                    F_r(mask) = nan;               
                else
                    F_r = zeros(ax.DEM_s1,ax.DEM_s2);
                    F_r(mask) = nan;
                end                
            end
            F_v = double(v_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1));
            F_v(mask) = nan;
            if ax.flags.flag_infiltration == 1
                ax.ax_i = app.UIAxes_4;
                bm_i = mapshow(ax.ax_i, A, RA, "AlphaData",0.35);hold(ax.ax_i, 'on');
                F_i = Maps.Hydro.I_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_i(mask) = nan;
                ax.monitor_i = pcolor(ax.ax_i,ax.y_grid,ax.x_grid,F_i);
                set(ax.monitor_i,'EdgeColor', 'none');
                shp_i = mapshow(ax.ax_i,S_p,'FaceColor','n');hold(ax.ax_i, 'on');
                set(ax.ax_i,'FontName','Garamond');
                ax.ax_i.XAxis.Exponent = 0; ax.ax_i.XAxis.TickLabelFormat = '%.0f';
                ax.ax_i.YAxis.Exponent = 0; ax.ax_i.YAxis.TickLabelFormat = '%.0f';
                ax.ax_i.YAxis.TickLabelRotation = 90;
                ax.ax_i.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
                colormap(ax.ax_i,WSE_RAS); hh = colorbar(ax.ax_i);  hh.TickDirection = 'out';

                ax.ax_C = app.UIAxes_3;
                bm_C = mapshow(ax.ax_C, A, RA, "AlphaData",0.35);hold(ax.ax_C, 'on');
                F_C = Maps.Hydro.C([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_C(mask) = nan;
                ax.monitor_C = pcolor(ax.ax_C,ax.y_grid,ax.x_grid,F_C);
                set(ax.monitor_C,'EdgeColor', 'none');
                shp_C = mapshow(ax.ax_C,S_p,'FaceColor','n');hold(ax.ax_C, 'on');
                set(ax.ax_C,'FontName','Garamond');
                ax.ax_C.XAxis.Exponent = 0; ax.ax_C.XAxis.TickLabelFormat = '%.0f';
                ax.ax_C.YAxis.Exponent = 0; ax.ax_C.YAxis.TickLabelFormat = '%.0f';
                ax.ax_C.YAxis.TickLabelRotation = 90;
                ax.ax_C.YAxis.TickValues = ax.ax_C.YAxis.TickValues(1:2:end);
                colormap(ax.ax_C,Spectrum);hh = colorbar(ax.ax_C);  hh.TickDirection = 'out';


                ax.ax_f = app.UIAxes_5;
                bm_i = mapshow(ax.ax_f, A, RA, "AlphaData",0.35);hold(ax.ax_f, 'on');
                F_f = Maps.Hydro.f([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_f(mask) = nan;
                ax.monitor_f = pcolor(ax.ax_f,ax.y_grid,ax.x_grid,F_f);
                set(ax.monitor_f,'EdgeColor', 'none');
                shp_f = mapshow(ax.ax_f,S_p,'FaceColor','n');hold(ax.ax_f, 'on');
                set(ax.ax_f,'FontName','Garamond');
                ax.ax_f.XAxis.Exponent = 0; ax.ax_f.XAxis.TickLabelFormat = '%.0f';
                ax.ax_f.YAxis.Exponent = 0; ax.ax_f.YAxis.TickLabelFormat = '%.0f';
                ax.ax_f.YAxis.TickLabelRotation = 90;
                ax.ax_f.YAxis.TickValues = ax.ax_f.YAxis.TickValues(1:2:end);
                colormap(ax.ax_f,Velocity_RAS);hh = colorbar(ax.ax_f);  hh.TickDirection = 'out';

            end
            if ax.flags.flag_reservoir == 1
                ax.ax_bc = app.UIAxes_12;
                ax.ax_list_bc = app.boundary_list;
                bm_bc = mapshow(ax.ax_bc, A, RA, "AlphaData",0.35);hold(ax.ax_bc, 'on');
                ax.monitor_bc = pcolor(ax.ax_bc,ax.y_grid,ax.x_grid,F_d);
                set(ax.monitor_bc,'EdgeColor', 'none');
                colormap(ax.ax_bc,Depths_RAS);hh = colorbar(ax.ax_bc);  hh.TickDirection = 'out';
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
                F_ETR(mask) = nan;
                ax.monitor_ETR = pcolor(ax.ax_ETR,ax.y_grid,ax.x_grid,F_ETR);
                set(ax.monitor_ETR,'EdgeColor', 'none');
                colormap(ax.ax_ETR,Spectrum);hh = colorbar(ax.ax_ETR);  hh.TickDirection = 'out';
                shp_ETR = mapshow(ax.ax_ETR,S_p,'FaceColor','n');hold(ax.ax_ETR, 'on');
                set(ax.ax_ETR,'FontName','Garamond');
                ax.ax_ETR.XAxis.Exponent = 0; ax.ax_ETR.XAxis.TickLabelFormat = '%.0f';
                ax.ax_ETR.YAxis.Exponent = 0; ax.ax_ETR.YAxis.TickLabelFormat = '%.0f';
                ax.ax_ETR.YAxis.TickLabelRotation = 90;
                ax.ax_ETR.YAxis.TickValues = ax.ax_ETR.YAxis.TickValues(1:2:end);               

                
                ax.ax_ETP = app.UIAxes_6;
                bm_etp = mapshow(ax.ax_ETP, A, RA, "AlphaData",0.35);hold(ax.ax_ETP, 'on');
                F_ETP = Maps.Hydro.ETP_save([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_ETP(mask) = nan;
                ax.monitor_ETP = pcolor(ax.ax_ETP,ax.y_grid,ax.x_grid,F_ETP);
                set(ax.monitor_ETP,'EdgeColor', 'none');
                colormap(ax.ax_ETP,Spectrum);hh = colorbar(ax.ax_ETP);  hh.TickDirection = 'out';
                shp_ETP = mapshow(ax.ax_ETP,S_p,'FaceColor','n');hold(ax.ax_ETP, 'on');
                set(ax.ax_ETP,'FontName','Garamond');
                ax.ax_ETP.XAxis.Exponent = 0; ax.ax_ETP.XAxis.TickLabelFormat = '%.0f';
                ax.ax_ETP.YAxis.Exponent = 0; ax.ax_ETP.YAxis.TickLabelFormat = '%.0f';
                ax.ax_ETP.YAxis.TickLabelRotation = 90;
                ax.ax_ETP.YAxis.TickValues = ax.ax_ETP.YAxis.TickValues(1:2:end);                 
            end

            % Water Quality
            if ax.flags.flag_waterquality == 1
                % Concentration
                ax.ax_Pol_Conc = app.UIAxes_8;
                bm_ax_Pol_Conc = mapshow(ax.ax_Pol_Conc, A, RA, "AlphaData",0.35);hold(ax.ax_Pol_Conc, 'on');
                F_Conc = Maps.WQ_States.Pol_Conc_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_Conc(mask) = nan;
                ax.monitor_Pol_Conc = pcolor(ax.ax_Pol_Conc,ax.y_grid,ax.x_grid,F_Conc);
                set(ax.monitor_Pol_Conc,'EdgeColor', 'none');
                colormap(ax.ax_Pol_Conc,Spectrum);hh = colorbar(ax.ax_Pol_Conc);  hh.TickDirection = 'out';
                shp_Pol_Conc = mapshow(ax.ax_Pol_Conc,S_p,'FaceColor','n');hold(ax.ax_Pol_Conc, 'on');
                set(ax.ax_Pol_Conc,'FontName','Garamond');
                ax.ax_Pol_Conc.XAxis.Exponent = 0; ax.ax_Pol_Conc.XAxis.TickLabelFormat = '%.0f';
                ax.ax_Pol_Conc.YAxis.Exponent = 0; ax.ax_Pol_Conc.YAxis.TickLabelFormat = '%.0f';
                ax.ax_Pol_Conc.YAxis.TickLabelRotation = 90;
                ax.ax_Pol_Conc.YAxis.TickValues = ax.ax_Pol_Conc.YAxis.TickValues(1:2:end);   

                % Load
                ax.ax_Load = app.UIAxes_9;
                bm_ax_Load = mapshow(ax.ax_Load, A, RA, "AlphaData",0.35);hold(ax.ax_Load, 'on');
                F_Load = Maps.WQ_States.Pol_Load_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_Load(mask) = nan;
                ax.monitor_Load = pcolor(ax.ax_Load,ax.y_grid,ax.x_grid,F_Load);
                set(ax.monitor_Load,'EdgeColor', 'none');
                colormap(ax.ax_Load,Spectrum);hh = colorbar(ax.ax_Load);  hh.TickDirection = 'out';
                shp_Load = mapshow(ax.ax_Load,S_p,'FaceColor','n');hold(ax.ax_Load, 'on');
                set(ax.ax_Load,'FontName','Garamond');
                ax.ax_Load.XAxis.Exponent = 0; ax.ax_Load.XAxis.TickLabelFormat = '%.0f';
                ax.ax_Load.YAxis.Exponent = 0; ax.ax_Load.YAxis.TickLabelFormat = '%.0f';
                ax.ax_Load.YAxis.TickLabelRotation = 90;
                ax.ax_Load.YAxis.TickValues = ax.ax_Load.YAxis.TickValues(1:2:end);      

                % Build-up
                ax.ax_Buildup = app.UIAxes_11;
                bm_ax_Buildup = mapshow(ax.ax_Buildup, A, RA, "AlphaData",0.35);hold(ax.ax_Buildup, 'on');
                F_Buildup = 1000*1/(Resolution^2)*Maps.WQ_States.Pol_Mass_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],1);
                F_Buildup(mask) = nan;
                ax.monitor_Buildup = pcolor(ax.ax_Buildup,ax.y_grid,ax.x_grid,F_Buildup);
                set(ax.monitor_Buildup,'EdgeColor', 'none');
                colormap(ax.ax_Buildup,Spectrum);hh = colorbar(ax.ax_Buildup);  hh.TickDirection = 'out';
                shp_Buildup = mapshow(ax.ax_Buildup,S_p,'FaceColor','n');hold(ax.ax_Buildup, 'on');
                set(ax.ax_Buildup,'FontName','Garamond');
                ax.ax_Buildup.XAxis.Exponent = 0; ax.ax_Buildup.XAxis.TickLabelFormat = '%.0f';
                ax.ax_Buildup.YAxis.Exponent = 0; ax.ax_Buildup.YAxis.TickLabelFormat = '%.0f';
                ax.ax_Buildup.YAxis.TickLabelRotation = 90;
                ax.ax_Buildup.YAxis.TickValues = ax.ax_Buildup.YAxis.TickValues(1:2:end);                  
            end

            ax.monitor_d = pcolor(ax.ax_d,ax.y_grid,ax.x_grid,F_d);
            set(ax.monitor_d,'EdgeColor', 'none');
            colormap(ax.ax_d,Depths_RAS); hh = colorbar(ax.ax_d); hh.TickDirection = "out";  hh.TickDirection = 'out';
            ax.monitor_r = pcolor(ax.ax_r,ax.y_grid,ax.x_grid,F_r); 
            set(ax.monitor_r,'EdgeColor', 'none'); 
            colormap(ax.ax_r,Spectrum);hh = colorbar(ax.ax_r);   hh.TickDirection = 'out';
            ax.monitor_v = pcolor(ax.ax_v,ax.y_grid,ax.x_grid,F_v); 
            set(ax.monitor_v,'EdgeColor', 'none');
            colormap(ax.ax_v,depth_ramp);hh = colorbar(ax.ax_v); 

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
            try
                for i = 1:numel(gauges.labels_observed_string)
                    char_vector_cell = [char_vector_cell; char(gauges.labels_observed_string{i})];
                end
            end
            ax.gauges = char_vector_cell;
            % Assign the modified cell array to ax_list.Items
            ax.ax_list.Items = char_vector_cell;               

            drawnow
        else
            ax.DEM_s1=size(DEM_raster.Z,1); ax.DEM_s2=size(DEM_raster.Z,2);
            % Update plots
            idx_g = Maps.Hydro.d([1:1:size(DEM_raster.Z,1)],[1:1:size(DEM_raster.Z,2)],layer);
            idx_g(idx_g == 0) = NaN;  
            idx_g(mask) = nan;
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
                if ax.flags.flag_rainfall == 1
                    idx_g = BC_States.delta_p_agg*ones(ax.DEM_s1,ax.DEM_s2); % ADD Rainfall here
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
                idx_g(idx_g == 0) = nan;
                idx_g(mask) = nan;
                set(ax.monitor_f, 'CData', idx_g);
                idx_g = Maps.Hydro.I_t([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
                idx_g(idx_g == 0) = NaN;  
                idx_g(mask) = nan;
                set(ax.monitor_i, 'CData', idx_g);
            end         

            % Water Quality
            if ax.flags.flag_waterquality == 1
                % Concentration                
                set(ax.monitor_Pol_Conc, 'CData', Maps.WQ_States.Pol_Conc_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
                idx_g = Maps.WQ_States.Pol_Conc_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
                idx_g(idx_g == 0) = nan;
                idx_g(mask) = nan;
                set(ax.monitor_Pol_Conc, 'CData', idx_g);                

                % Load

                set(ax.monitor_Load, 'CData', Maps.WQ_States.Pol_Load_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
                idx_g = Maps.WQ_States.Pol_Load_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
                idx_g(idx_g == 0) = nan;
                idx_g(mask) = nan;
                set(ax.monitor_Load, 'CData', idx_g);                

                % Build-up         
                set(ax.monitor_Buildup, 'CData', 1000*1/(Resolution^2)*Maps.WQ_States.Pol_Mass_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer));
                idx_g = 1000*(1/Resolution^2)*Maps.WQ_States.Pol_Mass_Map([1:1:ax.DEM_s1],[1:1:ax.DEM_s2],layer);
                idx_g(idx_g == 0) = nan;
                idx_g(mask) = nan;
                set(ax.monitor_Buildup, 'CData', idx_g);               
            end


            
            % set(ax.monitor_i, 'CData', idx_g, 'YData', ax.x_grid, 'XData', ax.y_grid)
            v_t(v_t == 0) = NaN;
            set(ax.monitor_v, 'CData', v_t);
            % set(ax.monitor_v, 'CData', v_t, 'YData', ax.x_grid, 'XData', ax.y_grid)


            set(ax.ax_date,'Value',datestr(ax.timer));
            set(ax.ax_system_output,'Value',strcat('Main model execution at', 32,num2str(ax.percentage),'%'));
            set(ax.ax_iter, 'Text', num2str(layer));

            drawnow
            pause(0.1)
        end
    catch e
        disp('Error occurred');
        disp(e.message);
        disp(e.stack(1));
    end
end




