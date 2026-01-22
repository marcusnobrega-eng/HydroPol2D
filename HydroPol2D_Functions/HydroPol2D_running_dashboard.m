function [ax] = HydroPol2D_running_dashboard(ax, Maps, v_t, DEM_raster, gauges, BC_States, time_step, Resolution, first_time, layer, C_a)
% HydroPol2D_RUNNING_DASHBOARD sets up or updates the dashboard for HydroPol2D.
%
%   Inputs:
%       ax          - structure with axes and other GUI handles
%       Maps        - structure containing hydrological maps data
%       v_t         - velocity matrix (time-dependent)
%       DEM_raster  - structure containing DEM data and spatial reference
%       gauges      - structure with gauge information
%       BC_States   - structure with boundary condition states (e.g., rainfall)
%       time_step   - simulation time step (in seconds)
%       Resolution  - spatial resolution used for scaling some maps
%       first_time  - flag indicating if this is the first run (1 for initialization)
%       layer       - current simulation layer or time index
%       C_a         - additional scaling constant (for rainfall)
%
%   Output:
%       ax          - updated structure with axes handles and monitor objects

% Create mask for invalid DEM points
mask = isnan(DEM_raster.Z);

web_mercator_crs = projcrs(3857);
try
    DEM_CRS = DEM_raster.georef.SpatialRef.ProjectedCRS;
    isWebMercator = isequal(DEM_CRS, web_mercator_crs);
catch
    warning("DEM_raster does not define a ProjectedCRS. Assuming non-Web Mercator.");
    isWebMercator = false;
end

% Load color ramps and colormaps for later use
[Spectrum, depth_ramp, terrain_ramp, blue_ramp, blues_2, pallete, ...
    Depth_RAS, Terrain_RAS, Velocity_RAS, WSE_RAS] = coloramps();
try
    if first_time == 1
        %% INITIALIZATION SECTION
        disp('Starting HydroPol2D_running_dashboard');
        % Create and show the dashboard application
        app = HydroPol2D_Monitor_base();
        app.UIFigure.Visible = 'on';
        disp('App created and set to visible');

        % Assign GUI handles to ax structure
        ax.ax_d = app.UIAxes;
        ax.ax_r = app.UIAxes2;
        ax.ax_date = app.DateTimeTextArea;
        ax.ax_iter = app.iter;
        ax.ax_v = app.UIAxes_2;
        ax.ax_system_output = app.SystemOutput;
        ax.ax_list = app.gauges_list;
        ax.app = app;
        app.flags = ax.flags; % pass flag settings to the app


        % Newer GUI handles
        ax.ax_GW = app.GWdepth;
        ax.ax_C = app.UIAxes_3;
        ax.ax_ETR = app.UIAxes_7;
        ax.ax_ETP = app.UIAxes_6;
        ax.ax_Abstraction = app.UIAxes_13;
        ax.ax_Snowpack = app.UIAxes_13;
        ax.ax_C = app.UIAxes_3;
        ax.ax_i = app.UIAxes_4;
        ax.ax_f = app.UIAxes_5;

        %% Creating a custom basemap
        basemapName = "openstreetmap";
        url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png";
        copyright = char(uint8(169));
        attribution = copyright + "OpenStreetMap contributors";
        addCustomBasemap(basemapName, url, "Attribution", attribution);

        % Get latitude and longitude limits from DEM spatial reference
        if isWebMercator
            [lat, lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS, ...
                DEM_raster.georef.SpatialRef.XWorldLimits, ...
                DEM_raster.georef.SpatialRef.YWorldLimits);
            latlim = [lat(1), lat(2)];
            lonlim = [lon(1), lon(2)];

            try
                [A, RA, attribA] = readBasemapImage(basemapName, latlim, lonlim);
            catch ME
                warning("Failed to load basemap")
                A = []; RA = [];
            end
        else
            A = []; RA = [];
            warning("DEM is not in Web Mercator. Basemap will be skipped.");
        end


        %% Create a reference shapefile from DEM
        % Create a binary mask and extract boundaries
        binaryMask = ~isnan(DEM_raster.Z);
        boundaries = bwboundaries(binaryMask);
        combinedX = [];
        combinedY = [];
        for k = 1:numel(boundaries)
            boundary = boundaries{k};
            X = boundary(:, 2);
            Y = boundary(:, 1);
            % Concatenate boundaries with NaN separators for plotting
            combinedX = [combinedX; X; NaN];
            combinedY = [combinedY; Y; NaN];
        end
        % Remove trailing NaN (optional)
        combinedX = combinedX(1:end-1);
        combinedY = combinedY(1:end-1);

        % Check CRS and build geostruct for shapefile if using Web Mercator
        web_mercator_crs = projcrs(3857);
        no_plot = 0;

        % Shapefile Strucutre
        S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
        S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); ...
            DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)];
        % Convert pixel indices to world coordinates
        S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1) + ...
            combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
        S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2) - ...
            combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';



        %% Setting up the spatial grid for monitors
        % Define grid coordinates based on DEM dimensions
        ax.DEM_s1 = size(DEM_raster.Z, 1);
        ax.DEM_s2 = size(DEM_raster.Z, 2);
        xbegin = 1;
        ybegin = 1;
        xend = ax.DEM_s1;
        yend = ax.DEM_s2;
        ax.x_grid = DEM_raster.georef.SpatialRef.YWorldLimits(2) - ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX * (xbegin:xend);
        ax.y_grid = DEM_raster.georef.SpatialRef.XWorldLimits(1) + ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX * (ybegin:yend);

        %% Plot Initialization on various monitors
        % Plot basemap on depth, rainfall, and velocity axes
        if isWebMercator && ~isempty(A) && ~isempty(RA)
            mapshow(ax.ax_d, A, RA, "AlphaData", 0.35); hold(ax.ax_d, 'on');
            mapshow(ax.ax_r, A, RA, "AlphaData", 0.35); hold(ax.ax_r, 'on');
            mapshow(ax.ax_v, A, RA, "AlphaData", 0.35); hold(ax.ax_v, 'on');
            mapshow(ax.ax_GW, A, RA, "AlphaData", 0.35); hold(ax.ax_GW, 'on');
            mapshow(ax.ax_C, A, RA, "AlphaData", 0.35); hold(ax.ax_C, 'on');
            mapshow(ax.ax_ETR, A, RA, "AlphaData", 0.35); hold(ax.ax_ETR, 'on');
            mapshow(ax.ax_ETP, A, RA, "AlphaData", 0.35); hold(ax.ax_ETP, 'on');
            mapshow(ax.ax_Abstraction, A, RA, "AlphaData", 0.35); hold(ax.ax_Abstraction, 'on');
            mapshow(ax.ax_Snowpack, A, RA, "AlphaData", 0.35); hold(ax.ax_Snowpack, 'on');
            mapshow(ax.ax_i, A, RA, "AlphaData", 0.35); hold(ax.ax_i, 'on');
            mapshow(ax.ax_f, A, RA, "AlphaData", 0.35); hold(ax.ax_f, 'on');
        end
        % Depth field from Maps.Hydro.d; filter out NaNs using mask
        F_d = Maps.Hydro.d(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
        F_d(mask) = nan;

        % Rainfall field: choose between spatial rainfall map and uniform BC_States rainfall
        if ax.flags.flag_spatial_rainfall == 1
            F_r = Maps.Hydro.spatial_rainfall_maps(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_r(mask) = nan;
        else
            if ax.flags.flag_rainfall == 1
                F_r = (BC_States.delta_p_agg / (time_step/60)) * ones(ax.DEM_s1, ax.DEM_s2);
                F_r(mask) = nan;
            else
                F_r = zeros(ax.DEM_s1, ax.DEM_s2);
                F_r(mask) = nan;
            end
        end

        % Velocity field (converted to double); filter out zeros using mask
        F_v = double(v_t(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
        F_v(mask) = nan;

        % Plot infiltration (and related fields) if flag is set
        if ax.flags.flag_infiltration == 1
            % Infiltration plot on UIAxes_4
            ax.ax_i = app.UIAxes_4;
            bm_i = mapshow(ax.ax_i, A, RA, "AlphaData", 0.35); hold(ax.ax_i, 'on');
            F_i = Maps.Hydro.I_t(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_i(mask) = nan;
            ax.monitor_i = pcolor(ax.ax_i, ax.y_grid, ax.x_grid, F_i);
            set(ax.monitor_i, 'EdgeColor', 'none');
            % Overlay shapefile reference
            shp_i = mapshow(ax.ax_i, S_p, 'FaceColor', 'n'); hold(ax.ax_i, 'on');
            % Format axis
            set(ax.ax_i, 'FontName', 'Montserrat');
            ax.ax_i.XAxis.Exponent = 0; ax.ax_i.XAxis.TickLabelFormat = '%.0f';
            ax.ax_i.YAxis.Exponent = 0; ax.ax_i.YAxis.TickLabelFormat = '%.0f';
            ax.ax_i.YAxis.TickLabelRotation = 90;
            ax.ax_i.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
            colormap(ax.ax_i, WSE_RAS);
            % colorbar(ax.ax_i, 'TickDirection', 'out');
            style_colorbar(ax.ax_i)

            % Similarly, plot field 'C' on UIAxes_3 and field 'f' on UIAxes_5
            ax.ax_C = app.UIAxes_3;
            bm_C = mapshow(ax.ax_C, A, RA, "AlphaData", 0.35); hold(ax.ax_C, 'on');
            F_C = Maps.Hydro.C(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_C(mask) = nan;
            ax.monitor_C = pcolor(ax.ax_C, ax.y_grid, ax.x_grid, F_C);
            set(ax.monitor_C, 'EdgeColor', 'none');
            shp_C = mapshow(ax.ax_C, S_p, 'FaceColor', 'n'); hold(ax.ax_C, 'on');
            set(ax.ax_C, 'FontName', 'Montserrat');
            ax.ax_C.XAxis.Exponent = 0; ax.ax_C.XAxis.TickLabelFormat = '%.0f';
            ax.ax_C.YAxis.Exponent = 0; ax.ax_C.YAxis.TickLabelFormat = '%.0f';
            ax.ax_C.YAxis.TickLabelRotation = 90;
            ax.ax_C.YAxis.TickValues = ax.ax_C.YAxis.TickValues(1:2:end);
            colormap(ax.ax_C, Spectrum);
            % colorbar(ax.ax_C, 'TickDirection', 'out');
            style_colorbar(ax.ax_C);

            ax.ax_f = app.UIAxes_5;
            bm_i = mapshow(ax.ax_f, A, RA, "AlphaData", 0.35); hold(ax.ax_f, 'on');
            F_f = Maps.Hydro.f(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_f(mask) = nan;
            ax.monitor_f = pcolor(ax.ax_f, ax.y_grid, ax.x_grid, F_f);
            set(ax.monitor_f, 'EdgeColor', 'none');
            shp_f = mapshow(ax.ax_f, S_p, 'FaceColor', 'n'); hold(ax.ax_f, 'on');
            set(ax.ax_f, 'FontName', 'Montserrat');
            ax.ax_f.XAxis.Exponent = 0; ax.ax_f.XAxis.TickLabelFormat = '%.0f';
            ax.ax_f.YAxis.Exponent = 0; ax.ax_f.YAxis.TickLabelFormat = '%.0f';
            ax.ax_f.YAxis.TickLabelRotation = 90;
            ax.ax_f.YAxis.TickValues = ax.ax_f.YAxis.TickValues(1:2:end);
            colormap(ax.ax_f, Velocity_RAS);
            % colorbar(ax.ax_f, 'TickDirection', 'out');
            style_colorbar(ax.ax_f)
        end

        % Plot baseflow if flag set
        if ax.flags.flag_baseflow == 1
            ax.ax_GW = app.GWdepth;
            bm_i = mapshow(ax.ax_GW, A, RA, "AlphaData", 0.35); hold(ax.ax_GW, 'on');
            GW_depth = Maps.Hydro.GWdepth_save(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            GW_depth(mask) = nan;
            ax.monitor_GW = pcolor(ax.ax_GW, ax.y_grid, ax.x_grid, GW_depth);
            set(ax.monitor_GW, 'EdgeColor', 'none');
            shp_f = mapshow(ax.ax_GW, S_p, 'FaceColor', 'n'); hold(ax.ax_GW, 'on');
            set(ax.ax_GW, 'FontName', 'Montserrat');
            ax.ax_GW.XAxis.Exponent = 0; ax.ax_GW.XAxis.TickLabelFormat = '%.0f';
            ax.ax_GW.YAxis.Exponent = 0; ax.ax_GW.YAxis.TickLabelFormat = '%.0f';
            ax.ax_GW.YAxis.TickLabelRotation = 90;
            ax.ax_GW.YAxis.TickValues = ax.ax_GW.YAxis.TickValues(1:2:end);
            colormap(ax.ax_GW, Velocity_RAS);
            style_colorbar(ax.ax_GW);
            % colorbar(ax.ax_GW, 'TickDirection', 'out');
        end

        % Plot reservoir if flag set
        if ax.flags.flag_reservoir == 1
            ax.ax_bc = app.UIAxes_12;
            ax.ax_list_bc = app.boundary_list;
            bm_bc = mapshow(ax.ax_bc, A, RA, "AlphaData", 0.35); hold(ax.ax_bc, 'on');
            ax.monitor_bc = pcolor(ax.ax_bc, ax.y_grid, ax.x_grid, F_d);
            set(ax.monitor_bc, 'EdgeColor', 'none');
            colormap(ax.ax_bc, Velocity_RAS);
            % colorbar(ax.ax_bc, 'TickDirection', 'out');
            style_colorbar(ax.ax_bc);
            shp_bc = mapshow(ax.ax_bc, S_p, 'FaceColor', 'n'); hold(ax.ax_bc, 'on');
            set(ax.ax_bc, 'FontName', 'Montserrat');
            ax.ax_bc.XAxis.Exponent = 0; ax.ax_bc.XAxis.TickLabelFormat = '%.0f';
            ax.ax_bc.YAxis.Exponent = 0; ax.ax_bc.YAxis.TickLabelFormat = '%.0f';
            ax.ax_bc.YAxis.TickLabelRotation = 90;
            ax.ax_bc.YAxis.TickValues = ax.ax_bc.YAxis.TickValues(1:2:end);
        end

        % Plot evaporation (ETP/ETR) if flag set
        if ax.flags.flag_ETP == 1
            % ETR Plot on UIAxes_7
            ax.ax_ETR = app.UIAxes_7;
            bm_etr = mapshow(ax.ax_ETR, A, RA, "AlphaData", 0.35); hold(ax.ax_ETR, 'on');
            F_ETR = Maps.Hydro.ETR_save(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_ETR(mask) = nan; F_ETR(F_ETR == 0) = nan;
            ax.monitor_ETR = pcolor(ax.ax_ETR, ax.y_grid, ax.x_grid, F_ETR);
            set(ax.monitor_ETR, 'EdgeColor', 'none');
            colormap(ax.ax_ETR, Spectrum);
            style_colorbar(ax.ax_ETR);
            % colorbar(ax.ax_ETR, 'TickDirection', 'out');
            shp_ETR = mapshow(ax.ax_ETR, S_p, 'FaceColor', 'n'); hold(ax.ax_ETR, 'on');
            set(ax.ax_ETR, 'FontName', 'Montserrat');
            ax.ax_ETR.XAxis.Exponent = 0; ax.ax_ETR.XAxis.TickLabelFormat = '%.0f';
            ax.ax_ETR.YAxis.Exponent = 0; ax.ax_ETR.YAxis.TickLabelFormat = '%.0f';
            ax.ax_ETR.YAxis.TickLabelRotation = 90;
            ax.ax_ETR.YAxis.TickValues = ax.ax_ETR.YAxis.TickValues(1:2:end);

            % ETP Plot on UIAxes_6
            ax.ax_ETP = app.UIAxes_6;
            bm_etp = mapshow(ax.ax_ETP, A, RA, "AlphaData", 0.35); hold(ax.ax_ETP, 'on');
            F_ETP = Maps.Hydro.ETP_save(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_ETP(mask) = nan;
            ax.monitor_ETP = pcolor(ax.ax_ETP, ax.y_grid, ax.x_grid, F_ETP);
            set(ax.monitor_ETP, 'EdgeColor', 'none');
            colormap(ax.ax_ETP, Spectrum);
            style_colorbar(ax.ax_ETP)
            % colorbar(ax.ax_ETP, 'TickDirection', 'out');
            shp_ETP = mapshow(ax.ax_ETP, S_p, 'FaceColor', 'n'); hold(ax.ax_ETP, 'on');
            set(ax.ax_ETP, 'FontName', 'Montserrat');
            ax.ax_ETP.XAxis.Exponent = 0; ax.ax_ETP.XAxis.TickLabelFormat = '%.0f';
            ax.ax_ETP.YAxis.Exponent = 0; ax.ax_ETP.YAxis.TickLabelFormat = '%.0f';
            ax.ax_ETP.YAxis.TickLabelRotation = 90;
            ax.ax_ETP.YAxis.TickValues = ax.ax_ETP.YAxis.TickValues(1:2:end);
        end

        % Interception
        if ax.flags.flag_abstraction == 1
            % Abstraction on UIAxes_13
            ax.ax_Abstraction = app.UIAxes_13;
            bm_Abstraction = mapshow(ax.ax_Abstraction, A, RA, "AlphaData", 0.35); hold(ax.ax_Abstraction, 'on');
            F_Abstraction = Maps.Hydro.Abstraction(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_Abstraction(mask) = nan;
            F_Abstraction(F_Abstraction <= 1e-3) = nan;
            ax.monitor_Abstraction = pcolor(ax.ax_Abstraction, ax.y_grid, ax.x_grid, F_Abstraction);
            set(ax.monitor_Abstraction, 'EdgeColor', 'none');
            colormap(ax.ax_Abstraction, Depth_RAS);
            style_colorbar(ax.ax_Abstraction);
            % colorbar(ax.ax_Abstraction, 'TickDirection', 'out');
            shp_Abstraction = mapshow(ax.ax_Abstraction, S_p, 'FaceColor', 'n'); hold(ax.ax_Abstraction, 'on');
            set(ax.ax_ETP, 'FontName', 'Montserrat');
            ax.ax_Abstraction.XAxis.Exponent = 0; ax.ax_Abstraction.XAxis.TickLabelFormat = '%.0f';
            ax.ax_Abstraction.YAxis.Exponent = 0; ax.ax_Abstraction.YAxis.TickLabelFormat = '%.0f';
            ax.ax_Abstraction.YAxis.TickLabelRotation = 90;
            ax.ax_Abstraction.YAxis.TickValues = ax.ax_Abstraction.YAxis.TickValues(1:2:end);
        end

        if ax.flags.flag_snow_modeling == 1
            % Snowpack on UIAxes_14
            ax.ax_Snowpack = app.UIAxes_14;
            bm_Snowpack = mapshow(ax.ax_Snowpack, A, RA, "AlphaData", 0.35); hold(ax.ax_Snowpack, 'on');
            F_Snowpack = Maps.Hydro.Snowpack(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_Snowpack(mask) = nan;
            F_Snowpack(F_Snowpack == 0) = nan;
            ax.monitor_Snowpack = pcolor(ax.ax_Snowpack, ax.y_grid, ax.x_grid, F_Snowpack);
            set(ax.monitor_Snowpack, 'EdgeColor', 'none');
            colormap(ax.ax_Snowpack, Spectrum);
            % colorbar(ax.ax_Snowpack, 'TickDirection', 'out');
            style_colorbar(ax.ax_Snowpack);
            shp_Snowpack = mapshow(ax.ax_Snowpack, S_p, 'FaceColor', 'n'); hold(ax.ax_Snowpack, 'on');
            set(ax.ax_ETP, 'FontName', 'Montserrat');
            ax.ax_Snowpack.XAxis.Exponent = 0; ax.ax_Snowpack.XAxis.TickLabelFormat = '%.0f';
            ax.ax_Snowpack.YAxis.Exponent = 0; ax.ax_Snowpack.YAxis.TickLabelFormat = '%.0f';
            ax.ax_Snowpack.YAxis.TickLabelRotation = 90;
            ax.ax_Snowpack.YAxis.TickValues = ax.ax_Snowpack.YAxis.TickValues(1:2:end);
        end

        %  Quality Plots
        if ax.flags.flag_waterquality == 1
            % Concentration Map on UIAxes_8
            ax.ax_Pol_Conc = app.UIAxes_8;
            bm_ax_Pol_Conc = mapshow(ax.ax_Pol_Conc, A, RA, "AlphaData", 0.35); hold(ax.ax_Pol_Conc, 'on');
            F_Conc = Maps.WQ_States.Pol_Conc_Map(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_Conc(mask) = nan;
            ax.monitor_Pol_Conc = pcolor(ax.ax_Pol_Conc, ax.y_grid, ax.x_grid, F_Conc);
            set(ax.monitor_Pol_Conc, 'EdgeColor', 'none');
            colormap(ax.ax_Pol_Conc, Spectrum);
            % colorbar(ax.ax_Pol_Conc, 'TickDirection', 'out');
            style_colorbar(ax.ax_Pol_Conc);
            shp_Pol_Conc = mapshow(ax.ax_Pol_Conc, S_p, 'FaceColor', 'n'); hold(ax.ax_Pol_Conc, 'on');
            set(ax.ax_Pol_Conc, 'FontName', 'Montserrat');
            ax.ax_Pol_Conc.XAxis.Exponent = 0; ax.ax_Pol_Conc.XAxis.TickLabelFormat = '%.0f';
            ax.ax_Pol_Conc.YAxis.Exponent = 0; ax.ax_Pol_Conc.YAxis.TickLabelFormat = '%.0f';
            ax.ax_Pol_Conc.YAxis.TickLabelRotation = 90;
            ax.ax_Pol_Conc.YAxis.TickValues = ax.ax_Pol_Conc.YAxis.TickValues(1:2:end);

            % Load Map on UIAxes_9
            ax.ax_Load = app.UIAxes_9;
            bm_ax_Load = mapshow(ax.ax_Load, A, RA, "AlphaData", 0.35); hold(ax.ax_Load, 'on');
            F_Load = Maps.WQ_States.Pol_Load_Map(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_Load(mask) = nan;
            ax.monitor_Load = pcolor(ax.ax_Load, ax.y_grid, ax.x_grid, F_Load);
            set(ax.monitor_Load, 'EdgeColor', 'none');
            colormap(ax.ax_Load, Spectrum);
            % colorbar(ax.ax_Load, 'TickDirection', 'out');
            style_colorbar(ax.ax_Load);
            shp_Load = mapshow(ax.ax_Load, S_p, 'FaceColor', 'n'); hold(ax.ax_Load, 'on');
            set(ax.ax_Load, 'FontName', 'Montserrat');
            ax.ax_Load.XAxis.Exponent = 0; ax.ax_Load.XAxis.TickLabelFormat = '%.0f';
            ax.ax_Load.YAxis.Exponent = 0; ax.ax_Load.YAxis.TickLabelFormat = '%.0f';
            ax.ax_Load.YAxis.TickLabelRotation = 90;
            ax.ax_Load.YAxis.TickValues = ax.ax_Load.YAxis.TickValues(1:2:end);

            % Build-up Map on UIAxes_11
            ax.ax_Buildup = app.UIAxes_11;
            bm_ax_Buildup = mapshow(ax.ax_Buildup, A, RA, "AlphaData", 0.35); hold(ax.ax_Buildup, 'on');
            F_Buildup = 1000 * 1/(Resolution^2) * Maps.WQ_States.Pol_Mass_Map(1:ax.DEM_s1, 1:ax.DEM_s2, 1);
            F_Buildup(mask) = nan;
            ax.monitor_Buildup = pcolor(ax.ax_Buildup, ax.y_grid, ax.x_grid, F_Buildup);
            set(ax.monitor_Buildup, 'EdgeColor', 'none');
            colormap(ax.ax_Buildup, Spectrum);
            % colorbar(ax.ax_Buildup, 'TickDirection', 'out');
            style_colorbar(ax.ax_Buildup);

            shp_Buildup = mapshow(ax.ax_Buildup, S_p, 'FaceColor', 'n'); hold(ax.ax_Buildup, 'on');
            set(ax.ax_Buildup, 'FontName', 'Montserrat');
            ax.ax_Buildup.XAxis.Exponent = 0; ax.ax_Buildup.XAxis.TickLabelFormat = '%.0f';
            ax.ax_Buildup.YAxis.Exponent = 0; ax.ax_Buildup.YAxis.TickLabelFormat = '%.0f';
            ax.ax_Buildup.YAxis.TickLabelRotation = 90;
            ax.ax_Buildup.YAxis.TickValues = ax.ax_Buildup.YAxis.TickValues(1:2:end);
        end

        % Finalize depth, rainfall and velocity plots on main axes
        ax.monitor_d = pcolor(ax.ax_d, ax.y_grid, ax.x_grid, F_d);
        set(ax.monitor_d, 'EdgeColor', 'none');
        colormap(ax.ax_d, Spectrum);
        style_colorbar(ax.ax_d)
        % colorbar(ax.ax_d, 'TickDirection', 'out');

        ax.monitor_r = pcolor(ax.ax_r, ax.y_grid, ax.x_grid, F_r);
        set(ax.monitor_r, 'EdgeColor', 'none');
        colormap(ax.ax_r, Spectrum);
        style_colorbar(ax.ax_r);
        % colorbar(ax.ax_r, 'TickDirection', 'out');

        ax.monitor_v = pcolor(ax.ax_v, ax.y_grid, ax.x_grid, F_v);
        set(ax.monitor_v, 'EdgeColor', 'none');
        colormap(ax.ax_v, Velocity_RAS);
        style_colorbar(ax.ax_v)
        % colorbar(ax.ax_v, 'TickDirection', 'out');

        % Overlay the shapefile reference on main axes
        mapshow(ax.ax_d, S_p, 'FaceColor', 'n'); hold(ax.ax_d, 'on');
        mapshow(ax.ax_r, S_p, 'FaceColor', 'n'); hold(ax.ax_d, 'on');
        mapshow(ax.ax_v, S_p, 'FaceColor', 'n'); hold(ax.ax_d, 'on');

        % Set initial status texts
        system_output = 'Initializing HydroPol2D Dashboard..';
        set(ax.ax_date, 'Value', num2str(1));
        set(ax.ax_system_output, 'Value', system_output);
        set(ax.ax_iter, 'Text', num2str(1));

        % Format main axes
        set(ax.ax_d, 'FontName', 'Montserrat');
        ax.ax_d.XAxis.Exponent = 0; ax.ax_d.XAxis.TickLabelFormat = '%.0f';
        ax.ax_d.YAxis.Exponent = 0; ax.ax_d.YAxis.TickLabelFormat = '%.0f';
        ax.ax_d.YAxis.TickLabelRotation = 90;
        ax.ax_d.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
        set(ax.ax_r, 'FontName', 'Montserrat');
        ax.ax_r.XAxis.Exponent = 0; ax.ax_r.XAxis.TickLabelFormat = '%.0f';
        ax.ax_r.YAxis.Exponent = 0; ax.ax_r.YAxis.TickLabelFormat = '%.0f';
        ax.ax_r.YAxis.TickLabelRotation = 90;
        ax.ax_r.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);
        set(ax.ax_v, 'FontName', 'Montserrat');
        ax.ax_v.XAxis.Exponent = 0; ax.ax_v.XAxis.TickLabelFormat = '%.0f';
        ax.ax_v.YAxis.Exponent = 0; ax.ax_v.YAxis.TickLabelFormat = '%.0f';
        ax.ax_v.YAxis.TickLabelRotation = 90;
        ax.ax_v.YAxis.TickValues = ax.ax_d.YAxis.TickValues(1:2:end);

        % Process gauge labels (UTM coordinates conversion if needed)
        char_vector_cell = {};
        try
            for i = 1:numel(gauges.easting_obs_gauges_absolute)
                char_vector_cell = [char_vector_cell; char(string(gauges.labels_observed_string{i}))];
            end
        end
        ax.gauges = char_vector_cell;
        % Update the gauge list in the GUI
        ax.ax_list.Items = char_vector_cell;

        % Ensure AxesHandles exists only once at initialization
        app.AxesHandles = containers.Map();
        
        % Create a vertical layout to hold gauge axes (optional if not using a panel layout)
        if isprop(app, 'GaugePanel')  % Use an existing App Designer panel
            parentContainer = app.GaugePanel;
        else
            parentContainer = app.UIFigure;  % Fallback if GaugePanel doesn't exist
        end
        
        for i = 1:numel(gauges.labels_observed_string)
            gname = strtrim(gauges.labels_observed_string{i});
            
            % Create a new uiaxes inside the chosen panel or UIFigure
            axGauge = uiaxes(parentContainer);
            axGauge.Position = [20, 400 - (i - 1) * 150, 350, 140];  % Adjust as needed
            title(axGauge, gname);
            xlabel(axGauge, 'Time');
            ylabel(axGauge, 'Depth');
        
            % Register the axes to the map
            app.AxesHandles(gname) = axGauge;
        end

        drawnow;
    else
        %% UPDATE SECTION
        % Update DEM dimensions
        ax.DEM_s1 = size(DEM_raster.Z, 1);
        ax.DEM_s2 = size(DEM_raster.Z, 2);

        % Update Depth Plot
        idx_g = Maps.Hydro.d(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
        idx_g(idx_g <= 0.001) = NaN;
        idx_g(mask) = nan;
        set(ax.monitor_d, 'CData', idx_g/1000);
        if ax.flags.flag_reservoir == 1
            set(ax.monitor_bc, 'CData', idx_g/1000);
        end

        % Update ETR plot if needed
        if ax.flags.flag_ETP == 1
            F_ETR = Maps.Hydro.ETR_save(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            F_ETR(F_ETR == 0) = nan;
            set(ax.monitor_ETR, 'CData', F_ETR);

            F_ETP = Maps.Hydro.ETP_save(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            F_ETP(F_ETP == 0) = nan;
            set(ax.monitor_ETP, 'CData', F_ETP);
        end

        % Update Abstraction plot if needed
        if ax.flags.flag_abstraction == 1
            F_Abstraction = Maps.Hydro.Abstraction(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            F_Abstraction(F_Abstraction == 0) = nan;
            set(ax.monitor_Abstraction, 'CData', F_Abstraction);
        end

        % Update Snowpack plot if needed
        if ax.flags.flag_snow_modeling == 1
            F_Snowpack = Maps.Hydro.Snowpack(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            F_Snowpack(F_Snowpack == 0) = nan;
            set(ax.monitor_Snowpack, 'CData', F_Snowpack);
        end

        % Update Baseflow plot if needed
        if ax.flags.flag_baseflow == 1
            set(ax.monitor_GW, 'CData', Maps.Hydro.GWdepth_save(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
        end

        % Update Rainfall Plot
        if ax.flags.flag_spatial_rainfall == 1
            idx_g = Maps.Hydro.spatial_rainfall_maps(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            idx_g(idx_g == 0) = NaN;
        else
            if ax.flags.flag_rainfall == 1
                idx_g = (BC_States.delta_p_agg .* C_a ./ (Resolution^2)) .* ones(ax.DEM_s1, ax.DEM_s2) / (time_step/60);
                idx_g(mask) = nan;
                idx_g(idx_g == 0) = nan;
            else
                idx_g = zeros(ax.DEM_s1, ax.DEM_s2);
                idx_g(mask) = nan;
                idx_g(idx_g == 0) = nan;
            end
        end
        set(ax.monitor_r, 'CData', gather(idx_g));

        % Update Infiltration and related plots if flag set
        if ax.flags.flag_infiltration == 1
            set(ax.monitor_i, 'CData', Maps.Hydro.I_t(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            set(ax.monitor_C, 'CData', Maps.Hydro.C(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_g = Maps.Hydro.f(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            idx_g(idx_g == 0) = nan;
            idx_g(mask) = nan;
            set(ax.monitor_f, 'CData', idx_g);
            idx_g = Maps.Hydro.I_t(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            idx_g(idx_g == 0) = NaN;
            idx_g(mask) = nan;
            set(ax.monitor_i, 'CData', idx_g);
        end

        % Update Water Quality plots if flag set
        if ax.flags.flag_waterquality == 1
            % Concentration
            idx_g = Maps.WQ_States.Pol_Conc_Map(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            idx_g(idx_g == 0) = nan;
            idx_g(mask) = nan;
            set(ax.monitor_Pol_Conc, 'CData', idx_g);

            % Load
            idx_g = Maps.WQ_States.Pol_Load_Map(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            idx_g(idx_g == 0) = nan;
            idx_g(mask) = nan;
            set(ax.monitor_Load, 'CData', idx_g);

            % Build-up
            idx_g = 1000 * (1/Resolution^2) * Maps.WQ_States.Pol_Mass_Map(1:ax.DEM_s1, 1:ax.DEM_s2, layer);
            idx_g(idx_g == 0) = nan;
            idx_g(mask) = nan;
            set(ax.monitor_Buildup, 'CData', idx_g);
        end

        % Update Velocity plot
        v_t(v_t == 0) = NaN;
        set(ax.monitor_v, 'CData', v_t);

        % Update system status text and iteration counter
        set(ax.ax_date, 'Value', datestr(ax.timer));
        % set(ax.ax_system_output, 'Value', strcat('Percentage Complete :', num2str(ax.percentage), '%'));
        set(ax.ax_system_output, 'Value', sprintf('Percentage Complete : %.2f %%', ax.percentage));
        set(ax.ax_iter, 'Text', num2str(layer));

        % === UNIFORM STYLING FOR ALL PLOTTED AXES ===
        standard_fields = {
            'ax_d','ax_r','ax_v','ax_C','ax_f','ax_i',...
            'ax_ETR','ax_ETP','ax_Abstraction','ax_Snowpack','ax_GW',...
            'ax_bc','ax_Pol_Conc','ax_Load','ax_Buildup'
        };
        
        for k = 1:numel(standard_fields)
            field = standard_fields{k};
            if isfield(ax, field) && ~isempty(ax.(field)) && isgraphics(ax.(field))
                style_axes(ax.(field));
                style_colorbar(ax.(field));
            end
        end


        drawnow;
        pause(0.1);
    end
catch e
    disp('Error occurred');
    disp(e.message);
    disp(e.stack(1));
end
end


function style_colorbar(ax_handle)
% STYLE_COLORBAR standardizes colorbar appearance for a given axis
    cb = colorbar(ax_handle, 'TickDirection', 'in');
    cb.Box = 'on';
    cb.TickLength = 0.02;
    cb.LineWidth = 1.5;
    cb.FontName = 'Montserrat';
    cb.FontSize = 12;
end

function style_axes(ax_handle, varargin)
% STYLE_AXES formats a UIAxes or standard axes for consistent style.
%
% Optional name-value pairs:
%   'XLabel', 'YLabel', 'Title', 'FontName', 'FontSize', 'Grid', 'TickDirection', 'LineWidth'

    p = inputParser;
    addParameter(p, 'XLabel', '', @ischar);
    addParameter(p, 'YLabel', '', @ischar);
    addParameter(p, 'Title', '', @ischar);
    addParameter(p, 'FontName', 'Montserrat', @ischar);
    addParameter(p, 'FontSize', 14, @isnumeric);
    addParameter(p, 'TickDirection', 'in', @ischar);
    addParameter(p, 'Grid', 'on', @(x) ismember(x, {'on', 'off'}));
    addParameter(p, 'LineWidth', 1.8, @isnumeric);
    parse(p, varargin{:});

    % Apply styles to axis
    set(ax_handle, ...
        'FontName', p.Results.FontName, ...
        'FontSize', p.Results.FontSize, ...
        'LineWidth', p.Results.LineWidth, ...
        'XGrid', p.Results.Grid, ...
        'YGrid', p.Results.Grid, ...
        'Box', 'on', ...
        'TickDir', p.Results.TickDirection);

    ax_handle.XAxis.TickLabelFormat = '%.0f';
    ax_handle.YAxis.TickLabelFormat = '%.0f';
    ax_handle.XAxis.Exponent = 0;
    ax_handle.YAxis.Exponent = 0;
    ax_handle.YAxis.TickLabelRotation = 90;

    % Apply optional labels
    if ~isempty(p.Results.XLabel)
        xlabel(ax_handle, p.Results.XLabel);
    end
    if ~isempty(p.Results.YLabel)
        ylabel(ax_handle, p.Results.YLabel);
    end
    if ~isempty(p.Results.Title)
        title(ax_handle, p.Results.Title);
    end
end
