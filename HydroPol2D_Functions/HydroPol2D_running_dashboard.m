function [ax] = HydroPol2D_running_dashboard(ax, Maps, v_t, DEM_raster, gauges, BC_States, time_step, Resolution, first_time, layer, C_a)
% HydroPol2D_RUNNING_DASHBOARD sets up or updates the dashboard for HydroPol2D.
%
% FULL HPC-SAFE VERSION (Sherlock / headless / no-internet / UIHTML blocked):
%   ✅ Detects whether UIFigures can truly be shown (DISPLAY + JVM + ShowFigureWindows)
%   ✅ Skips dashboard entirely when headless (never breaks simulation)
%   ✅ HARD-DISABLES any uihtml components inside the App (prevents "Not Found")
%   ✅ Never attempts basemap downloads unless explicitly enabled AND safe
%   ✅ Gathers gpuArray data before plotting where needed
%   ✅ Keeps your original plotting logic / axes styling
%
% Inputs:
%   ax, Maps, v_t, DEM_raster, gauges, BC_States, time_step, Resolution,
%   first_time, layer, C_a
%
% Output:
%   ax - updated structure with axes handles and monitor objects

% =========================
%  0) BASIC SAFETY / MASKS
% =========================
mask = isnan(DEM_raster.Z);

% =========================
%  1) HPC / HEADLESS CHECK  (ROBUST)
% =========================
hasFigureWindows = logical(feature('ShowFigureWindows'));
hasDisplay = ~isempty(getenv('DISPLAY')) || ~isempty(getenv('WAYLAND_DISPLAY'));  % Linux / some systems
hasJVM = usejava('jvm');

% Consider headless only if we truly cannot show figure windows
isHeadless = ~(hasFigureWindows && hasDisplay && hasJVM);

% =========================
%  1.25) FORCE NO-INTERNET DEFAULT (IMPORTANT)
% =========================
% Default: do NOT use basemaps (prevents ANY tile download attempts).
% You can re-enable locally by setting ax.flags.flag_basemap = 1.
useBasemap = false;

if isfield(ax,'flags') && isfield(ax.flags,'flag_basemap')
    useBasemap = (ax.flags.flag_basemap == 1) && ~isHeadless;
end

% If headless and first_time init, do NOT create UI.
if first_time == 1 && isHeadless
    % Dashboard skipped gracefully; simulation can proceed.
    if isfield(ax,'ax_system_output')
        try
            set(ax.ax_system_output,'Value','Dashboard skipped (headless/HPC mode).');
        catch
        end
    end
    return;
end

% =========================
%  1.5) GAUGE LABELS (CPU SAFE)
% =========================
gaugeLabels = {};
if isfield(ax,'gaugeLabels') && ~isempty(ax.gaugeLabels)
    gaugeLabels = ax.gaugeLabels;
elseif ~isempty(gauges) && isfield(gauges,'labels_observed_string') && ~isempty(gauges.labels_observed_string)
    gaugeLabels = gauges.labels_observed_string;
end

if isempty(gaugeLabels) && isfield(ax,'app') && isprop(ax.app,'gauges_data')
    try
        gaugeLabels = ax.app.gauges_data.Properties.VariableNames;
    catch
    end
end

try
    gaugeLabels = cellstr(string(gaugeLabels));
catch
    gaugeLabels = {};
end

% =========================
%  2) CRS CHECK
% =========================
web_mercator_crs = projcrs(3857);
try
    DEM_CRS = DEM_raster.georef.SpatialRef.ProjectedCRS;
    isWebMercator = isequal(DEM_CRS, web_mercator_crs);
catch
    warning("DEM_raster does not define a ProjectedCRS. Assuming non-Web Mercator.");
    isWebMercator = false;
end

% =========================
%  3) COLORRAMPS
% =========================
[Spectrum, ~, ~, ~, ~, ~, Depth_RAS, ~, Velocity_RAS, WSE_RAS] = coloramps(); %#ok<ASGLU>

% =========================
%  4) MAIN TRY/CATCH
% =========================
try
    if first_time == 1
        %% ================================================================
        %  INITIALIZATION SECTION
        % ================================================================
        disp('Starting HydroPol2D_running_dashboard');

        % ---- Create and show the dashboard application ----
        app = HydroPol2D_Monitor_base();

        % ---- HARD FIX: kill any UIHTML panels that try to load URLs ----
        % This is the most common reason for the big "Not Found" screen.
        try
            htmlObjs = findall(app.UIFigure,'Type','uihtml');
            for ii = 1:numel(htmlObjs)
                delete(htmlObjs(ii));
            end
        catch
        end

        % Show the app
        try
            app.UIFigure.Visible = 'on';
        catch
        end
        disp('App created and set to visible');

        % ---- Assign GUI handles ----
        ax.ax_d = app.UIAxes;
        ax.ax_r = app.UIAxes2;
        ax.ax_date = app.DateTimeTextArea;
        ax.ax_iter = app.iter;
        ax.ax_v = app.UIAxes_2;
        ax.ax_system_output = app.SystemOutput;
        ax.ax_list = app.gauges_list;
        ax.app = app;
        if isfield(ax,'flags')
            app.flags = ax.flags; % pass flag settings to the app
        end

        % Newer GUI handles
        ax.ax_GW = app.GWdepth;
        ax.ax_C = app.UIAxes_3;
        ax.ax_ETR = app.UIAxes_7;
        ax.ax_ETP = app.UIAxes_6;
        ax.ax_Abstraction = app.UIAxes_13;
        ax.ax_i = app.UIAxes_4;
        ax.ax_f = app.UIAxes_5;

        % ================================================================
        %  Basemap (STRICTLY OPTIONAL; OFF by default)
        % ================================================================
        A = []; RA = [];

        if useBasemap && isWebMercator
            basemapName = "openstreetmap";
            url = "https://a.tile.openstreetmap.org/${z}/${x}/${y}.png"; % FIXED
            copyright = char(uint8(169));
            attribution = copyright + "OpenStreetMap contributors";

            try
                addCustomBasemap(basemapName, url, "Attribution", attribution);
            catch
            end

            try
                [lat, lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS, ...
                    DEM_raster.georef.SpatialRef.XWorldLimits, ...
                    DEM_raster.georef.SpatialRef.YWorldLimits);
                latlim = [lat(1), lat(2)];
                lonlim = [lon(1), lon(2)];

                [A, RA] = readBasemapImage(basemapName, latlim, lonlim);
            catch ME
                warning(ME.message);
                A = []; RA = [];
            end
        else
            if ~isWebMercator
                warning("DEM is not in Web Mercator. Basemap will be skipped.");
            end
        end

        %% ================================================================
        %  Create a reference polygon (domain outline) from DEM
        % ================================================================
        binaryMask = ~isnan(DEM_raster.Z);
        boundaries = bwboundaries(binaryMask);

        combinedX = [];
        combinedY = [];
        for kk = 1:numel(boundaries)
            boundary = boundaries{kk};
            X = boundary(:, 2);
            Y = boundary(:, 1);
            combinedX = [combinedX; X; NaN]; %#ok<AGROW>
            combinedY = [combinedY; Y; NaN]; %#ok<AGROW>
        end
        if ~isempty(combinedX)
            combinedX = combinedX(1:end-1);
            combinedY = combinedY(1:end-1);
        end

        S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
        S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); ...
            DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)];

        S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1) + ...
            combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';

        S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2) - ...
            combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';

        %% ================================================================
        %  Spatial grid for monitors
        % ================================================================
        ax.DEM_s1 = size(DEM_raster.Z, 1);
        ax.DEM_s2 = size(DEM_raster.Z, 2);

        xbegin = 1; ybegin = 1;
        xend = ax.DEM_s1;
        yend = ax.DEM_s2;

        ax.x_grid = DEM_raster.georef.SpatialRef.YWorldLimits(2) - ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX * (xbegin:xend);

        ax.y_grid = DEM_raster.georef.SpatialRef.XWorldLimits(1) + ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX * (ybegin:yend);

        %% ================================================================
        %  Optional basemap underlay
        % ================================================================
        if ~isempty(A) && ~isempty(RA)
            mapshow(ax.ax_d, A, RA, "AlphaData", 0.35); hold(ax.ax_d, 'on');
            mapshow(ax.ax_r, A, RA, "AlphaData", 0.35); hold(ax.ax_r, 'on');
            mapshow(ax.ax_v, A, RA, "AlphaData", 0.35); hold(ax.ax_v, 'on');
            mapshow(ax.ax_GW, A, RA, "AlphaData", 0.35); hold(ax.ax_GW, 'on');
            mapshow(ax.ax_C, A, RA, "AlphaData", 0.35); hold(ax.ax_C, 'on');
            mapshow(ax.ax_ETR, A, RA, "AlphaData", 0.35); hold(ax.ax_ETR, 'on');
            mapshow(ax.ax_ETP, A, RA, "AlphaData", 0.35); hold(ax.ax_ETP, 'on');
            mapshow(ax.ax_Abstraction, A, RA, "AlphaData", 0.35); hold(ax.ax_Abstraction, 'on');
            mapshow(ax.ax_i, A, RA, "AlphaData", 0.35); hold(ax.ax_i, 'on');
            mapshow(ax.ax_f, A, RA, "AlphaData", 0.35); hold(ax.ax_f, 'on');
        end

        %% ================================================================
        %  Initial Fields (GPU-safe)
        % ================================================================
        F_d = localGather(Maps.Hydro.d(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
        F_d(mask) = nan;

        if isfield(ax,'flags') && isfield(ax.flags,'flag_spatial_rainfall') && ax.flags.flag_spatial_rainfall == 1
            F_r = localGather(Maps.Hydro.spatial_rainfall_maps(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
            F_r(mask) = nan;
        else
            if isfield(ax,'flags') && isfield(ax.flags,'flag_rainfall') && ax.flags.flag_rainfall == 1
                F_r = localGather((BC_States.delta_p_agg / (time_step/60)) * ones(ax.DEM_s1, ax.DEM_s2));
                F_r(mask) = nan;
            else
                F_r = zeros(ax.DEM_s1, ax.DEM_s2);
                F_r(mask) = nan;
            end
        end

        F_v = localGather(double(v_t(1:ax.DEM_s1, 1:ax.DEM_s2, 1)));
        F_v(mask) = nan;

        %% ================================================================
        %  Infiltration plots
        % ================================================================
        if ax.flags.flag_infiltration == 1
            F_i = localGather(Maps.Hydro.I_t(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_i(mask) = nan;
            ax.monitor_i = pcolor(ax.ax_i, ax.y_grid, ax.x_grid, F_i); set(ax.monitor_i,'EdgeColor','none');
            mapshow(ax.ax_i, S_p, 'FaceColor', 'n'); hold(ax.ax_i,'on');
            colormap(ax.ax_i, WSE_RAS); style_colorbar(ax.ax_i);

            F_C = localGather(Maps.Hydro.C(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_C(mask) = nan;
            ax.monitor_C = pcolor(ax.ax_C, ax.y_grid, ax.x_grid, F_C); set(ax.monitor_C,'EdgeColor','none');
            mapshow(ax.ax_C, S_p, 'FaceColor', 'n'); hold(ax.ax_C,'on');
            colormap(ax.ax_C, Spectrum); style_colorbar(ax.ax_C);

            F_f = localGather(Maps.Hydro.f(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_f(mask) = nan;
            ax.monitor_f = pcolor(ax.ax_f, ax.y_grid, ax.x_grid, F_f); set(ax.monitor_f,'EdgeColor','none');
            mapshow(ax.ax_f, S_p, 'FaceColor', 'n'); hold(ax.ax_f,'on');
            colormap(ax.ax_f, Velocity_RAS); style_colorbar(ax.ax_f);
        end

        %% ================================================================
        %  Baseflow plot
        % ================================================================
        if ax.flags.flag_baseflow == 1
            GW_depth = localGather(Maps.Hydro.GWdepth_save(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
            GW_depth(mask) = nan;
            ax.monitor_GW = pcolor(ax.ax_GW, ax.y_grid, ax.x_grid, GW_depth);
            set(ax.monitor_GW, 'EdgeColor', 'none');
            mapshow(ax.ax_GW, S_p, 'FaceColor', 'n'); hold(ax.ax_GW, 'on');
            colormap(ax.ax_GW, Velocity_RAS);
            style_colorbar(ax.ax_GW);
        end

        %% ================================================================
        %  Reservoir plot (kept)
        % ================================================================
        if isfield(ax.flags,'flag_reservoir') && ax.flags.flag_reservoir == 1
            ax.ax_bc = app.UIAxes_12;
            ax.ax_list_bc = app.boundary_list;
            ax.monitor_bc = pcolor(ax.ax_bc, ax.y_grid, ax.x_grid, F_d);
            set(ax.monitor_bc, 'EdgeColor', 'none');
            colormap(ax.ax_bc, Velocity_RAS);
            style_colorbar(ax.ax_bc);
            mapshow(ax.ax_bc, S_p, 'FaceColor', 'n'); hold(ax.ax_bc, 'on');
        end

        %% ================================================================
        %  ETP/ETR plots
        % ================================================================
        if ax.flags.flag_ETP == 1
            F_ETR = localGather(Maps.Hydro.ETR_save(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_ETR(mask)=nan; F_ETR(F_ETR==0)=nan;
            ax.monitor_ETR = pcolor(ax.ax_ETR, ax.y_grid, ax.x_grid, F_ETR); set(ax.monitor_ETR,'EdgeColor','none');
            colormap(ax.ax_ETR, Spectrum); style_colorbar(ax.ax_ETR);
            mapshow(ax.ax_ETR, S_p, 'FaceColor', 'n'); hold(ax.ax_ETR,'on');

            F_ETP = localGather(Maps.Hydro.ETP_save(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_ETP(mask)=nan;
            ax.monitor_ETP = pcolor(ax.ax_ETP, ax.y_grid, ax.x_grid, F_ETP); set(ax.monitor_ETP,'EdgeColor','none');
            colormap(ax.ax_ETP, Spectrum); style_colorbar(ax.ax_ETP);
            mapshow(ax.ax_ETP, S_p, 'FaceColor', 'n'); hold(ax.ax_ETP,'on');
        end

        %% ================================================================
        %  Abstraction plot
        % ================================================================
        if isfield(ax.flags,'flag_abstraction') && ax.flags.flag_abstraction == 1
            F_Abstraction = localGather(Maps.Hydro.Abstraction(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
            F_Abstraction(mask) = nan;
            F_Abstraction(F_Abstraction <= 1e-3) = nan;
            ax.monitor_Abstraction = pcolor(ax.ax_Abstraction, ax.y_grid, ax.x_grid, F_Abstraction);
            set(ax.monitor_Abstraction, 'EdgeColor', 'none');
            colormap(ax.ax_Abstraction, Depth_RAS);
            style_colorbar(ax.ax_Abstraction);
            mapshow(ax.ax_Abstraction, S_p, 'FaceColor', 'n'); hold(ax.ax_Abstraction, 'on');
        end

        %% ================================================================
        %  Snowpack plot
        % ================================================================
        if isfield(ax.flags,'flag_snow_modeling') && ax.flags.flag_snow_modeling == 1
            ax.ax_Snowpack = app.UIAxes_14;
            F_Snowpack = localGather(Maps.Hydro.Snowpack(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
            F_Snowpack(mask) = nan; F_Snowpack(F_Snowpack==0)=nan;
            ax.monitor_Snowpack = pcolor(ax.ax_Snowpack, ax.y_grid, ax.x_grid, F_Snowpack);
            set(ax.monitor_Snowpack,'EdgeColor','none');
            colormap(ax.ax_Snowpack, Spectrum);
            style_colorbar(ax.ax_Snowpack);
            mapshow(ax.ax_Snowpack, S_p, 'FaceColor','n'); hold(ax.ax_Snowpack,'on');
        end

        %% ================================================================
        %  Water Quality plots
        % ================================================================
        if isfield(ax.flags,'flag_waterquality') && ax.flags.flag_waterquality == 1
            ax.ax_Pol_Conc = app.UIAxes_8;
            F_Conc = localGather(Maps.WQ_States.Pol_Conc_Map(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_Conc(mask)=nan;
            ax.monitor_Pol_Conc = pcolor(ax.ax_Pol_Conc, ax.y_grid, ax.x_grid, F_Conc);
            set(ax.monitor_Pol_Conc,'EdgeColor','none');
            colormap(ax.ax_Pol_Conc,Spectrum); style_colorbar(ax.ax_Pol_Conc);
            mapshow(ax.ax_Pol_Conc,S_p,'FaceColor','n'); hold(ax.ax_Pol_Conc,'on');

            ax.ax_Load = app.UIAxes_9;
            F_Load = localGather(Maps.WQ_States.Pol_Load_Map(1:ax.DEM_s1, 1:ax.DEM_s2, 1)); F_Load(mask)=nan;
            ax.monitor_Load = pcolor(ax.ax_Load, ax.y_grid, ax.x_grid, F_Load);
            set(ax.monitor_Load,'EdgeColor','none');
            colormap(ax.ax_Load,Spectrum); style_colorbar(ax.ax_Load);
            mapshow(ax.ax_Load,S_p,'FaceColor','n'); hold(ax.ax_Load,'on');

            ax.ax_Buildup = app.UIAxes_11;
            F_Buildup = localGather(1000*(1/(Resolution^2))*Maps.WQ_States.Pol_Mass_Map(1:ax.DEM_s1, 1:ax.DEM_s2, 1));
            F_Buildup(mask)=nan;
            ax.monitor_Buildup = pcolor(ax.ax_Buildup, ax.y_grid, ax.x_grid, F_Buildup);
            set(ax.monitor_Buildup,'EdgeColor','none');
            colormap(ax.ax_Buildup,Spectrum); style_colorbar(ax.ax_Buildup);
            mapshow(ax.ax_Buildup,S_p,'FaceColor','n'); hold(ax.ax_Buildup,'on');
        end

        %% ================================================================
        %  Main Depth/Rain/Velocity monitors
        % ================================================================
        ax.monitor_d = pcolor(ax.ax_d, ax.y_grid, ax.x_grid, F_d);
        set(ax.monitor_d, 'EdgeColor', 'none');
        colormap(ax.ax_d, Spectrum); style_colorbar(ax.ax_d);
        mapshow(ax.ax_d, S_p, 'FaceColor', 'n'); hold(ax.ax_d, 'on');

        ax.monitor_r = pcolor(ax.ax_r, ax.y_grid, ax.x_grid, F_r);
        set(ax.monitor_r, 'EdgeColor', 'none');
        colormap(ax.ax_r, Spectrum); style_colorbar(ax.ax_r);
        mapshow(ax.ax_r, S_p, 'FaceColor', 'n'); hold(ax.ax_r, 'on');

        ax.monitor_v = pcolor(ax.ax_v, ax.y_grid, ax.x_grid, F_v);
        set(ax.monitor_v, 'EdgeColor', 'none');
        colormap(ax.ax_v, Velocity_RAS); style_colorbar(ax.ax_v);
        mapshow(ax.ax_v, S_p, 'FaceColor', 'n'); hold(ax.ax_v, 'on');

        %% ================================================================
        %  Initial status texts
        % ================================================================
        system_output = 'Initializing HydroPol2D Dashboard..';
        try, set(ax.ax_date, 'Value', num2str(1)); catch, end
        try, set(ax.ax_system_output, 'Value', system_output); catch, end
        try, set(ax.ax_iter, 'Text', num2str(1)); catch, end

        %% ================================================================
        %  Gauge list population
        % ================================================================
        ax.gaugeLabels = gaugeLabels;
        char_vector_cell = {};
        try
            for i = 1:numel(gaugeLabels)
                char_vector_cell = [char_vector_cell; char(string(gaugeLabels{i}))]; %#ok<AGROW>
            end
        catch
        end
        ax.gauges = char_vector_cell;
        try
            ax.ax_list.Items = char_vector_cell;
        catch
        end

        drawnow;

    else
        %% ================================================================
        %  UPDATE SECTION
        % ================================================================
        ax.DEM_s1 = size(DEM_raster.Z, 1);
        ax.DEM_s2 = size(DEM_raster.Z, 2);

        % ---- Depth ----
        idx_g = localGather(Maps.Hydro.d(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
        idx_g(idx_g <= 0.001) = NaN;
        idx_g(mask) = nan;
        set(ax.monitor_d, 'CData', idx_g/1000);
        if isfield(ax.flags,'flag_reservoir') && ax.flags.flag_reservoir == 1 && isfield(ax,'monitor_bc')
            set(ax.monitor_bc, 'CData', idx_g/1000);
        end

        % ---- ETR/ETP ----
        if ax.flags.flag_ETP == 1
            F_ETR = localGather(Maps.Hydro.ETR_save(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            F_ETR(F_ETR == 0) = nan;
            set(ax.monitor_ETR, 'CData', F_ETR);

            F_ETP = localGather(Maps.Hydro.ETP_save(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            F_ETP(F_ETP == 0) = nan;
            set(ax.monitor_ETP, 'CData', F_ETP);
        end

        % ---- Abstraction ----
        if isfield(ax.flags,'flag_abstraction') && ax.flags.flag_abstraction == 1
            F_Abstraction = localGather(Maps.Hydro.Abstraction(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            F_Abstraction(F_Abstraction == 0) = nan;
            set(ax.monitor_Abstraction, 'CData', F_Abstraction);
        end

        % ---- Snowpack ----
        if isfield(ax.flags,'flag_snow_modeling') && ax.flags.flag_snow_modeling == 1
            F_Snowpack = localGather(Maps.Hydro.Snowpack(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            F_Snowpack(F_Snowpack == 0) = nan;
            set(ax.monitor_Snowpack, 'CData', F_Snowpack);
        end

        % ---- Baseflow ----
        if ax.flags.flag_baseflow == 1
            set(ax.monitor_GW, 'CData', localGather(Maps.Hydro.GWdepth_save(1:ax.DEM_s1, 1:ax.DEM_s2, layer)));
        end

        % ---- Rainfall ----
        if isfield(ax.flags,'flag_spatial_rainfall') && ax.flags.flag_spatial_rainfall == 1
            idx_r = localGather(Maps.Hydro.spatial_rainfall_maps(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_r(idx_r == 0) = NaN;
        else
            if isfield(ax.flags,'flag_rainfall') && ax.flags.flag_rainfall == 1
                idx_r = localGather((BC_States.delta_p_agg .* C_a ./ (Resolution^2)) .* ones(ax.DEM_s1, ax.DEM_s2) / (time_step/60));
                idx_r(mask) = nan;
                idx_r(idx_r == 0) = nan;
            else
                idx_r = zeros(ax.DEM_s1, ax.DEM_s2);
                idx_r(mask) = nan;
                idx_r(idx_r == 0) = nan;
            end
        end
        set(ax.monitor_r, 'CData', idx_r);

        % ---- Infiltration ----
        if ax.flags.flag_infiltration == 1
            idx_i = localGather(Maps.Hydro.I_t(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_i(idx_i == 0) = NaN;
            idx_i(mask) = nan;
            set(ax.monitor_i, 'CData', idx_i);

            set(ax.monitor_C, 'CData', localGather(Maps.Hydro.C(1:ax.DEM_s1, 1:ax.DEM_s2, layer)));

            idx_f = localGather(Maps.Hydro.f(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_f(idx_f == 0) = nan;
            idx_f(mask) = nan;
            set(ax.monitor_f, 'CData', idx_f);
        end

        % ---- Water Quality ----
        if isfield(ax.flags,'flag_waterquality') && ax.flags.flag_waterquality == 1
            idx_c = localGather(Maps.WQ_States.Pol_Conc_Map(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_c(idx_c == 0) = nan; idx_c(mask) = nan;
            set(ax.monitor_Pol_Conc, 'CData', idx_c);

            idx_l = localGather(Maps.WQ_States.Pol_Load_Map(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_l(idx_l == 0) = nan; idx_l(mask) = nan;
            set(ax.monitor_Load, 'CData', idx_l);

            idx_b = localGather(1000 * (1/Resolution^2) * Maps.WQ_States.Pol_Mass_Map(1:ax.DEM_s1, 1:ax.DEM_s2, layer));
            idx_b(idx_b == 0) = nan; idx_b(mask) = nan;
            set(ax.monitor_Buildup, 'CData', idx_b);
        end

        % ---- Velocity ----
        vnow = localGather(v_t);
        vnow(vnow == 0) = NaN;
        set(ax.monitor_v, 'CData', vnow);

        % ---- Status ----
        try, set(ax.ax_date, 'Value', datestr(ax.timer)); catch, end
        try, set(ax.ax_system_output, 'Value', sprintf('Percentage Complete : %.2f %%', ax.percentage)); catch, end
        try, set(ax.ax_iter, 'Text', num2str(layer)); catch, end

        % ---- Styling ----
        standard_fields = {
            'ax_d','ax_r','ax_v','ax_C','ax_f','ax_i',...
            'ax_ETR','ax_ETP','ax_Abstraction','ax_Snowpack','ax_GW',...
            'ax_bc','ax_Pol_Conc','ax_Load','ax_Buildup'
        };
        for kk = 1:numel(standard_fields)
            field = standard_fields{kk};
            if isfield(ax, field) && ~isempty(ax.(field)) && isgraphics(ax.(field))
                style_axes(ax.(field));
                style_colorbar(ax.(field));
            end
        end

        drawnow;
        pause(0.05);
    end

catch e
    disp('Error occurred in HydroPol2D_running_dashboard');
    disp(e.message);
    if ~isempty(e.stack)
        disp(e.stack(1));
    end
end
end

% =========================================================================
% Helpers
% =========================================================================
function A = localGather(A)
% Gather only if gpuArray; otherwise return as-is.
if isa(A,'gpuArray')
    A = gather(A);
end
end

function style_colorbar(ax_handle)
cb = colorbar(ax_handle, 'TickDirection', 'in');
cb.Box = 'on';
cb.TickLength = 0.02;
cb.LineWidth = 1.5;
cb.FontName = 'Montserrat';
cb.FontSize = 12;
end

function style_axes(ax_handle, varargin)
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

if ~isempty(p.Results.XLabel); xlabel(ax_handle, p.Results.XLabel); end
if ~isempty(p.Results.YLabel); ylabel(ax_handle, p.Results.YLabel); end
if ~isempty(p.Results.Title);  title(ax_handle, p.Results.Title);  end
end