function HydroPol2D_running_dashboard_plot(Maps, v_t, DEM_raster, gauges, BC_States, time_step, Resolution, C_a, k)
% =========================================================================
% HydroPol2D_running_dashboard_plot
% =========================================================================
% Optimized runtime dashboard for HydroPol2D using a standard MATLAB figure.
%
% Main features
% -------------
% 1) Creates figure/axes/objects only once
% 2) Updates only CData/text on subsequent calls
% 3) Uses original HydroPol2D coloramps()
% 4) Uses LaTeX + Garamond + outward ticks + thicker colorbars
% 5) Preserves domain outline
% 6) HPC/headless safe
%
% Inputs
% ------
% Maps       : model state/output structure
% v_t        : velocity map (2D or 3D, CPU or gpuArray)
% DEM_raster : DEM raster structure with DEM_raster.Z and georef
% gauges     : gauge structure (optional)
% BC_States  : BC structure
% time_step  : current model time step [s]
% Resolution : cell size [m]
% C_a        : area factor used in rainfall conversion
% k          : current layer / iteration index
% =========================================================================

    persistent DASH

    % =====================================================================
    % 0) HEADLESS SAFETY
    % =====================================================================
    hasFigureWindows = logical(feature('ShowFigureWindows'));
    hasDisplay = ~isempty(getenv('DISPLAY')) || ~isempty(getenv('WAYLAND_DISPLAY'));
    hasJVM = usejava('jvm');
    isHeadless = ~(hasFigureWindows && hasDisplay && hasJVM);

    if isHeadless
        return;
    end

    % =====================================================================
    % 1) STYLE / LABEL CUSTOMIZATION
    % =====================================================================
    S = local_style();

    % =====================================================================
    % 2) BASIC GEOMETRY
    % =====================================================================
    DEM = localGather(DEM_raster.Z);
    DEM = double(DEM);

    if isempty(DEM) || ~ismatrix(DEM)
        warning('HydroPol2D_running_dashboard_plot: invalid DEM_raster.Z.');
        return;
    end

    [nr, nc] = size(DEM);
    mask = isnan(DEM);

    % =====================================================================
    % 3) INITIALIZE / REBUILD ONLY IF NEEDED
    % =====================================================================
    needRebuild = false;

    if isempty(DASH) || ~isstruct(DASH)
        needRebuild = true;
    elseif ~isfield(DASH, 'fig') || ~isgraphics(DASH.fig)
        needRebuild = true;
    elseif ~isfield(DASH, 'nr') || ~isfield(DASH, 'nc')
        needRebuild = true;
    elseif DASH.nr ~= nr || DASH.nc ~= nc
        needRebuild = true;
    end

    if needRebuild
        DASH = local_initialize_dashboard(DEM_raster, DEM, mask, Resolution, S);
    end

    % =====================================================================
    % 4) LOAD ORIGINAL COLORAMPS
    % =====================================================================
    try
        [Spectrum, ~, ~, ~, ~, ~, Depth_RAS, ~, Velocity_RAS, WSE_RAS] = coloramps(); %#ok<ASGLU>
    catch
        warning('coloramps() not found. Falling back to MATLAB colormaps.');
        Spectrum     = parula(256);
        Depth_RAS    = turbo(256);
        Velocity_RAS = turbo(256);
        WSE_RAS      = parula(256);
    end

    % =====================================================================
    % 5) EXTRACT / PREPARE FIELDS
    % =====================================================================

    % ---- Water depth [model mm -> plot m]
    F_d = local_get_map(Maps, {'Hydro','d'}, k);
    if ~isempty(F_d)
        F_d = double(F_d);
        F_d(F_d <= 1e-3) = NaN;
        F_d(mask) = NaN;
        F_d = F_d ./ 1000;
    end

    % ---- Rainfall [mm h^-1]
    F_r = [];
    spatialRain = local_get_map(Maps, {'Hydro','spatial_rainfall_maps'}, k);
    if ~isempty(spatialRain)
        F_r = double(spatialRain);
        F_r(mask) = NaN;
        F_r(F_r == 0) = NaN;
    else
        if isfield(BC_States, 'delta_p_agg') && ~isempty(BC_States.delta_p_agg)
            try
                F_r = double((BC_States.delta_p_agg .* C_a ./ (Resolution^2)) .* ones(nr, nc) / max(time_step/60, eps));
            catch
                F_r = double((BC_States.delta_p_agg) .* ones(nr, nc) / max(time_step/60, eps));
            end
            F_r(mask) = NaN;
            F_r(F_r == 0) = NaN;
        end
    end

    % ---- Velocity [m s^-1]
    F_v = local_get_array(v_t, k);
    if ~isempty(F_v)
        F_v = double(F_v);
        F_v(mask) = NaN;
        F_v(F_v == 0) = NaN;
    end

    % ---- Infiltration storage I_t [mm]
    F_i = local_get_map(Maps, {'Hydro','I_t'}, k);
    if ~isempty(F_i)
        F_i = double(F_i);
        F_i(mask) = NaN;
        F_i(F_i == 0) = NaN;
    end

    % ---- State C [mm]  <-- adjust here if your internal unit differs
    F_C = local_get_map(Maps, {'Hydro','C'}, k);
    if ~isempty(F_C)
        F_C = double(F_C);
        F_C(mask) = NaN;
        F_C(F_C == 0) = NaN;
    end

    % ---- Infiltration rate f [mm h^-1]
    F_f = local_get_map(Maps, {'Hydro','f'}, k);
    if ~isempty(F_f)
        F_f = double(F_f);
        F_f(mask) = NaN;
        F_f(F_f == 0) = NaN;
    end

    % ---- Groundwater depth [m]
    F_GW = local_get_map(Maps, {'Hydro','GWdepth_save'}, k);
    if ~isempty(F_GW)
        F_GW = double(F_GW);
        F_GW(mask) = NaN;
        F_GW(F_GW == 0) = NaN;
    end

    % ---- Actual ET [mm h^-1]
    F_ETR = local_get_map(Maps, {'Hydro','ETR_save'}, k);
    if ~isempty(F_ETR)
        F_ETR = double(F_ETR);
        F_ETR(mask) = NaN;
        F_ETR(F_ETR == 0) = NaN;
    end

    % ---- Potential ET [mm h^-1]
    F_ETP = local_get_map(Maps, {'Hydro','ETP_save'}, k);
    if ~isempty(F_ETP)
        F_ETP = double(F_ETP);
        F_ETP(mask) = NaN;
        F_ETP(F_ETP == 0) = NaN;
    end

    % ---- Optional footer-only fields
    F_Abs = local_get_map(Maps, {'Hydro','Abstraction'}, k);
    if ~isempty(F_Abs)
        F_Abs = double(F_Abs);
        F_Abs(mask) = NaN;
        F_Abs(F_Abs == 0) = NaN;
    end

    F_Snow = local_get_map(Maps, {'Hydro','Snowpack'}, k);
    if ~isempty(F_Snow)
        F_Snow = double(F_Snow);
        F_Snow(mask) = NaN;
        F_Snow(F_Snow == 0) = NaN;
    end

    % =====================================================================
    % 6) HEADER
    % =====================================================================
    DASH.mainTitle.String = sprintf(['HydroPol2D Runtime Dashboard --- Iteration / Layer = %d --- ' ...
                                     '$\\Delta t = %.3f\\ \\mathrm{s}$'], ...
                                     k, double(time_step));

    % =====================================================================
    % 7) UPDATE PANELS (USING ORIGINAL COLORAMP LOGIC)
    % =====================================================================
    local_update_panel(DASH.P.depth, F_d,   Spectrum,     S.labels.depth,      DASH.mask);
    local_update_panel(DASH.P.rain,  F_r,   WSE_RAS,      S.labels.rain,       DASH.mask);
    local_update_panel(DASH.P.vel,   F_v,   Velocity_RAS, S.labels.vel,        DASH.mask);

    local_update_panel(DASH.P.It,    F_i,   WSE_RAS,      S.labels.It,         DASH.mask);
    local_update_panel(DASH.P.f,     F_f,   Velocity_RAS, S.labels.f,          DASH.mask);
    local_update_panel(DASH.P.C,     F_C,   Spectrum,     S.labels.C,          DASH.mask);

    local_update_panel(DASH.P.GW,    F_GW,  Depth_RAS,    S.labels.GW,         DASH.mask);
    local_update_panel(DASH.P.ETR,   F_ETR, Spectrum,     S.labels.ETR,        DASH.mask);
    local_update_panel(DASH.P.ETP,   F_ETP, Spectrum,     S.labels.ETP,        DASH.mask);

    % =====================================================================
    % 8) FOOTER
    % =====================================================================
    DASH.footer.String = local_build_footer_text(gauges, F_Abs, F_Snow);

    drawnow limitrate nocallbacks;

    pause(0.05);

end

% =========================================================================
% DASHBOARD INITIALIZATION
% =========================================================================
function DASH = local_initialize_dashboard(DEM_raster, DEM, mask, Resolution, S)

    [nr, nc] = size(DEM);

    % ---------------------------------------------------------------------
    % Coordinates using original logic as much as possible
    % ---------------------------------------------------------------------
    [x_plot, y_plot, x_grid, y_grid, xlab, ylab] = local_build_coordinates(DEM_raster, Resolution, nr, nc);

    % ---------------------------------------------------------------------
    % Domain outline
    % ---------------------------------------------------------------------
    [outlineX, outlineY] = local_build_domain_outline(DEM_raster, DEM);

    % ---------------------------------------------------------------------
    % Figure
    % ---------------------------------------------------------------------
    fig = figure( ...
        'Name', 'HydroPol2D Runtime Dashboard', ...
        'NumberTitle', 'off', ...
        'Color', S.figureColor, ...
        'Units', 'pixels', ...
        'Position', S.figurePosition, ...
        'Renderer', S.renderer);

    tl = tiledlayout(fig, 3, 3, ...
        'TileSpacing', S.tileSpacing, ...
        'Padding', S.padding);

    mainTitle = title(tl, '', ...
        'Interpreter', 'latex', ...
        'FontName', S.fontName, ...
        'FontSize', S.titleFontSize, ...
        'FontWeight', S.titleFontWeight, ...
        'Color', S.titleColor);

    % ---------------------------------------------------------------------
    % Create panels
    % ---------------------------------------------------------------------
    P = struct();

    P.depth = local_create_panel(nexttile(tl,1), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.depth, xlab, ylab, outlineX, outlineY, S);

    P.rain  = local_create_panel(nexttile(tl,2), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.rain, xlab, ylab, outlineX, outlineY, S);

    P.vel   = local_create_panel(nexttile(tl,3), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.vel, xlab, ylab, outlineX, outlineY, S);

    P.It    = local_create_panel(nexttile(tl,4), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.It, xlab, ylab, outlineX, outlineY, S);

    P.f     = local_create_panel(nexttile(tl,5), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.f, xlab, ylab, outlineX, outlineY, S);

    P.C     = local_create_panel(nexttile(tl,6), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.C, xlab, ylab, outlineX, outlineY, S);

    P.GW    = local_create_panel(nexttile(tl,7), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.GW, xlab, ylab, outlineX, outlineY, S);

    P.ETR   = local_create_panel(nexttile(tl,8), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.ETR, xlab, ylab, outlineX, outlineY, S);

    P.ETP   = local_create_panel(nexttile(tl,9), x_plot, y_plot, x_grid, y_grid, nr, nc, ...
        S.labels.ETP, xlab, ylab, outlineX, outlineY, S);

    % ---------------------------------------------------------------------
    % Footer
    % ---------------------------------------------------------------------
    footer = annotation(fig, 'textbox', S.footerPosition, ...
        'String', '', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Interpreter', 'latex', ...
        'FontName', S.fontName, ...
        'FontSize', S.footerFontSize, ...
        'Color', S.footerColor);

    DASH = struct();
    DASH.fig       = fig;
    DASH.tl        = tl;
    DASH.mainTitle = mainTitle;
    DASH.footer    = footer;
    DASH.P         = P;

    DASH.nr     = nr;
    DASH.nc     = nc;
    DASH.mask   = mask;
    DASH.DEM    = DEM;
    DASH.x_plot = x_plot;
    DASH.y_plot = y_plot;
    DASH.x_grid = x_grid;
    DASH.y_grid = y_grid;
end

% =========================================================================
% CREATE ONE PANEL
% =========================================================================
function P = local_create_panel(axh, x_plot, y_plot, x_grid, y_grid, nr, nc, ttl, xlab, ylab, outlineX, outlineY, S)

    % Use surface instead of pcolor:
    % - keeps full matrix size
    % - faster updates
    % - avoids pcolor dropping last row/column
    [X, Y] = meshgrid(y_grid, x_grid);

    surfHandle = surface(axh, X, Y, zeros(nr,nc), nan(nr,nc), ...
        'EdgeColor', 'none', ...
        'FaceColor', 'flat');

    view(axh, 2);
    axis(axh, 'tight');
    axis(axh, 'equal');
    hold(axh, 'on');

    % Domain outline
    if ~isempty(outlineX)
        plot(axh, outlineX, outlineY, 'k-', 'LineWidth', S.domainLineWidth);
    end

    % Missing text
    txt = text(axh, 0.5, 0.5, 'Not available', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Interpreter', 'latex', ...
        'FontName', S.fontName, ...
        'FontSize', S.missingFontSize, ...
        'FontWeight', S.missingFontWeight, ...
        'Color', S.missingColor, ...
        'Visible', 'off');

    % Axes labels and title
    title(axh, ttl, ...
        'Interpreter', 'latex', ...
        'FontName', S.fontName, ...
        'FontSize', S.panelTitleFontSize, ...
        'FontWeight', S.panelTitleFontWeight);

    xlabel(axh, xlab, ...
        'Interpreter', 'latex', ...
        'FontName', S.fontName, ...
        'FontSize', S.labelFontSize);

    ylabel(axh, ylab, ...
        'Interpreter', 'latex', ...
        'FontName', S.fontName, ...
        'FontSize', S.labelFontSize);

    set(axh, ...
        'FontName', S.fontName, ...
        'FontSize', S.axisFontSize, ...
        'LineWidth', S.axisLineWidth, ...
        'Box', 'on', ...
        'TickDir', S.tickDir, ...
        'TickLength', S.tickLength, ...
        'Layer', 'top', ...
        'XColor', S.axisColor, ...
        'YColor', S.axisColor, ...
        'Color', S.axesColor);

    % Use original integer-looking labels
    try
        axh.XAxis.TickLabelFormat = '%.0f';
        axh.YAxis.TickLabelFormat = '%.0f';
        axh.XAxis.Exponent = 0;
        axh.YAxis.Exponent = 0;
        axh.YAxis.TickLabelRotation = 90;
    catch
    end

    cb = colorbar(axh, 'TickDirection', S.tickDir);
    cb.Box = 'on';
    cb.TickLength = S.colorbarTickLength;
    cb.LineWidth = S.colorbarLineWidth;
    cb.FontName = S.fontName;
    cb.FontSize = S.colorbarFontSize;
    cb.Color = S.axisColor;

    P = struct();
    P.ax   = axh;
    P.surf = surfHandle;
    P.cb   = cb;
    P.txt  = txt;
end

% =========================================================================
% UPDATE ONE PANEL
% =========================================================================
function local_update_panel(P, Z, cmap, ttl, mask)

    title(P.ax, ttl, 'Interpreter', 'latex');

    if isempty(Z)
        set(P.surf, 'CData', nan(size(mask)), 'AlphaData', zeros(size(mask)));
        colormap(P.ax, gray(16));
        caxis(P.ax, [0 1]);
        P.cb.Visible = 'off';
        P.txt.Visible = 'on';
        return;
    end

    Z(mask) = NaN;
    finiteVals = Z(isfinite(Z));

    if isempty(finiteVals)
        set(P.surf, 'CData', nan(size(Z)), 'AlphaData', zeros(size(Z)));
        colormap(P.ax, gray(16));
        caxis(P.ax, [0 1]);
        P.cb.Visible = 'off';
        P.txt.Visible = 'on';
        return;
    end

    set(P.surf, 'CData', Z, 'AlphaData', double(~isnan(Z)));
    colormap(P.ax, cmap);

    zmin = min(finiteVals(:));
    zmax = max(finiteVals(:));

    if zmin == zmax
        dz = max(abs(zmin)*1e-6, eps);
        zmin = zmin - dz;
        zmax = zmax + dz;
    end

    caxis(P.ax, [zmin zmax]);
    P.cb.Visible = 'on';
    P.txt.Visible = 'off';
end

% =========================================================================
% STYLE / LABELS
% =========================================================================
function S = local_style()

    S = struct();

    % Figure
    S.figureColor    = 'w';
    S.figurePosition = [20 20 1000 500];
    S.renderer       = 'painters';

    % Layout
    S.tileSpacing = 'compact';
    S.padding     = 'compact';

    % Typography
    S.fontName           = 'Garamond';
    S.titleFontSize      = 12;
    S.panelTitleFontSize = 12;
    S.labelFontSize      = 12;
    S.axisFontSize       = 12;
    S.colorbarFontSize   = 12;
    S.footerFontSize     = 12;
    S.missingFontSize    = 12;

    S.titleFontWeight      = 'bold';
    S.panelTitleFontWeight = 'bold';
    S.missingFontWeight    = 'bold';

    % Colors
    S.titleColor   = [0 0 0];
    S.axisColor    = [0 0 0];
    S.footerColor  = [0 0 0];
    S.missingColor = [0.25 0.25 0.25];
    S.axesColor    = 'w';

    % Lines / ticks
    S.axisLineWidth      = 1.5;
    S.tickDir            = 'out';
    S.tickLength         = [0.012 0.012];
    S.colorbarLineWidth  = 1.5;
    S.colorbarTickLength = 0.02;
    S.domainLineWidth    = 1.1;

    % Footer
    S.footerPosition = [0.005 0.002 0.99 0.032];

    % Panel labels
    % Edit these here if any unit needs to be changed.
    S.labels.depth = 'Water Depth [$m$]';
    S.labels.rain  = 'Rainfall Intensity [$mm\,h^{-1}$]';
    S.labels.vel   = 'Velocity [$m\,s^{-1}$]';

    S.labels.It    = 'Infiltration Storage $I_t$ [$mm$]';
    S.labels.f     = 'Infiltration Rate $f$ [$mm\,h^{-1}$]';
    S.labels.C     = 'State $C$ [$mm/h$]';

    S.labels.GW    = 'Groundwater Depth [$m$]';
    S.labels.ETR   = 'Actual ET ($ETR$) [$mm\,day^{-1}$]';
    S.labels.ETP   = 'Potential ET ($ETP$) [$mm\,day^{-1}$]';
end

% =========================================================================
% GATHER GPU IF NEEDED
% =========================================================================
function A = localGather(A)
    if isa(A, 'gpuArray')
        A = gather(A);
    end
end

% =========================================================================
% GET 2D/3D ARRAY
% =========================================================================
function A = local_get_array(Ain, k)

    A = localGather(Ain);

    if isempty(A)
        A = [];
        return;
    end

    if ndims(A) == 2
        return;
    elseif ndims(A) >= 3
        kk = min(max(1, k), size(A,3));
        A = A(:,:,kk);
    else
        A = [];
    end
end

% =========================================================================
% GET NESTED FIELD MAP
% =========================================================================
function A = local_get_map(S, pathCell, k)

    A = [];

    try
        ref = S;
        for ii = 1:numel(pathCell)
            if isstruct(ref) && isfield(ref, pathCell{ii})
                ref = ref.(pathCell{ii});
            else
                return;
            end
        end

        ref = localGather(ref);

        if isempty(ref)
            return;
        end

        if ndims(ref) == 2
            A = ref;
        elseif ndims(ref) >= 3
            kk = min(max(1, k), size(ref,3));
            A = ref(:,:,kk);
        end
    catch
        A = [];
    end
end

% =========================================================================
% BUILD COORDINATES
% =========================================================================
function [x_plot, y_plot, x_grid, y_grid, xlab, ylab] = local_build_coordinates(DEM_raster, Resolution, nr, nc)

    % Fallback
    x_plot = (0:nc-1) * Resolution / 1000;
    y_plot = ((nr-1):-1:0) * Resolution / 1000;

    % For surface/pcolor-like plotting, preserve original axis logic
    x_grid = ((nr):-1:1) * Resolution / 1000;
    y_grid = (1:nc) * Resolution / 1000;

    xlab = '$x\ [km]$';
    ylab = '$y\ [km]$';

    try
        if isfield(DEM_raster,'georef') && isfield(DEM_raster.georef,'SpatialRef')
            R = DEM_raster.georef.SpatialRef;

            if isprop(R,'XWorldLimits') && isprop(R,'YWorldLimits') && isprop(R,'CellExtentInWorldX')
                dx = double(R.CellExtentInWorldX);

                % Similar to original dashboard
                x_grid = double(R.YWorldLimits(2) - dx*(1:nr));
                y_grid = double(R.XWorldLimits(1) + dx*(1:nc));

                % convert to km when projected in meters
                if max(abs(x_grid)) > 1000 || max(abs(y_grid)) > 1000
                    x_grid = x_grid / 1000;
                    y_grid = y_grid / 1000;
                    xlab = '$x\ [km]$';
                    ylab = '$y\ [km]$';
                else
                    xlab = '$x$';
                    ylab = '$y$';
                end

                % imagesc-style helpers if ever needed
                x_plot = y_grid;
                y_plot = x_grid;
            end
        end
    catch
    end
end

% =========================================================================
% BUILD DOMAIN OUTLINE
% =========================================================================
function [outlineX, outlineY] = local_build_domain_outline(DEM_raster, DEM)

    outlineX = [];
    outlineY = [];

    try
        binaryMask = ~isnan(DEM);
        boundaries = bwboundaries(binaryMask);

        if isempty(boundaries)
            return;
        end

        combinedX = [];
        combinedY = [];

        for kk = 1:numel(boundaries)
            boundary = boundaries{kk};
            X = boundary(:,2);
            Y = boundary(:,1);

            combinedX = [combinedX; X; NaN]; %#ok<AGROW>
            combinedY = [combinedY; Y; NaN]; %#ok<AGROW>
        end

        if ~isempty(combinedX)
            combinedX = combinedX(1:end-1);
            combinedY = combinedY(1:end-1);
        end

        if isfield(DEM_raster,'georef') && isfield(DEM_raster.georef,'SpatialRef')
            R = DEM_raster.georef.SpatialRef;

            if isprop(R,'XWorldLimits') && isprop(R,'YWorldLimits') && isprop(R,'CellExtentInWorldX')
                dx = double(R.CellExtentInWorldX);

                outlineX = (double(R.XWorldLimits(1)) + combinedX * dx - dx/2)';
                outlineY = (double(R.YWorldLimits(2)) - combinedY * dx + dx/2)';

                if max(abs(outlineX)) > 1000 || max(abs(outlineY)) > 1000
                    outlineX = outlineX / 1000;
                    outlineY = outlineY / 1000;
                end
            end
        end
    catch
        outlineX = [];
        outlineY = [];
    end
end

% =========================================================================
% FOOTER
% =========================================================================
function txt = local_build_footer_text(gauges, F_Abs, F_Snow)

    ng = 0;
    try
        if ~isempty(gauges)
            if isstruct(gauges)
                if isfield(gauges, 'coordinates')
                    ng = size(gauges.coordinates,1);
                elseif isfield(gauges, 'x')
                    ng = numel(gauges.x);
                elseif isfield(gauges, 'labels_observed_string')
                    ng = numel(gauges.labels_observed_string);
                end
            end
        end
    catch
        ng = 0;
    end

    hasAbs  = ~isempty(F_Abs);
    hasSnow = ~isempty(F_Snow);

    txt = sprintf(['Gauges detected: %d \\quad | \\quad Optional fields available ' ...
                   '\\rightarrow \\quad Abstraction: %s \\quad | \\quad Snowpack: %s'], ...
                   ng, local_yesno(hasAbs), local_yesno(hasSnow));
end

function s = local_yesno(tf)
    if tf
        s = 'yes';
    else
        s = 'no';
    end
end