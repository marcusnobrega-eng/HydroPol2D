% HydroPol2D Snapshot
% Developer: Marcus Nobrega
% Updated version following the same style as the original code
% Main updates:
%   - NO LaTeX interpreters
%   - x and y shown in km
%   - safer color limits
%   - publication-style ticks/colorbars
%   - same subplot logic and visual structure as original script

if k > 1
    close all
    clf

    % ---------------------------------------------------------------------
    % Coordinates in km
    % ---------------------------------------------------------------------
    x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution * (1:size(DEM_raster.Z,2));
    y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution * (1:size(DEM_raster.Z,1));

    x_grid = x_grid ./ 1000;   % km
    y_grid = y_grid ./ 1000;   % km

    [X_grid,Y_grid] = meshgrid(x_grid,y_grid);

    set(gcf,'units','inches','position',[2,0,10,8])
    sgtitle(['Elapsed Time = ' num2str(t/60/24,'%.3f') ' days'], ...
        'Interpreter','none','FontName','Garamond','FontSize',13,'FontWeight','bold');

    % ---------------------------------------------------------------------
    % Creating the custom basemap
    % NOTE:
    % Since axes are now in km, the old mapshow basemap is NOT plotted
    % directly underneath the data because it would be in projected coords.
    % To keep this script stable and publication-ready, we keep the boundary
    % extraction logic but skip the actual basemap underlay in km coordinates.
    % ---------------------------------------------------------------------
    web_mercator_crs = projcrs(3857);

    no_plot = 1;
    S_p = [];

    if DEM_raster.georef.SpatialRef.ProjectedCRS.Name == web_mercator_crs.Name
        binaryMask = ~isnan(DEM_raster.Z);
        boundaries = bwboundaries(binaryMask);

        combinedX = [];
        combinedY = [];

        for kk = 1:numel(boundaries)
            boundary = boundaries{kk};
            X = boundary(:, 2);
            Y = boundary(:, 1);
            combinedX = [combinedX; X; NaN];
            combinedY = [combinedY; Y; NaN];
        end

        if ~isempty(combinedX)
            combinedX = combinedX(1:end-1);
            combinedY = combinedY(1:end-1);
        end

        S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
        S_p.BoundingBox = [ ...
            DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); ...
            DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)];

        S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + ...
            combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)' ./ 1000;

        S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - ...
            combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + ...
            DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)' ./ 1000;

        no_plot = 0;
    end

    % =====================================================================
    % ------------------------------ ETP -----------------------------------
    % =====================================================================
    t_title = 'ETP';
    ax1 = subplot(3,2,1);
    hold on
    axis tight; grid on; box on;

    if no_plot == 0
        try
            plot(S_p.X, S_p.Y, 'k', 'LineWidth', 0.8);
        catch
        end
    end

    z = gather(Hydro_States.ETP);
    z(idx_nan) = nan;
    F = z;

    map = surf_plot(t, 'ETP', '$\mathrm{mm day^-1', F, X_grid, Y_grid);
    set(map,'LineStyle','none');

    title(t_title,'Interpreter','none','FontSize',12,'FontName','Garamond')
    xlabel('x (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    ylabel('y (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    format_axes_publication(gca)

    % =====================================================================
    % ----------------------- Infiltration Rate ----------------------------
    % =====================================================================
    ax2 = subplot(3,2,2);
    hold on
    axis tight; grid on; box on;

    if no_plot == 0
        try
            plot(S_p.X, S_p.Y, 'k', 'LineWidth', 0.8);
        catch
        end
    end

    t_title = 'Infiltration Rate';

    if flags.flag_infiltration == 1
        z = gather(Hydro_States.f);
    else
        z = zeros(size(idx_nan));
    end

    z(idx_nan) = nan;
    z(z < 0) = nan;
    z(isinf(z)) = nan;
    F = z;

    map = surf_plot(t, 'f', 'mm h^-1', F, X_grid, Y_grid);
    set(map,'LineStyle','none');

    title(t_title,'Interpreter','none','FontSize',12,'FontName','Garamond')
    xlabel('x (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    ylabel('y (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    format_axes_publication(gca)

    % =====================================================================
    % ----------------------------- Depths ---------------------------------
    % =====================================================================
    ax3 = subplot(3,2,3);
    hold on
    axis tight; grid on; box on;

    if no_plot == 0
        try
            plot(S_p.X, S_p.Y, 'k', 'LineWidth', 0.8);
        catch
        end
    end

    t_title = 'Depths';

    z = gather(depths.d_t/1000);   % m
    z(idx_nan) = nan;
    F = z;

    map = surf_plot(t, 'd', 'm', F, X_grid, Y_grid);
    set(map,'LineStyle','none');

    title(t_title,'Interpreter','none','FontSize',12,'FontName','Garamond')
    xlabel('x (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    ylabel('y (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    format_axes_publication(gca)

    % =====================================================================
    % ---------------------------- Velocity --------------------------------
    % =====================================================================
    ax4 = subplot(3,2,4);
    hold on
    axis tight; grid on; box on;

    if no_plot == 0
        try
            plot(S_p.X, S_p.Y, 'k', 'LineWidth', 0.8);
        catch
        end
    end

    t_title = 'Velocity Raster';

    z = gather(velocities.velocity_raster);
    z(idx_nan) = nan;
    F = z;

    map = surf_plot(t, 'v', 'm s^-1', F, X_grid, Y_grid);
    set(map,'LineStyle','none');

    title(t_title,'Interpreter','none','FontSize',12,'FontName','Garamond')
    xlabel('x (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    ylabel('y (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    format_axes_publication(gca)

    % =====================================================================
    % ------------------------ Rainfall Intensity --------------------------
    % =====================================================================
    ax5 = subplot(3,2,5);
    hold on
    axis tight; grid on; box on;

    if no_plot == 0
        try
            plot(S_p.X, S_p.Y, 'k', 'LineWidth', 0.8);
        catch
        end
    end

    t_title = 'Rainfall Intensity';

    if flags.flag_rainfall == 1
        if flags.flag_alternated_blocks == 1 || flags.flag_huff == 1 || flags.flag_spatial_rainfall == 0
            z = gather(BC_States.delta_p_agg) * ones(size(idx_nan)) / (time_step/60);
        else
            z = gather(BC_States.delta_p_agg / (time_step/60));
        end
    else
        z = zeros(size(idx_nan));
    end

    z(idx_nan) = nan;
    z(~isfinite(z)) = nan;
    F = z;

    map = surf_plot(t, 'i', 'mm h^-1', F, X_grid, Y_grid);
    set(map,'LineStyle','none');

    title(t_title,'Interpreter','none','FontSize',12,'FontName','Garamond')
    xlabel('x (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    ylabel('y (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    format_axes_publication(gca)

    % =====================================================================
    % ------------------- Cumulative Infiltration -------------------------
    % =====================================================================
    ax6 = subplot(3,2,6);
    hold on
    axis tight; grid on; box on;

    if no_plot == 0
        try
            plot(S_p.X, S_p.Y, 'k', 'LineWidth', 0.8);
        catch
        end
    end

    t_title = 'Cumulative Infiltration';

    z = gather(Soil_Properties.I_t);
    z(idx_nan) = nan;
    z(~isfinite(z)) = nan;
    F = z;

    map = surf_plot(t, 'I_t', 'mm', F, X_grid, Y_grid);
    set(map,'LineStyle','none');

    title(t_title,'Interpreter','none','FontSize',12,'FontName','Garamond')
    xlabel('x (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    ylabel('y (km)','Interpreter','none','FontSize',11,'FontName','Garamond')
    format_axes_publication(gca)

end


% =========================================================================
% LOCAL FUNCTION: axes formatting
% =========================================================================
function format_axes_publication(ax)

    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;

    xtickformat(ax,'%.1f');
    ytickformat(ax,'%.1f');

    set(ax, 'FontName', 'Garamond', 'FontSize', 12)
    set(ax, 'TickLength', [0.02 0.01]);
    set(ax, 'TickDir', 'out')
    set(ax, 'LineWidth', 1.0)

    grid(ax,'on')
    box(ax,'on')
    axis(ax,'tight')
end