% HydroPol2D Snapshot
% Developer: Marcus Nobrega
% Date 8/14/2023
% Goal: Plot Snapshot of the Results

if k > 1
    close all
    clf

    x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,2)]; y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,1)];
    filename = 'Input_Maps';
    set(gcf,'units','inches','position',[2,0,10,8])
    sgtitle(strcat('Elapsed Time =~',num2str(t/60/24),' days'),'interpreter','latex');

    %% Creating the custom basemap
    web_mercator_crs = projcrs(3857);

    if DEM_raster.georef.SpatialRef.ProjectedCRS.Name == web_mercator_crs.Name
        basemapName = "openstreetmap";
        url = "a.tile.openstreetmap.org/${z}/${x}/${y}.png";
        copyright = char(uint8(169));
        attribution = copyright + "OpenStreetMap contributors";
        addCustomBasemap(basemapName,url,"Attribution",attribution)
        % Getting lat and lon from the study area

        [lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
        latlim = [lat(1) lat(2)];
        lonlim = [lon(1) lon(2)];
        % Retriving the basemap image
        try
            [A,RA,attribA] = readBasemapImage(basemapName,latlim,lonlim);
        catch ME            
        end

        %% Creating a Shapefile from the DEM as reference

        % Creating a binary mask
        binaryMask = ~isnan(DEM_raster.Z);
        boundaries = bwboundaries(binaryMask);
        % Pre-allocate arrays to store combined X and Y coordinates
        combinedX = [];
        combinedY = [];

        % Combine all boundary coordinates into a single array
        for kk = 1:numel(boundaries)
            boundary = boundaries{kk};
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
        if DEM_raster.georef.SpatialRef.ProjectedCRS.Name ~= web_mercator_crs.Name;
            no_plot = 1;
        else
            S_p = struct('Geometry', 'Polygon', 'BoundingBox', [], 'X', [], 'Y', [], 'fid', 1, 'DN', 0);
            S_p.BoundingBox = [DEM_raster.georef.SpatialRef.XWorldLimits(1,1), DEM_raster.georef.SpatialRef.YWorldLimits(1,1); DEM_raster.georef.SpatialRef.XWorldLimits(1,2), DEM_raster.georef.SpatialRef.YWorldLimits(1,2)]; % Calculate bounding box for each polygon
            S_p.X = (DEM_raster.georef.SpatialRef.XWorldLimits(1)  + combinedX * DEM_raster.georef.SpatialRef.CellExtentInWorldX - DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
            S_p.Y = (DEM_raster.georef.SpatialRef.YWorldLimits(2)  - combinedY * DEM_raster.georef.SpatialRef.CellExtentInWorldX + DEM_raster.georef.SpatialRef.CellExtentInWorldX/2)';
        end
    else
        no_plot = 1;


    end

    % -------- ETP  -------- %
    t_title = 'ETP';
    ax1 = subplot(3,2,1);
    if no_plot==0;
        try
            ax1 = mapshow(A,RA,"AlphaData",0.45);hold on;
            ax1 = mapshow(S_p,'FaceColor','n'); hold on;
        catch ME
        end
    end
    axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    z = gather(Hydro_States.ETP); z(idx_nan) = nan;
    idx = z < 0;
    z(idx) = nan;
    idx = isinf(z);
    z(idx) = nan;
    xmax = size(z,2);
    xend = xmax;
    ymax = size(z,1);
    yend = ymax;
    h_min = min(min(z));
    F = z;
    zmax = max(max(z(~isnan(z))));
    if isempty(zmax) || isinf(zmax) || zmax == 0
        zmax = 0.1;
    end
    zmin = min(min(z(~isnan(z))));
    map = surf(x_grid,y_grid,F);
    set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    title((t_title),'Interpreter','Latex','FontSize',12)
    view(0,90)
    if h_min == zmax
        zmax = 2*h_min;
    end
    colormap(jet)
    hold on
    kk = colorbar ;
    caxis([zmin zmax]);
    ylabel(kk,'$\mathrm{ETP}$~($\mathrm{mm~day{-1}}$)','Interpreter','Latex','FontSize',12)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
    ax = ancestor(ax1, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');

    zlabel ('$\mathrm{ETP}$~($\mathrm{mm~day{-1}}$)','Interpreter','Latex','FontSize',12)
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
    % ---------- h_0 --------------- %

    ax2 = subplot(3,2,2);
    if no_plot==0;
        try
            ax2 = mapshow(A,RA,"AlphaData",0.45);hold on;
            ax2 = mapshow(S_p,'FaceColor','n'); hold on;
        catch ME

        end
    end
    t_title = 'Infiltration Rate';
    axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    if flags.flag_infiltration == 1
        z = gather(Hydro_States.f); z(idx_nan) = nan;
    else
        z = zeros(size(idx_nan)); z(idx_nan) = nan;
    end
    idx = z < 0;
    z(idx) = nan;
    idx = isinf(z);
    z(idx) = nan;
    xmax = size(z,2);
    xend = xmax;
    ymax = size(z,1);
    yend = ymax;
    h_min = min(min(z));
    F = z;
    zmax = max(max(z(~isnan(z))));
    if isempty(zmax) || isinf(zmax) || zmax == 0
        zmax = 0.1;
    end
    map = surf(x_grid,y_grid,F);
    set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    title((t_title),'Interpreter','Latex','FontSize',12)
    view(0,90)
    if h_min == zmax
        zmax = 2*h_min;
    end
    caxis([h_min zmax]);
    colormap(jet)
    hold on
    kk = colorbar ;
    ylabel(kk,'$f~(\mathrm{mm~h^{-1}})$','Interpreter','Latex','FontSize',12)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
    ax = ancestor(ax2, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');

    zlabel ('$f$ ($\mathrm{mm~h^{-1}})$)','Interpreter','Latex','FontSize',12)
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')

    % ----------  Depths ------------- %
    ax3 = subplot(3,2,3);
    if no_plot==0;
        try
            ax3 = mapshow(A,RA,"AlphaData",0.45);hold on;
            ax3 = mapshow(S_p,'FaceColor','n'); hold on;
        catch ME
        end
    end
    t_title = 'Depths';
    axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    z = gather(depths.d_t/1000); z(idx_nan) = nan;
    idx = z < 0;
    z(idx) = nan;
    idx = isinf(z);
    z(idx) = nan;
    xmax = size(z,2);
    xend = xmax;
    ymax = size(z,1);
    yend = ymax;
    h_min = min(min(z));
    F = z;
    zmax = max(max(z(~isnan(z))));
    if isempty(zmax) || isinf(zmax) || zmax == 0
        zmax = 0.1;
    end
    map = surf(x_grid,y_grid,F);
    set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    title((t_title),'Interpreter','Latex','FontSize',12)
    view(0,90)
    if h_min == zmax
        zmax = 2*h_min;
    end
    caxis([h_min zmax]);
    colormap(jet)
    hold on
    kk = colorbar ;
    ylabel(kk,'$d$ ($\mathrm{m})$','Interpreter','Latex','FontSize',12)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
    ax = ancestor(ax3, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');

    zlabel ('$d$ ($\mathrm{m})$)','Interpreter','Latex','FontSize',12)
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')

    % ----------  Velocity ------------- %
    ax4 = subplot(3,2,4);
    if no_plot==0;
        try
        ax4 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax4 = mapshow(S_p,'FaceColor','n'); hold on;
        catch ME
        end
    end
    t_title = 'Velocity Raster';
    axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    z = gather((velocities.velocity_raster)); z(idx_nan) = nan;
    idx = z < 0;
    z(idx) = nan;
    idx = isinf(z);
    z(idx) = nan;
    xmax = size(z,2);
    xend = xmax;
    ymax = size(z,1);
    yend = ymax;
    h_min = min(min(z));
    F = z;
    zmax = max(max(z(~isnan(z))));
    if isempty(zmax) || isinf(zmax) || zmax == 0
        zmax = 0.1;
    end
    map = surf(x_grid,y_grid,F);
    set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    title((t_title),'Interpreter','Latex','FontSize',12)
    view(0,90)
    if h_min == zmax
        zmax = 2*h_min;
    end
    caxis([h_min zmax]);
    colormap(jet)
    hold on
    kk = colorbar ;
    ylabel(kk,'$v$ ($\mathrm{m \cdot s^{-1}})$','Interpreter','Latex','FontSize',12)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
    ax = ancestor(ax4, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');

    zlabel ('$\Delta \theta$ ($\mathrm{cm^3.cm^{-3}}$)','Interpreter','Latex','FontSize',12)
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')

    % ---------- i --------------- %

    ax5 = subplot(3,2,5);
    if no_plot==0;
        try
        ax2 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax2 = mapshow(S_p,'FaceColor','n'); hold on;
        catch ME
        end
    end
    t_title = 'Rainfall Intensity';
    axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    if flags.flag_rainfall == 1
        z = gather(BC_States.delta_p_agg/(time_step/60)); z(idx_nan) = nan;
    else
        z = zeros(size(idx_nan)); z(idx_nan) = nan;
    end
    idx = z < 0;
    z(idx) = nan;
    idx = isinf(z);
    z(idx) = nan;
    xmax = size(z,2);
    xend = xmax;
    ymax = size(z,1);
    yend = ymax;
    h_min = min(min(z));
    F = z;
    zmax = max(max(z(~isnan(z))));
    if isempty(zmax) || isinf(zmax) || zmax == 0
        zmax = 0.1;
    end
    map = surf(x_grid,y_grid,F);
    set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    title((t_title),'Interpreter','Latex','FontSize',12)
    view(0,90)
    if h_min == zmax
        zmax = 2*h_min;
    end
    caxis([h_min zmax]);
    colormap(jet)
    hold on
    kk = colorbar ;
    ylabel(kk,'$i~(\mathrm{mm~h^{-1}})$','Interpreter','Latex','FontSize',12)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
    ax = ancestor(ax2, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');

    zlabel ('$i$ ($\mathrm{mm~h^{-1}})$)','Interpreter','Latex','FontSize',12)
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')

    % ---------- i --------------- %

    ax5 = subplot(3,2,6);
    if no_plot==0;
        try
        ax2 = mapshow(A,RA,"AlphaData",0.45);hold on;
        ax2 = mapshow(S_p,'FaceColor','n'); hold on;
        catch ME
        end
    end
    t_title = 'Cumulative Infiltration';
    axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    z = gather(Soil_Properties.I_t); z(idx_nan) = nan;
    idx = z < 0;
    z(idx) = nan;
    idx = isinf(z);
    z(idx) = nan;
    xmax = size(z,2);
    xend = xmax;
    ymax = size(z,1);
    yend = ymax;
    h_min = min(min(z));
    F = z;
    zmax = max(max(z(~isnan(z))));
    if isempty(zmax) || isinf(zmax) || zmax == 0
        zmax = 0.1;
    end
    map = surf(x_grid,y_grid,F);
    set(map,'LineStyle','none'); axis tight; grid on; box on; % this ensures that getframe() returns a consistent size; axis tight; grid on; box on; % this ensures that getframe() returns a consistent size
    title((t_title),'Interpreter','Latex','FontSize',12)
    view(0,90)
    if h_min == zmax
        zmax = 2*h_min;
    end
    caxis([h_min zmax]);
    colormap(jet)
    hold on
    kk = colorbar ;
    ylabel(kk,'$I_t~(\mathrm{mm})$','Interpreter','Latex','FontSize',12)
    xlabel(' x (m) ','Interpreter','Latex','FontSize',12)
    ylabel ('y (m) ','Interpreter','Latex','FontSize',12)
    ax = ancestor(ax2, 'axes');
    ax.XAxis.Exponent = 0;xtickformat('%.0f');
    ax.YAxis.Exponent = 0;ytickformat('%.0f');

    zlabel ('$I_t~(\mathrm{mm})$)','Interpreter','Latex','FontSize',12)
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)
    set(gca, 'TickLength', [0.02 0.01]);
    set(gca,'Tickdir','out')
end
