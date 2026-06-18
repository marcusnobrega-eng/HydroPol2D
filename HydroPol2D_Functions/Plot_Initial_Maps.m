%% ========================================================================
%  HydroPol2D - Input / Parameter Raster Plotter (FULLY UPDATED)
%
%  Main improvements:
%    1) NaN values rendered transparent
%    2) Smarter axis units (m or km)
%    3) Better figure dimensions for export
%    4) Better categorical plotting
%    5) Smarter colormap selection by variable type
%    6) Fixed DEM stream overlay bug
%    7) More robust coordinate and overlay handling
%    8) Safer handling of masks / gauge labels / reservoir positions
%
%  Run AFTER HydroPol2D preprocessing, in the same workspace.
% ========================================================================

%% ========================= USER CONTROLS ================================

% -------------------------
% Export folder
% -------------------------
export_folder = Paths.FigPDF;   % leave [] or "" to choose manually

% -------------------------
% Figure visibility
% -------------------------
show_figures = true;   % true/false

% -------------------------
% Styling
% -------------------------
title_fs     = 12;   % subplot title fontsize
label_fs     = 12;   % axis label fontsize
tick_fs      = 12;    % axes tick fontsize
text_fs      = 12;   % fallback text fontsize
cb_fs        = 12;    % colorbar tick fontsize
suptitle_fs  = 14;   % figure supertitle fontsize

font_name    = 'Helvetica';
interp_latex = 'latex';

% -------------------------
% Colorbar controls
% -------------------------
cb_linewidth = 1.5;
cb_tickdir   = 'in';     % 'out' or 'in'
cb_location  = 'eastoutside';

% -------------------------
% Axes controls
% -------------------------
axes_linewidth = 1.5;
axes_tickdir   = 'out';
axis_equal_xy  = true;

% -------------------------
% Export controls
% -------------------------
pdf_resolution = 300;
paper_bg       = 'white';

% -------------------------
% Layout controls
% -------------------------
tile_spacing = 'compact';
tile_padding = 'compact';

% -------------------------
% Figure size controls [inches]
% -------------------------
figsize_2x2 = [6.5, 5];
figsize_2x3 = [10, 5];

% -------------------------
% Overlay controls
% -------------------------
stream_linewidth = 0.75;
stream_color     = [0.08 0.08 0.08];

gauge_markersize = 5;
gauge_linewidth  = 0.9;
reservoir_markersize = 6;

label_gauges           = true;   % set false to suppress labels
max_gauge_labels       = 8;      % if more than this, labels are thinned
gauge_label_every_n    = 2;      % used when too many gauges
gauge_text_offset_frac = 0.01;   % fraction of map range

% -------------------------
% Coordinate unit behavior
% -------------------------
use_km_threshold_m = 5000;  % if x or y extent exceeds this, plot in km

% -------------------------
% Stream plotting behavior
% -------------------------
prefer_stream_xy_overlay = true;   % preferred if STREAMobj2XY exists
fallback_stream_contour  = true;   % fallback to contour of STREAMMASK

%% ======================== CHECK EXPORT FOLDER ===========================

if isempty(export_folder) || strlength(string(export_folder)) == 0
    export_folder = uigetdir(pwd, 'Choose folder to export HydroPol2D raster figures');
    if isequal(export_folder, 0)
        error('No export folder selected. Script stopped.');
    end
end

if ~exist(export_folder, 'dir')
    mkdir(export_folder);
end

%% ======================== CHECK REQUIRED VARIABLES ======================

required_vars = {'DEM_raster','LULC_raster','SOIL_raster'};
for ii = 1:numel(required_vars)
    if ~exist(required_vars{ii}, 'var') && ...
            ~evalin('base', sprintf('exist(''%s'',''var'')', required_vars{ii}))
        error('Required variable "%s" was not found in the workspace.', required_vars{ii});
    end
end

DEM_raster  = getVarFromWorkspace('DEM_raster', []);
LULC_raster = getVarFromWorkspace('LULC_raster', []);
SOIL_raster = getVarFromWorkspace('SOIL_raster', []);

flags              = getVarFromWorkspace('flags', []);
LAI_raster         = getVarFromWorkspace('LAI_raster', []);
Albedo_raster      = getVarFromWorkspace('Albedo_raster', []);
DTB_raster         = getVarFromWorkspace('DTB_raster', []);
Wshed_Properties   = getVarFromWorkspace('Wshed_Properties', struct());
Soil_Properties    = getVarFromWorkspace('Soil_Properties', struct());
LULC_Properties    = getVarFromWorkspace('LULC_Properties', struct());
WQ_States          = getVarFromWorkspace('WQ_States', struct());
gauges             = getVarFromWorkspace('gauges', struct());
Reservoir_Data     = getVarFromWorkspace('Reservoir_Data', struct());
outlet_index       = getVarFromWorkspace('outlet_index', []);
Subgrid_Properties = getVarFromWorkspace('Subgrid_Properties', struct());
GIS_data           = getVarFromWorkspace('GIS_data', struct());
LULC_name          = getVarFromWorkspace('LULC_name', []);
SOIL_name          = getVarFromWorkspace('SOIL_name', []);
LULC_index         = getVarFromWorkspace('LULC_index', []);
SOIL_index         = getVarFromWorkspace('SOIL_index', []);

%% ======================== LOAD CUSTOM COLORMAPS =========================

try
    [Spectrum, depth_ramp, terrain_ramp, blue_ramp, blues_2, pallete, ...
        Depth_RAS, Terrain_RAS, Velocity_RAS, WSE_RAS] = coloramps();
catch
    warning('Could not call coloramps(). Falling back to MATLAB native colormaps.');
    Spectrum     = turbo(255);
    depth_ramp   = parula(255);
    terrain_ramp = parula(256);
    blue_ramp    = winter(100);
    blues_2      = parula(256);
    pallete      = struct();
    Depth_RAS    = parula(255);
    Terrain_RAS  = parula(255);
    Velocity_RAS = turbo(255);
    WSE_RAS      = turbo(255);
end

% --- Build semantic colormaps ---
CM = buildSemanticColormaps(Spectrum, depth_ramp, terrain_ramp, blue_ramp, ...
    blues_2, Depth_RAS, Terrain_RAS, Velocity_RAS, WSE_RAS);

%% ======================== BASIC MAPS ====================================

DEM  = double(DEM_raster.Z);
LULC = double(LULC_raster.Z);
SOIL = double(SOIL_raster.Z);

domain_mask = ~isnan(DEM);

if isprop(DEM_raster,'cellsize')
    dx = double(DEM_raster.cellsize);
elseif isfield(DEM_raster,'cellsize')
    dx = double(DEM_raster.cellsize);
elseif isfield(Wshed_Properties,'Resolution') && ~isempty(Wshed_Properties.Resolution)
    dx = double(Wshed_Properties.Resolution);
else
    dx = 1;
end

[nr, nc] = size(DEM);

% Prefer actual map coordinates if available
[xv_m, yv_m] = buildMapCoordinates(DEM_raster, GIS_data, dx, nr, nc);

% Smart unit conversion
coord_meta = chooseCoordinateUnits(xv_m, yv_m, use_km_threshold_m);
xv = coord_meta.xv_plot;
yv = coord_meta.yv_plot;

%% ======================== DEM-DERIVED MAPS ==============================

SLP = [];
try
    SLPobj = arcslope(DEM_raster);
    SLP = double(SLPobj.Z);
catch
    warning('Could not compute slope with arcslope. Using gradient fallback.');
    [gx, gy] = gradient(DEM, dx, dx);
    SLP = atand(sqrt(gx.^2 + gy.^2));
    SLP(~domain_mask) = nan;
end

FLOWACC    = [];
STREAMMASK = [];
S          = [];
FD         = [];

try
    FD = FLOWobj(DEM_raster);
    FA = flowacc(FD);
    FLOWACC = double(FA.Z);

    minarea_km2 = getNestedFieldSafe(GIS_data, 'min_area', []);
    if isempty(minarea_km2)
        if isfield(Wshed_Properties,'fac_area') && ~isempty(Wshed_Properties.fac_area)
            minarea_km2 = max(1, 0.01 * nanmax(Wshed_Properties.fac_area(:)));
        else
            minarea_km2 = 1;
        end
    end

    area_cells = max(1, round((minarea_km2 * 1e6) / (dx^2)));
    S = STREAMobj(FD, 'minarea', area_cells);

    if ~isempty(S)
        try
            STREAMMASK = false(size(DEM));
            STREAMMASK(S.IXgrid) = true;
        catch
            STREAMMASK = [];
        end
    end
catch
    warning('Could not compute flow accumulation / streams. These panels may be blank.');
end

%% ======================== PARAMETER MAPS ================================

LAI_map     = getRasterZ(LAI_raster);
ALB_map     = getRasterZ(Albedo_raster);
Manning_map = getNestedFieldSafe(LULC_Properties, 'roughness', []);

ksat_map      = getNestedFieldSafe(Soil_Properties, 'ksat', []);
ksat_gw_map   = getNestedFieldSafe(Soil_Properties, 'ksat_gw', []);
theta_sat_map = getNestedFieldSafe(Soil_Properties, 'theta_sat', []);
theta_r_map   = getNestedFieldSafe(Soil_Properties, 'theta_r', []);
theta_i_map   = getNestedFieldSafe(Soil_Properties, 'theta_i', []);
alpha_map     = getNestedFieldSafe(Soil_Properties, 'alpha_vg', []);
nvg_map       = getNestedFieldSafe(Soil_Properties, 'n_vg', []);
Sy_map        = getNestedFieldSafe(Soil_Properties, 'Sy', []);
DTB_map       = getNestedFieldSafe(Soil_Properties, 'Soil_Depth', []);

if isempty(DTB_map)
    DTB_map = getRasterZ(DTB_raster);
end

bedrock_map = [];
if ~isempty(DTB_map) && isequal(size(DTB_map), size(DEM))
    bedrock_map = DEM - DTB_map;
    bedrock_map(~domain_mask) = nan;
end

RiverWidth_map = getNestedFieldSafe(Wshed_Properties, 'River_Width', []);
RiverDepth_map = getNestedFieldSafe(Wshed_Properties, 'River_Depth', []);
inlet_mask     = getNestedFieldSafe(Wshed_Properties, 'inflow_mask', []);
perimeter_mask = getNestedFieldSafe(Wshed_Properties, 'perimeter', []);

C1_map = getNestedFieldSafe(LULC_Properties, 'C_1', []);
C2_map = getNestedFieldSafe(LULC_Properties, 'C_2', []);
C3_map = getNestedFieldSafe(LULC_Properties, 'C_3', []);
C4_map = getNestedFieldSafe(LULC_Properties, 'C_4', []);
B0_map = getNestedFieldSafe(WQ_States, 'B_0', []);

%% ======================== FIGURE 1 ======================================
fig1 = createStyledFigure(show_figures, paper_bg, figsize_2x2);
tl = tiledlayout(fig1, 2, 2, 'TileSpacing', tile_spacing, 'Padding', tile_padding);

nexttile;
plotRasterPanel(DEM, xv, yv, domain_mask, CM.terrain, ...
    'DEM', '$z$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(SLP, xv, yv, domain_mask, CM.slope, ...
    'Slope', '$S$ [m/m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, [0 maxFinite(SLP)]);

nexttile;
if ~isempty(FLOWACC)
    % Convert number of cells to contributing area [km^2]
    cell_area_km2 = (dx^2) / 1e6;

    FLOWACC_km2 = double(FLOWACC) * cell_area_km2;
    FLOWACC_km2(FLOWACC_km2 <= 0) = nan;
    FLOWACC_km2(~domain_mask) = nan;

    % Plot log10(area)
    FLOWACC_log = log10(FLOWACC_km2);

    plotRasterPanel(FLOWACC_log, xv, yv, domain_mask, CM.flowacc, ...
        'Flow accumulation', '$\log_{10}(A\,[\mathrm{km}^2])$ [-]', true, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, ...
        cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);
else
    plotMissingPanel('Flow accumulation not available', text_fs, font_name, interp_latex);
end

nexttile;
plotRasterPanel(DEM, xv, yv, domain_mask, CM.terrain, ...
    'Stream network over DEM', '$z$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);
hold on;
plotStreamOverlayTopo(S, STREAMMASK, coord_meta, stream_color, stream_linewidth, ...
    prefer_stream_xy_overlay, fallback_stream_contour);
hold off;

title(tl, 'Terrain and Hydrology', ...
    'Interpreter', interp_latex, 'FontSize', suptitle_fs, ...
    'FontName', font_name, 'FontWeight', 'bold');

addCRSLabel(fig1, DEM_raster, coord_meta, font_name);
exportFigurePDF(fig1, export_folder, 'Figure_1_Terrain_Hydrology.pdf', pdf_resolution);

%% ======================== FIGURE 2 ======================================
fig2 = createStyledFigure(show_figures, paper_bg, figsize_2x2);
tl = tiledlayout(fig2, 2, 2, 'TileSpacing', tile_spacing, 'Padding', tile_padding);

nexttile;
plotCategoricalPanel(LULC, xv, yv, domain_mask, ...
    'LULC', '', LULC_index, LULC_name, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex);

nexttile;
plotRasterPanel(LAI_map, xv, yv, domain_mask, CM.vegetation, ...
    'LAI', 'LAI [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, [0 maxFinite(LAI_map)]);

nexttile;
plotRasterPanel(ALB_map, xv, yv, domain_mask, CM.albedo, ...
    'Albedo', '$\alpha$ [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(Manning_map, xv, yv, domain_mask, CM.roughness, ...
    'Manning roughness', '$n$ [$\mathrm{s \cdot m^{-1/3}}$]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

title(tl, 'Surface Properties', ...
    'Interpreter', interp_latex, 'FontSize', suptitle_fs, ...
    'FontName', font_name, 'FontWeight', 'bold');

addCRSLabel(fig2, DEM_raster, coord_meta, font_name);
exportFigurePDF(fig2, export_folder, 'Figure_2_Surface_Properties.pdf', pdf_resolution);

%% ======================== FIGURE 3A =====================================
fig3a = createStyledFigure(show_figures, paper_bg, figsize_2x3);
tl = tiledlayout(fig3a, 2, 3, 'TileSpacing', tile_spacing, 'Padding', tile_padding);

nexttile;
plotCategoricalPanel(SOIL, xv, yv, domain_mask, ...
    'Soil classes', 'Soil class', SOIL_index, SOIL_name, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex);

nexttile;
plotRasterPanel(ksat_map, xv, yv, domain_mask, CM.conductivity, ...
    '$K_{sat}$', '$K_{sat}$ [mm h$^{-1}$]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(ksat_gw_map, xv, yv, domain_mask, CM.conductivity, ...
    '$K_{sat,gw}$', '$K_{sat,gw}$ [mm h$^{-1}$]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(theta_sat_map, xv, yv, domain_mask, CM.moisture, ...
    '$\theta_{sat}$', '$\theta_{sat}$ [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(theta_r_map, xv, yv, domain_mask, CM.moisture, ...
    '$\theta_r$', '$\theta_r$ [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(theta_i_map, xv, yv, domain_mask, CM.moisture, ...
    '$\theta_i$', '$\theta_i$ [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

%title(tl, 'Soil Hydraulic Properties', ...
%    'Interpreter', interp_latex, 'FontSize', suptitle_fs, ...
%    'FontName', font_name, 'FontWeight', 'bold');

addCRSLabel(fig3a, DEM_raster, coord_meta, font_name);
exportFigurePDF(fig3a, export_folder, 'Figure_3A_Soil_Hydraulic_Properties.pdf', pdf_resolution);

%% ======================== FIGURE 3B =====================================
fig3b = createStyledFigure(show_figures, paper_bg, figsize_2x3);
tl = tiledlayout(fig3b, 2, 3, 'TileSpacing', tile_spacing, 'Padding', tile_padding);

nexttile;
plotCategoricalPanel(SOIL, xv, yv, domain_mask, ...
    'Soil classes', 'Soil class', SOIL_index, SOIL_name, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex);

nexttile;
plotRasterPanel(alpha_map, xv, yv, domain_mask, CM.retention, ...
    '$\alpha_{vg}$', '$\alpha_{vg}$ [m$^{-1}$]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(nvg_map, xv, yv, domain_mask, CM.retention, ...
    '$n_{vg}$', '$n_{vg}$ [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(Sy_map, xv, yv, domain_mask, CM.storage, ...
    '$S_y$', '$S_y$ [-]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(DTB_map, xv, yv, domain_mask, CM.depth, ...
    'DTB / Soil depth', 'z [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(bedrock_map, xv, yv, domain_mask, CM.bedrock, ...
    'Bedrock elevation', '$z_b$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

%title(tl, 'Soil Retention and Storage Properties', ...
%    'Interpreter', interp_latex, 'FontSize', suptitle_fs, ...
%    'FontName', font_name, 'FontWeight', 'bold');

addCRSLabel(fig3b, DEM_raster, coord_meta, font_name);
exportFigurePDF(fig3b, export_folder, 'Figure_3B_Soil_Retention_Storage_Properties.pdf', pdf_resolution);

%% ======================== FIGURE 4 ======================================
fig4 = createStyledFigure(show_figures, paper_bg, figsize_2x3);
tl = tiledlayout(fig4, 2, 3, 'TileSpacing', tile_spacing, 'Padding', tile_padding);

nexttile;
plotRasterPanel(RiverWidth_map, xv, yv, domain_mask, CM.channel, ...
    'River width', '$B$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(RiverDepth_map, xv, yv, domain_mask, CM.depth, ...
    'River depth', '$H$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotRasterPanel(bedrock_map, xv, yv, domain_mask, CM.bedrock, ...
    'Bedrock elevation', '$z_b$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

nexttile;
plotMaskPanel(outlet_index, xv, yv, domain_mask, 'Outlet mask', coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex);

nexttile;
plotMaskPanel(inlet_mask, xv, yv, domain_mask, 'Inlet mask', coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex);

nexttile;
plotRasterPanel(DEM, xv, yv, domain_mask, CM.terrain, ...
    'Gauges / reservoirs', '$z$ [m]', true, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);
hold on;

plotGaugeAndReservoirOverlay(gauges, Reservoir_Data, coord_meta, ...
    gauge_markersize, gauge_linewidth, reservoir_markersize, ...
    label_gauges, max_gauge_labels, gauge_label_every_n, ...
    gauge_text_offset_frac, text_fs, font_name);
hold off;

% title(tl, 'Channel and Boundary Structure', ...
%     'Interpreter', interp_latex, 'FontSize', suptitle_fs, ...
%    'FontName', font_name, 'FontWeight', 'bold');

addCRSLabel(fig4, DEM_raster, coord_meta, font_name);
exportFigurePDF(fig4, export_folder, 'Figure_4_Channel_Boundary_Structure.pdf', pdf_resolution);

%% ======================== FIGURE 4B =====================================
fig4b = createStyledFigure(show_figures, paper_bg, [11.0, 8.5]);

ax4b = axes(fig4b, ...
    'Units', 'normalized', ...
    'Position', [0.08 0.10 0.78 0.80]);

plotRasterPanel(DEM, xv, yv, domain_mask, CM.terrain, ...
    'Gauges and reservoirs', '$z$ [m]', true, coord_meta, ...
    title_fs+1, label_fs, tick_fs, cb_fs, ...
    cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

hold(ax4b, 'on');

% --- Rivers ---
plotStreamOverlayTopo(S, STREAMMASK, coord_meta, ...
    stream_color, stream_linewidth+1, ...
    prefer_stream_xy_overlay, fallback_stream_contour);

% --- Gauges and reservoirs ---
plotGaugeAndReservoirOverlay(gauges, Reservoir_Data, coord_meta, ...
    gauge_markersize+1, gauge_linewidth, reservoir_markersize+1, ...
    true, inf, 1, ...
    gauge_text_offset_frac, text_fs, font_name);

hold(ax4b, 'off');

sgtitle(fig4b, 'Observation Gauges and Reservoirs', ...
    'Interpreter', interp_latex, ...
    'FontSize', suptitle_fs, ...
    'FontName', font_name, ...
    'FontWeight', 'bold');

addCRSLabel(fig4b, DEM_raster, coord_meta, font_name);
exportFigurePDF(fig4b, export_folder, 'Figure_4B_Gauges_Reservoirs.pdf', pdf_resolution);

%% ======================== FIGURE 5 ======================================
hasWQ = false;
if ~isempty(flags)
    try
        hasWQ = isfield(flags,'flag_waterquality') && (flags.flag_waterquality == 1);
    catch
        hasWQ = false;
    end
end

hasAnyWQMap = ~isempty(C1_map) || ~isempty(C2_map) || ~isempty(C3_map) || ...
    ~isempty(C4_map) || ~isempty(B0_map);

if hasWQ || hasAnyWQMap
    fig5 = createStyledFigure(show_figures, paper_bg, figsize_2x3);
    tl = tiledlayout(fig5, 2, 3, 'TileSpacing', tile_spacing, 'Padding', tile_padding);

    nexttile;
    plotCategoricalPanel(LULC, xv, yv, domain_mask, ...
        'LULC', '', LULC_index, LULC_name, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex);

    nexttile;
    plotRasterPanel(C1_map, xv, yv, domain_mask, CM.wq, ...
        '$C_1$', '$C_1$ [-]', true, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, ...
        cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

    nexttile;
    plotRasterPanel(C2_map, xv, yv, domain_mask, CM.wq, ...
        '$C_2$', '$C_2$ [-]', true, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, ...
        cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

    nexttile;
    plotRasterPanel(C3_map, xv, yv, domain_mask, CM.wq, ...
        '$C_3$', '$C_3$ [-]', true, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, ...
        cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

    nexttile;
    plotRasterPanel(C4_map, xv, yv, domain_mask, CM.wq, ...
        '$C_4$', '$C_4$ [-]', true, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, ...
        cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

    nexttile;
    plotRasterPanel(B0_map, xv, yv, domain_mask, CM.wqmass, ...
        '$B_0$', '$B_0$ [kg cell$^{-1}$]', true, coord_meta, ...
        title_fs, label_fs, tick_fs, cb_fs, ...
        cb_linewidth, cb_tickdir, cb_location, ...
        axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex, []);

    %title(tl, 'Water Quality Parameters', ...
    %    'Interpreter', interp_latex, 'FontSize', suptitle_fs, ...
    %   'FontName', font_name, 'FontWeight', 'bold');

    addCRSLabel(fig5, DEM_raster, coord_meta, font_name);
    exportFigurePDF(fig5, export_folder, 'Figure_5_Water_Quality_Parameters.pdf', pdf_resolution);
end

disp(' ');
disp('==============================================================');
disp('HydroPol2D input data and parametrization finished.');
disp(['PDF figures exported to: ' char(string(export_folder))]);
disp('==============================================================');

close all

%% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function val = getVarFromWorkspace(name, defaultVal)
if evalin('caller', sprintf('exist(''%s'',''var'')', name))
    val = evalin('caller', name);
elseif evalin('base', sprintf('exist(''%s'',''var'')', name))
    val = evalin('base', name);
else
    val = defaultVal;
end
end

function fig = createStyledFigure(show_figures, paper_bg, figsize_inches)

if show_figures
    vis = 'on';
else
    vis = 'off';
end

fig = figure( ...
    'Color', paper_bg, ...
    'Visible', vis, ...
    'WindowStyle', 'normal', ...
    'Units', 'inches', ...
    'Position', [1 1 figsize_inches(1) figsize_inches(2)], ...
    'Resize', 'off', ...
    'PaperUnits', 'inches', ...
    'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 figsize_inches(1) figsize_inches(2)], ...
    'PaperSize', [figsize_inches(1) figsize_inches(2)] );

drawnow;
set(fig, 'Units', 'inches');
pos = get(fig, 'Position');
pos(3) = figsize_inches(1);
pos(4) = figsize_inches(2);
set(fig, 'Position', pos);
drawnow;
end

function exportFigurePDF(fig, outfolder, fname, dpi_res)
drawnow;
exportgraphics(fig, fullfile(outfolder, fname), ...
    'ContentType', 'vector', ...
    'Resolution', dpi_res, ...
    'BackgroundColor', 'white');
end

function Z = getRasterZ(R)
Z = [];
if isempty(R), return; end
try
    if isobject(R) && isprop(R,'Z')
        Z = double(R.Z);
    elseif isstruct(R) && isfield(R,'Z')
        Z = double(R.Z);
    end
catch
    Z = [];
end
end

function out = getNestedFieldSafe(S, fieldname, defaultVal)
out = defaultVal;
try
    if isstruct(S) && isfield(S, fieldname)
        out = S.(fieldname);
    end
catch
    out = defaultVal;
end
end

function [xv_m, yv_m] = buildMapCoordinates(DEM_raster, GIS_data, dx, nr, nc)
% Prefer GIS_data corners if available
if isstruct(GIS_data) && isfield(GIS_data,'xulcorner') && isfield(GIS_data,'yulcorner') ...
        && ~isempty(GIS_data.xulcorner) && ~isempty(GIS_data.yulcorner)
    xv_m = double(GIS_data.xulcorner) + (0:nc-1) * dx;
    yv_m = double(GIS_data.yulcorner) - (0:nr-1) * dx;
    return;
end

% Fallback local coordinates
xv_m = (0:nc-1) * dx;
yv_m = (0:nr-1) * dx;
end

function coord_meta = chooseCoordinateUnits(xv_m, yv_m, threshold_m)
xrange = abs(max(xv_m) - min(xv_m));
yrange = abs(max(yv_m) - min(yv_m));

use_km = max(xrange, yrange) >= threshold_m;

if use_km
    scale = 1000;
    unit_lbl = 'km';
else
    scale = 1;
    unit_lbl = 'm';
end

coord_meta.use_km   = use_km;
coord_meta.scale    = scale;
coord_meta.unit_lbl = unit_lbl;
coord_meta.xv_m     = xv_m(:).';
coord_meta.yv_m     = yv_m(:).';
coord_meta.xv_plot  = coord_meta.xv_m ./ scale;
coord_meta.yv_plot  = coord_meta.yv_m ./ scale;
end

function CM = buildSemanticColormaps(Spectrum, depth_ramp, terrain_ramp, blue_ramp, ...
    blues_2, Depth_RAS, Terrain_RAS, Velocity_RAS, WSE_RAS)

CM = struct();

% Topography / elevations
CM.terrain      = Terrain_RAS;
CM.bedrock      = Terrain_RAS;

% Slopes / energetic gradients
CM.slope        = Spectrum;

% Contributing area / hydraulic importance
CM.flowacc      = flipud(gray(256));

% Vegetation / ecological intensity
CM.vegetation   = winter;

% Albedo / reflective surface property
CM.albedo       = gray;

% Roughness / drag
CM.roughness    = Velocity_RAS;

% Conductivity / transmissivity
CM.conductivity = Velocity_RAS;

% Water-content-like variables
CM.moisture     = parula;
CM.storage      = parula;

% Retention curve parameters
CM.retention    = Spectrum;

% Depth-like variables
CM.depth        = Depth_RAS;

% Channel geometry
CM.channel      = WSE_RAS;

% Water quality
CM.wq           = Spectrum;
CM.wqmass       = WSE_RAS;
end

function plotRasterPanel(Z, xv, yv, mask, cmap_in, ttl, cblab, addcb, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, ...
    cb_location, axes_linewidth, axes_tickdir, axis_equal_xy, ...
    font_name, interp_latex, clim_in)

if isempty(Z) || ~isequal(size(Z), size(mask))
    plotMissingPanel([ttl ' not available'], tick_fs+2, font_name, interp_latex);
    return;
end

Z = double(Z);
Z(~mask) = nan;

ax = gca;
cla(ax);
ax.Color = 'white';

h = imagesc(xv, yv, Z, 'Parent', ax);
set(ax, 'YDir', 'normal');
set(h, 'AlphaData', ~isnan(Z));  % transparent NaNs

if ~isempty(cmap_in)
    colormap(ax, cmap_in);
else
    colormap(ax, parula(256));
end

if ~isempty(clim_in) && isnumeric(clim_in) && numel(clim_in)==2 && all(isfinite(clim_in))
    clim_fixed = makeSafeCLim(clim_in(1), clim_in(2));
    clim(ax, clim_fixed);
else
    finiteVals = Z(isfinite(Z));
    if ~isempty(finiteVals)
        zmin = min(finiteVals);
        zmax = max(finiteVals);
        clim_fixed = makeSafeCLim(zmin, zmax);
        clim(ax, clim_fixed);
    end
end

axis(ax, 'tight');
if axis_equal_xy, axis(ax, 'equal'); end
box(ax, 'on');

ax.LineWidth = axes_linewidth;
ax.TickDir   = axes_tickdir;
ax.FontSize  = tick_fs;
ax.FontName  = font_name;
ax.TickLabelInterpreter = 'none';
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;

title(ax, ttl, 'Interpreter', interp_latex, 'FontSize', title_fs, ...
    'FontName', font_name, 'FontWeight', 'bold');

xlabel(ax, sprintf('$x$ [%s]', coord_meta.unit_lbl), ...
    'Interpreter', interp_latex, 'FontSize', label_fs, 'FontName', font_name);
ylabel(ax, sprintf('$y$ [%s]', coord_meta.unit_lbl), ...
    'Interpreter', interp_latex, 'FontSize', label_fs, 'FontName', font_name);

if addcb
    cb = colorbar(ax, cb_location);
    cb.TickLabelInterpreter = interp_latex;
    cb.FontSize = cb_fs;
    cb.LineWidth = cb_linewidth;
    cb.TickDirection = cb_tickdir;
    cb.Label.String = cblab;
    cb.Label.Interpreter = interp_latex;
    cb.Label.FontSize = label_fs;
    cb.Label.FontName = font_name;
end
end

function plotCategoricalPanel(Z, xv, yv, mask, ttl, cblab, class_index, class_names, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, ...
    cb_location, axes_linewidth, axes_tickdir, axis_equal_xy, ...
    font_name, interp_latex)

if isempty(Z) || ~isequal(size(Z), size(mask))
    plotMissingPanel([ttl ' not available'], tick_fs+2, font_name, interp_latex);
    return;
end

Z = double(Z);
Z(~mask) = nan;

% ---------------------------------------------------------
% Read class_index in the SAME order used in preprocessing
% ---------------------------------------------------------
ordered_codes = [];
if ~isempty(class_index)
    try
        if istable(class_index)
            ordered_codes = table2array(class_index(:,1));
        elseif iscell(class_index)
            ordered_codes = cell2mat(class_index(:,1));
        else
            ordered_codes = class_index(:);
        end
        ordered_codes = double(ordered_codes(:));
    catch
        ordered_codes = [];
    end
end

% Classes actually present in raster
raster_codes = unique(Z(isfinite(Z)));

% ---------------------------------------------------------
% Build plotting order:
% 1) follow class_index order
% 2) keep only codes present in raster
% 3) append any extra raster codes not listed in class_index
% ---------------------------------------------------------
valid = [];

if ~isempty(ordered_codes)
    for i = 1:numel(ordered_codes)
        if any(raster_codes == ordered_codes(i))
            valid(end+1,1) = ordered_codes(i); %#ok<AGROW>
        end
    end

    extras = setdiff(raster_codes, valid, 'stable');
    valid  = [valid; extras(:)];
else
    valid = raster_codes(:);
end

ncat = numel(valid);

if ncat == 0
    plotMissingPanel([ttl ' not available'], tick_fs+2, font_name, interp_latex);
    return;
end

% ---------------------------------------------------------
% Re-map actual class codes -> 1:ncat for plotting only
% ---------------------------------------------------------
Zcat = nan(size(Z));
for i = 1:ncat
    Zcat(Z == valid(i)) = i;
end

cmap = categoricalColormap(ncat);

ax = gca;
cla(ax);
ax.Color = 'white';

h = imagesc(xv, yv, Zcat, 'Parent', ax);
set(ax, 'YDir', 'normal');
set(h, 'AlphaData', ~isnan(Zcat));

colormap(ax, cmap);
clim(ax, [0.5 ncat+0.5]);

axis(ax, 'tight');
if axis_equal_xy, axis(ax, 'equal'); end
box(ax, 'on');

ax.LineWidth = axes_linewidth;
ax.TickDir   = axes_tickdir;
ax.FontSize  = tick_fs;
ax.FontName  = font_name;
ax.TickLabelInterpreter = 'none';
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;

title(ax, ttl, 'Interpreter', interp_latex, 'FontSize', title_fs, ...
    'FontName', font_name, 'FontWeight', 'bold');

xlabel(ax, sprintf('$x$ [%s]', coord_meta.unit_lbl), ...
    'Interpreter', interp_latex, 'FontSize', label_fs, 'FontName', font_name);
ylabel(ax, sprintf('$y$ [%s]', coord_meta.unit_lbl), ...
    'Interpreter', interp_latex, 'FontSize', label_fs, 'FontName', font_name);

% ---------------------------------------------------------
% Build tick labels from class_index -> class_names mapping
% ---------------------------------------------------------
tick_labels = strings(ncat,1);
for i = 1:ncat
    tick_labels(i) = string(valid(i)); % fallback
end

if ~isempty(class_index) && ~isempty(class_names)
    try
        % indices
        if istable(class_index)
            idx_vals = table2array(class_index(:,1));
        elseif iscell(class_index)
            idx_vals = cell2mat(class_index(:,1));
        else
            idx_vals = class_index(:);
        end
        idx_vals = double(idx_vals(:));

        % names
        if istable(class_names)
            name_vals = string(table2cell(class_names(:,1)));
        elseif iscell(class_names)
            name_vals = string(class_names(:,1));
        elseif isstring(class_names)
            name_vals = class_names(:);
        else
            name_vals = string(class_names);
        end
        name_vals = name_vals(:);

        nmap = min(numel(idx_vals), numel(name_vals));
        idx_vals  = idx_vals(1:nmap);
        name_vals = name_vals(1:nmap);

        for i = 1:ncat
            hit = find(idx_vals == valid(i), 1, 'first');
            if ~isempty(hit)
                tick_labels(i) = name_vals(hit);
            end
        end
    catch
        % keep fallback numeric labels
    end
end

cb = colorbar(ax, cb_location);
cb.TickLabelInterpreter = 'none';
cb.FontSize = 8;
cb.LineWidth = cb_linewidth;
cb.TickDirection = cb_tickdir;
cb.Ticks = 1:ncat;
cb.TickLabels = cellstr(tick_labels);

cb.Label.String = cblab;
cb.Label.Interpreter = interp_latex;
cb.Label.FontSize = label_fs;
cb.Label.FontName = font_name;
end

function cmap = categoricalColormap(ncat)
base = [
    0.1216    0.4667    0.7059
    1.0000    0.4980    0.0549
    0.1725    0.6275    0.1725
    0.8392    0.1529    0.1569
    0.5804    0.4039    0.7412
    0.5490    0.3373    0.2941
    0.8902    0.4667    0.7608
    0.4980    0.4980    0.4980
    0.7373    0.7412    0.1333
    0.0902    0.7451    0.8118
    0.6500    0.8078    0.8902
    0.9843    0.6039    0.6000
    ];

if ncat <= size(base,1)
    cmap = base(1:ncat,:);
else
    cmap = lines(ncat);
end
end

function plotMaskPanel(Z, xv, yv, mask, ttl, coord_meta, ...
    title_fs, label_fs, tick_fs, cb_fs, cb_linewidth, cb_tickdir, cb_location, ...
    axes_linewidth, axes_tickdir, axis_equal_xy, font_name, interp_latex)

if isempty(Z) || ~isequal(size(Z), size(mask))
    plotMissingPanel([ttl ' not available'], tick_fs+2, font_name, interp_latex);
    return;
end

Z = double(Z);
Z(~mask) = nan;
Z = double(Z > 0);

cmap = [0.95 0.95 0.95; 0.80 0.10 0.10];

ax = gca;
cla(ax);
ax.Color = 'white';

h = imagesc(xv, yv, Z, 'Parent', ax);
set(ax, 'YDir', 'normal');
set(h, 'AlphaData', ~isnan(Z));

colormap(ax, cmap);
clim(ax, [0 1]);

axis(ax, 'tight');
if axis_equal_xy, axis(ax, 'equal'); end
box(ax, 'on');

ax.LineWidth = axes_linewidth;
ax.TickDir   = axes_tickdir;
ax.FontSize  = tick_fs;
ax.FontName  = font_name;
ax.TickLabelInterpreter = 'none';
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;

title(ax, ttl, 'Interpreter', interp_latex, 'FontSize', title_fs, ...
    'FontName', font_name, 'FontWeight', 'bold');

xlabel(ax, sprintf('$x$ [%s]', coord_meta.unit_lbl), ...
    'Interpreter', interp_latex, 'FontSize', label_fs, 'FontName', font_name);
ylabel(ax, sprintf('$y$ [%s]', coord_meta.unit_lbl), ...
    'Interpreter', interp_latex, 'FontSize', label_fs, 'FontName', font_name);

cb = colorbar(ax, cb_location);
cb.TickLabelInterpreter = interp_latex;
cb.FontSize = cb_fs;
cb.LineWidth = cb_linewidth;
cb.TickDirection = cb_tickdir;
cb.Ticks = [0 1];
cb.TickLabels = {'0','1'};
cb.Label.String = 'Mask [-]';
cb.Label.Interpreter = interp_latex;
cb.Label.FontSize = label_fs;
cb.Label.FontName = font_name;
end

function plotMissingPanel(msg, text_fs, font_name, interp_latex)
axis off;
text(0.5, 0.5, msg, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Interpreter', interp_latex, ...
    'FontSize', text_fs, 'FontName', font_name);
end

function plotGaugeAndReservoirOverlay(gauges, Reservoir_Data, coord_meta, ...
    gauge_markersize, gauge_linewidth, reservoir_markersize, ...
    label_gauges, max_gauge_labels, gauge_label_every_n, ...
    gauge_text_offset_frac, text_fs, font_name)

xrange = max(coord_meta.xv_plot) - min(coord_meta.xv_plot);
yrange = max(coord_meta.yv_plot) - min(coord_meta.yv_plot);
dx_txt = gauge_text_offset_frac * max(xrange, eps);
dy_txt = gauge_text_offset_frac * max(yrange, eps);

% Gauges
try
    if isfield(gauges,'easting_obs_gauges_absolute') && isfield(gauges,'northing_obs_gauges_absolute') ...
            && ~isempty(gauges.easting_obs_gauges_absolute)

        xg = double(gauges.easting_obs_gauges_absolute) ./ coord_meta.scale;
        yg = double(gauges.northing_obs_gauges_absolute) ./ coord_meta.scale;

    elseif isfield(gauges,'easting_obs_gauges') && isfield(gauges,'northing_obs_gauges') ...
            && ~isempty(gauges.easting_obs_gauges)

        cellsize_plot = mean(diff(coord_meta.xv_plot));
        if isempty(cellsize_plot) || ~isfinite(cellsize_plot)
            cellsize_plot = 1;
        end

        xg = coord_meta.xv_plot(1) + (double(gauges.easting_obs_gauges)-1) * cellsize_plot;
        yg = coord_meta.yv_plot(1) + (double(gauges.northing_obs_gauges)-1) * cellsize_plot;
    else
        xg = [];
        yg = [];
    end

    if ~isempty(xg)
        plot(xg, yg, 'ko', ...
            'MarkerFaceColor', [1.00 0.95 0.10], ...
            'MarkerSize', gauge_markersize, ...
            'LineWidth', gauge_linewidth);

        if label_gauges && isfield(gauges,'labels_observed_string') && ~isempty(gauges.labels_observed_string)
            nlab = min(numel(gauges.labels_observed_string), numel(xg));

            if isinf(max_gauge_labels)
                idxlab = 1:nlab;
            elseif nlab > max_gauge_labels
                idxlab = 1:gauge_label_every_n:nlab;
            else
                idxlab = 1:nlab;
            end

            for i = idxlab
                txt = char(string(gauges.labels_observed_string{i}));
                text(xg(i)+dx_txt, yg(i)+dy_txt, txt, ...
                    'Interpreter', 'none', ...
                    'FontSize', max(text_fs-2,8), ...
                    'FontName', font_name, ...
                    'Color', 'k', ...
                    'BackgroundColor', [1 1 1], ...
                    'Margin', 1, ...
                    'Clipping', 'on');
            end
        end
    end
catch
end

% Reservoir
try
    if isfield(Reservoir_Data,'x_us') && isfield(Reservoir_Data,'y_us') ...
            && ~isempty(Reservoir_Data.x_us)
        xr = double(Reservoir_Data.x_us) ./ coord_meta.scale;
        yr = double(Reservoir_Data.y_us) ./ coord_meta.scale;
        plot(xr, yr, 's', ...
            'Color', [0.85 0.10 0.10], ...
            'MarkerFaceColor', [0.85 0.10 0.10], ...
            'MarkerSize', reservoir_markersize, ...
            'LineWidth', 1.0);

    elseif isfield(Reservoir_Data,'x_index') && isfield(Reservoir_Data,'y_index') ...
            && ~isempty(Reservoir_Data.x_index)
        cellsize_plot = mean(diff(coord_meta.xv_plot));
        if isempty(cellsize_plot) || ~isfinite(cellsize_plot)
            cellsize_plot = 1;
        end

        xr = coord_meta.xv_plot(1) + (double(Reservoir_Data.x_index)-1) * cellsize_plot;
        yr = coord_meta.yv_plot(1) + (double(Reservoir_Data.y_index)-1) * cellsize_plot;

        plot(xr, yr, 's', ...
            'Color', [0.85 0.10 0.10], ...
            'MarkerFaceColor', [0.85 0.10 0.10], ...
            'MarkerSize', reservoir_markersize, ...
            'LineWidth', 1.0);
    end
catch
end
end

function plotStreamOverlayTopo(S, STREAMMASK, coord_meta, stream_color, stream_linewidth, ...
    prefer_stream_xy_overlay, fallback_stream_contour)

% Preferred approach: use STREAMobj2XY because it plots real stream lines
if prefer_stream_xy_overlay && ~isempty(S)
    try
        [xs, ys] = STREAMobj2XY(S);
        for k = 1:numel(xs)
            plot(xs{k} ./ coord_meta.scale, ys{k} ./ coord_meta.scale, ...
                '-', 'Color', stream_color, 'LineWidth', stream_linewidth);
        end
        return;
    catch
        % continue to fallback
    end
end

% Fallback: contour the stream mask
if fallback_stream_contour && ~isempty(STREAMMASK)
    try
        contour(coord_meta.xv_plot, coord_meta.yv_plot, double(STREAMMASK), [1 1], ...
        'Color', stream_color, ...
        'LineWidth', stream_linewidth, ...
        'LineStyle', '-');
        return;
    catch
    end
end
end

function m = maxFinite(Z)
m = [];
if isempty(Z), return; end
zf = Z(isfinite(Z));
if isempty(zf), return; end
m = max(zf);
end

function crs_string = detectCRSString(DEM_raster, GIS_data)

crs_string = '';

% ---------------------------------------------------
% Attempt 1 — Modern MATLAB geospatial raster objects
% ---------------------------------------------------
try
    if isprop(DEM_raster,'SpatialRef') && ~isempty(DEM_raster.SpatialRef)
        R = DEM_raster.SpatialRef;

        if isprop(R,'ProjectedCRS') && ~isempty(R.ProjectedCRS)

            crs_name = R.ProjectedCRS.Name;

            if isprop(R.ProjectedCRS,'EPSGCode') && ~isempty(R.ProjectedCRS.EPSGCode)
                crs_string = sprintf('%s (EPSG:%d)', crs_name, R.ProjectedCRS.EPSGCode);
            else
                crs_string = crs_name;
            end

            return
        end
    end
catch
end

% ---------------------------------------------------
% Attempt 2 — TopoToolbox GRIDobj georef
% ---------------------------------------------------
try
    if isprop(DEM_raster,'georef') && isfield(DEM_raster.georef,'ProjectedCRS')
        crs_name = DEM_raster.georef.ProjectedCRS.Name;

        if isfield(DEM_raster.georef.ProjectedCRS,'EPSGCode')
            epsg = DEM_raster.georef.ProjectedCRS.EPSGCode;
            crs_string = sprintf('%s (EPSG:%d)', crs_name, epsg);
        else
            crs_string = crs_name;
        end

        return
    end
catch
end

% ---------------------------------------------------
% Attempt 3 — GIS_data metadata
% ---------------------------------------------------
try
    if isfield(GIS_data,'projection') && ~isempty(GIS_data.projection)
        crs_string = char(GIS_data.projection);
        return
    end
end

% ---------------------------------------------------
% Final fallback
% ---------------------------------------------------
crs_string = 'Unknown projected coordinate system';

end

function addCRSLabel(fig, DEM_raster, coord_meta, font_name)

% Extract CRS directly from DEM raster
try
    crs_name = DEM_raster.georef.SpatialRef.ProjectedCRS.Name;
catch
    crs_name = 'Unavailable';
end

% Unit text
if coord_meta.use_km
    unit_text = 'Axes in km';
else
    unit_text = 'Axes in m';
end

label_text = sprintf('CRS: %s', crs_name);

annotation(fig,'textbox', ...
    [0.55 0.952 0.43 0.035], ...
    'String', label_text, ...
    'EdgeColor','none', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',8, ...
    'FontName',font_name, ...
    'Interpreter','none');

end

function clim_out = makeSafeCLim(zmin, zmax)
% Ensures color limits are always valid for MATLAB clim().
% Handles constant maps by expanding limits around the constant value.

% fallback if something weird comes in
if ~isfinite(zmin) || ~isfinite(zmax)
    clim_out = [0 1];
    return;
end

% normal case
if zmax > zmin
    clim_out = [zmin zmax];
    return;
end

% constant-value case
z0 = zmin;

if z0 > 0
    % user-requested fallback: 90% / 110%
    clim_out = [0.9*z0, 1.1*z0];

elseif z0 < 0
    % for negative constant values, keep lower < upper
    clim_out = [1.1*z0, 0.9*z0];

else
    % z0 == 0, so 90% and 110% are still zero
    clim_out = [-1e-6, 1e-6];
end

% absolute safety check
if ~(isfinite(clim_out(1)) && isfinite(clim_out(2))) || clim_out(2) <= clim_out(1)
    delta = max(1e-6, eps(max(abs(z0),1)));
    clim_out = [z0 - delta, z0 + delta];
end
end
