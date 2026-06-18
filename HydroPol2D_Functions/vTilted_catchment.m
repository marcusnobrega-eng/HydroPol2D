%% ================================================================
%  SYNTHETIC V-TILTED CATCHMENT GENERATOR
%  HydroPol2D testing
%
%  Updated version:
%  - No NaN cells are included inside the active domain.
%  - Boundary cells are excluded from all zone masks.
%  - GeoTIFFs are written without NaN values.
%  - A domain mask is still created for clarity.
% ================================================================

clear; clc; close all;

%% ================= USER PARAMETERS =================

outDir = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/v-Tilted/Static';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

dx = 20;
lengthY = 1000;

leftWidth    = 800;
channelWidth = 20;
rightWidth   = 800;

leftSlope  = 0.05;
rightSlope = 0.05;
longSlope  = 0.02;

z_outlet = 0;
channelIncision = 0;

DTB_left = 1;
DTB_channel = 1;
DTB_right = 1;

Albedo_left = 0.5;
Albedo_channel = 0.5;
Albedo_right = 0.5;

LAI_left = 0.5;
LAI_channel = 0.5;
LAI_right = 0.5;

InitialSM_left = 50;
InitialSM_channel = 50;
InitialSM_right = 50;

GWdepth_left = 0.5;
GWdepth_channel = 0.5;
GWdepth_right = 0.5;

%% ================= DOMAIN =================

widthTotal = leftWidth + channelWidth + rightWidth;

nx = round(widthTotal / dx);
ny = round(lengthY / dx);

[X, Y] = meshgrid(1:nx, 1:ny);

x = (X-1) * dx;
y = (Y-1) * dx;

% Active model domain: no NaNs anywhere
activeMask = true(ny, nx);

%% ================= ZONES =================

x_left_end    = leftWidth;
x_channel_end = leftWidth + channelWidth;

leftMask    = activeMask & x < x_left_end;
channelMask = activeMask & x >= x_left_end & x <= x_channel_end;
rightMask   = activeMask & x > x_channel_end;

%% ================= DEM =================

DEM = z_outlet + longSlope .* y;

channelBase = DEM - channelIncision;

DEM(leftMask) = channelBase(leftMask) + ...
    leftSlope .* (x_left_end - x(leftMask));

DEM(channelMask) = channelBase(channelMask);

DEM(rightMask) = channelBase(rightMask) + ...
    rightSlope .* (x(rightMask) - x_channel_end);

%% ================= LULC / SOIL =================

LULC = zeros(ny,nx);
SOIL = zeros(ny,nx);

LULC(leftMask)    = 1;
LULC(channelMask) = 2;
LULC(rightMask)   = 3;

SOIL(leftMask)    = 1;
SOIL(channelMask) = 2;
SOIL(rightMask)   = 3;

%% ================= OTHER RASTERS =================

DTB        = zeros(ny,nx);
Albedo     = zeros(ny,nx);
LAI        = zeros(ny,nx);
Initial_SM = zeros(ny,nx);
GW_depth   = zeros(ny,nx);

DTB(leftMask)    = DTB_left;
DTB(channelMask) = DTB_channel;
DTB(rightMask)   = DTB_right;

Albedo(leftMask)    = Albedo_left;
Albedo(channelMask) = Albedo_channel;
Albedo(rightMask)   = Albedo_right;

LAI(leftMask)    = LAI_left;
LAI(channelMask) = LAI_channel;
LAI(rightMask)   = LAI_right;

Initial_SM(leftMask)    = InitialSM_left;
Initial_SM(channelMask) = InitialSM_channel;
Initial_SM(rightMask)   = InitialSM_right;

GW_depth(leftMask)    = GWdepth_left;
GW_depth(channelMask) = GWdepth_channel;
GW_depth(rightMask)   = GWdepth_right;

GW_table = DEM - DTB + GW_depth;

%% ================= SAFETY CHECK: REMOVE / BLOCK NaNs =================

rasterNames = {'DEM','LULC','SOIL','DTB','Albedo','LAI','Initial_SM','GW_depth','GW_table'};

for k = 1:numel(rasterNames)
    A = eval(rasterNames{k});

    if any(isnan(A(:)))
        error('%s still contains NaN values.', rasterNames{k});
    end

    if any(isinf(A(:)))
        error('%s contains Inf values.', rasterNames{k});
    end
end

%% ================= GEOREFERENCE =================

x0 = 0;
y0 = lengthY;

R = maprefcells( ...
    [x0, x0 + nx*dx], ...
    [y0 - ny*dx, y0], ...
    [ny, nx], ...
    'ColumnsStartFrom','north');

epsgCode = 3857;

%% ================= WRITE GEOTIFF FILES =================

geotiffwrite(fullfile(outDir,'DEM.tif'),        single(DEM),        R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'LULC.tif'),       single(LULC),       R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'SOIL.tif'),       single(SOIL),       R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'DTB.tif'),        single(DTB),        R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'Albedo.tif'),     single(Albedo),     R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'LAI.tif'),        single(LAI),        R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'Initial_SM.tif'), single(Initial_SM), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'GW_table.tif'),   single(GW_table),   R, 'CoordRefSysCode', epsgCode);

fprintf('Rasters saved to:\n%s\n', outDir);

%% ================= QUICK PLOT =================

figure('Color','w','Position',[100 100 1600 850]);

subplot(2,4,1)
imagesc(DEM); axis image off; colorbar;
title('DEM [m]')

subplot(2,4,2)
imagesc(LULC); axis image off; colorbar;
title('LULC: 1 left, 2 channel, 3 right')

subplot(2,4,3)
imagesc(SOIL); axis image off; colorbar;
title('SOIL: 1 left, 2 channel, 3 right')

subplot(2,4,4)
imagesc(DTB); axis image off; colorbar;
title('DTB [m]')

subplot(2,4,5)
imagesc(Albedo); axis image off; colorbar;
title('Albedo')

subplot(2,4,6)
imagesc(LAI); axis image off; colorbar;
title('LAI')

subplot(2,4,7)
imagesc(Initial_SM); axis image off; colorbar;
title('Initial\_SM [mm]')

subplot(2,4,8)
imagesc(GW_table); axis image off; colorbar;
title('GW\_table [m]')

sgtitle('Synthetic V-Tilted Catchment Inputs');

exportgraphics(gcf, fullfile(outDir,'overview.png'), 'Resolution',300);

fprintf('Overview figure saved to:\n%s\n', fullfile(outDir,'overview.png'));