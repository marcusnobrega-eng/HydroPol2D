%% Phase 1 V-tilted synthetic catchment generator
% Generates static rasters for the shared Phase 1 validation domain.
%
% Run from the repository root:
%   run('HydroPol2D_Model/Validation/Phase1_VTilted_Catchment/generate_phase1_vtilted_domain.m')

clear; clc;

case_dir = fileparts(mfilename('fullpath'));
static_dir = fullfile(case_dir, 'Static');
config_file = fullfile(case_dir, 'Config', 'Domain_Config.csv');

if ~exist(static_dir, 'dir')
    mkdir(static_dir);
end

cfg = readtable(config_file, 'TextType', 'string');
getv = @(name) cfg.value(cfg.parameter == name);
getn = @(name) str2double(getv(name));

dx = getn('dx');
lengthY = getn('length_y');
leftWidth = getn('left_width');
channelWidth = getn('channel_width');
rightWidth = getn('right_width');
leftSlope = getn('left_slope');
rightSlope = getn('right_slope');
longSlope = getn('longitudinal_slope');
z_outlet = getn('z_outlet');
channelIncision = getn('channel_incision');

DTB_left = getn('dtb_left');
DTB_channel = getn('dtb_channel');
DTB_right = getn('dtb_right');

Albedo_left = getn('albedo_left');
Albedo_channel = getn('albedo_channel');
Albedo_right = getn('albedo_right');

LAI_left = getn('lai_left');
LAI_channel = getn('lai_channel');
LAI_right = getn('lai_right');

InitialSM_left = getn('initial_sm_left');
InitialSM_channel = getn('initial_sm_channel');
InitialSM_right = getn('initial_sm_right');

GWdepth_left = getn('gw_depth_left');
GWdepth_channel = getn('gw_depth_channel');
GWdepth_right = getn('gw_depth_right');
epsgCode = getn('epsg');

widthTotal = leftWidth + channelWidth + rightWidth;
nx = round(widthTotal / dx);
ny = round(lengthY / dx);

[X, Y] = meshgrid(1:nx, 1:ny);
x = (X - 1) * dx;
y = (Y - 1) * dx;

activeMask = true(ny, nx);

x_left_end = leftWidth;
x_channel_end = leftWidth + channelWidth;

leftMask = activeMask & x < x_left_end;
channelMask = activeMask & x >= x_left_end & x <= x_channel_end;
rightMask = activeMask & x > x_channel_end;

DEM = z_outlet + longSlope .* y;
channelBase = DEM - channelIncision;

DEM(leftMask) = channelBase(leftMask) + leftSlope .* (x_left_end - x(leftMask));
DEM(channelMask) = channelBase(channelMask);
DEM(rightMask) = channelBase(rightMask) + rightSlope .* (x(rightMask) - x_channel_end);

LULC = zeros(ny, nx);
SOIL = zeros(ny, nx);
Zone_ID = zeros(ny, nx);

LULC(leftMask) = 1;
LULC(channelMask) = 2;
LULC(rightMask) = 3;

SOIL(leftMask) = 1;
SOIL(channelMask) = 2;
SOIL(rightMask) = 3;

Zone_ID(leftMask) = 1;
Zone_ID(channelMask) = 2;
Zone_ID(rightMask) = 3;

DTB = zeros(ny, nx);
Albedo = zeros(ny, nx);
LAI = zeros(ny, nx);
Initial_SM = zeros(ny, nx);
GW_depth = zeros(ny, nx);

DTB(leftMask) = DTB_left;
DTB(channelMask) = DTB_channel;
DTB(rightMask) = DTB_right;

Albedo(leftMask) = Albedo_left;
Albedo(channelMask) = Albedo_channel;
Albedo(rightMask) = Albedo_right;

LAI(leftMask) = LAI_left;
LAI(channelMask) = LAI_channel;
LAI(rightMask) = LAI_right;

Initial_SM(leftMask) = InitialSM_left;
Initial_SM(channelMask) = InitialSM_channel;
Initial_SM(rightMask) = InitialSM_right;

GW_depth(leftMask) = GWdepth_left;
GW_depth(channelMask) = GWdepth_channel;
GW_depth(rightMask) = GWdepth_right;

GW_table = DEM - DTB + GW_depth;

rasterNames = {'DEM','LULC','SOIL','DTB','Albedo','LAI','Initial_SM','GW_depth','GW_table','Zone_ID'};
for k = 1:numel(rasterNames)
    A = eval(rasterNames{k});
    if any(isnan(A(:)))
        error('%s contains NaN values.', rasterNames{k});
    end
    if any(isinf(A(:)))
        error('%s contains Inf values.', rasterNames{k});
    end
end

x0 = 0;
y0 = lengthY;
R = maprefcells([x0, x0 + nx * dx], [y0 - ny * dx, y0], [ny, nx], ...
    'ColumnsStartFrom', 'north');

geotiffwrite(fullfile(static_dir, 'DEM.tif'), single(DEM), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'LULC.tif'), single(LULC), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'SOIL.tif'), single(SOIL), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'DTB.tif'), single(DTB), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'Albedo.tif'), single(Albedo), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'LAI.tif'), single(LAI), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'Initial_SM.tif'), single(Initial_SM), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'GW_depth.tif'), single(GW_depth), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'GW_table.tif'), single(GW_table), R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(static_dir, 'Zone_ID.tif'), single(Zone_ID), R, 'CoordRefSysCode', epsgCode);

summary = table( ...
    ["left_hillslope"; "channel_strip"; "right_hillslope"], ...
    [sum(leftMask(:)); sum(channelMask(:)); sum(rightMask(:))], ...
    [LAI_left; LAI_channel; LAI_right], ...
    [InitialSM_left; InitialSM_channel; InitialSM_right], ...
    [DTB_left; DTB_channel; DTB_right], ...
    'VariableNames', {'zone','cell_count','lai','initial_sm_mm','dtb_m'});
writetable(summary, fullfile(case_dir, 'Outputs', 'Validation', 'VTilted_Domain_Summary.csv'));

figure('Color', 'w', 'Position', [100 100 1600 850]);
subplot(2,4,1); imagesc(DEM); axis image off; colorbar; title('DEM [m]');
subplot(2,4,2); imagesc(Zone_ID); axis image off; colorbar; title('Zone ID');
subplot(2,4,3); imagesc(LULC); axis image off; colorbar; title('LULC');
subplot(2,4,4); imagesc(SOIL); axis image off; colorbar; title('SOIL');
subplot(2,4,5); imagesc(DTB); axis image off; colorbar; title('DTB [m]');
subplot(2,4,6); imagesc(Albedo); axis image off; colorbar; title('Albedo');
subplot(2,4,7); imagesc(LAI); axis image off; colorbar; title('LAI');
subplot(2,4,8); imagesc(GW_table); axis image off; colorbar; title('GW table [m]');
sgtitle('Phase 1 V-tilted Catchment Inputs');
exportgraphics(gcf, fullfile(static_dir, 'overview.png'), 'Resolution', 300);

fprintf('Phase 1 v-tilted catchment written to:\n%s\n', case_dir);
