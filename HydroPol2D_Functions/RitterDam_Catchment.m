%% ================================================================
%  RITTER DAM-BREAK BENCHMARK GENERATOR
%  HydroPol2D local inertial validation
%
%  Physics:
%   - Flat bed
%   - Initial reservoir depth h0 on the left
%   - Dry bed on the right
%   - Negligible friction
%   - Compare depth profile h(x,t) against Ritter analytical solution
%
%  HydroPol2D required outputs:
%   - DEM.tif
%   - LULC.tif
%   - SOIL.tif
%   - DTB.tif
%   - Albedo.tif
%   - LAI.tif
%   - GW_table.tif
%   - Warmup_Depth.tif              [mm]
%   - Initial_Soil_Moisture.tif     [mm]
% ================================================================

clear; clc; close all;

%% ================= USER PARAMETERS =================

outDir = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/RitterDam/Static';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

dx = 5;             % grid resolution [m]

Lx = 1000;          % domain length [m]
Ly = 100;           % domain width [m]

h0 = 1.0;           % initial upstream water depth [m]
x_dam = 250;        % dam-break location [m]

z0 = 0;             % flat bed elevation [m]

% Use very small roughness to approximate frictionless Ritter solution.
% Do NOT use exactly zero if any model equation divides by roughness.
n_manning = 1e-6; %#ok<NASGU>

%% Analytical comparison times
compare_times_sec = [10 20 30 40 50 60];

g = 9.81;
c0 = sqrt(g*h0);

%% Uniform properties
DTB_val = 1;
Albedo_val = 0.5;
LAI_val = 0.5;
InitialSoilMoisture_val = 0;  % [mm]
GWdepth_val = 1;

%% ================= DOMAIN =================

nx = round(Lx / dx);
ny = round(Ly / dx);

[X, Y] = meshgrid(1:nx, 1:ny);

x = (X-1) * dx;
y = (Y-1) * dx; %#ok<NASGU>

%% ================= DEM =================

DEM = z0 * ones(ny,nx);
DEM(:,end) = z0 - 1e-8;
%% ================= INITIAL WATER DEPTH / WARMUP DEPTH =================

Initial_Depth_m = zeros(ny,nx);       % [m]
Initial_Depth_m(x <= x_dam) = h0;     % reservoir on left side

Initial_WSE = DEM + Initial_Depth_m;  % [m], exported only for diagnostic use

% HydroPol2D warmup depth must be in m
Warmup_Depth = Initial_Depth_m; % [m]

%% ================= UNIFORM RASTERS =================

LULC = ones(ny,nx);
SOIL = ones(ny,nx);

DTB        = DTB_val    * ones(ny,nx);
Albedo     = Albedo_val * ones(ny,nx);
LAI        = LAI_val    * ones(ny,nx);
GW_depth   = GWdepth_val * ones(ny,nx);

% HydroPol2D soil moisture raster
Initial_Soil_Moisture = InitialSoilMoisture_val * ones(ny,nx); % [mm]

GW_table = DEM - DTB + GW_depth;

%% ================= SAFETY CHECK =================

rasterNames = {'DEM','LULC','SOIL','DTB','Albedo','LAI', ...
               'GW_table','Warmup_Depth','Initial_Soil_Moisture', ...
               'Initial_Depth_m','Initial_WSE'};

for k = 1:numel(rasterNames)
    A = eval(rasterNames{k});
    if any(isnan(A(:)))
        error('%s contains NaN values.', rasterNames{k});
    end
    if any(isinf(A(:)))
        error('%s contains Inf values.', rasterNames{k});
    end
end

%% ================= GEOREFERENCE =================

x0 = 0;
y0 = Ly;

R = maprefcells( ...
    [x0, x0 + nx*dx], ...
    [y0 - ny*dx, y0], ...
    [ny, nx], ...
    'ColumnsStartFrom','north');

epsgCode = 3857;

%% ================= WRITE GEOTIFF FILES =================

geotiffwrite(fullfile(outDir,'DEM.tif'),                   single(DEM),                   R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'LULC.tif'),                  single(LULC),                  R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'SOIL.tif'),                  single(SOIL),                  R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'DTB.tif'),                   single(DTB),                   R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'Albedo.tif'),                single(Albedo),                R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'LAI.tif'),                   single(LAI),                   R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'GW_table.tif'),              single(GW_table),              R, 'CoordRefSysCode', epsgCode);

% HydroPol2D initial-condition rasters
geotiffwrite(fullfile(outDir,'Warmup_Depth.tif'),          single(Warmup_Depth),          R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'Initial_Soil_Moisture.tif'), single(Initial_Soil_Moisture), R, 'CoordRefSysCode', epsgCode);

% Optional diagnostic rasters
geotiffwrite(fullfile(outDir,'Initial_Depth_m.tif'),       single(Initial_Depth_m),       R, 'CoordRefSysCode', epsgCode);
geotiffwrite(fullfile(outDir,'Initial_WSE.tif'),           single(Initial_WSE),           R, 'CoordRefSysCode', epsgCode);

fprintf('Ritter rasters saved to:\n%s\n', outDir);
fprintf('Warmup depth saved as [mm]:\n%s\n', fullfile(outDir,'Warmup_Depth.tif'));
fprintf('Initial soil moisture saved as [mm]:\n%s\n', fullfile(outDir,'Initial_Soil_Moisture.tif'));

%% ================= ANALYTICAL RITTER SOLUTION =================

x_profile = ((0:nx-1)' * dx);
x_rel = x_profile - x_dam;

Analytical = table();
Analytical.x_m = x_profile;

for it = 1:numel(compare_times_sec)
    t = compare_times_sec(it);

    h = zeros(size(x_profile));
    u = zeros(size(x_profile));

    left_region = x_rel < -c0*t;
    fan_region  = x_rel >= -c0*t & x_rel <= 2*c0*t;

    h(left_region) = h0;
    u(left_region) = 0;

    h(fan_region) = (4/(9*g)) .* ...
        (c0 - x_rel(fan_region)./(2*t)).^2;

    u(fan_region) = (2/3) .* ...
        (c0 + x_rel(fan_region)./t);

    h(h < 0) = 0;

    Analytical.(sprintf('h_t%03ds_m',t))  = h;
    Analytical.(sprintf('u_t%03ds_ms',t)) = u;
end

writetable(Analytical, fullfile(outDir,'Ritter_Analytical_Profiles.csv'));

fprintf('Analytical Ritter profiles saved to:\n%s\n', ...
    fullfile(outDir,'Ritter_Analytical_Profiles.csv'));

%% ================= QUICK PLOT =================

figure('Color','w','Position',[100 100 1500 750]);

subplot(2,2,1)
imagesc(DEM); axis image off; colorbar;
title('DEM [m]')

subplot(2,2,2)
imagesc(Warmup_Depth); axis image off; colorbar;
title('Warmup\_Depth [mm]')

subplot(2,2,3)
plot(x_profile, Initial_Depth_m(round(ny/2),:)', 'k-', 'LineWidth', 2);
grid on;
xlabel('x [m]');
ylabel('Initial depth [m]');
title('Initial dam-break condition')

subplot(2,2,4)
hold on;
for it = 1:numel(compare_times_sec)
    t = compare_times_sec(it);
    plot(x_profile, Analytical.(sprintf('h_t%03ds_m',t)), 'LineWidth', 1.5);
end
grid on;
xlabel('x [m]');
ylabel('Depth h [m]');
title('Analytical Ritter depth profiles')
legend(compose('t = %d s', compare_times_sec), 'Location','best');

sgtitle('Ritter Dam-Break Benchmark');

exportgraphics(gcf, fullfile(outDir,'overview_ritter.png'), 'Resolution',300);

fprintf('Overview figure saved to:\n%s\n', fullfile(outDir,'overview_ritter.png'));