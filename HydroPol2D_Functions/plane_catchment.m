%% ================================================================
%  SYNTHETIC TILTED PLANE CATCHMENT + ANALYTICAL FLOW BENCHMARK
%  HydroPol2D testing
%
%  Purpose:
%   - Uniform slope
%   - Uniform roughness
%   - No channel
%   - Analytical kinematic-wave outlet hydrograph
% ================================================================

clear; clc; close all;

%% ================= USER PARAMETERS =================

outDir = '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Plane_Infiltration/Static';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

dx = 20;           % grid resolution [m]

Lx = 2000;         % length in flow direction [m]
Ly = 1000;         % width [m]

S0 = 0.01;         % slope in flow direction [-]
n_manning = 0.015;  % Manning roughness [s/m^(1/3)]

z_outlet = 0;      % outlet elevation

%% Rainfall benchmark
rain_mm_h = 100;        % rainfall excess [mm/h]
rain_duration_min = 240;
sim_duration_min  = 240;
dt_benchmark_sec  = 10*60;

%% Uniform properties
DTB_val = 1;
Albedo_val = 0.5;
LAI_val = 0.5;
InitialSM_val = 0;
GWdepth_val = 0;

%% ================= DOMAIN =================

nx = round(Ly / dx);
ny = round(Lx / dx);

[X, Y] = meshgrid(1:nx, 1:ny);

x = (X-1) * dx;
y = (Y-1) * dx;

%% ================= DEM =================

% Flow direction: TOP -> BOTTOM (y increasing)
DEM = z_outlet + S0 * y;

%% ================= UNIFORM FIELDS =================

LULC = ones(ny,nx);
SOIL = ones(ny,nx);

DTB        = DTB_val        * ones(ny,nx);
Albedo     = Albedo_val     * ones(ny,nx);
LAI        = LAI_val        * ones(ny,nx);
Initial_SM = InitialSM_val  * ones(ny,nx);
GW_depth   = GWdepth_val    * ones(ny,nx);

GW_table = DEM - DTB + GW_depth;

%% ================= ANALYTICAL KINEMATIC-WAVE BENCHMARK =================
% Unit-width Manning relation:
%   q = alpha h^(5/3)
% where:
%   alpha = sqrt(S0)/n
%
% For constant rainfall excess r:
%   rising limb: q_out(t) = alpha * (r*t)^(5/3)
%   equilibrium: q_eq = r * Lx
%
% Total discharge:
%   Q_out = q_out * Ly

r = rain_mm_h / 1000 / 3600;       % [m/s]
rain_duration_sec = rain_duration_min * 60;
sim_duration_sec  = sim_duration_min  * 60;

alpha = sqrt(S0) / n_manning;

% Time to equilibrium / time of concentration [s]
t_eq_sec = (Lx / (alpha * r^(2/3)))^(3/5);

% Steady unit-width discharge [m2/s]
q_eq = r * Lx;

% Total steady discharge [m3/s]
Q_eq = q_eq * Ly;

t_sec = (0:dt_benchmark_sec:sim_duration_sec)';
Q_analytic = zeros(size(t_sec));

for ii = 1:numel(t_sec)
    tt = t_sec(ii);

    if tt <= rain_duration_sec
        % Rising limb during rainfall
        if tt < t_eq_sec
            q_unit = alpha * (r * tt)^(5/3);
        else
            q_unit = q_eq;
        end

    else
        % Simple recession approximation for this first benchmark:
        % after rainfall stops, no new rainfall enters.
        % For now, keep NaN because exact recession requires characteristic solution.
        q_unit = NaN;
    end

    Q_analytic(ii) = q_unit * Ly;
end

Benchmark = table();
Benchmark.time_sec = t_sec;
Benchmark.time_min = t_sec / 60;
Benchmark.Q_analytic_m3s = Q_analytic;

writetable(Benchmark, fullfile(outDir,'Analytical_KinematicWave_Hydrograph.csv'));

fprintf('\n================ ANALYTICAL BENCHMARK ================\n');
fprintf('Rainfall excess        = %.3f mm/h\n', rain_mm_h);
fprintf('Plane length Lx        = %.2f m\n', Lx);
fprintf('Plane width Ly         = %.2f m\n', Ly);
fprintf('Slope S0               = %.4f\n', S0);
fprintf('Manning n              = %.4f s/m^(1/3)\n', n_manning);
fprintf('Time to equilibrium    = %.2f min\n', t_eq_sec/60);
fprintf('Steady outlet flow Qeq = %.6f m3/s\n', Q_eq);
fprintf('Benchmark CSV saved to:\n%s\n', fullfile(outDir,'Analytical_KinematicWave_Hydrograph.csv'));
fprintf('======================================================\n\n');

%% ================= SAFETY CHECK =================

rasterNames = {'DEM','LULC','SOIL','DTB','Albedo','LAI','Initial_SM','GW_depth','GW_table'};

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
y0 = Lx;

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

fprintf('Plane rasters saved to:\n%s\n', outDir);

%% ================= QUICK PLOT =================

figure('Color','w','Position',[100 100 1500 800])

subplot(2,2,1)
imagesc(DEM); axis image off; colorbar
title('DEM [m]')

subplot(2,2,2)
imagesc(y); axis image off; colorbar
title('Flow direction coordinate y [m]')

subplot(2,2,3)
imagesc(LULC); axis image off; colorbar
title('Uniform LULC')

subplot(2,2,4)
plot(Benchmark.time_min, Benchmark.Q_analytic_m3s, 'LineWidth', 2)
grid on
xlabel('Time [min]')
ylabel('Outlet discharge [m^3/s]')
title(sprintf('Analytical rising limb: Q_{eq}=%.3f m^3/s, t_{eq}=%.1f min', Q_eq, t_eq_sec/60))

sgtitle('Tilted Plane Benchmark + Analytical Outlet Flow');

exportgraphics(gcf, fullfile(outDir,'overview_plane_with_flow.png'), 'Resolution',300);

fprintf('Overview figure saved to:\n%s\n', fullfile(outDir,'overview_plane_with_flow.png'));