
folder_path = 'C:\Users\marcu\Documents\GitHub\HydroPol2D\Synthetic_Rasters';


synthetic_rasters(folder_path, ...
    3000, 3000, ...     % Domain 3 km x 3 km
    100, 10);            % Main resolution: 30 m, subgrid DEM: 10 m


visualize_synthetic_groundwater_rasters(folder_path)

visualize_synthetic_forcings(folder_path);

function synthetic_rasters(outputFolder, domainWidth_m, domainHeight_m, resMain_m, resSubgrid_m)
% =========================================================================
% 🛠️ Generate synthetic rasters for groundwater model (Marcus style)
% =========================================================================
% INPUTS:
%   outputFolder     - path to save rasters (e.g., 'D:\GW_Inputs\')
%   domainWidth_m    - domain width in meters (e.g., 3000)
%   domainHeight_m   - domain height in meters (e.g., 3000)
%   resMain_m        - main resolution (e.g., 30 m)
%   resSubgrid_m     - subgrid DEM resolution (e.g., 10 m)
% =========================================================================

% Compute grid dimensions (rounded to ensure consistency with maprefcells)
nCols = round(domainWidth_m / resMain_m);
nRows = round(domainHeight_m / resMain_m);

fprintf('📏 Domain = %d x %d m (%d rows x %d cols)\n', ...
        domainWidth_m, domainHeight_m, nRows, nCols);

% Create folder if needed
if ~exist(outputFolder, 'dir'); mkdir(outputFolder); end
R = maprefcells([0, domainHeight_m], [0, domainWidth_m], [nRows, nCols]);

% === 1. DEM - slope from west (100 m) to east (50 m) ===
DEM = repmat(linspace(100, 97, nCols), nRows, 1);
write_geotiff(fullfile(outputFolder, 'DEM.tif'), DEM, R);

% === 2. Categorical: LULC and SOIL ===
C1 = ones(nRows, nCols, 'uint8');
write_geotiff(fullfile(outputFolder, 'LULC.tif'), C1, R);
write_geotiff(fullfile(outputFolder, 'SOIL_new.tif'), C1, R);

% === 3. Vegetation/Surface Properties ===
write_geotiff(fullfile(outputFolder, 'Albedo.tif'), 0.15 * ones(nRows, nCols), R);
write_geotiff(fullfile(outputFolder, 'LAI.tif'),    2.0 * ones(nRows, nCols), R);
write_geotiff(fullfile(outputFolder, 'NDVI.tif'),   0.6 * ones(nRows, nCols), R);

% === 4. River properties ===
write_geotiff(fullfile(outputFolder, 'River_Width.tif'), 10 * ones(nRows, nCols), R);
write_geotiff(fullfile(outputFolder, 'River_Depth.tif'), 2 * ones(nRows, nCols), R);

% === 5. Depth to aquifer (DTB) ===
write_geotiff(fullfile(outputFolder, 'DTB.tif'), 5 * ones(nRows, nCols), R);

% === 6. Subgrid DEM ===
sgRows = domainHeight_m / resSubgrid_m;
sgCols = domainWidth_m / resSubgrid_m;
R_sub = maprefcells([0, domainHeight_m], [0, domainWidth_m], [sgRows, sgCols]);
Sub_DEM = repmat(linspace(100, 50, sgCols), sgRows, 1);
write_geotiff(fullfile(outputFolder, 'Subgrid_DEM.tif'), Sub_DEM, R_sub);

% === 7. Placeholder layers ===
for name = ["B1_Raster", "B2_Raster", "W1_Raster", "W2_Raster"]
    write_geotiff(fullfile(outputFolder, name + ".tif"), ones(nRows, nCols), R);
end

fprintf('✅ Synthetic rasters saved in: %s\n', outputFolder);

%% === 8. Generate masks and synthetic time series based on simulation time ===
fprintf("🚧 Generating masks and synthetic time series based on simulation time...\n")

% === Configurable parameters ===
sim_duration_hr = 365*24;   % total simulation duration in hours (modifiable)
dt_hr = 24;              % time step in hours (modifiable)
nTime = sim_duration_hr / dt_hr;
t = (0:dt_hr:(sim_duration_hr - dt_hr))';  % time vector in hours

% === Systematic synthetic points for each boundary condition ===
% We distribute 3 points across the grid using nRows and nCols

points_rainfall = [ ...
    round(nRows * 0.2), round(nCols * 0.2); ...
    round(nRows * 0.5), round(nCols * 0.5); ...
    round(nRows * 0.8), round(nCols * 0.8)];

points_inflow = [ ...
    round(nRows * 0.3), round(nCols * 0.7); ...
    round(nRows * 0.5), round(nCols * 0.2); ...
    round(nRows * 0.7), round(nCols * 0.5)];

points_stage = [ ...
    round(nRows * 0.2), round(nCols * 0.8); ...
    round(nRows * 0.6), round(nCols * 0.4); ...
    round(nRows * 0.9), round(nCols * 0.1)];

% === Initialize mask rasters ===
Rainfall_Mask = zeros(nRows, nCols);
Inflow_Mask   = zeros(nRows, nCols);
Stage_Mask    = zeros(nRows, nCols);

% Assign rainfall points
for i = 1:size(points_rainfall,1)
    r = points_rainfall(i,1);
    c = points_rainfall(i,2);
    Rainfall_Mask(r,c) = 1;
end

% Assign inflow points
for i = 1:size(points_inflow,1)
    r = points_inflow(i,1);
    c = points_inflow(i,2);
    Inflow_Mask(r,c) = 1;
end

% Assign stage points
for i = 1:size(points_stage,1)
    r = points_stage(i,1);
    c = points_stage(i,2);
    Stage_Mask(r,c) = 1;
end


% === Save mask rasters as GeoTIFF ===
write_geotiff(fullfile(outputFolder, 'Rainfall_Station_Mask.tif'), Rainfall_Mask, R);
write_geotiff(fullfile(outputFolder, 'Inflow_Mask.tif'), Inflow_Mask, R);
write_geotiff(fullfile(outputFolder, 'Stage_Mask.tif'), Stage_Mask, R);

% === Create synthetic signals ===
rain_signal   = 10 + 5 * sin(2 * pi * t / sim_duration_hr);   % mm/h rainfall
inflow_signal = 1 + 10 * sin(2 * pi * t / sim_duration_hr);  % m³/s inflow
stage_signal  = 0 + 0.5 * sin(2 * pi * t / sim_duration_hr); % m water level

% === Convert row/col to X/Y ===
[Xgrid, Ygrid] = worldGrid(R);

% === Export rainfall input ===
Rainfall_Table = table((1:3)', points_rainfall(:,1), points_rainfall(:,2), ...
    zeros(3,1), zeros(3,1), 'VariableNames', {'ID','row','col','X','Y'});
for i = 1:3
    [x, y] = rowcol_to_xy(Xgrid, Ygrid, points_rainfall(i,1), points_rainfall(i,2));
    Rainfall_Table.X(i) = x; Rainfall_Table.Y(i) = y;
end
for i = 1:nTime
    Rainfall_Table.(['T' num2str(i,'%03d')]) = rain_signal(i) * ones(3,1);
end
writetable(Rainfall_Table, fullfile(outputFolder, 'Rainfall_Spatial_Input.xlsx'));


% === Export inflow input ===
Inflow_Table = table(points_inflow(:,1), points_inflow(:,2), ...
    zeros(3,1), zeros(3,1), 'VariableNames', {'row','col','X','Y'});
for i = 1:3
    [x, y] = rowcol_to_xy(Xgrid, Ygrid, points_inflow(i,1), points_inflow(i,2));
    Inflow_Table.X(i) = x; Inflow_Table.Y(i) = y;
end
for i = 1:nTime
    Inflow_Table.(['Q' num2str(i,'%03d')]) = inflow_signal(i) * ones(3,1);
end
writetable(Inflow_Table, fullfile(outputFolder, 'Inflow_Input.xlsx'));


% === Export stage input ===
Stage_Table = table(points_stage(:,1), points_stage(:,2), ...
    zeros(3,1), zeros(3,1), 'VariableNames', {'row','col','X','Y'});
for i = 1:3
    [x, y] = rowcol_to_xy(Xgrid, Ygrid, points_stage(i,1), points_stage(i,2));
    Stage_Table.X(i) = x; Stage_Table.Y(i) = y;
end
for i = 1:nTime
    Stage_Table.(['H' num2str(i,'%03d')]) = stage_signal(i) * ones(3,1);
end
writetable(Stage_Table, fullfile(outputFolder, 'Stage_Input.xlsx'));


fprintf("✅ Synthetic series created with %d time steps (%.2f hrs, Δt = %.2f hrs)\n", ...
    nTime, sim_duration_hr, dt_hr);

end

function visualize_synthetic_groundwater_rasters(folderPath)
% =========================================================================
% 🖼️ Visualize all .tif rasters in 3D or 2D from the given folder
% =========================================================================

% Get list of all GeoTIFF files in folder
files = dir(fullfile(folderPath, '*.tif'));
nFiles = length(files);

% Set up plot layout
nCols = 4;
nRows = ceil(nFiles / nCols);

figure('Color', 'w', 'Position', [100, 100, 1400, 800]);
tiledlayout(nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:nFiles
    file = files(k).name;
    [Z, R] = readgeoraster(fullfile(folderPath, file));

    % Handle NaNs and types
    Z = double(Z);
    Z(Z == -9999 | Z == -9999.0) = NaN;

    % Create x-y mesh
    [x, y] = meshgrid(...
        R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2)-R.CellExtentInWorldX, ...
        R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1)+R.CellExtentInWorldY ...
    );

    nexttile;
    if max(Z(:)) - min(Z(:)) > 1e-2  % Use 3D surf for continuous variables
        surf(x, y, Z, 'EdgeColor', 'none');
        view(35, 55);
        shading interp;
        camlight; lighting gouraud;
    else
        imagesc(Z);
        axis image;
    end
    title(strrep(file, '_', '\_'), 'Interpreter', 'tex');
    colorbar;
end

sgtitle('🧭 Synthetic Groundwater Input Layers', 'FontSize', 16, 'FontWeight', 'bold');

end

function [X, Y] = worldGrid(R)
% Helper to generate X and Y coordinate grids from referencing matrix
x = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2)-R.CellExtentInWorldX;
y = R.YWorldLimits(2):-R.CellExtentInWorldY:R.YWorldLimits(1)+R.CellExtentInWorldY;
[X, Y] = meshgrid(x, y);
end

function write_geotiff(filename, data, R)
% =========================================================================
% 📝 Helper to save 2D array as GeoTIFF using map reference
% INPUTS:
%   filename - full path to output file (e.g., 'DEM.tif')
%   data     - 2D matrix to export
%   R        - maprefcells spatial referencing object
% =========================================================================

geotiffwrite(filename, data, R, 'CoordRefSysCode', 3857);  % 
end

function [x, y] = rowcol_to_xy(Xgrid, Ygrid, r, c)
    x = Xgrid(r, c);
    y = Ygrid(r, c);
end


function visualize_synthetic_forcings(folderPath)
% =========================================================================
% 🌧️📈 Visualize rainfall mask, inflow hydrographs, and stage hydrographs
% =========================================================================

fprintf('📊 Visualizing synthetic rainfall, inflow, and stage series...\n');

% === Load raster masks ===
Rainfall_Mask = readgeoraster(fullfile(folderPath, 'Rainfall_Station_Mask.tif'));
Inflow_Mask   = readgeoraster(fullfile(folderPath, 'Inflow_Mask.tif'));
Stage_Mask    = readgeoraster(fullfile(folderPath, 'Stage_Mask.tif'));

% === Plot the masks ===
figure('Color','w','Position',[100 100 1400 400]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% Rainfall
nexttile;
imagesc(Rainfall_Mask);
axis image;
title('Rainfall Station Mask');
colorbar;

% Inflow
nexttile;
imagesc(Inflow_Mask);
axis image;
title('Inflow Gauge Mask');
colorbar;

% Stage
nexttile;
imagesc(Stage_Mask);
axis image;
title('Stage Gauge Mask');
colorbar;

sgtitle('🗺️ Boundary Condition Locations','FontSize',14,'FontWeight','bold');

% === Plot Time Series ===
% Rainfall
Rain_tbl = readtable(fullfile(folderPath, 'Rainfall_Spatial_Input.xlsx'));
timeCols = startsWith(Rain_tbl.Properties.VariableNames, 'T');
timeVec  = 1:sum(timeCols);  % assuming Δt = 1 hr

figure('Color','w','Position',[100 100 1400 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
for i = 1:height(Rain_tbl)
    nexttile;
    plot(timeVec, Rain_tbl{i, timeCols}, '-o', 'LineWidth', 1.5);
    title(sprintf('Rainfall @ Gauge %d', Rain_tbl.ID(i)));
    xlabel('Time [hr]'); ylabel('Rainfall [mm/h]');
    grid on;
end
sgtitle('🌧️ Synthetic Rainfall Time Series','FontSize',14);

% Inflow
Inflow_tbl = readtable(fullfile(folderPath, 'Inflow_Input.xlsx'));
timeCols = startsWith(Inflow_tbl.Properties.VariableNames, 'Q');
timeVec = 1:sum(timeCols);

figure('Color','w','Position',[100 100 1400 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
for i = 1:height(Inflow_tbl)
    nexttile;
    plot(timeVec, Inflow_tbl{i, timeCols}, '-s', 'LineWidth', 1.5);
    title(sprintf('Inflow @ Gauge %d', i));
    xlabel('Time [hr]'); ylabel('Flow [m³/s]');
    grid on;
end
sgtitle('💧 Synthetic Inflow Hydrographs','FontSize',14);

% Stage
Stage_tbl = readtable(fullfile(folderPath, 'Stage_Input.xlsx'));
timeCols = startsWith(Stage_tbl.Properties.VariableNames, 'H');
timeVec = 1:sum(timeCols);

figure('Color','w','Position',[100 100 1400 300]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
for i = 1:height(Stage_tbl)
    nexttile;
    plot(timeVec, Stage_tbl{i, timeCols}, '-^', 'LineWidth', 1.5);
    title(sprintf('Stage @ Gauge %d', i));
    xlabel('Time [hr]'); ylabel('Stage [m]');
    grid on;
end
sgtitle('📶 Synthetic Stage Hydrographs','FontSize',14);
end

