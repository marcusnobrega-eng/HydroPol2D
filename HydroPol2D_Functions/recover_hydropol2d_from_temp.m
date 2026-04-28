function recover_hydropol2d_from_temp(tempDir, demFile, outputRoot, varargin)
%% ========================================================================
% recover_hydropol2d_from_temp.m
%
% Recover HydroPol2D post-processing products from Temporary_Files/
% save_map_hydro_*.mat files.
%
% MAIN FEATURES
% -------------------------------------------------------------------------
% 1) Wrapped as a function with explicit inputs:
%       tempDir, demFile, outputRoot
% 2) Keeps DEM input for WSE reconstruction and georeferenced exports
% 3) Infers the useful model domain from the saved depth maps:
%       validMask = any(isfinite(depth) & depth>threshold, over time)
% 4) All plan-view plots/animations use transparency outside validMask,
%    so invalid areas have NO color.
% 5) Static rasters are also masked outside validMask.
% 6) Slightly improved animation quality.
%
% REQUIRED INPUTS
% -------------------------------------------------------------------------
% tempDir    : folder containing save_map_hydro_*.mat
% demFile    : DEM GeoTIFF used in the simulation
% outputRoot : output folder for recovered products
%
% OPTIONAL NAME-VALUE INPUTS
% -------------------------------------------------------------------------
% 'date_begin'             : datetime(2025,6,1,0,0,0)
% 'date_end_intended'      : datetime(...) only for sanity check
% 'map_save_timestep'      : hours(3)
% 'flag_elapsed_time'      : false
% 'depthThreshold_m'       : 0.01
% 'deletePreviousOutputs'  : true
% 'progressEveryNFrames'   : 25
%
% EXAMPLE
% -------------------------------------------------------------------------
% recover_hydropol2d_from_temp( ...
%     '/path/to/Temporary_Files', ...
%     '/path/to/DEM.tif', ...
%     '/path/to/Recovered_PostProcessing', ...
%     'date_begin', datetime(2025,6,1,0,0,0), ...
%     'map_save_timestep', hours(3), ...
%     'flag_elapsed_time', false);
%
% ========================================================================

%% ========================================================================
% INPUT PARSING
% ========================================================================

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'tempDir', @(x) ischar(x) || isstring(x));
addRequired(p, 'demFile', @(x) ischar(x) || isstring(x));
addRequired(p, 'outputRoot', @(x) ischar(x) || isstring(x));

addParameter(p, 'date_begin', datetime(2025,6,1,0,0,0), @(x) isdatetime(x) && isscalar(x));
addParameter(p, 'date_end_intended', datetime(2025,6,30,0,0,0), @(x) isdatetime(x) && isscalar(x));
addParameter(p, 'map_save_timestep', hours(3), @(x) isduration(x) && isscalar(x));
addParameter(p, 'flag_elapsed_time', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'depthThreshold_m', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'deletePreviousOutputs', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'progressEveryNFrames', 25, @(x) isnumeric(x) && isscalar(x) && x >= 1);

parse(p, tempDir, demFile, outputRoot, varargin{:});
OPT = p.Results;

tempDir   = char(OPT.tempDir);
demFile   = char(OPT.demFile);
outputRoot = char(OPT.outputRoot);

date_begin            = OPT.date_begin;
date_end_intended     = OPT.date_end_intended;
map_save_timestep     = OPT.map_save_timestep;
flag_elapsed_time     = logical(OPT.flag_elapsed_time);
depthThreshold_m      = OPT.depthThreshold_m;
deletePreviousOutputs = logical(OPT.deletePreviousOutputs);
progressEveryNFrames  = OPT.progressEveryNFrames;

%% ========================================================================
% BUILD OPTIONS
% ========================================================================

BUILD = struct();
BUILD.timeVaryingDepthRasters = true;
BUILD.timeVaryingWQRasters    = true;
BUILD.staticRasters           = true;
BUILD.depthVideo              = true;
BUILD.wseDepthVideo           = true;
BUILD.wqConcGif               = false;   % GIFs are slower
BUILD.wqMassGif               = false;   % GIFs are slower
BUILD.snowVideo               = true;

%% ========================================================================
% VIDEO / GIF SETTINGS
% ========================================================================

VIDEO = struct();
VIDEO.visible                = false;    % faster offscreen
VIDEO.fps                    = 4;
VIDEO.targetHeightPx         = 900;      % slightly higher quality
VIDEO.aviQuality             = 95;       % slightly higher quality
VIDEO.convertToMp4           = false;
VIDEO.mp4Crf                 = 20;       % slightly better quality
VIDEO.mp4Preset              = 'medium';
VIDEO.deleteAviAfterMp4      = false;
VIDEO.aviProfile             = 'Motion JPEG AVI';
VIDEO.forceConstFrameSize    = true;
VIDEO.interpreter            = 'latex';
VIDEO.fontName               = 'Helvetica';
VIDEO.fontSize               = 12;
VIDEO.titleFontSize          = 13;

%% ========================================================================
% OUTPUT FOLDERS
% ========================================================================

OUT = struct();
OUT.ROOT          = string(outputRoot);
OUT.RastersDepth  = fullfile(OUT.ROOT, 'Rasters_Water_Depths');
OUT.RastersWQ     = fullfile(OUT.ROOT, 'Rasters_WQ');
OUT.RastersStatic = fullfile(OUT.ROOT, 'Rasters_Static');
OUT.Anim          = fullfile(OUT.ROOT, 'GIFs_MP4');

mk(OUT.ROOT);
mk(OUT.RastersDepth);
mk(OUT.RastersWQ);
mk(OUT.RastersStatic);
mk(OUT.Anim);

if deletePreviousOutputs
    deleteMatchingFiles(OUT.RastersDepth,  {'*.tif','*.tiff','*.tif.aux'});
    deleteMatchingFiles(OUT.RastersWQ,     {'*.tif','*.tiff','*.tif.aux'});
    deleteMatchingFiles(OUT.RastersStatic, {'*.tif','*.tiff','*.tif.aux'});
    deleteMatchingFiles(OUT.Anim,          {'*.avi','*.mp4','*.gif'});
end

fprintf('\n============================================================\n');
fprintf('HydroPol2D TEMP RECOVERY\n');
fprintf('Temporary_Files : %s\n', tempDir);
fprintf('DEM file        : %s\n', demFile);
fprintf('Output folder   : %s\n', OUT.ROOT);
fprintf('Start datetime  : %s\n', datestr(date_begin,'yyyy-mm-dd HH:MM:SS'));
fprintf('Intended end    : %s\n', datestr(date_end_intended,'yyyy-mm-dd HH:MM:SS'));
fprintf('Save timestep   : %s\n', char(string(map_save_timestep)));
fprintf('Elapsed naming  : %d\n', flag_elapsed_time);
fprintf('============================================================\n\n');

%% ========================================================================
% READ DEM
% ========================================================================

if ~isfile(demFile)
    error('DEM file not found: %s', demFile);
end

DEM = struct();
DEM.Z = [];
DEM.SpatialRef = [];
DEM.GeoKeyDirectoryTag = [];

try
    [DEM.Z, DEM.SpatialRef] = readgeoraster(demFile);
    DEM.Z = double(DEM.Z);
catch ME
    error('Failed to read DEM file %s\n%s', demFile, ME.message);
end

try
    infoDEM = geotiffinfo(demFile);
    if isfield(infoDEM, 'GeoTIFFTags') && isfield(infoDEM.GeoTIFFTags, 'GeoKeyDirectoryTag')
        DEM.GeoKeyDirectoryTag = infoDEM.GeoTIFFTags.GeoKeyDirectoryTag;
    end
catch
    DEM.GeoKeyDirectoryTag = [];
end

refInfo = struct('SpatialRef', DEM.SpatialRef, ...
                 'GeoKeyDirectoryTag', DEM.GeoKeyDirectoryTag, ...
                 'RefRaster', []);

elevation2D = DEM.Z;

%% ========================================================================
% FIND CHUNKS
% ========================================================================

files = dir(fullfile(tempDir, 'save_map_hydro_*.mat'));
if isempty(files)
    error('No save_map_hydro_*.mat files found in %s', tempDir);
end

ids = nan(numel(files),1);
for i = 1:numel(files)
    tok = regexp(files(i).name, 'save_map_hydro_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tok)
        error('Unexpected file name: %s', files(i).name);
    end
    ids(i) = str2double(tok{1});
end

[~, ord] = sort(ids);
files = files(ord);

fprintf('Found %d chunk files.\n', numel(files));

%% ========================================================================
% PRE-SCAN CHUNKS (PASS 1)
% ========================================================================

totalFrames   = 0;
chunkFrames   = zeros(numel(files),1);

hasDepth      = false;
hasWQConc     = false;
hasWQMass     = false;
hasSnow       = false;
hasElevation  = false;
hasRisk       = false;

maxDepth_m  = [];
maxWQConc   = [];
maxSnow     = [];
finalWQMass = [];

minDepthPositive_m = inf;
minWQConcPositive  = inf;
minSnowPositive    = inf;

% Useful-domain mask inferred from any valid depth over time
validMask = [];

for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name));
    Maps = getMapsStruct(S);
    if isempty(Maps)
        error('Could not find Maps structure in %s', files(i).name);
    end

    if isfield(Maps,'Hydro') && isfield(Maps.Hydro,'d')
        d = gatherIfNeeded(Maps.Hydro.d);
        chunkFrames(i) = size(d,3);
        totalFrames = totalFrames + chunkFrames(i);
        hasDepth = true;

        d_m = d / 1000;
        d_m(~isfinite(d_m)) = nan;

        localValidMask = any(isfinite(d_m) & (d_m > depthThreshold_m), 3);
        if isempty(validMask)
            validMask = localValidMask;
        else
            validMask = validMask | localValidMask;
        end

        d_m(d_m <= depthThreshold_m) = nan;

        localMax = max(d_m, [], 3, 'omitnan');
        if isempty(maxDepth_m)
            maxDepth_m = localMax;
        else
            maxDepth_m = max(maxDepth_m, localMax);
        end

        localMin = min(d_m(:), [], 'omitnan');
        if isfinite(localMin)
            minDepthPositive_m = min(minDepthPositive_m, localMin);
        end
    end

    if isfield(Maps,'WQ_States') && isfield(Maps.WQ_States,'Pol_Conc_Map')
        hasWQConc = true;
        c = gatherIfNeeded(Maps.WQ_States.Pol_Conc_Map);
        c(~isfinite(c)) = nan;
        c(c < 0) = nan;

        localMax = max(c, [], 3, 'omitnan');
        if isempty(maxWQConc)
            maxWQConc = localMax;
        else
            maxWQConc = max(maxWQConc, localMax);
        end

        localMin = min(c(:), [], 'omitnan');
        if isfinite(localMin)
            minWQConcPositive = min(minWQConcPositive, localMin);
        end
    end

    if isfield(Maps,'WQ_States') && isfield(Maps.WQ_States,'Pol_mass_map')
        hasWQMass = true;
        m = gatherIfNeeded(Maps.WQ_States.Pol_mass_map);
        if i == numel(files)
            finalWQMass = m(:,:,end);
        end
    end

    if isfield(Maps,'Hydro') && isfield(Maps.Hydro,'Snowpack')
        hasSnow = true;
        s = gatherIfNeeded(Maps.Hydro.Snowpack);
        s(~isfinite(s)) = nan;
        s(s <= 0) = nan;

        localMax = max(s, [], 3, 'omitnan');
        if isempty(maxSnow)
            maxSnow = localMax;
        else
            maxSnow = max(maxSnow, localMax);
        end

        localMin = min(s(:), [], 'omitnan');
        if isfinite(localMin)
            minSnowPositive = min(minSnowPositive, localMin);
        end
    end

    if isfield(Maps,'Hydro') && isfield(Maps.Hydro,'risk')
        hasRisk = true;
    end
end

if ~hasDepth
    error('No Maps.Hydro.d found in the temporary files.');
end

if isempty(validMask)
    error('Could not infer a valid model mask from Maps.Hydro.d.');
end

%% ========================================================================
% RESAMPLE DEM TO MATCH RECOVERED MAP GRID
% ========================================================================

[nyMap, nxMap] = size(maxDepth_m);
[nyDEM, nxDEM] = size(DEM.Z);

fprintf('DEM size           = [%d %d]\n', nyDEM, nxDEM);
fprintf('Recovered map size = [%d %d]\n', nyMap, nxMap);

if isempty(DEM.SpatialRef)
    error('DEM spatial reference could not be read from %s', demFile);
end

if isequal(size(DEM.Z), size(maxDepth_m))
    fprintf('DEM already matches recovered map grid.\n');
    elevation2D = DEM.Z;
    hasElevation = true;
else
    fprintf('Resampling DEM to recovered map grid...\n');
    [elevation2D, newRef] = resampleDEMToMatchMaps(DEM.Z, DEM.SpatialRef, size(maxDepth_m));

    DEM.Z = elevation2D;
    DEM.SpatialRef = newRef;
    refInfo.SpatialRef = newRef;
    refInfo.RefRaster  = maxDepth_m;

    fprintf('DEM resampled successfully.\n');
    fprintf('Resampled DEM size = [%d %d]\n', size(elevation2D,1), size(elevation2D,2));

    hasElevation = true;
end

if ~isequal(size(elevation2D), size(maxDepth_m))
    error(['DEM resampling failed to match recovered depth maps.\n' ...
           'DEM size after resample = [%d %d]\nRecovered size = [%d %d]'], ...
           size(elevation2D,1), size(elevation2D,2), size(maxDepth_m,1), size(maxDepth_m,2));
end

refInfo.SpatialRef = DEM.SpatialRef;
refInfo.RefRaster  = maxDepth_m;

% Mask DEM itself outside useful domain so derived WSE products do not show
% spurious color there
elevation2D(~validMask) = nan;

fprintf('Total recovered frames : %d\n', totalFrames);
fprintf('Depth maps found       : %d\n', hasDepth);
fprintf('Elevation found        : %d\n', hasElevation);
fprintf('WQ concentration found : %d\n', hasWQConc);
fprintf('WQ mass found          : %d\n', hasWQMass);
fprintf('Snowpack found         : %d\n', hasSnow);
fprintf('Risk found             : %d\n\n', hasRisk);

%% ========================================================================
% REBUILD TIME AXIS
% ========================================================================

if ~isduration(map_save_timestep)
    error('map_save_timestep must be a duration, e.g. hours(1) or minutes(30).');
end

time_dt = date_begin + (0:totalFrames-1).' * map_save_timestep;
time_minutes = minutes(time_dt - time_dt(1));
time_minutes = time_minutes(:);
date_end_actual = time_dt(end);

fprintf('Actual recovered end   : %s\n', datestr(date_end_actual,'yyyy-mm-dd HH:MM:SS'));

if date_end_actual < date_end_intended
    fprintf(['Note: recovered simulation ended before intended end date.\n' ...
             '      This is expected when the run failed early.\n\n']);
elseif date_end_actual > date_end_intended
    fprintf(['Note: recovered simulation extends beyond intended end date.\n' ...
             '      Check that the chosen save timestep is correct.\n\n']);
else
    fprintf('Recovered end matches intended end date.\n\n');
end

%% ========================================================================
% COLORMAPS
% ========================================================================

try
    [Spectrum,Depth_Purple,~,~,~,~,Velocity_RAS,~,~,WSE_RAS] = coloramps(); %#ok<ASGLU>
catch
    Spectrum      = parula(256);
    Depth_Purple  = turbo(256);
    Velocity_RAS  = turbo(256);
    WSE_RAS       = turbo(256);
end

%% ========================================================================
% BUILD XY COORDINATES FROM DEM GRID
% ========================================================================

[ny, nx] = size(maxDepth_m);
SR = DEM.SpatialRef;

try
    x_grid = linspace(SR.XWorldLimits(1) + SR.CellExtentInWorldX/2, ...
                      SR.XWorldLimits(2) - SR.CellExtentInWorldX/2, nx);

    y_grid = linspace(SR.YWorldLimits(2) - SR.CellExtentInWorldY/2, ...
                      SR.YWorldLimits(1) + SR.CellExtentInWorldY/2, ny);
catch
    x_grid = 1:nx;
    y_grid = ny:-1:1;
end

xlimVec = [min(x_grid) max(x_grid)];
ylimVec = [min(y_grid) max(y_grid)];

%% ========================================================================
% GLOBAL PLOT SCALES
% ========================================================================

maxDepth_m(~validMask) = nan;

zmaxDepth = max(maxDepth_m(:), [], 'omitnan');
if isempty(zmaxDepth) || ~isfinite(zmaxDepth) || zmaxDepth <= 0
    zmaxDepth = 1;
end

if ~isfinite(minDepthPositive_m)
    minDepthPositive_m = depthThreshold_m;
end

zminDepth = min(depthThreshold_m, minDepthPositive_m);
if zminDepth >= zmaxDepth
    zminDepth = 0;
end

if hasElevation
    zminWSE = min(elevation2D(:), [], 'omitnan');
    zmaxWSE = max((maxDepth_m + elevation2D), [], 'all', 'omitnan');
    if ~isfinite(zminWSE), zminWSE = 0; end
    if ~isfinite(zmaxWSE) || zmaxWSE <= zminWSE
        zmaxWSE = zminWSE + 1;
    end
end

if hasWQConc
    maxWQConc(~validMask) = nan;
    zmaxWQ = max(maxWQConc(:), [], 'omitnan');
    if ~isfinite(zmaxWQ) || zmaxWQ <= 0
        zmaxWQ = 1;
    end
    zminWQ = minWQConcPositive;
    if ~isfinite(zminWQ) || zminWQ >= zmaxWQ
        zminWQ = 0;
    end
end

if hasSnow
    maxSnow(~validMask) = nan;
    zmaxSnow = max(maxSnow(:), [], 'omitnan');
    if ~isfinite(zmaxSnow) || zmaxSnow <= 0
        zmaxSnow = 1;
    end
    zminSnow = minSnowPositive;
    if ~isfinite(zminSnow) || zminSnow >= zmaxSnow
        zminSnow = 0;
    end
end

%% ========================================================================
% EXPORT STATIC RASTERS
% ========================================================================

if BUILD.staticRasters
    fprintf('Exporting recoverable static rasters...\n');

    writeRecoveredRaster(fullfile(OUT.RastersStatic, 'Maximum_Depths.tif'), maxDepth_m, refInfo);

    if hasElevation
        maxWSE = maxDepth_m + elevation2D;
        maxWSE(~validMask) = nan;
        writeRecoveredRaster(fullfile(OUT.RastersStatic, 'Max_Water_Surface_Elevation.tif'), maxWSE, refInfo);
        writeRecoveredRaster(fullfile(OUT.RastersStatic, 'Resampled_DEM_Model_Grid.tif'), elevation2D, refInfo);
    end

    if hasWQConc && ~isempty(maxWQConc)
        writeRecoveredRaster(fullfile(OUT.RastersStatic, 'Maximum_Pol_Conc_min.tif'), maxWQConc, refInfo);
    end

    if hasWQMass && ~isempty(finalWQMass)
        finalWQMass = double(finalWQMass);
        finalWQMass(~validMask) = nan;
        writeRecoveredRaster(fullfile(OUT.RastersStatic, 'Final_Mass_Of_Pollutant.tif'), finalWQMass, refInfo);
    end

    if hasSnow && ~isempty(maxSnow)
        writeRecoveredRaster(fullfile(OUT.RastersStatic, 'Maximum_Snowpack.tif'), maxSnow, refInfo);
    end

    fprintf('Static rasters complete.\n\n');
end

%% ========================================================================
% PREPARE ANIMATION WRITERS / FIGURES
% ========================================================================

ANIM = struct();

if BUILD.depthVideo && hasDepth
    fprintf('Preparing Depths video...\n');
    [ANIM.depth.video, ANIM.depth.aviPath, ANIM.depth.mp4Path, ANIM.depth.fig] = ...
        startGlobalVideo(OUT.Anim, 'Depths', VIDEO);
    ANIM.depth.ax = axes('Parent', ANIM.depth.fig);
    ANIM.depth.im = setupRasterAxes(ANIM.depth.ax, x_grid, y_grid, nan(ny,nx), ...
        Depth_Purple, [zminDepth zmaxDepth], 'Depth [m]', VIDEO, xlimVec, ylimVec, validMask);
    ANIM.depth.title = title(ANIM.depth.ax, '', 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.titleFontSize);
    ANIM.depth.targetH = [];
    ANIM.depth.targetW = [];
end

if BUILD.wseDepthVideo && hasElevation && hasDepth
    fprintf('Preparing WSE_Depths video...\n');
    [ANIM.wse.video, ANIM.wse.aviPath, ANIM.wse.mp4Path, ANIM.wse.fig] = ...
        startGlobalVideo(OUT.Anim, 'WSE_Depths', VIDEO);
    ANIM.wse.ax1 = subplot(2,1,1, 'Parent', ANIM.wse.fig);
    ANIM.wse.ax2 = subplot(2,1,2, 'Parent', ANIM.wse.fig);
    ANIM.wse.im1 = setupRasterAxes(ANIM.wse.ax1, x_grid, y_grid, nan(ny,nx), ...
        WSE_RAS, [zminWSE zmaxWSE], 'WSE [m]', VIDEO, xlimVec, ylimVec, validMask);
    ANIM.wse.im2 = setupRasterAxes(ANIM.wse.ax2, x_grid, y_grid, nan(ny,nx), ...
        Depth_Purple, [zminDepth zmaxDepth], 'Depth [m]', VIDEO, xlimVec, ylimVec, validMask);
    ANIM.wse.title1 = title(ANIM.wse.ax1, '', 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.titleFontSize);
    ANIM.wse.title2 = title(ANIM.wse.ax2, '', 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.titleFontSize);
    ANIM.wse.targetH = [];
    ANIM.wse.targetW = [];
end

if BUILD.wqConcGif && hasWQConc
    fprintf('Preparing Pollutant_Concentration.gif ...\n');
    ANIM.wqGif.gifPath = fullfile(OUT.Anim, 'Pollutant_Concentration.gif');
    safeDelete(ANIM.wqGif.gifPath);
    ANIM.wqGif.fig = figure('Visible', onOff(VIDEO.visible), 'Color', 'w', ...
        'Units', 'pixels', 'Position', [60 60 950 700], 'Renderer', 'opengl');
    ANIM.wqGif.ax = axes('Parent', ANIM.wqGif.fig);
    ANIM.wqGif.im = setupRasterAxes(ANIM.wqGif.ax, x_grid, y_grid, nan(ny,nx), ...
        Spectrum, [zminWQ zmaxWQ], 'Concentration (mg/L)', VIDEO, xlimVec, ylimVec, validMask);
    ANIM.wqGif.title = title(ANIM.wqGif.ax, '', 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.titleFontSize);
    ANIM.wqGif.firstFrame = true;
end

if BUILD.wqMassGif && hasWQMass
    fprintf('Preparing Mass_of_pollutant.gif ...\n');
    ANIM.massGif.gifPath = fullfile(OUT.Anim, 'Mass_of_pollutant.gif');
    safeDelete(ANIM.massGif.gifPath);
    ANIM.massGif.fig = figure('Visible', onOff(VIDEO.visible), 'Color', 'w', ...
        'Units', 'pixels', 'Position', [60 60 950 700], 'Renderer', 'opengl');
    ANIM.massGif.ax = axes('Parent', ANIM.massGif.fig);

    [massLogMin, massLogMax] = computeGlobalLogMassLimits(files, validMask);
    ANIM.massGif.cLim = [massLogMin massLogMax];
    ANIM.massGif.im = setupRasterAxes(ANIM.massGif.ax, x_grid, y_grid, nan(ny,nx), ...
        Velocity_RAS, ANIM.massGif.cLim, 'Log mass', VIDEO, xlimVec, ylimVec, validMask);
    ANIM.massGif.title = title(ANIM.massGif.ax, '', 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.titleFontSize);
    ANIM.massGif.firstFrame = true;
end

if BUILD.snowVideo && hasSnow
    fprintf('Preparing Snowpack video...\n');
    [ANIM.snow.video, ANIM.snow.aviPath, ANIM.snow.mp4Path, ANIM.snow.fig] = ...
        startGlobalVideo(OUT.Anim, 'Snowpack', VIDEO);
    ANIM.snow.ax = axes('Parent', ANIM.snow.fig);
    ANIM.snow.im = setupRasterAxes(ANIM.snow.ax, x_grid, y_grid, nan(ny,nx), ...
        Depth_Purple, [zminSnow zmaxSnow], 'Snowpack [mm]', VIDEO, xlimVec, ylimVec, validMask);
    ANIM.snow.title = title(ANIM.snow.ax, '', 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.titleFontSize);
    ANIM.snow.targetH = [];
    ANIM.snow.targetW = [];
end

%% ========================================================================
% MAIN EXPORT LOOP (PASS 2)
% ========================================================================

fprintf('Starting main export loop...\n');

globalIdx = 0;

for iFile = 1:numel(files)
    fpath = fullfile(files(iFile).folder, files(iFile).name);
    S = load(fpath);
    Maps = getMapsStruct(S);
    if isempty(Maps)
        error('Could not find Maps structure in %s', files(iFile).name);
    end

    d_all = [];
    c_all = [];
    m_all = [];
    s_all = [];

    if isfield(Maps,'Hydro') && isfield(Maps.Hydro,'d')
        d_all = gatherIfNeeded(Maps.Hydro.d);
    end
    if hasWQConc && isfield(Maps,'WQ_States') && isfield(Maps.WQ_States,'Pol_Conc_Map')
        c_all = gatherIfNeeded(Maps.WQ_States.Pol_Conc_Map);
    end
    if hasWQMass && isfield(Maps,'WQ_States') && isfield(Maps.WQ_States,'Pol_mass_map')
        m_all = gatherIfNeeded(Maps.WQ_States.Pol_mass_map);
    end
    if hasSnow && isfield(Maps,'Hydro') && isfield(Maps.Hydro,'Snowpack')
        s_all = gatherIfNeeded(Maps.Hydro.Snowpack);
    end

    nLocal = chunkFrames(iFile);

    for j = 1:nLocal
        globalIdx = globalIdx + 1;
        titleStr = makeFrameTitle(globalIdx, time_dt, flag_elapsed_time, time_minutes);

        depth_m = [];
        conc    = [];
        logMass = [];
        snow    = [];
        wse     = [];

        if ~isempty(d_all)
            depth_m = d_all(:,:,j) / 1000;
            depth_m(~isfinite(depth_m)) = nan;
            depth_m(depth_m <= depthThreshold_m) = nan;
            depth_m(~validMask) = nan;

            if hasElevation
                wse = elevation2D + depth_m;
                wse(~validMask) = nan;
            end
        end

        if ~isempty(c_all)
            conc = c_all(:,:,j);
            conc(~isfinite(conc)) = nan;
            conc(conc < 0) = nan;
            conc(~validMask) = nan;
        end

        if ~isempty(m_all)
            logMass = double(m_all(:,:,j));
            logMass(logMass <= 0) = nan;
            logMass(~validMask) = nan;
            logMass = log(logMass);
        end

        if ~isempty(s_all)
            snow = s_all(:,:,j);
            snow(~isfinite(snow)) = nan;
            snow(snow <= 0) = nan;
            snow(~validMask) = nan;
        end

        %% Time-varying rasters
        if BUILD.timeVaryingDepthRasters && ~isempty(depth_m)
            if flag_elapsed_time
                hh = time_minutes(globalIdx) / 60;
                timeLabelDepth = sprintf('t_%s_h', stripTrailingZeros(hh));
            else
                timeLabelDepth = datestr(time_dt(globalIdx), 'yyyy_mm_dd_HH_MM_ss');
            end

            outName = fullfile(OUT.RastersDepth, sprintf('Flood_Depths_%s.tif', timeLabelDepth));
            writeRecoveredRaster(outName, depth_m, refInfo);
        end

        if BUILD.timeVaryingWQRasters && ~isempty(conc)
            if flag_elapsed_time
                outNameWQ = fullfile(OUT.RastersWQ, sprintf('Pollutant_Concentration_%smin.tif', ...
                    stripTrailingZeros(time_minutes(globalIdx))));
            else
                outNameWQ = fullfile(OUT.RastersWQ, sprintf('Pollutant_Concentration_%s.tif', ...
                    datestr(time_dt(globalIdx), 'yyyy_mm_dd_HH_MM_ss')));
            end
            writeRecoveredRaster(outNameWQ, conc, refInfo);
        end

        %% Depth video
        if isfield(ANIM, 'depth') && ~isempty(depth_m)
            set(ANIM.depth.im, 'CData', depth_m, 'AlphaData', double(~isnan(depth_m)));
            set(ANIM.depth.title, 'String', titleStr);
            drawnow limitrate nocallbacks;
            img = getframeImage(ANIM.depth.fig, VIDEO, ANIM.depth.targetH, ANIM.depth.targetW);
            if isempty(ANIM.depth.targetH)
                [ANIM.depth.targetH, ANIM.depth.targetW, ~] = size(img);
            end
            writeVideo(ANIM.depth.video, img);
        end

        %% WSE + Depth video
        if isfield(ANIM, 'wse') && ~isempty(depth_m) && ~isempty(wse)
            set(ANIM.wse.im1, 'CData', wse,     'AlphaData', double(~isnan(wse)));
            set(ANIM.wse.im2, 'CData', depth_m, 'AlphaData', double(~isnan(depth_m)));
            set(ANIM.wse.title1, 'String', titleStr);
            set(ANIM.wse.title2, 'String', titleStr);
            drawnow limitrate nocallbacks;
            img = getframeImage(ANIM.wse.fig, VIDEO, ANIM.wse.targetH, ANIM.wse.targetW);
            if isempty(ANIM.wse.targetH)
                [ANIM.wse.targetH, ANIM.wse.targetW, ~] = size(img);
            end
            writeVideo(ANIM.wse.video, img);
        end

        %% WQ concentration GIF
        if isfield(ANIM, 'wqGif') && ~isempty(conc)
            set(ANIM.wqGif.im, 'CData', conc, 'AlphaData', double(~isnan(conc)));
            set(ANIM.wqGif.title, 'String', titleStr);
            drawnow limitrate nocallbacks;
            ANIM.wqGif.firstFrame = appendGifFrame(ANIM.wqGif.fig, ANIM.wqGif.gifPath, VIDEO.fps, ANIM.wqGif.firstFrame);
        end

        %% Mass GIF
        if isfield(ANIM, 'massGif') && ~isempty(logMass)
            set(ANIM.massGif.im, 'CData', logMass, 'AlphaData', double(~isnan(logMass)));
            set(ANIM.massGif.title, 'String', titleStr);
            drawnow limitrate nocallbacks;
            ANIM.massGif.firstFrame = appendGifFrame(ANIM.massGif.fig, ANIM.massGif.gifPath, VIDEO.fps, ANIM.massGif.firstFrame);
        end

        %% Snow video
        if isfield(ANIM, 'snow') && ~isempty(snow)
            set(ANIM.snow.im, 'CData', snow, 'AlphaData', double(~isnan(snow)));
            set(ANIM.snow.title, 'String', titleStr);
            drawnow limitrate nocallbacks;
            img = getframeImage(ANIM.snow.fig, VIDEO, ANIM.snow.targetH, ANIM.snow.targetW);
            if isempty(ANIM.snow.targetH)
                [ANIM.snow.targetH, ANIM.snow.targetW, ~] = size(img);
            end
            writeVideo(ANIM.snow.video, img);
        end

        if mod(globalIdx, progressEveryNFrames) == 0 || globalIdx == totalFrames
            fprintf('  Frame %d / %d processed\n', globalIdx, totalFrames);
        end
    end

    clear S Maps d_all c_all m_all s_all
end

fprintf('Main export loop complete.\n\n');

%% ========================================================================
% CLOSE WRITERS / FIGURES
% ========================================================================

if isfield(ANIM, 'depth')
    finishGlobalVideo(ANIM.depth.video, ANIM.depth.aviPath, ANIM.depth.mp4Path, VIDEO);
    close(ANIM.depth.fig);
    fprintf('Depths video complete.\n');
end

if isfield(ANIM, 'wse')
    finishGlobalVideo(ANIM.wse.video, ANIM.wse.aviPath, ANIM.wse.mp4Path, VIDEO);
    close(ANIM.wse.fig);
    fprintf('WSE_Depths video complete.\n');
end

if isfield(ANIM, 'wqGif')
    close(ANIM.wqGif.fig);
    fprintf('Pollutant_Concentration.gif complete.\n');
end

if isfield(ANIM, 'massGif')
    close(ANIM.massGif.fig);
    fprintf('Mass_of_pollutant.gif complete.\n');
end

if isfield(ANIM, 'snow')
    finishGlobalVideo(ANIM.snow.video, ANIM.snow.aviPath, ANIM.snow.mp4Path, VIDEO);
    close(ANIM.snow.fig);
    fprintf('Snowpack video complete.\n');
end

fprintf('\n============================================================\n');
fprintf('RECOVERY COMPLETE.\n');
fprintf('Recovered outputs written to:\n%s\n', OUT.ROOT);
fprintf('Actual recovered end time:\n%s\n', datestr(date_end_actual,'yyyy-mm-dd HH:MM:SS'));
fprintf('============================================================\n');

end

%% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function mk(folderPath)
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
end

function deleteMatchingFiles(folderPath, patterns)
    if ~exist(folderPath, 'dir')
        return;
    end
    if ischar(patterns) || isstring(patterns)
        patterns = cellstr(patterns);
    end
    for k = 1:numel(patterns)
        ff = dir(fullfile(folderPath, patterns{k}));
        for i = 1:numel(ff)
            try
                delete(fullfile(ff(i).folder, ff(i).name));
            catch
            end
        end
    end
end

function safeDelete(fpath)
    if exist(fpath, 'file')
        try
            delete(fpath);
        catch
        end
    end
end

function out = onOff(tf)
    if tf
        out = 'on';
    else
        out = 'off';
    end
end

function A = gatherIfNeeded(A)
    try
        if isa(A, 'gpuArray')
            A = gather(A);
        end
    catch
    end
    A = double(A);
end

function Maps = getMapsStruct(S)
    Maps = [];
    if isfield(S, 'Maps')
        Maps = S.Maps;
        return;
    end

    f = fieldnames(S);
    for i = 1:numel(f)
        val = S.(f{i});
        if isstruct(val)
            if isfield(val, 'Hydro') || isfield(val, 'WQ_States')
                Maps = val;
                return;
            end
        end
    end
end

function [Zout, Rout] = resampleDEMToMatchMaps(Zin, Rin, outSize)
% Resample DEM to target raster size while preserving same world limits.
%
% outSize = [nRows nCols]

    nRows = outSize(1);
    nCols = outSize(2);

    xIn = linspace(Rin.XWorldLimits(1) + Rin.CellExtentInWorldX/2, ...
                   Rin.XWorldLimits(2) - Rin.CellExtentInWorldX/2, size(Zin,2));

    yIn = linspace(Rin.YWorldLimits(2) - Rin.CellExtentInWorldY/2, ...
                   Rin.YWorldLimits(1) + Rin.CellExtentInWorldY/2, size(Zin,1));

    xOut = linspace(Rin.XWorldLimits(1) + (diff(Rin.XWorldLimits)/nCols)/2, ...
                    Rin.XWorldLimits(2) - (diff(Rin.XWorldLimits)/nCols)/2, nCols);

    yOut = linspace(Rin.YWorldLimits(2) - (diff(Rin.YWorldLimits)/nRows)/2, ...
                    Rin.YWorldLimits(1) + (diff(Rin.YWorldLimits)/nRows)/2, nRows);

    [Xq, Yq] = meshgrid(xOut, yOut);
    [X,  Y ] = meshgrid(xIn, yIn);

    Zout = interp2(X, Y, double(Zin), Xq, Yq, 'linear', nan);

    try
        Rout = maprefcells(Rin.XWorldLimits, Rin.YWorldLimits, [nRows nCols], ...
            'ColumnsStartFrom', Rin.ColumnsStartFrom, ...
            'RowsStartFrom', Rin.RowsStartFrom);
    catch
        Rout = maprefcells(Rin.XWorldLimits, Rin.YWorldLimits, [nRows nCols]);
    end
end

function writeRecoveredRaster(outName, A, refInfo)
% Write GeoTIFF when possible. Otherwise fallback to TIFF.
% Input array already expected to be masked outside valid domain.

    A = double(A);

    try
        if ~isempty(refInfo) && isfield(refInfo,'SpatialRef') && ~isempty(refInfo.SpatialRef)
            if isfield(refInfo, 'GeoKeyDirectoryTag') && ~isempty(refInfo.GeoKeyDirectoryTag)
                geotiffwrite(outName, A, refInfo.SpatialRef, ...
                    'GeoKeyDirectoryTag', refInfo.GeoKeyDirectoryTag);
            else
                geotiffwrite(outName, A, refInfo.SpatialRef);
            end
        else
            imwrite(single(A), outName, 'tif');
        end
    catch
        [pth, nam, ~] = fileparts(outName);
        fallbackName = fullfile(pth, [nam '.tif']);
        try
            imwrite(single(A), fallbackName, 'tif');
        catch ME
            warning('Failed writing raster %s\n%s', outName, ME.message);
        end
    end
end

function hImg = setupRasterAxes(ax, x, y, C, cmap, climVals, cbarLabel, VIDEO, xlimVec, ylimVec, validMask)
% Create a masked image plot.
% Invalid cells are transparent, so they show no color.

    imagesc(ax, x, y, C);
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    set(ax, 'Layer', 'top');
    set(ax, 'FontName', VIDEO.fontName, 'FontSize', VIDEO.fontSize);
    xlim(ax, xlimVec);
    ylim(ax, ylimVec);
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    colormap(ax, cmap);
    caxis(ax, climVals);

    hImg = findobj(ax, 'Type', 'image');
    if isempty(hImg)
        hImg = image(ax, 'XData', x, 'YData', y, 'CData', C);
        set(ax, 'YDir', 'normal');
    else
        hImg = hImg(1);
    end

    alpha0 = double(validMask);
    alpha0(~isfinite(C)) = 0;
    set(hImg, 'AlphaData', alpha0);

    cb = colorbar(ax);
    ylabel(cb, cbarLabel, 'Interpreter', VIDEO.interpreter, 'FontSize', VIDEO.fontSize);
    box(ax, 'on');
end

function [vObj, aviPath, mp4Path, fig] = startGlobalVideo(outFolder, baseName, VIDEO)
    aviPath = fullfile(outFolder, [baseName '.avi']);
    mp4Path = fullfile(outFolder, [baseName '.mp4']);

    safeDelete(aviPath);
    if VIDEO.convertToMp4
        safeDelete(mp4Path);
    end

    fig = figure('Visible', onOff(VIDEO.visible), ...
                 'Color', 'w', ...
                 'Units', 'pixels', ...
                 'Position', [50 50 1100 850], ...
                 'Renderer', 'opengl');

    vObj = VideoWriter(aviPath, VIDEO.aviProfile);
    vObj.FrameRate = VIDEO.fps;
    try
        vObj.Quality = VIDEO.aviQuality;
    catch
    end
    open(vObj);
end

function finishGlobalVideo(vObj, aviPath, mp4Path, VIDEO)
    try
        close(vObj);
    catch
    end

    if VIDEO.convertToMp4
        try
            convertAviToMp4_ffmpeg(aviPath, mp4Path, VIDEO.mp4Crf, VIDEO.mp4Preset);
            if VIDEO.deleteAviAfterMp4 && isfile(mp4Path)
                safeDelete(aviPath);
            end
        catch ME
            warning('Could not convert AVI to MP4 for %s\n%s', aviPath, ME.message);
        end
    end
end

function convertAviToMp4_ffmpeg(aviPath, mp4Path, crf, preset)
% Requires ffmpeg in system path.

    if ~isfile(aviPath)
        error('AVI not found: %s', aviPath);
    end

    cmd = sprintf('ffmpeg -y -i "%s" -c:v libx264 -pix_fmt yuv420p -crf %d -preset %s "%s"', ...
                  aviPath, crf, preset, mp4Path);

    [status, msg] = system(cmd);
    if status ~= 0
        error('ffmpeg failed:\n%s', msg);
    end
end

function img = getframeImage(fig, VIDEO, targetH, targetW)
% Capture frame and optionally normalize frame size.

    drawnow;
    fr = getframe(fig);
    img = fr.cdata;

    if VIDEO.forceConstFrameSize
        if ~isempty(targetH) && ~isempty(targetW)
            if size(img,1) ~= targetH || size(img,2) ~= targetW
                img = imresize(img, [targetH targetW]);
            end
        else
            scale = VIDEO.targetHeightPx / size(img,1);
            targetH = round(size(img,1) * scale);
            targetW = round(size(img,2) * scale);
            img = imresize(img, [targetH targetW]);
        end
    end
end

function firstFrame = appendGifFrame(fig, gifPath, fps, firstFrame)
    fr = getframe(fig);
    [A, map] = rgb2ind(frame2im(fr), 256);

    delayTime = 1 / max(fps, 1);

    if firstFrame
        imwrite(A, map, gifPath, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
        firstFrame = false;
    else
        imwrite(A, map, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end

function s = makeFrameTitle(globalIdx, time_dt, flag_elapsed_time, time_minutes)
    if flag_elapsed_time
        hh = time_minutes(globalIdx)/60;
        s = sprintf('Elapsed time = %s h', stripTrailingZeros(hh));
    else
        s = datestr(time_dt(globalIdx), 'yyyy-mm-dd HH:MM:SS');
    end
end

function out = stripTrailingZeros(x)
    if abs(x - round(x)) < 1e-12
        out = sprintf('%d', round(x));
    else
        out = regexprep(sprintf('%.6f', x), '0+$', '');
        out = regexprep(out, '\.$', '');
    end
end

function [logMin, logMax] = computeGlobalLogMassLimits(files, validMask)
    logMin = inf;
    logMax = -inf;

    for i = 1:numel(files)
        S = load(fullfile(files(i).folder, files(i).name));
        Maps = getMapsStruct(S);
        if isempty(Maps)
            continue;
        end
        if isfield(Maps,'WQ_States') && isfield(Maps.WQ_States,'Pol_mass_map')
            M = gatherIfNeeded(Maps.WQ_States.Pol_mass_map);
            M = double(M);
            M(M <= 0) = nan;

            for k = 1:size(M,3)
                A = M(:,:,k);
                A(~validMask) = nan;
                A = log(A);
                amin = min(A(:), [], 'omitnan');
                amax = max(A(:), [], 'omitnan');
                if isfinite(amin), logMin = min(logMin, amin); end
                if isfinite(amax), logMax = max(logMax, amax); end
            end
        end
    end

    if ~isfinite(logMin) || ~isfinite(logMax) || logMin >= logMax
        logMin = 0;
        logMax = 1;
    end
end