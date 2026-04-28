%% ========================================================================
% input_paths_bypass_script.m
% ========================================================================
% HydroPol2D | External file-path definition for full bypass workflow
%
% PURPOSE
% -------------------------------------------------------------------------
% This script defines ONLY the locations of external files and folders used
% by HydroPol2D when running in bypass mode.
%
% It is intentionally separated from input_data_bypass_script.m:
%
%   1) input_paths_bypass_script.m  -> where files are
%   2) input_data_bypass_script.m   -> how those files are converted into
%                                      HydroPol2D input structures
%
% This file should be called BEFORE preprocessing reads DEM/LULC/SOIL.
%
% OUTPUT
% -------------------------------------------------------------------------
% This script must create one struct:
%
%   InputPaths
%
% containing:
%   - tool paths
%   - static raster paths
%   - optional raster paths
%   - forcing folders
%   - hydrograph CSV paths
%   - forcing file lists already sorted in time
%
% IMPORTANT
% -------------------------------------------------------------------------
% For raster forcing folders, filenames should ideally contain timestamps in
% a sortable format such as:
%
%   rain_2025_06_01_00_00.tif
%   rain_2025_06_01_03_00.tif
%
% or
%
%   GPM_180min_2025_06_01_00_00.tif
%
% If filenames follow a YYYY_MM_DD_HH_mm pattern, this script can sort them
% robustly from oldest to newest or newest to oldest.
% ========================================================================

clear InputPaths
InputPaths = struct();

%% ========================================================================
% 1) GENERAL SETTINGS
% ========================================================================

% -------------------------------------------------------------------------
% Root folder for this case
% -------------------------------------------------------------------------
InputPaths.case_root = ...
    '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Pune_Clipped';

% -------------------------------------------------------------------------
% Forcing file order
%
% Options:
%   'ascend'  -> oldest to newest  (recommended for simulation)
%   'descend' -> newest to oldest
% -------------------------------------------------------------------------
InputPaths.forcing_sort_order = 'ascend';

% -------------------------------------------------------------------------
% Accepted raster extensions
% -------------------------------------------------------------------------
InputPaths.raster_extensions = {'.tif', '.tiff'};

%% ========================================================================
% 2) TOOL PATHS
% ========================================================================

InputPaths.topo_path = ...
    '/oak/stanford/groups/gorelick/HydroPol2D/Dependencies/topotoolbox';

InputPaths.hydropol2d_tools = ...
    '/oak/stanford/groups/gorelick/HydroPol2D/HydroPol2D_Functions';

%% ========================================================================
% 3) STATIC RASTER PATHS
% ========================================================================

static_root = fullfile(InputPaths.case_root, 'Static');

InputPaths.DEM_path   = fullfile(static_root, 'DEM.tif');
InputPaths.LULC_path  = fullfile(static_root, 'LULC.tif');
InputPaths.SOIL_path  = fullfile(static_root, 'SOIL.tif');

% Optional static rasters
InputPaths.DTB_path            = fullfile(static_root, 'DTB.tif');
InputPaths.LAI_path            = fullfile(static_root, 'LAI.tif');
InputPaths.Albedo_path         = fullfile(static_root, 'Albedo.tif');
InputPaths.Subgrid_DEM_path    = fullfile(static_root, 'Subgrid_DEM.tif');
InputPaths.RiverWidths_path    = fullfile(static_root, 'RiverWidths.tif');
InputPaths.RiverDepths_path    = fullfile(static_root, 'RiverDepths.tif');

% Optional warmup / initial-condition rasters
InputPaths.Warmup_Depth_path          = fullfile(static_root, 'Warmup_Depth.tif');
InputPaths.Initial_Buildup_path       = fullfile(static_root, 'Initial_Buildup.tif');
InputPaths.Initial_Soil_Moisture_path = fullfile(static_root, 'Initial_Soil_Moisture.tif');

% Optional water-quality rasters
InputPaths.B1_path = fullfile(static_root, 'B1.tif');
InputPaths.B2_path = fullfile(static_root, 'B2.tif');
InputPaths.W1_path = fullfile(static_root, 'W1.tif');
InputPaths.W2_path = fullfile(static_root, 'W2.tif');

%% ========================================================================
% 4) FORCING FOLDERS
% ========================================================================

forcing_root = fullfile(InputPaths.case_root, 'Forcing');

InputPaths.Rainfall_Rasters_Folder      = fullfile(forcing_root, 'Rainfall');
InputPaths.Transpiration_Rasters_Folder = fullfile(forcing_root, 'Transpiration');
InputPaths.Evaporation_Rasters_Folder   = fullfile(forcing_root, 'Evaporation');

InputPaths.Inflow_Hydrograph_CSV = fullfile(forcing_root, 'Inflow', 'inflow_hydrograph.csv');
InputPaths.Stage_Hydrograph_CSV  = fullfile(forcing_root, 'Stage',  'stage_hydrograph.csv');
InputPaths.Observed_Gauges_CSV   = fullfile(forcing_root, 'Observed_Gauges', 'observed_gauges.csv');

%% ========================================================================
% 5) DISCOVER RASTER FORCING FILES
% ========================================================================

InputPaths.Rainfall_Raster_Files = list_and_sort_rasters( ...
    InputPaths.Rainfall_Rasters_Folder, ...
    InputPaths.raster_extensions, ...
    InputPaths.forcing_sort_order);

InputPaths.Transpiration_Raster_Files = list_and_sort_rasters( ...
    InputPaths.Transpiration_Rasters_Folder, ...
    InputPaths.raster_extensions, ...
    InputPaths.forcing_sort_order);

InputPaths.Evaporation_Raster_Files = list_and_sort_rasters( ...
    InputPaths.Evaporation_Rasters_Folder, ...
    InputPaths.raster_extensions, ...
    InputPaths.forcing_sort_order);

%% ========================================================================
% 6) OPTIONAL: DISPLAY SUMMARY
% ========================================================================

fprintf('\n============================================================\n');
fprintf('HydroPol2D | input_paths_bypass_script loaded successfully\n');
fprintf('============================================================\n');
fprintf('Case root: %s\n', InputPaths.case_root);
fprintf('Sort order: %s\n', InputPaths.forcing_sort_order);
fprintf('DEM:   %s\n', InputPaths.DEM_path);
fprintf('LULC:  %s\n', InputPaths.LULC_path);
fprintf('SOIL:  %s\n', InputPaths.SOIL_path);
fprintf('Rainfall folder:      %s\n', InputPaths.Rainfall_Rasters_Folder);
fprintf('Transpiration folder: %s\n', InputPaths.Transpiration_Rasters_Folder);
fprintf('Evaporation folder:   %s\n', InputPaths.Evaporation_Rasters_Folder);
fprintf('Inflow CSV:         %s\n', InputPaths.Inflow_Hydrograph_CSV);
fprintf('Stage CSV:          %s\n', InputPaths.Stage_Hydrograph_CSV);
fprintf('Observed Gauges CSV:%s\n', InputPaths.Observed_Gauges_CSV);
fprintf('Rainfall rasters found:      %d\n', numel(InputPaths.Rainfall_Raster_Files));
fprintf('Transpiration rasters found: %d\n', numel(InputPaths.Transpiration_Raster_Files));
fprintf('Evaporation rasters found:   %d\n', numel(InputPaths.Evaporation_Raster_Files));
fprintf('============================================================\n\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function file_list = list_and_sort_rasters(folder_path, extensions, sort_order)
%LIST_AND_SORT_RASTERS
% Returns a cell array of full file paths sorted chronologically when
% possible. If no timestamp can be parsed for all files, falls back to
% filename sorting.
%
% PERFORMANCE NOTES
% -------------------------------------------------------------------------
% This version is optimized for large forcing folders on network filesystems:
%
%   1) avoids repeated struct concatenation during dir() calls
%   2) removes directory entries early
%   3) preallocates datetime arrays
%   4) sorts using parsed datetimes only when available for all files
%   5) otherwise falls back to alphabetical filename sorting
%
% INPUTS
% -------------------------------------------------------------------------
% folder_path : folder containing raster files
% extensions  : cell array of extensions, e.g. {'.tif','.tiff'}
% sort_order  : 'ascend' or 'descend'
%
% OUTPUT
% -------------------------------------------------------------------------
% file_list   : cell array of full file paths
%
% EXAMPLES
% -------------------------------------------------------------------------
% file_list = list_and_sort_rasters('/path/to/Rainfall');
% file_list = list_and_sort_rasters('/path/to/Rainfall', {'.tif','.tiff'});
% file_list = list_and_sort_rasters('/path/to/Rainfall', {'.tif'}, 'descend');

    file_list = {};

    % ---------------------------------------------------------------------
    % Defaults
    % ---------------------------------------------------------------------
    if nargin < 2 || isempty(extensions)
        extensions = {'.tif', '.tiff'};
    end

    if nargin < 3 || isempty(sort_order)
        sort_order = 'ascend';
    end

    % ---------------------------------------------------------------------
    % Input validation
    % ---------------------------------------------------------------------
    if ~(ischar(folder_path) || isstring(folder_path))
        return;
    end

    folder_path = char(folder_path);
    sort_order  = char(sort_order);

    if ~exist(folder_path, 'dir')
        warning('Forcing folder not found: %s', folder_path);
        return;
    end

    if ~iscell(extensions)
        extensions = cellstr(extensions);
    end

    % Make sure extensions are char row vectors and start with '.'
    for i = 1:numel(extensions)
        extensions{i} = char(extensions{i});
        if ~startsWith(extensions{i}, '.')
            extensions{i} = ['.' extensions{i}];
        end
    end

    % ---------------------------------------------------------------------
    % Directory listing without repeated concatenation
    % ---------------------------------------------------------------------
    dir_parts = cell(numel(extensions), 1);

    for i = 1:numel(extensions)
        ext = extensions{i};
        dir_parts{i} = dir(fullfile(folder_path, ['*' ext]));
    end

    non_empty = ~cellfun(@isempty, dir_parts);

    if any(non_empty)
        files_all = vertcat(dir_parts{non_empty});
    else
        files_all = struct( ...
            'name',   {}, ...
            'folder', {}, ...
            'date',   {}, ...
            'bytes',  {}, ...
            'isdir',  {}, ...
            'datenum',{});
    end

    % ---------------------------------------------------------------------
    % Exit early if nothing found
    % ---------------------------------------------------------------------
    if isempty(files_all)
        return;
    end

    % Remove directories if any appear
    files_all = files_all(~[files_all.isdir]);

    if isempty(files_all)
        return;
    end

    % ---------------------------------------------------------------------
    % Build names and full paths
    % ---------------------------------------------------------------------
    n_files = numel(files_all);

    names = {files_all.name}';
    fulls = fullfile({files_all.folder}', names);

    % ---------------------------------------------------------------------
    % Parse timestamps from filenames
    % ---------------------------------------------------------------------
    t = NaT(n_files, 1);

    for i = 1:n_files
        t(i) = parse_datetime_from_filename(names{i});
    end

    % ---------------------------------------------------------------------
    % Sort files
    % ---------------------------------------------------------------------
    if all(~isnat(t))
        [~, idx] = sort(t);
    else
        [~, idx] = sort(lower(string(names)));
    end

    if strcmpi(sort_order, 'descend')
        idx = flipud(idx(:));
    else
        idx = idx(:);
    end

    % ---------------------------------------------------------------------
    % Return sorted full paths
    % ---------------------------------------------------------------------
    file_list = fulls(idx);
end

function dt = parse_datetime_from_filename(fname)
%PARSE_DATETIME_FROM_FILENAME
% Attempts to extract a datetime from filenames like:
%   rain_2025_06_01_00_00.tif
%   GPM_180min_2025_06_01_03_00.tif
%   2025-06-01-03-00.tif
%
% Returns NaT if parsing fails.

    dt = NaT;

    if ~(ischar(fname) || isstring(fname))
        return;
    end

    fname = char(fname);

    % Pattern 1: YYYY_MM_DD_HH_mm
    token = regexp(fname, '(\d{4})[_-](\d{2})[_-](\d{2})[_-](\d{2})[_-](\d{2})', ...
        'tokens', 'once');

    if isempty(token)
        % Pattern 2: YYYY_MM_DD only
        token = regexp(fname, '(\d{4})[_-](\d{2})[_-](\d{2})', ...
            'tokens', 'once');

        if isempty(token)
            return;
        else
            yyyy = str2double(token{1});
            mm   = str2double(token{2});
            dd   = str2double(token{3});
            HH   = 0;
            MM   = 0;
        end
    else
        yyyy = str2double(token{1});
        mm   = str2double(token{2});
        dd   = str2double(token{3});
        HH   = str2double(token{4});
        MM   = str2double(token{5});
    end

    try
        dt = datetime(yyyy, mm, dd, HH, MM, 0);
    catch
        dt = NaT;
    end
end