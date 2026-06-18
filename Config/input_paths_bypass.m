%% ========================================================================
% input_paths_bypass_script.m
% ========================================================================
%
% HydroPol2D | External file-path definition for full bypass workflow
%
% PURPOSE
% -------------------------------------------------------------------------
% This function defines the locations of all external files and folders used
% by HydroPol2D when running in bypass mode.
%
% DESIGN
% -------------------------------------------------------------------------
% This version supports:
%   1) default case-root based paths
%   2) optional user overrides through the input struct "Overrides"
%
% This allows maximum flexibility while preserving a clean default layout.
%
% USAGE
% -------------------------------------------------------------------------
% Example 1 — use all defaults:
%
%   InputPaths = input_paths_bypass(topo_path, hydropol2d_tools);
%
% Example 2 — override only rainfall folder:
%
%   Overrides = struct();
%   Overrides.Rainfall_Rasters_Folder = ...
%       '/oak/stanford/groups/gorelick/Marcus/MyRainfallFolder';
%
%   InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides);
%
% Example 3 — override multiple files:
%
%   Overrides = struct();
%   Overrides.DEM_path = '/path/to/MyDEM.tif';
%   Overrides.LULC_path = '/path/to/MyLULC.tif';
%   Overrides.Rainfall_Rasters_Folder = '/path/to/Rainfall';
%   Overrides.Inflow_Hydrograph_CSV = '/path/to/inflow.csv';
%
%   InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides);
%
% OUTPUT
% -------------------------------------------------------------------------
% Returns one struct:
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
% sortable formats such as:
%
%   rain_2025_06_01_00_00.tif
%   GPM_180min_2025_06_01_03_00.tif
%   2025-06-01-03-00.tif
%
% If filenames follow YYYY_MM_DD_HH_mm or YYYY-MM-DD-HH-mm, this function
% will sort them chronologically.
% ========================================================================

function InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides)

InputPaths = struct();

% -------------------------------------------------------------------------
% Validate required inputs
% -------------------------------------------------------------------------
if nargin < 2
    error('input_paths_bypass requires at least two inputs: topo_path and hydropol2d_tools.');
end

if nargin < 3 || isempty(Overrides)
    Overrides = struct();
end

if ~(ischar(topo_path) || isstring(topo_path))
    error('topo_path must be a character array or string.');
end

if ~(ischar(hydropol2d_tools) || isstring(hydropol2d_tools))
    error('hydropol2d_tools must be a character array or string.');
end

if ~isstruct(Overrides)
    error('Overrides must be a struct.');
end

topo_path = char(topo_path);
hydropol2d_tools = char(hydropol2d_tools);

if ~exist(topo_path, 'dir')
    error('TopoToolbox folder not found:\n  %s', topo_path);
end

if ~exist(hydropol2d_tools, 'dir')
    error('HydroPol2D tools folder not found:\n  %s', hydropol2d_tools);
end

%% ========================================================================
% 1) GENERAL SETTINGS
% ========================================================================

% -------------------------------------------------------------------------
% Root folder for this case (auto-detected from this script location)
% -------------------------------------------------------------------------
this_file = mfilename('fullpath');
config_folder = fileparts(this_file);
case_root = fileparts(config_folder);

InputPaths.case_root = case_root;

% -------------------------------------------------------------------------
% Forcing file order
%
% Default:
%   'ascend'  -> oldest to newest
%
% Override field:
%   Overrides.forcing_sort_order = 'descend';
% -------------------------------------------------------------------------
InputPaths.forcing_sort_order = get_override( ...
    Overrides, 'forcing_sort_order', 'ascend');

% -------------------------------------------------------------------------
% Accepted raster extensions
%
% Default:
%   {'.tif', '.tiff'}
%
% Override field:
%   Overrides.raster_extensions = {'.tif','.tiff','.asc'};
% -------------------------------------------------------------------------
InputPaths.raster_extensions = get_override( ...
    Overrides, 'raster_extensions', {'.tif', '.tiff'});

%% ========================================================================
% 2) TOOL PATHS
% ========================================================================

InputPaths.topo_path = topo_path;
InputPaths.hydropol2d_tools = hydropol2d_tools;

%% ========================================================================
% 3) DEFAULT ROOTS
% ========================================================================

static_root  = fullfile(InputPaths.case_root, 'Static');
forcing_root = fullfile(InputPaths.case_root, 'Forcing');

%% ========================================================================
% 4) RESOLVE STATIC RASTER PATHS
% ========================================================================

InputPaths.DEM_path  = get_override(Overrides, 'DEM_path',  fullfile(static_root, 'DEM.tif'));
InputPaths.LULC_path = get_override(Overrides, 'LULC_path', fullfile(static_root, 'LULC.tif'));
InputPaths.SOIL_path = get_override(Overrides, 'SOIL_path', fullfile(static_root, 'SOIL.tif'));

% Optional static rasters
InputPaths.DTB_path = get_override( ...
    Overrides, 'DTB_path', fullfile(static_root, 'DTB.tif'));

InputPaths.GW_table_path = get_override( ...
    Overrides, 'GW_table_path', fullfile(static_root, 'GW_table.tif'));

InputPaths.LAI_path = get_override( ...
    Overrides, 'LAI_path', fullfile(static_root, 'LAI.tif'));

InputPaths.Albedo_path = get_override( ...
    Overrides, 'Albedo_path', fullfile(static_root, 'Albedo.tif'));

InputPaths.Subgrid_DEM_path = get_override(Overrides, 'Subgrid_DEM_path', fullfile(static_root, 'Subgrid_DEM.tif'));
InputPaths.RiverWidths_path = get_override(Overrides, 'RiverWidths_path', fullfile(static_root, 'RiverWidths.tif'));
InputPaths.RiverDepths_path = get_override(Overrides, 'RiverDepths_path', fullfile(static_root, 'RiverDepths.tif'));

% Optional warmup / initial-condition rasters
InputPaths.Warmup_Depth_path = get_override( ...
    Overrides, 'Warmup_Depth_path', fullfile(static_root, 'Warmup_Depth.tif'));

InputPaths.Initial_Buildup_path = get_override( ...
    Overrides, 'Initial_Buildup_path', fullfile(static_root, 'Initial_Buildup.tif'));

InputPaths.Initial_Soil_Moisture_path = get_override( ...
    Overrides, 'Initial_Soil_Moisture_path', fullfile(static_root, 'Initial_Soil_Moisture.tif'));

% Optional water-quality rasters
InputPaths.B1_path = get_override(Overrides, 'B1_path', fullfile(static_root, 'B1.tif'));
InputPaths.B2_path = get_override(Overrides, 'B2_path', fullfile(static_root, 'B2.tif'));
InputPaths.W1_path = get_override(Overrides, 'W1_path', fullfile(static_root, 'W1.tif'));
InputPaths.W2_path = get_override(Overrides, 'W2_path', fullfile(static_root, 'W2.tif'));


%% ========================================================================
% 5) RESOLVE FORCING FOLDERS
% ========================================================================

InputPaths.Rainfall_Rasters_Folder = get_override( ...
    Overrides, 'Rainfall_Rasters_Folder', fullfile(forcing_root, 'Rainfall'));

InputPaths.Transpiration_Rasters_Folder = get_override( ...
    Overrides, 'Transpiration_Rasters_Folder', fullfile(forcing_root, 'Transpiration'));

InputPaths.Evaporation_Rasters_Folder = get_override( ...
    Overrides, 'Evaporation_Rasters_Folder', fullfile(forcing_root, 'Evaporation'));

%% ========================================================================
% 6) RESOLVE CSV / TABULAR INPUTS
% ========================================================================

InputPaths.Inflow_Hydrograph_CSV = get_override( ...
    Overrides, 'Inflow_Hydrograph_CSV', fullfile(forcing_root, 'Inflow', 'inflow_hydrograph.csv'));

InputPaths.Outlet_Cells_CSV = get_override( ...
    Overrides, 'Outlet_Cells_CSV', fullfile(forcing_root, 'Outlet', 'outlet_cells.csv'));

InputPaths.Stage_Hydrograph_CSV = get_override( ...
    Overrides, 'Stage_Hydrograph_CSV', fullfile(forcing_root, 'Stage', 'stage_hydrograph.csv'));

InputPaths.Observed_Gauges_CSV = get_override( ...
    Overrides, 'Observed_Gauges_CSV', fullfile(forcing_root, 'Observed_Gauges', 'observed_gauges.csv'));

% Non-raster meteorological forcing for INTERNAL ETP calculation
% Used when:
%   flags.flag_ETP = 1
%   flags.flag_input_ETP_map = 0
InputPaths.ETP_input_spreadsheet = get_override( ...
    Overrides, 'ETP_input_spreadsheet', ...
    fullfile(forcing_root, 'Evapotranspiration', 'ETP_input_data.xlsx'));

% Non-raster lumped rainfall forcing file
% Used when:
%   flags.flag_rainfall = 1
%   flags.flag_spatial_rainfall = 0
%   flags.flag_input_rainfall_map = 0
%   flags.flag_satellite_rainfall = 0
%   flags.flag_real_time_satellite_rainfall = 0
%   flags.flag_alternated_blocks = 0
%   flags.flag_huff = 0
InputPaths.Rainfall_Timeseries_File = get_override( ...
    Overrides, 'Rainfall_Timeseries_File', ...
    fullfile(forcing_root, 'Rainfall', 'Rainfall_Intensity_Data.xlsx'));

%% ========================================================================
% 7) DISCOVER RASTER FORCING FILES
% ========================================================================
% IMPORTANT
% -------------------------------------------------------------------------
% Rainfall rasters are NOT scanned here.
%
% Reason:
%   Long rainfall archives (e.g. 30-min data for decades) may contain a
%   very large number of files, making directory listing too slow on
%   network storage.
%
% Instead, rainfall file paths are generated later in
% input_data_bypass_script.m from:
%   - InputData_Bypass.general.date_begin
%   - InputData_Bypass.general.date_end
%   - InputData_Bypass.general.dt_rainfall_maps_min
%   - InputData_Bypass.general.rainfall_filename_example
%
% Therefore, here we only keep the rainfall folder path and leave the file
% list empty.
InputPaths.Rainfall_Raster_Files = {};

% Transpiration / evaporation folders are usually much smaller, so they can
% still be scanned directly.
InputPaths.Transpiration_Raster_Files = list_and_sort_rasters( ...
    InputPaths.Transpiration_Rasters_Folder, ...
    InputPaths.raster_extensions, ...
    InputPaths.forcing_sort_order);

InputPaths.Evaporation_Raster_Files = list_and_sort_rasters( ...
    InputPaths.Evaporation_Rasters_Folder, ...
    InputPaths.raster_extensions, ...
    InputPaths.forcing_sort_order);

%% ========================================================================
% 8) DISPLAY SUMMARY
% ========================================================================

fprintf('\n============================================================\n');
fprintf('HydroPol2D | input_paths_bypass loaded successfully\n');
fprintf('============================================================\n');
fprintf('Case root: %s\n', InputPaths.case_root);
fprintf('Sort order: %s\n', string(InputPaths.forcing_sort_order));

fprintf('\n--- REQUIRED STATIC RASTERS ---\n');
print_path_status('DEM',  InputPaths.DEM_path,  'file');
print_path_status('LULC', InputPaths.LULC_path, 'file');
print_path_status('SOIL', InputPaths.SOIL_path, 'file');

fprintf('\n--- OPTIONAL STATIC RASTERS ---\n');
print_path_status('DTB',         InputPaths.DTB_path,         'file');
print_path_status('LAI',         InputPaths.LAI_path,         'file');
print_path_status('Albedo',      InputPaths.Albedo_path,      'file');
print_path_status('Subgrid DEM', InputPaths.Subgrid_DEM_path, 'file');
print_path_status('RiverWidths', InputPaths.RiverWidths_path, 'file');
print_path_status('RiverDepths', InputPaths.RiverDepths_path, 'file');

fprintf('\n--- WARMUP / INITIAL CONDITIONS ---\n');
print_path_status('Warmup Depth',          InputPaths.Warmup_Depth_path,          'file');
print_path_status('Initial Buildup',       InputPaths.Initial_Buildup_path,       'file');
print_path_status('Initial Soil Moisture', InputPaths.Initial_Soil_Moisture_path, 'file');

fprintf('\n--- WATER QUALITY RASTERS ---\n');
print_path_status('B1', InputPaths.B1_path, 'file');
print_path_status('B2', InputPaths.B2_path, 'file');
print_path_status('W1', InputPaths.W1_path, 'file');
print_path_status('W2', InputPaths.W2_path, 'file');

fprintf('\n--- FORCING FOLDERS ---\n');
print_path_status('Rainfall folder',      InputPaths.Rainfall_Rasters_Folder,      'dir');
print_path_status('Transpiration folder', InputPaths.Transpiration_Rasters_Folder, 'dir');
print_path_status('Evaporation folder',   InputPaths.Evaporation_Rasters_Folder,   'dir');

fprintf('\n--- CSV / TABULAR INPUTS ---\n');
print_path_status('Inflow CSV',          InputPaths.Inflow_Hydrograph_CSV, 'file');
print_path_status('Outlet cells CSV',    InputPaths.Outlet_Cells_CSV, 'file');
print_path_status('Stage CSV',           InputPaths.Stage_Hydrograph_CSV,  'file');
print_path_status('Observed Gauges CSV', InputPaths.Observed_Gauges_CSV,   'file');
print_path_status('ETP spreadsheet',     InputPaths.ETP_input_spreadsheet, 'file');
print_path_status('Rainfall timeseries', InputPaths.Rainfall_Timeseries_File, 'file');

fprintf('\n--- RASTER FORCING COUNTS ---\n');
fprintf('Rainfall rasters found:      %d\n', numel(InputPaths.Rainfall_Raster_Files));
fprintf('Transpiration rasters found: %d\n', numel(InputPaths.Transpiration_Raster_Files));
fprintf('Evaporation rasters found:   %d\n', numel(InputPaths.Evaporation_Raster_Files));
fprintf('============================================================\n\n');

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function value = get_override(S, field_name, default_value)
%GET_OVERRIDE
% Returns S.(field_name) if it exists and is not empty.
% Otherwise returns default_value.

    if isstruct(S) && isfield(S, field_name) && ~isempty(S.(field_name))
        value = S.(field_name);
    else
        value = default_value;
    end

    % Convert strings to char for path-like entries
    if isstring(value) && isscalar(value)
        value = char(value);
    end
end

function print_path_status(label, path_value, kind)
%PRINT_PATH_STATUS
% Prints the resolved path and whether it exists.

    if nargin < 3
        kind = 'file';
    end

    exists_flag = false;

    if strcmpi(kind, 'dir')
        exists_flag = exist(path_value, 'dir') == 7;
    else
        exists_flag = exist(path_value, 'file') == 2;
    end

    fprintf('%-22s: %s\n', label, path_value);
    fprintf('%-22s  exists = %d\n', '', exists_flag);
end

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

    % Pattern 1: YYYY_MM_DD_HH_mm or YYYY-MM-DD-HH-mm
    token = regexp(fname, ...
        '(\d{4})[_-](\d{2})[_-](\d{2})[_-](\d{2})[_-](\d{2})', ...
        'tokens', 'once');

    if isempty(token)
        % Pattern 2: YYYY_MM_DD or YYYY-MM-DD
        token = regexp(fname, ...
            '(\d{4})[_-](\d{2})[_-](\d{2})', ...
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
