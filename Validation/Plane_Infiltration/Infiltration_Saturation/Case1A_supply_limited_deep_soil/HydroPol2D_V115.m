%% ========================================================================
% HydroPol2D | Main Run Script
% Developer: Marcus Nobrega, Ph.D.
% Main launcher for HydroPol2D
% ----------------------------- Version 15.0 ------------------------------
% Last official model update: 07/18/2024
%
% PURPOSE
%   This script is the main entry point to run a complete HydroPol2D
%   simulation. It performs:
%       1) Pre-processing
%       2) Main numerical simulation
%       3) Post-processing
%
% FOR USERS
%   Edit ONLY the section called:
%       "USER INPUTS (EDIT ONLY THIS SECTION)"
%
%   Do NOT modify the remaining sections unless you are developing the
%   model itself.
%
% WORKFLOW
%   This script calls the following saved scripts:
%       - HydroPol2D_preprocessing : reads inputs and prepares all variables
%       - HydroPol2D_Main_While    : runs the main solver
%       - post_processing          : exports figures/tables/maps/animations
%
% IMPORTANT
%   If "clean_output_folder = true", the selected export folder will be
%   cleaned before the new simulation starts.
%
% BYPASS MODE
%   In bypass mode this launcher uses:
%
%       InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides)
%
%   together with:
%
%       input_data_bypass_script.m
%
%   where:
%       - input_paths_bypass.m resolves all file/folder locations
%       - Overrides is an optional struct with user-defined path overrides
%       - input_data_bypass_script.m defines the input-data content
% ========================================================================

%% Clean MATLAB environment
clear; clc;

%% ========================================================================
% USER INPUTS (EDIT ONLY THIS SECTION)
% ========================================================================

% -------------------------------------------------------------------------
% [MODE SELECTOR]
%   'excel'  -> legacy spreadsheet workflow
%   'bypass' -> Config/input_paths_bypass.m(topo_path, hydropol2d_tools, Overrides)
%               + Config/input_data_bypass_script.m
% -------------------------------------------------------------------------
run_mode = 'bypass';   % 'excel' or 'bypass'

% -------------------------------------------------------------------------
% Output / run behavior flags
% -------------------------------------------------------------------------
clean_output_folder = true;   % true = wipe Outputs folder before run
run_postprocessing  = true;   % run post_processing step
enable_logging      = true;   % write log file

%% ========================================================================
% SECTION A — EXCEL MODE INPUTS
% Only used when run_mode = 'excel'
% ========================================================================

% Main HydroPol2D Excel input file
input_excel_file = ...
    '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India_Full/Input_Data_Sheets/General_Data.xlsx';

% Optional folder containing input spreadsheets
add_input_sheets_to_path = true;
input_sheets_folder      = 'Input_Data_Sheets';

%% ========================================================================
% SECTION B — BYPASS MODE INPUTS
% Only used when run_mode = 'bypass'
% ========================================================================

% -------------------------------------------------------------------------
% Path to the function that builds InputPaths for this case
% -------------------------------------------------------------------------
input_paths_function = 'Config/input_paths_bypass.m';

% -------------------------------------------------------------------------
% Path to the case-specific MATLAB bypass input-data script
% -------------------------------------------------------------------------
input_data_bypass_script_path = 'Config/input_data_bypass_script.m';

% -------------------------------------------------------------------------
% Global tool folders (machine-specific, not case-specific)
% -------------------------------------------------------------------------
topo_path_user = ...
    '/oak/stanford/groups/gorelick/HydroPol2D/Topotoolbox/topotoolbox-master';

hydropol2d_tools_user = ...
    '/oak/stanford/groups/gorelick/HydroPol2D/HydroPol2D_Functions';

% -------------------------------------------------------------------------
% Optional path/file overrides for bypass mode
%
% DEFAULT CASE STRUCTURE (NO OVERRIDES)
% -------------------------------------------------------------------------
% If Overrides is empty, the model assumes a standard case folder layout:
%
%   CaseFolder/
%       Static/
%           DEM.tif
%           LULC.tif
%           SOIL.tif
%           DTB.tif
%           LAI.tif
%           Albedo.tif
%           Subgrid_DEM.tif
%           RiverWidths.tif
%           RiverDepths.tif
%           Warmup_Depth.tif
%           Initial_Buildup.tif
%           Initial_Soil_Moisture.tif
%           B1.tif, B2.tif, W1.tif, W2.tif
%
%       Forcing/
%           Rainfall/               -> rainfall rasters (*.tif)
%           Transpiration/          -> transpiration rasters (*.tif)
%           Evaporation/            -> evaporation rasters (*.tif)
%
%           Inflow/
%               inflow_hydrograph.csv
%
%           Stage/
%               stage_hydrograph.csv
%
%           Observed_Gauges/
%               observed_gauges.csv
%
%           Evapotranspiration/
%               ETP_input_data.xlsx   (used for INTERNAL ETP calculation)
%
% OUTPUT STRUCTURE (DEFAULT)
% -------------------------------------------------------------------------
% If no override is provided, results are written to:
%
%   ./Outputs/
%       Modeling_Results/
%       Temporary_Files/
%       Logs/
%
% PURPOSE OF OVERRIDES
% -------------------------------------------------------------------------
% The Overrides struct allows the user to:
%   - relocate any input file or folder
%   - use external datasets (e.g., different rainfall directory)
%   - redirect output results to a custom location
%
% Overrides are passed to:
%
%   InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides)
%
% and take precedence over the default structure above.
%
% HOW TO USE
% -------------------------------------------------------------------------
% Uncomment only the fields you want to override.
%
% EXAMPLES
%   Overrides.DEM_path = '/oak/.../DEM.tif';
%   Overrides.Rainfall_Rasters_Folder = '/oak/.../Rainfall';
%   Overrides.ETP_input_spreadsheet = '/oak/.../ETP_input_data.xlsx';
%   Overrides.export_root_dir = '/oak/.../Outputs';
%
% IMPORTANT NOTES
% -------------------------------------------------------------------------
% 1) If a field is NOT defined in Overrides → default folder structure is used
% 2) Raster folders must contain consistent GeoTIFFs (.tif)
% 3) Time series files (CSV/XLSX) must match expected format
% 4) Output override ONLY affects where results are written (not inputs)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
Overrides = struct();

% ===== Example overrides (leave commented unless needed) ==================
% Overrides.DEM_path                     = '/oak/stanford/groups/gorelick/Marcus/Case/Static/DEM.tif';
% Overrides.LULC_path                    = '/oak/stanford/groups/gorelick/Marcus/Case/Static/LULC.tif';
% Overrides.SOIL_path                    = '/oak/stanford/groups/gorelick/Marcus/Case/Static/SOIL.tif';
% Overrides.DTB_path                     = '/oak/stanford/groups/gorelick/Marcus/Case/Static/DTB.tif';
% Overrides.LAI_path                     = '/oak/stanford/groups/gorelick/Marcus/Case/Static/LAI.tif';
% Overrides.Albedo_path                  = '/oak/stanford/groups/gorelick/Marcus/Case/Static/Albedo.tif';
% Overrides.Subgrid_DEM_path             = '/oak/stanford/groups/gorelick/Marcus/Case/Static/Subgrid_DEM.tif';
% Overrides.RiverWidths_path             = '/oak/stanford/groups/gorelick/Marcus/Case/Static/RiverWidths.tif';
% Overrides.RiverDepths_path             = '/oak/stanford/groups/gorelick/Marcus/Case/Static/RiverDepths.tif';
%
% Overrides.Warmup_Depth_path            = '/oak/stanford/groups/gorelick/Marcus/Case/Static/Warmup_Depth.tif';
% Overrides.Initial_Buildup_path         = '/oak/stanford/groups/gorelick/Marcus/Case/Static/Initial_Buildup.tif';
% Overrides.Initial_Soil_Moisture_path   = '/oak/stanford/groups/gorelick/Marcus/Case/Static/Initial_Soil_Moisture.tif';
%
% Overrides.B1_path                      = '/oak/stanford/groups/gorelick/Marcus/Case/Static/B1.tif';
% Overrides.B2_path                      = '/oak/stanford/groups/gorelick/Marcus/Case/Static/B2.tif';
% Overrides.W1_path                      = '/oak/stanford/groups/gorelick/Marcus/Case/Static/W1.tif';
% Overrides.W2_path                      = '/oak/stanford/groups/gorelick/Marcus/Case/Static/W2.tif';

%
Overrides.Rainfall_Rasters_Folder      = '/oak/stanford/groups/gorelick/Marcus/India/IMERG_30min_mmhr';
% Overrides.Transpiration_Rasters_Folder = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Transpiration';
% Overrides.Evaporation_Rasters_Folder   = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Evaporation';
%
% Overrides.Inflow_Hydrograph_CSV        = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Inflow/inflow_hydrograph.csv';
% Overrides.Stage_Hydrograph_CSV         = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Stage/stage_hydrograph.csv';
% Overrides.Observed_Gauges_CSV          = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Observed_Gauges/observed_gauges.csv';
% Overrides.ETP_input_spreadsheet        = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Evapotranspiration/ETP_input_data.xlsx';
% Overrides.Rainfall_Timeseries_File     = '/oak/stanford/groups/gorelick/Marcus/Case/Forcing/Rainfall/Rainfall_Intensity_Data.xlsx';

% % Overrides.export_root_dir             = '/oak/stanford/groups/gorelick/Marcus/Case/Outputs';
%
% Overrides.forcing_sort_order           = 'ascend';     % or 'descend'
% Overrides.raster_extensions            = {'.tif','.tiff'};
% ========================================================================

%% ========================================================================
% SECTION C — COMMON RUN / OUTPUT OPTIONS
% Used in both modes
% ========================================================================

% Export root directory for this run
% export_root_dir = ...
%     '/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/IndiaFull/2500m_CPU/Modeling_Results';
if ~exist('export_root_dir', 'var') || isempty(export_root_dir)
    export_root_dir = fullfile(pwd, 'Outputs');
end

% -------------------------------------------------------------------------
% Optional output-folder override from bypass mode
% -------------------------------------------------------------------------
% If provided through Overrides, it takes precedence over the local
% export_root_dir defined above.
if strcmpi(run_mode,'bypass') && isfield(Overrides,'export_root_dir') && ...
        ~isempty(Overrides.export_root_dir)
    export_root_dir = char(Overrides.export_root_dir);
end

if ~isfolder(export_root_dir)
    mkdir(export_root_dir);
end

%% ========================================================================
% INTERNAL MODEL LOGIC (DO NOT EDIT BELOW THIS LINE UNLESS DEVELOPING)
% ========================================================================

fprintf('\n============================================================\n');
fprintf('               HydroPol2D | Simulation Launcher             \n');
fprintf('============================================================\n');
fprintf('Run mode         : %s\n', run_mode);

if strcmpi(run_mode,'excel')
    fprintf('Input Excel file : %s\n', input_excel_file);
else
    fprintf('Input paths func : %s\n', input_paths_function);
    fprintf('Input data script: %s\n', input_data_bypass_script_path);
    fprintf('TopoToolbox path : %s\n', topo_path_user);
    fprintf('HP2D tools path  : %s\n', hydropol2d_tools_user);
    fprintf('Overrides fields : %d\n', numel(fieldnames(Overrides)));
end

fprintf('Export root dir  : %s\n', export_root_dir);
fprintf('Clean outputs    : %d\n', clean_output_folder);
fprintf('Run post-process : %d\n', run_postprocessing);
fprintf('Logging enabled  : %d\n', enable_logging);
fprintf('============================================================\n\n');

%% ------------------------------------------------------------------------
% Validate main user inputs
% -------------------------------------------------------------------------
if ~(ischar(run_mode) || isstring(run_mode))
    error('The variable "run_mode" must be a character array or string.');
end

if ~(ischar(export_root_dir) || isstring(export_root_dir))
    error('The variable "export_root_dir" must be a character array or string.');
end

run_mode        = char(run_mode);
export_root_dir = char(export_root_dir);

if strcmpi(run_mode,'excel')

    if ~(ischar(input_excel_file) || isstring(input_excel_file))
        error('The variable "input_excel_file" must be a character array or string.');
    end

    input_excel_file = char(input_excel_file);

    if ~exist(input_excel_file, 'file')
        error('Input Excel file not found:\n  %s', input_excel_file);
    end

elseif strcmpi(run_mode,'bypass')

    if ~(ischar(input_paths_function) || isstring(input_paths_function))
        error('The variable "input_paths_function" must be a character array or string.');
    end
    input_paths_function = char(input_paths_function);

    if ~exist(input_paths_function, 'file')
        error('Input paths bypass function not found:\n  %s', input_paths_function);
    end

    if ~(ischar(input_data_bypass_script_path) || isstring(input_data_bypass_script_path))
        error('The variable "input_data_bypass_script_path" must be a character array or string.');
    end
    input_data_bypass_script_path = char(input_data_bypass_script_path);

    if ~exist(input_data_bypass_script_path, 'file')
        error('Bypass input-data script not found:\n  %s', input_data_bypass_script_path);
    end

    if ~(ischar(topo_path_user) || isstring(topo_path_user))
        error('The variable "topo_path_user" must be a character array or string.');
    end

    if ~(ischar(hydropol2d_tools_user) || isstring(hydropol2d_tools_user))
        error('The variable "hydropol2d_tools_user" must be a character array or string.');
    end

    topo_path_user = char(topo_path_user);
    hydropol2d_tools_user = char(hydropol2d_tools_user);

    if ~exist(topo_path_user, 'dir')
        error('TopoToolbox folder not found:\n  %s', topo_path_user);
    end

    if ~exist(hydropol2d_tools_user, 'dir')
        error('HydroPol2D tools folder not found:\n  %s', hydropol2d_tools_user);
    end

    if ~isstruct(Overrides)
        error('In bypass mode, "Overrides" must be a struct.');
    end

else
    error('run_mode must be either "excel" or "bypass".');
end

%% ------------------------------------------------------------------------
% Optional: make spreadsheet folder visible to MATLAB
% Only relevant in Excel mode
% -------------------------------------------------------------------------
if strcmpi(run_mode,'excel') && add_input_sheets_to_path
    if exist(input_sheets_folder, 'dir')
        addpath(genpath(input_sheets_folder));
        fprintf('Added input sheets folder to path:\n  %s\n\n', input_sheets_folder);
    else
        warning(['HydroPol2D:InputSheetsFolderNotFound\n' ...
                 'Requested input_sheets_folder was not found:\n  %s\n' ...
                 'Continuing without adding it to the MATLAB path.'], input_sheets_folder);
    end
end

%% ------------------------------------------------------------------------
% Load configuration source
% -------------------------------------------------------------------------
if strcmpi(run_mode,'excel')

    fprintf('Reading "General_Data" sheet from input Excel file...\n');
    generalDataSheet = readcell(input_excel_file, 'Sheet', 'General_Data');
    fprintf('General_Data loaded successfully.\n\n');

    fprintf('Reading required tool paths from General_Data...\n');

    hydropol2d_tools = string(xlget(generalDataSheet, 'hydropol2d_tools'));
    topo_path        = string(xlget(generalDataSheet, 'topo_path'));

    model_folder = input_excel_file;
    GD = generalDataSheet;

    use_inputpaths_bypass = 0;
    use_inputdata_bypass  = 0;

elseif strcmpi(run_mode,'bypass')
    
    fprintf('Running input_paths_bypass function...\n');

    % Add the config folder so MATLAB can find the bypass function by name
    config_folder = fileparts(input_paths_function);
    if ~isempty(config_folder)
        addpath(config_folder);
    end

    % Extract the function name from the .m file path
    [~, func_name, ~] = fileparts(input_paths_function);

    % ---------------------------------------------------------------------
    % NEW FUNCTION CALL WITH OVERRIDES
    % ---------------------------------------------------------------------
    % This is the key update:
    %   InputPaths = input_paths_bypass(topo_path, hydropol2d_tools, Overrides)
    % ---------------------------------------------------------------------
    InputPaths = feval(func_name, topo_path_user, hydropol2d_tools_user, Overrides);

    if ~exist('InputPaths','var') || ~isstruct(InputPaths)
        error('Bypass paths function must return a struct named InputPaths.');
    end

    if ~isfield(InputPaths,'hydropol2d_tools')
        error('InputPaths.hydropol2d_tools is missing.');
    end

    if ~isfield(InputPaths,'topo_path')
        error('InputPaths.topo_path is missing.');
    end

    hydropol2d_tools = string(InputPaths.hydropol2d_tools);
    topo_path        = string(InputPaths.topo_path);

    model_folder = '';
    GD = [];

    use_inputpaths_bypass = 1;
    use_inputdata_bypass  = 1;

    fprintf('Bypass paths resolved successfully.\n\n');

else
    error('Unsupported run mode.');
end

%% ------------------------------------------------------------------------
% Validate and add required tool folders
% -------------------------------------------------------------------------
if strlength(strtrim(hydropol2d_tools)) == 0
    error('"hydropol2d_tools" is empty.');
end

if strlength(strtrim(topo_path)) == 0
    error('"topo_path" is empty.');
end

if ~exist(char(hydropol2d_tools), 'dir')
    error('HydroPol2D tools folder not found:\n  %s', char(hydropol2d_tools));
end

if ~exist(char(topo_path), 'dir')
    error('Topography/tools folder not found:\n  %s', char(topo_path));
end

addpath(genpath(char(hydropol2d_tools)));
addpath(genpath(char(topo_path)));

fprintf('Added HydroPol2D tools to path:\n  %s\n', char(hydropol2d_tools));
fprintf('Added topo/tools folder to path:\n  %s\n\n', char(topo_path));

%% ------------------------------------------------------------------------
% Initialize organized output folder tree
% -------------------------------------------------------------------------
fprintf('Initializing output folder tree...\n');

Paths = init_results_tree(export_root_dir, clean_output_folder);

% Keep compatibility with legacy scripts that expect these variables
ExportRootDir = Paths.Root;
resultsDir    = Paths.Results;

% Make variables visible to scripts that run in the base workspace
assignin('base','Paths',Paths);
assignin('base','ExportRootDir',ExportRootDir);
assignin('base','resultsDir',resultsDir);
assignin('base','model_folder',model_folder);
assignin('base','GD',GD);
assignin('base','use_inputpaths_bypass',use_inputpaths_bypass);
assignin('base','use_inputdata_bypass',use_inputdata_bypass);

if strcmpi(run_mode,'bypass')
    assignin('base','InputPaths',InputPaths);
    assignin('base','Overrides',Overrides);
    assignin('base','input_data_bypass_script_path',input_data_bypass_script_path);
    assignin('base','input_paths_function',input_paths_function);
end

fprintf('\n================ HydroPol2D Output Tree ================\n');
fprintf('Export root:           %s\n', Paths.Root);
fprintf('Modeling_Results:      %s\n', Paths.Results);
fprintf('Temporary_Files:       %s\n', Paths.Temp);
fprintf('Figures_PDF:           %s\n', Paths.FigPDF);
fprintf('Figures_FIG:           %s\n', Paths.FigFIG);
fprintf('Tables_CSV:            %s\n', Paths.Tables);
fprintf('Rasters_Water_Depths:  %s\n', Paths.RastersWD);
fprintf('Rasters_WSE:           %s\n', Paths.RastersWSE);
fprintf('Rasters_Static:        %s\n', Paths.RastersStatic);
fprintf('Rasters_WQ:            %s\n', Paths.WQMaps);
fprintf('Rasters_Human_Risk:    %s\n', Paths.HRMaps);
fprintf('GIFs_MP4:              %s\n', Paths.Anim);
fprintf('Shapefiles:            %s\n', Paths.Shapes);
fprintf('Logs:                  %s\n', Paths.Logs);
fprintf('========================================================\n\n');

%% ------------------------------------------------------------------------
% Start logging if requested
% -------------------------------------------------------------------------
if enable_logging
    try
        if ~isempty(Paths.Logs)
            diary('off');
            diary(fullfile(Paths.Logs, 'run_log.txt'));
            diary('on');
            fprintf('Run logging enabled:\n  %s\n\n', fullfile(Paths.Logs, 'run_log.txt'));
        end
    catch ME
        warning('Could not start diary logging')
    end
end

%% ------------------------------------------------------------------------
% Run HydroPol2D workflow
% -------------------------------------------------------------------------

% 1) Pre-processing
fprintf('------------------------------------------------------------\n');
fprintf('STEP 1/3 | Running HydroPol2D_preprocessing\n');
fprintf('This step reads inputs, loads rasters/tables, and prepares\n');
fprintf('all variables needed by the numerical solver.\n');
fprintf('------------------------------------------------------------\n\n');

HydroPol2D_preprocessing;

% 2) Main solver
fprintf('\n------------------------------------------------------------\n');
fprintf('STEP 2/3 | Running HydroPol2D_Main_While\n');
fprintf('This step executes the main HydroPol2D numerical simulation.\n');
fprintf('------------------------------------------------------------\n\n');

HydroPol2D_Main_While

% 3) Post-processing
if run_postprocessing
    fprintf('\n------------------------------------------------------------\n');
    fprintf('STEP 3/3 | Running post_processing\n');
    fprintf('This step exports figures, tables, rasters, animations,\n');
    fprintf('and other organized outputs.\n');
    fprintf('------------------------------------------------------------\n\n');

    close all
    post_processing;
else
    fprintf('\nPost-processing was skipped because "run_postprocessing = false".\n');
end

%% ------------------------------------------------------------------------
% Finish logging
% -------------------------------------------------------------------------
if enable_logging
    try
        diary('off');
    catch
    end
end

fprintf('\n============================================================\n');
fprintf('HydroPol2D run completed successfully.\n');
fprintf('Results saved under:\n  %s\n', Paths.Root);
fprintf('============================================================\n\n');

%% ========================================================================
% LOCAL HELPER FUNCTIONS
% ========================================================================

function v = xlget(GD, key)
%XLGET Find a label (key) anywhere in a readcell() grid and return the cell
% to the right. Search is case-insensitive.
%
% INPUTS
%   GD  : cell array returned by readcell(...)
%   key : label to search for
%
% OUTPUT
%   v   : value stored in the cell immediately to the right of the key
%
% EXAMPLE
%   topo_path = xlget(GD, 'topo_path');

    S = strings(size(GD));
    for r = 1:size(GD,1)
        for c = 1:size(GD,2)
            if ischar(GD{r,c}) || isstring(GD{r,c})
                S(r,c) = string(GD{r,c});
            end
        end
    end

    [rr, cc] = find(strcmpi(strtrim(S), key), 1, 'first');

    if isempty(rr)
        error("General_Data: key '%s' not found.", key);
    end

    if cc == size(GD,2)
        error("General_Data: key '%s' found in last column; no value exists to the right.", key);
    end

    v = GD{rr, cc+1};
end

function Paths = init_results_tree(exportRootDir, cleanOutputFolder)
%INIT_RESULTS_TREE Create the organized output folder tree expected by
% HydroPol2D post-processing.
%
% INPUTS
%   exportRootDir      : root directory for this run
%   cleanOutputFolder  : logical flag
%       true  -> delete existing contents inside exportRootDir
%       false -> preserve existing contents
%
% OUTPUT
%   Paths : structure with all output folders
%
% FOLDER STRUCTURE CREATED
%   exportRootDir/
%       Modeling_Results/
%           Figures_PDF/
%           Figures_FIG/
%           Tables_CSV/
%           Rasters_Water_Depths/
%           Rasters_WSE/
%           Rasters_Static/
%           Rasters_WQ/
%           Rasters_Human_Risk/
%           GIFs_MP4/
%           Shapefiles/
%       Temporary_Files/
%       Logs/

    arguments
        exportRootDir (1,:) char
        cleanOutputFolder (1,1) logical = true
    end

    % ---------------------------------------------------------------------
    % Prepare export root
    % ---------------------------------------------------------------------
    if ~exist(exportRootDir,'dir')
        mkdir(exportRootDir);
    else
        if cleanOutputFolder && ~is_dir_empty(exportRootDir)
            warning(['HydroPol2D:ExportRootNotEmpty\n' ...
                     'Export root already exists and contains files:\n  %s\n' ...
                     'Deleting all existing contents before exporting new results.'], exportRootDir);
            delete_dir_contents(exportRootDir);
        elseif ~cleanOutputFolder && ~is_dir_empty(exportRootDir)
            warning(['HydroPol2D:ExportRootNotEmpty\n' ...
                     'Export root already exists and contains files:\n  %s\n' ...
                     'Existing contents will be preserved because clean_output_folder = false.'], exportRootDir);
        end
    end

    % ---------------------------------------------------------------------
    % Main folders
    % ---------------------------------------------------------------------
    Paths = struct();

    Paths.Root    = exportRootDir;
    Paths.Results = fullfile(exportRootDir, 'Modeling_Results');
    Paths.Temp    = fullfile(exportRootDir, 'Temporary_Files');
    Paths.Logs    = fullfile(exportRootDir, 'Logs');

    % Subfolders inside Modeling_Results
    Paths.FigPDF        = fullfile(Paths.Results, 'Figures_PDF');
    Paths.FigFIG        = fullfile(Paths.Results, 'Figures_FIG');
    Paths.Tables        = fullfile(Paths.Results, 'Tables_CSV');

    Paths.RastersWD     = fullfile(Paths.Results, 'Rasters_Water_Depths');
    Paths.RastersWSE    = fullfile(Paths.Results, 'Rasters_WSE');
    Paths.RastersStatic = fullfile(Paths.Results, 'Rasters_Static');

    Paths.WQMaps        = fullfile(Paths.Results, 'Rasters_WQ');
    Paths.HRMaps        = fullfile(Paths.Results, 'Rasters_Human_Risk');

    Paths.Anim          = fullfile(Paths.Results, 'GIFs_MP4');
    Paths.Shapes        = fullfile(Paths.Results, 'Shapefiles');

    % ---------------------------------------------------------------------
    % Create all folders
    % ---------------------------------------------------------------------
    mkdir_if_missing(Paths.Root);
    mkdir_if_missing(Paths.Results);
    mkdir_if_missing(Paths.Temp);
    mkdir_if_missing(Paths.Logs);

    mkdir_if_missing(Paths.FigPDF);
    mkdir_if_missing(Paths.FigFIG);
    mkdir_if_missing(Paths.Tables);

    mkdir_if_missing(Paths.RastersWD);
    mkdir_if_missing(Paths.RastersWSE);
    mkdir_if_missing(Paths.RastersStatic);

    mkdir_if_missing(Paths.WQMaps);
    mkdir_if_missing(Paths.HRMaps);

    mkdir_if_missing(Paths.Anim);
    mkdir_if_missing(Paths.Shapes);
end

function tf = is_dir_empty(d)
%IS_DIR_EMPTY Return true if directory contains no files/folders besides
% "." and "..".

    L = dir(d);
    names = {L.name};
    tf = all(ismember(names, {'.','..'}));
end

function delete_dir_contents(d)
%DELETE_DIR_CONTENTS Delete all files and folders inside directory d,
% while keeping directory d itself.

    L = dir(d);
    for i = 1:numel(L)
        name = L(i).name;

        if strcmp(name,'.') || strcmp(name,'..')
            continue;
        end

        target = fullfile(d, name);

        if L(i).isdir
            rmdir(target, 's');
        else
            delete(target);
        end
    end
end

function mkdir_if_missing(d)
%MKDIR_IF_MISSING Create directory if it does not already exist.

    if ~exist(d,'dir')
        mkdir(d);
    end
end