
%% ========================================================================
%  HydroPol2D Preprocessing v3
%  - Excel mode OR bypass mode
% ========================================================================

% =========================================
% REAL-TIME RUN TRACKING
% =========================================
run_start_datetime = datetime('now');
run_start_str = char(datetime(run_start_datetime, ...
    'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));

if ~exist('use_inputpaths_bypass','var') || isempty(use_inputpaths_bypass)
    use_inputpaths_bypass = 0;
end

% --- Inputs ---
% model_folder: full path to your Excel file

% ============================
% OUTPUT FOLDERS (ROBUST)
% ============================
results_dir  = resultsDir;
temp_dir     = fullfile(Paths.Temp);
if ~exist(results_dir,'dir'); mkdir(results_dir); end
if ~exist(temp_dir,'dir');    mkdir(temp_dir);    end

% keep legacy names
folderName   = results_dir;
folderName_2 = temp_dir;

disp(['Outputs will be saved to: ' results_dir]);

% variable to specify the size of matrices to store out of memory
saver_memory_maps = 12;

% =========================
% CONFIG SOURCE
% =========================
if use_inputpaths_bypass == 1

    % ------------------------------------------------------------
    % BYPASS MODE
    % Assumes InputPaths already exists in workspace
    % and flags are created later by input_data_bypass_script
    % or already exist from wrapper.
    % ------------------------------------------------------------
    if ~exist('InputPaths','var') || ~isstruct(InputPaths)
        error(['use_inputpaths_bypass = 1 but InputPaths was not found. ' ...
               'Call input_paths_bypass(...) before preprocessing.']);
    end

    if ~exist('flags','var') || ~isstruct(flags)
    
        % ------------------------------------------------------------
        % In bypass mode, load the bypass input definition EARLY so that
        % flags are already available for optional raster loading,
        % resampling, and other preprocessing decisions.
        % ------------------------------------------------------------
        if exist('input_data_bypass_script.m','file') ~= 2
            error(['use_inputpaths_bypass = 1 but input_data_bypass_script.m ' ...
                   'was not found in the current MATLAB path or working folder.']);
        end
    
        % NOTE:
        % We only want to preload InputData_Bypass and flags here.
        % The full legacy-variable translation will still happen later
        % through the normal call to input_data_script.
        if ~exist('input_data_bypass_script_path','var') || isempty(input_data_bypass_script_path)
            error(['use_inputpaths_bypass = 1 but input_data_bypass_script_path was not ' ...
                   'provided by the main wrapper.']);
        end
        
        if exist(input_data_bypass_script_path,'file') ~= 2
            error(['use_inputpaths_bypass = 1 but the bypass input script was not found:\n  %s'], ...
                  input_data_bypass_script_path);
        end
        
        clear InputData_Bypass
        run(input_data_bypass_script_path);
        
        if ~exist('InputData_Bypass','var') || ~isstruct(InputData_Bypass)
            error('input_data_bypass_script.m must create a struct named InputData_Bypass.');
        end
        
        if ~isfield(InputData_Bypass,'flags') || ~isstruct(InputData_Bypass.flags)
            error('InputData_Bypass.flags is missing or invalid in input_data_bypass_script.m.');
        end
        
        flags = InputData_Bypass.flags;
    end

    % No General_Data table in bypass mode
    GD = [];

    % Tool paths
    topo_path        = get_inputpaths_field(InputPaths,'topo_path','');
    hydropol2d_tools = get_inputpaths_field(InputPaths,'hydropol2d_tools','');

    % Required rasters
    DEM_path  = get_inputpaths_field(InputPaths,'DEM_path','');
    LULC_path = get_inputpaths_field(InputPaths,'LULC_path','');
    SOIL_path = get_inputpaths_field(InputPaths,'SOIL_path','');

    % Optional rasters
    Warmup_Depth_path          = get_inputpaths_field(InputPaths,'Warmup_Depth_path','');
    Initial_Buildup_path       = get_inputpaths_field(InputPaths,'Initial_Buildup_path','');
    Initial_Soil_Moisture_path = get_inputpaths_field(InputPaths,'Initial_Soil_Moisture_path','');

    Albedo_path = get_inputpaths_field(InputPaths,'Albedo_path','');
    LAI_path    = get_inputpaths_field(InputPaths,'LAI_path','');
    DTB_path    = get_inputpaths_field(InputPaths,'DTB_path','');
    GW_table_path = get_inputpaths_field(InputPaths,'GW_table_path','');

    B1_path = get_inputpaths_field(InputPaths,'B1_path','');
    B2_path = get_inputpaths_field(InputPaths,'B2_path','');
    W1_path = get_inputpaths_field(InputPaths,'W1_path','');
    W2_path = get_inputpaths_field(InputPaths,'W2_path','');

    Subgrid_DEM_path = get_inputpaths_field(InputPaths,'Subgrid_DEM_path','');

    RiverWidths_path = get_inputpaths_field(InputPaths,'RiverWidths_path','');
    RiverDepths_path = get_inputpaths_field(InputPaths,'RiverDepths_path','');

else

    % ------------------------------------------------------------
    % EXCEL MODE (legacy)
    % ------------------------------------------------------------
    GD = readcell(model_folder,'Sheet','General_Data');

    FlagsCell = readcell(model_folder,'Sheet','Flags');
    flags = read_flags_table(FlagsCell);

    topo_path        = xlgetstr(GD,'topo_path',"");
    hydropol2d_tools = xlgetstr(GD,'hydropol2d_tools',"");

    DEM_path  = xlgetstr(GD,'DEM_path',"");
    LULC_path = xlgetstr(GD,'LULC_path',"");
    SOIL_path = xlgetstr(GD,'SOIL_path',"");

    Warmup_Depth_path          = xlgetstr(GD,'Warmup_Depth_path',"");
    Initial_Buildup_path       = xlgetstr(GD,'Initial_Buildup_path',"");
    Initial_Soil_Moisture_path = xlgetstr(GD,'Initial_Soil_Moisture_path',"");

    Albedo_path = xlgetstr(GD,'Albedo_path',"");
    LAI_path    = xlgetstr(GD,'LAI_path',"");
    DTB_path    = xlgetstr(GD,'DTB_path',"");
    
    % Optional groundwater table raster path.
    % If the Excel field does not exist, assume no raster was provided.
    GW_table_path = "";
    try
        GW_table_path = xlgetstr(GD,'GW_table_path',"");
    catch
        GW_table_path = "";
    end

    B1_path = xlgetstr(GD,'B1_path',"");
    B2_path = xlgetstr(GD,'B2_path',"");
    W1_path = xlgetstr(GD,'W1_path',"");
    W2_path = xlgetstr(GD,'W2_path',"");

    Subgrid_DEM_path = xlgetstr(GD,'Subgrid_DEM_path',"");

    RiverWidths_path = xlgetstr(GD,'RiverWidths_path',"");
    RiverDepths_path = xlgetstr(GD,'RiverDepths_path',"");

end

% =========================
% ADD TOOL PATHS
% =========================
if strlength(string(topo_path)) > 0 && isfolder(char(topo_path))
    addpath(genpath(char(topo_path)));
end
if strlength(string(hydropol2d_tools)) > 0 && isfolder(char(hydropol2d_tools))
    addpath(genpath(char(hydropol2d_tools)));
end

% =========================
% REQUIRED RASTER CHECK
% =========================
assert(isfile(DEM_path),  'DEM_path not found or empty.');
assert(isfile(LULC_path), 'LULC_path not found or empty.');
assert(isfile(SOIL_path), 'SOIL_path not found or empty.');

% =========================
% READ BASE RASTERS
% =========================
DEM_raster  = GRIDobj(char(DEM_path));
LULC_raster = GRIDobj(char(LULC_path));
SOIL_raster = GRIDobj(char(SOIL_path));

% ---- Ensure CRS is present (copy from DEM geotiff if missing) ----
% DEM_raster  = ensure_projected_crs_from_geotiff(DEM_raster,  char(DEM_path));
% LULC_raster = ensure_projected_crs_like(        LULC_raster, DEM_raster);
% SOIL_raster = ensure_projected_crs_like(        SOIL_raster, DEM_raster);

% ---- Crop NaN border consistently (and re-attach CRS) ----
crs_save = DEM_raster.georef.SpatialRef.ProjectedCRS;
DEM_raster  = crop(DEM_raster);
LULC_raster = crop(LULC_raster);
SOIL_raster = crop(SOIL_raster);
DEM_raster.georef.SpatialRef.ProjectedCRS  = crs_save;
LULC_raster.georef.SpatialRef.ProjectedCRS = crs_save;
SOIL_raster.georef.SpatialRef.ProjectedCRS = crs_save;

% =========================
% OPTIONAL RASTERS (read safely)
% =========================
DTB_raster      = [];
GW_table_raster = [];
LAI_raster      = [];
Albedo_raster   = [];
widths_raster   = [];
depths_raster   = [];

if use_inputpaths_bypass == 1

    % ------------------------------------------------------------
    % BYPASS MODE
    % At this stage flags may not yet be fully built by input_data_script.
    % Therefore, load optional rasters based on file existence only.
    % Later model logic will decide whether to actually use them.
    % ------------------------------------------------------------
    if isfile(DTB_path)
        DTB_raster = GRIDobj(char(DTB_path));
    end
    
    if isfile(GW_table_path)
        GW_table_raster = GRIDobj(char(GW_table_path));
    end
    
    if isfile(LAI_path)
        LAI_raster = GRIDobj(char(LAI_path));
    end

    if isfile(Albedo_path)
        Albedo_raster = GRIDobj(char(Albedo_path));
    end

    if isfile(RiverWidths_path) && isfile(RiverDepths_path)
        widths_raster = GRIDobj(char(RiverWidths_path));
        depths_raster = GRIDobj(char(RiverDepths_path));
    end

else

    % ------------------------------------------------------------
    % EXCEL MODE
    % Preserve original flag-driven behavior
    % ------------------------------------------------------------
    if (flags.flag_baseflow == 1 || flags.flag_groundwater_modeling == 1) && isfile(DTB_path)
        DTB_raster = GRIDobj(char(DTB_path));
    end
    
    % Optional initial groundwater table raster.
    % If GW_table_path is missing, empty, or not a file, GW_table_raster remains [].
    if (flags.flag_baseflow == 1 || flags.flag_groundwater_modeling == 1) && ...
            strlength(strtrim(string(GW_table_path))) > 0 && ...
            isfile(char(GW_table_path))
    
        GW_table_raster = GRIDobj(char(GW_table_path));
    
    end
    
    if (flags.flag_abstraction == 1) && isfile(LAI_path)
        LAI_raster = GRIDobj(char(LAI_path));
    end

    if (flags.flag_spatial_albedo == 1) && isfile(Albedo_path)
        Albedo_raster = GRIDobj(char(Albedo_path));
    end

    if (flags.flag_subgrid == 1) && (flags.flag_river_rasters == 1) ...
            && isfile(RiverWidths_path) && isfile(RiverDepths_path)
        widths_raster = GRIDobj(char(RiverWidths_path));
        depths_raster = GRIDobj(char(RiverDepths_path));
    end

end


%% Load General Data Input
if use_inputpaths_bypass == 1
    if ~exist('GIS_data','var') || ~isstruct(GIS_data)
        GIS_data = struct();
    end

    if exist('InputData_Bypass','var') && isstruct(InputData_Bypass) && ...
       isfield(InputData_Bypass,'general') && isstruct(InputData_Bypass.general) && ...
       isfield(InputData_Bypass.general,'resolution_resample')
        GIS_data.resolution_resample = InputData_Bypass.general.resolution_resample;
    else
        GIS_data.resolution_resample = NaN;
    end
else
    GIS_data.resolution_resample = xlnum(GD,'resolution_resample');
end

% =========================
% RESAMPLE OPTION (apply ONCE)
% =========================
if flags.flag_resample == 1
    resolution = GIS_data.resolution_resample;

    DEM0 = DEM_raster;   % preserve original CRS/reference metadata

    DEM_raster = resample(DEM_raster, resolution, 'bilinear');

    % Reattach CRS if resample dropped it
    DEM_raster.georef.SpatialRef.ProjectedCRS = DEM0.georef.SpatialRef.ProjectedCRS;

    LULC_raster = resample(LULC_raster, DEM_raster, 'nearest');
    LULC_raster.Z = round(LULC_raster.Z);

    SOIL_raster = resample(SOIL_raster, DEM_raster, 'nearest');
    SOIL_raster.Z = round(SOIL_raster.Z);

    if ~isempty(DTB_raster),       DTB_raster      = resample(DTB_raster,      DEM_raster, 'bilinear'); end
    if ~isempty(GW_table_raster),  GW_table_raster = resample(GW_table_raster, DEM_raster, 'bilinear'); end
    if ~isempty(LAI_raster),       LAI_raster      = resample(LAI_raster,      DEM_raster, 'bilinear'); end
    if ~isempty(Albedo_raster),    Albedo_raster   = resample(Albedo_raster,   DEM_raster, 'bilinear'); end
    if ~isempty(widths_raster), widths_raster = resample(widths_raster, DEM_raster, 'nearest');  end
    if ~isempty(depths_raster), depths_raster = resample(depths_raster, DEM_raster, 'nearest');  end
end

% =========================
% HARD ALIGNMENT CHECK (one place only)
% =========================
[LULC_raster, SOIL_raster, DTB_raster, LAI_raster, Albedo_raster, widths_raster, depths_raster] = ...
    align_all_to_dem(DEM_raster, LULC_raster, SOIL_raster, DTB_raster, LAI_raster, Albedo_raster, widths_raster, depths_raster);

if ~isempty(GW_table_raster)
    GW_table_raster = resample(GW_table_raster, DEM_raster, 'bilinear');
end

% =========================
% Convert to matrices, domain mask
% =========================
min_dem_value = -200; % your choice

DEM  = double(DEM_raster.Z);
LULC = double(LULC_raster.Z);
SOIL = double(SOIL_raster.Z);

% Validity checks
neg_DEM  = DEM  < min_dem_value;
neg_LULC = LULC < 0;
neg_SOIL = SOIL < 0;

invalid_DEM  = isinf(DEM)  | isnan(DEM)  | neg_DEM;
invalid_LULC = isinf(LULC) | isnan(LULC) | neg_LULC;
invalid_SOIL = isinf(SOIL) | isnan(SOIL) | neg_SOIL;

% Base invalid mask: ONLY DEM and LULC define exclusion from domain
idx_nan = invalid_DEM | invalid_LULC;

% Cells inside domain but with missing soil --> assign soil class 0
idx_soil_nodata_inside = ~idx_nan & invalid_SOIL;
SOIL(idx_soil_nodata_inside) = 0;

% Apply NaN only to truly invalid cells
DEM(idx_nan)  = nan;
LULC(idx_nan) = nan;
SOIL(idx_nan) = nan;

DEM_raster.Z  = DEM;
LULC_raster.Z = LULC;
SOIL_raster.Z = SOIL;

if ~isempty(DTB_raster),       DTB_raster.Z(idx_nan) = nan; end
if ~isempty(GW_table_raster),  GW_table_raster.Z(idx_nan) = nan; end
if ~isempty(LAI_raster),       LAI_raster.Z(idx_nan) = nan; end
if ~isempty(Albedo_raster),    Albedo_raster.Z(idx_nan) = nan; end
if ~isempty(widths_raster),    widths_raster.Z(idx_nan) = nan; end
if ~isempty(depths_raster),    depths_raster.Z(idx_nan) = nan; end

GIS_data.xulcorner = DEM_raster.refmat(3,1);
GIS_data.yulcorner = DEM_raster.refmat(3,2);
Wshed_Properties.Resolution = DEM_raster.cellsize;
Wshed_Properties.cell_area  = Wshed_Properties.Resolution^2;

[ny,nx] = size(DEM);

input_data_script;  % Load general data, soil, and LULC parameters.

% ------------------------------------------------------------
% Groundwater flag consistency
% ------------------------------------------------------------
if ~isfield(flags, 'flag_groundwater_modeling')
    flags.flag_groundwater_modeling = 0;
end
if ~isfield(flags, 'flag_groundwater_async') || ~isfinite(double(flags.flag_groundwater_async))
    flags.flag_groundwater_async = 1;
end
if ~isfield(flags, 'flag_capillary_rise') || ~isfinite(double(flags.flag_capillary_rise))
    flags.flag_capillary_rise = 1;
end
if ~isfield(flags, 'groundwater_target_dt_min') || ~isfinite(double(flags.groundwater_target_dt_min))
    flags.groundwater_target_dt_min = 1440;
end
if ~isfield(flags, 'groundwater_min_dt_min') || ~isfinite(double(flags.groundwater_min_dt_min))
    flags.groundwater_min_dt_min = 1;
end
if ~isfield(flags, 'groundwater_max_head_change_m') || ~isfinite(double(flags.groundwater_max_head_change_m))
    flags.groundwater_max_head_change_m = 0.25;
end
if ~isfield(flags, 'groundwater_courant') || ~isfinite(double(flags.groundwater_courant))
    flags.groundwater_courant = 0.25;
end

if flags.flag_baseflow == 1 && flags.flag_groundwater_modeling == 0
    error(['Invalid flag combination: flag_baseflow = 1 requires ', ...
           'flag_groundwater_modeling = 1.']);
end


% ========================================================================
% From here, keep your existing pipeline:
% - input_data_script
% - gauges (but use the block-reader)
% - rainfall/ETP handling
% - subgrid creation
% - outlets, streams, etc.
% ========================================================================

% --- Optional raster paths (only required if corresponding flag is ON) ---
if use_inputpaths_bypass ~= 1
    Warmup_Depth_path          = xlpath(GD,'Warmup_Depth_path',          flags.flag_warmup == 1);
    Initial_Buildup_path       = xlpath(GD,'Initial_Buildup_path',       flags.flag_initial_buildup == 1);
    Initial_Soil_Moisture_path = xlpath(GD,'Initial_Soil_Moisture_path', flags.flag_warmup == 1);
    
    Albedo_path = xlpath(GD,'Albedo_path', flags.flag_spatial_albedo == 1);
    
    LAI_path    = xlpath(GD,'LAI_path',    flags.flag_abstraction == 1);     % or whatever your LAI flag name is
    % NDVI_path = xlpath(GD,'NDVI_path',   flags.flag_NDVI == 1);    % if you have it
    
    RiverWidths_path = xlpath(GD,'RiverWidths_path', flags.flag_subgrid == 1 && flags.flag_river_rasters == 1);
    RiverDepths_path = xlpath(GD,'RiverDepths_path', flags.flag_subgrid == 1 && flags.flag_river_rasters == 1);
    
    DTB_path = xlpath(GD,'DTB_path', flags.flag_baseflow == 1 || flags.flag_groundwater_modeling == 1);
    
    B1_path = xlpath(GD,'B1_path', flags.flag_waterquality == 1 && flags.flag_WQ_Rasters == 1);
    B2_path = xlpath(GD,'B2_path', flags.flag_waterquality == 1 && flags.flag_WQ_Rasters == 1);
    W1_path = xlpath(GD,'W1_path', flags.flag_waterquality == 1 && flags.flag_WQ_Rasters == 1);
    W2_path = xlpath(GD,'W2_path', flags.flag_waterquality == 1 && flags.flag_WQ_Rasters == 1);
    
    Subgrid_DEM_path = xlpath(GD,'Subgrid_DEM_path', flags.flag_subgrid == 1);
end

% Rasters
fname_LULC = LULC_path; fname_DEM = DEM_path;
fname_SOIL = SOIL_path;
min_dem_value = -200; % min value that a dem can have


% %% Checking Extent Problem
% crs_save = DEM_raster.georef.SpatialRef.ProjectedCRS;
% % croping all nan rows and columns
% DEM_raster = crop(DEM_raster);
% LULC_raster = crop(LULC_raster);
% SOIL_raster = crop(SOIL_raster);
% 
% baseflow_check = ~isempty(DTB_path);
% if baseflow_check == 1
%     DTB_raster = crop(DTB_raster);
% end
% %
% DEM_raster.georef.SpatialRef.ProjectedCRS = crs_save;
% LULC_raster.georef.SpatialRef.ProjectedCRS = crs_save;
% SOIL_raster.georef.SpatialRef.ProjectedCRS = crs_save;
% 
% % if min(min(DEM_raster.Z)) <= 0
% %     error('Please make sure that you actually have negative or 0 values in the DEM. Otherwise, treat non-value points as NaN or -9999.')
% % end
% if sum(size(DEM_raster.Z)) == sum(size(LULC_raster.Z)) && sum(size(DEM_raster.Z)) == sum(size(SOIL_raster.Z))
% else
%     if sum(size(DEM_raster.Z)) > sum(size(LULC_raster.Z)) || sum(size(DEM_raster.Z)) > sum(size(SOIL_raster.Z)) % DEM is larger
%         raster_resample = DEM_raster;
%         % Resample other two rasters
%         % ---- Constraint at LULC Raster
%         LULC_raster = resample(LULC_raster,raster_resample,'nearest');
%         % LULC_raster.Z = round(LULC_raster.Z,'nearest'); % Only Integers
%         % ---- Constraint at SOIL Raster
%         SOIL_raster = resample(SOIL_raster,raster_resample,'nearest');
%         % SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
%     end
% 
%     if sum(size(SOIL_raster.Z)) > sum(size(DEM_raster.Z)) || sum(size(SOIL_raster.Z)) > sum(size(LULC_raster.Z))  % SOIL is larger
%         raster_resample = SOIL_raster;
%         % Resample other two rasters
%         LULC_raster = resample(LULC_raster,raster_resample,'nearest');
%         % ---- Constraint at LULC Raster
%         LULC_raster = resample(LULC_raster,raster_resample,'nearest');
%         % LULC_raster.Z = round(LULC_raster.Z); % Only Integers
%         DEM_raster = resample(DEM_raster,raster_resample,'bilinear');
%     end
% 
%     if sum(size(LULC_raster.Z)) > sum(size(DEM_raster.Z)) || sum(size(LULC_raster.Z)) > sum(size(SOIL_raster.Z))  % SOIL is larger
%         raster_resample = LULC_raster;
%         % Resample other two rasters
%         % ---- Constraint at SOIL Raster
%         SOIL_raster = resample(SOIL_raster,raster_resample,'nearest');
%         % SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
%         DEM_raster = resample(DEM_raster,raster_resample,'bilinear');
%     end
% 
%     if baseflow_check == 1 % Case of baseflow simulation
%         DTB_raster = resample(DTB_raster,raster_resample,'bilinear');
%     end
% 
%     try
%         Albedo_raster = resample(Albedo_raster,raster_resample,'bilinear');
%     end
%     try
%         LAI_raster = resample(LAI_raster,raster_resample,'bilinear');
%     end
% 
% end
% % Checking if there are nan cells within the study area to avoid numerical
% % instability
% if any(~isnan(DEM_raster.Z).*isnan(LULC_raster.Z)) == 1
%     imagesc(~isnan(DEM_raster.Z).*isnan(LULC_raster.Z));
%     error('Please, check your DEM and LULC rasters. There are cells with no information which will produce numerical instability')
% elseif any(~isnan(DEM_raster.Z).*isnan(SOIL_raster.Z)) == 1
%     imagesc(~isnan(DEM_raster.Z).*isnan(SOIL_raster.Z));
%     error('Please, check your DEM and SOIL rasters. There are cells with no information which will produce numerical instability')
% end
% 
% % Raster Extent
% GIS_data.xulcorner = DEM_raster.refmat(3,1); % Up Left Corner
% GIS_data.yulcorner = DEM_raster.refmat(3,2);
% % - Extent is already solved, we can login the input data
% Wshed_Properties.Resolution = DEM_raster.cellsize; % m


% %% Resampling Maps
% % In case we want to resample the DEM 
% if flags.flag_resample == 1
%     resolution = GIS_data.resolution_resample; % m
%     % DEM
%     DEM_raster = resample(DEM_raster,resolution,'bilinear');
%     % LULC
%     LULC_raster = resample(LULC_raster,resolution);
%     LULC_raster.Z = round(LULC_raster.Z);
%     % SOIL
%     SOIL_raster = resample(SOIL_raster,resolution);
%     SOIL_raster.Z = round(SOIL_raster.Z);
% 
%     % Extent Problem
%     if sum(size(DEM_raster.Z)) > sum(size(LULC_raster.Z)) && sum(size(DEM_raster.Z)) >= sum(size(SOIL_raster.Z)) % DEM is larger
%         raster_resample = DEM_raster;
%         % Resample other two rasters
%         % ---- Constraint at LULC Raster
%         LULC_raster = resample(LULC_raster,raster_resample);
%         LULC_raster.Z = round(LULC_raster.Z); % Only Integers
%         % ---- Constraint at SOIL Raster
%         SOIL_raster = resample(SOIL_raster,raster_resample);
%         SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
%     end
% 
%     if sum(size(SOIL_raster.Z)) > sum(size(DEM_raster.Z)) && sum(size(SOIL_raster.Z)) >= sum(size(LULC_raster.Z))  % SOIL is larger
%         raster_resample = SOIL_raster;
%         % Resample other two rasters
%         LULC_raster = resample(LULC_raster,raster_resample);
%         % ---- Constraint at LULC Raster
%         LULC_raster = resample(LULC_raster,raster_resample);
%         LULC_raster.Z = round(LULC_raster.Z); % Only Integers
%         DEM_raster = resample(DEM_raster,raster_resample);
%     end
% 
%     if sum(size(LULC_raster.Z)) > sum(size(DEM_raster.Z)) && sum(size(LULC_raster.Z)) >= sum(size(SOIL_raster.Z))  % SOIL is larger
%         raster_resample = DEM_raster;
%         % Resample other two rasters
%         % ---- Constraint at SOIL Raster
%         SOIL_raster = resample(SOIL_raster,raster_resample);
%         SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
%         DEM_raster = resample(DEM_raster,raster_resample);
%     end
% 
%     % Checking CRS
%     try
%         if isempty(DEM_raster.georef.SpatialRef.ProjectedCRS)
%             % Enter with CRS code for the raster
%             %     prompt = "What is the CRS code for your raster data? ";
%             %     code = input(prompt);
%             %     projected_data = projcrs(code);
%             %     gtiffinfo.SpatialRef.ProjectedCRS = projected_data;
%             %     DEM_raster.georef.SpatialRef  = gtiffinfo.SpatialRef;
%             %     LULC_raster.georef.SpatialRef  = gtiffinfo.SpatialRef;
%             %     SOIL_raster.georef.SpatialRef  = gtiffinfo.SpatialRef;
%             [~,R] = readgeoraster(fname_DEM);
%             proj = R.ProjectedCRS;
%             DEM_raster.georef.SpatialRef.ProjectedCRS  = proj;
%             LULC_raster.georef.SpatialRef.ProjectedCRS  = proj;
%             SOIL_raster.georef.SpatialRef.ProjectedCRS  = proj;
%         end
%     catch me
%         if isempty(DEM_raster.georef)
%             % some code?
%         end
%     end
% 
%     % Raster Extent
%     GIS_data.xulcorner = DEM_raster.refmat(3,1); % Up Left Corner
%     GIS_data.yulcorner = DEM_raster.refmat(3,2);
%     % - Extent is already solved, we can login the input data
%     Wshed_Properties.Resolution = DEM_raster.cellsize; % m
% end
% 
% if baseflow_check == 1 && flags.flag_resample == 1% Resample depths if required
%     DTB_raster = resample(DTB_raster,resolution,'bilinear');
% end
% try
%     if LAI_raster.cellsize > 0
%         LAI_raster = resample(LAI_raster,resolution,'bilinear');
%     end
% end


%% Observed Gauges (name-based from General_Data table block)
if flags.flag_obs_gauges == 1

    if use_inputdata_bypass == 1
        if ~exist('Obs','var') || isempty(Obs)
            error(['flag_obs_gauges = 1 and bypass mode is active, but Obs was not ' ...
                   'created in input_data_bypass_script.']);
        end
    else
        Obs = read_block_table(GD, 'Gauge', {'Gauge','Easting (m)','Northing (m)','Label Name'});
    end

    obs_gauges = Obs.("Gauge");
    gauges.num_obs_gauges = sum(~isnan(obs_gauges));

    Obs = Obs(1:gauges.num_obs_gauges,:);
    
    % -------------------------------------------------------------------------
    % Read observed gauge coordinates with flexible column-name handling
    % -------------------------------------------------------------------------
    obsVars = Obs.Properties.VariableNames;
    
    % Possible names after readtable sanitization
    possibleEastingNames = { ...
        'Easting (m)', ...
        'Easting_m', ...
        'easting_m', ...
        'easting__m_', ...
        'easting_m_', ...
        'Easting_m_', ...
        'easting' ...
    };
    
    possibleNorthingNames = { ...
        'Northing (m)', ...
        'Northing_m', ...
        'northing_m', ...
        'northing__m_', ...
        'northing_m_', ...
        'Northing_m_', ...
        'northing' ...
    };
    
    % Find matching easting column
    idxE = find(ismember(obsVars, possibleEastingNames), 1);
    
    % Find matching northing column
    idxN = find(ismember(obsVars, possibleNorthingNames), 1);
    
    % Stop clearly if columns are not found
    if isempty(idxE)
        error('Could not find easting column in Obs. Available columns are: %s', ...
            strjoin(obsVars, ', '));
    end
    
    if isempty(idxN)
        error('Could not find northing column in Obs. Available columns are: %s', ...
            strjoin(obsVars, ', '));
    end
    
    % Assign coordinates
    gauges.easting_obs_gauges_absolute  = Obs.(obsVars{idxE});
    gauges.northing_obs_gauges_absolute = Obs.(obsVars{idxN});

    gauges.x_coord_gauges = gauges.easting_obs_gauges_absolute;
    gauges.y_coord_gauges = gauges.northing_obs_gauges_absolute;

    gauges.easting_obs_gauges  = round((-GIS_data.xulcorner + gauges.x_coord_gauges)/Wshed_Properties.Resolution);
    gauges.northing_obs_gauges = round((GIS_data.yulcorner - gauges.y_coord_gauges)/Wshed_Properties.Resolution);

    gauges.x_index_gauges = gauges.easting_obs_gauges;
    gauges.y_index_gauges = gauges.northing_obs_gauges;

    gauges.labels_observed_string = cellstr(string(Obs.("Label Name")));

    % Labels
    % -------------------------------------------------------------------------
    % Read observed gauge labels with flexible column-name handling
    % -------------------------------------------------------------------------
    obsVars = Obs.Properties.VariableNames;
    
    possibleLabelNames = { ...
        'Label Name', ...
        'Label_Name', ...
        'LabelName', ...
        'Name', ...
        'name', ...
        'Gauge_Name', ...
        'GaugeName', ...
        'label_name', ...
        'labelname' ...
    };
    
    idxLabel = find(ismember(obsVars, possibleLabelNames), 1);
    
    if isempty(idxLabel)
        error('Could not find label/name column in Obs. Available columns are: %s', ...
            strjoin(obsVars, ', '));
    end
    
    gauges.labels_observed_string = cellstr(string(Obs.(obsVars{idxLabel})));

    % OPTIONAL morphometric params (only if you add these columns to the table)
    % if ismember("alfa_1", Obs.Properties.VariableNames)
    %     GIS_data.alfa_1 = Obs.alfa_1;
    %     GIS_data.alfa_2 = Obs.alfa_2;
    %     GIS_data.beta_1 = Obs.beta_1;
    %     GIS_data.beta_2 = Obs.beta_2;
    %     River_Manning    = Obs.River_Manning;
    %     Lateral_Groundwater_Flux = Obs.Lateral_Groundwater_Flux;
    % end
end
%% Snow Parameters
if flags.flag_snow_modeling == 1

    % ------------------------------------------------------------
    % Initialize snow state rasters if they do not already exist
    % ------------------------------------------------------------
    if ~exist('Snow_Properties','var') || ~isstruct(Snow_Properties)
        Snow_Properties = struct();
    end

    if ~isfield(Snow_Properties,'H_snow_t') || isempty(Snow_Properties.H_snow_t)
        Snow_Properties.H_snow_t = zeros(size(DEM_raster.Z)); % [m]
        Snow_Properties.H_snow_t(isnan(DEM_raster.Z)) = nan;
    end

    if ~isfield(Snow_Properties,'SWE_t') || isempty(Snow_Properties.SWE_t)
        Snow_Properties.SWE_t = Snow_Properties.H_snow_t;
    end

    % ------------------------------------------------------------
    % Assign defaults ONLY if fields were not provided already
    % ------------------------------------------------------------
    if ~isfield(Snow_Properties,'alpha'),           Snow_Properties.alpha = 0.8; end
    if ~isfield(Snow_Properties,'epsilon'),         Snow_Properties.epsilon = 0.98; end
    if ~isfield(Snow_Properties,'C_e'),             Snow_Properties.C_e = 0.001; end
    if ~isfield(Snow_Properties,'DDF'),             Snow_Properties.DDF = 2; end
    if ~isfield(Snow_Properties,'T_thresh'),        Snow_Properties.T_thresh = 0; end
    if ~isfield(Snow_Properties,'rho_snow_init'),   Snow_Properties.rho_snow_init = 100; end
    if ~isfield(Snow_Properties,'rho_max'),         Snow_Properties.rho_max = 400; end
    if ~isfield(Snow_Properties,'k_t'),             Snow_Properties.k_t = 0.1; end
    if ~isfield(Snow_Properties,'k_swe'),           Snow_Properties.k_swe = 0.001; end
    if ~isfield(Snow_Properties,'k_D'),             Snow_Properties.k_D = 0.02; end
    if ~isfield(Snow_Properties,'snow_fraction_a'), Snow_Properties.snow_fraction_a = 0.2; end

end

%% ----- Transforming Raster into Matrix with Values ----- %

% LULC = double(LULC_raster.Z);
% DEM = double(DEM_raster.Z);
% SOIL = double(SOIL_raster.Z);

% % temporal modification, only to consider negative values in DEM due to
% % coastal negative values
% if flags.flag_input_rainfall_map
%     inf_nan_MAPS = isinf(DEM) + isnan(DEM); % Logical array
% else
%     neg_DEM = DEM < min_dem_value;
%     neg_LULC = LULC < 0;
%     neg_SOIL = SOIL < 0;
%     inf_nan_MAPS = isinf(DEM) + isnan(DEM) + neg_DEM + isnan(LULC) + isnan(SOIL) + neg_LULC + neg_SOIL + isinf(LULC) + isinf(SOIL); % Logical array
% end
% idx = inf_nan_MAPS > 0;
% 
% % Rebuilding Rasters to the Lowest Extent
% LULC_raster.Z = LULC;% Land Use and Land Cover Classification
% DEM_raster.Z = DEM; % Digital Elevation Model
% SOIL_raster.Z = SOIL; % Soil Map

% Replacing Values with Issues
% DEM_raster.Z(idx) = nan;
% LULC_raster.Z(idx) = nan;
% SOIL_raster.Z(idx) = nan;

LULC = double(LULC_raster.Z);
DEM = double(DEM_raster.Z);
SOIL = double(SOIL_raster.Z);

dem = DEM; % Further used in the elevation data
imp = LULC;
soil = SOIL;

Wshed_Properties.cell_area = Wshed_Properties.Resolution^2; % cell area in square meters
zzz = size(dem); % Dimensions of DEM matrix

% ---- DEM Dimensions --- %
[ny,nx] = size(dem);

if flags.flag_obs_gauges == 1
    gauges.easting_obs_gauges = round((-GIS_data.xulcorner + gauges.easting_obs_gauges_absolute)/Wshed_Properties.Resolution);
    gauges.northing_obs_gauges = round((GIS_data.yulcorner - gauges.northing_obs_gauges_absolute)/Wshed_Properties.Resolution);
end

% --- Converting coordinates to local coordinates in pixels
if flags.flag_reservoir == 1
    Reservoir_Data.x_index = round((-GIS_data.xulcorner + Reservoir_Data.x_us)/Wshed_Properties.Resolution);
    Reservoir_Data.y_index = round((GIS_data.yulcorner - Reservoir_Data.y_us)/Wshed_Properties.Resolution);
    Reservoir_Data.x_ds1_index = round((-GIS_data.xulcorner + Reservoir_Data.x_ds1)/Wshed_Properties.Resolution);
    Reservoir_Data.x_ds2_index = round((-GIS_data.xulcorner + Reservoir_Data.x_ds2)/Wshed_Properties.Resolution);
    Reservoir_Data.y_ds1_index = round((GIS_data.yulcorner - Reservoir_Data.y_ds1)/Wshed_Properties.Resolution);
    Reservoir_Data.y_ds2_index = round((GIS_data.yulcorner - Reservoir_Data.y_ds2)/Wshed_Properties.Resolution);
else
    Reservoir_Data.index = [];
    Reservoir_Data.x_index = [];
    Reservoir_Data.y_index = [];
    Reservoir_Data.k1 = [];
    Reservoir_Data.h1 = [];
    Reservoir_Data.k2 = [];
    Reservoir_Data.x_ds1_index = [];
    Reservoir_Data.y_ds1_index = [];
    Reservoir_Data.k3 = [];
    Reservoir_Data.h2 = [];
    Reservoir_Data.k4 = [];
    Reservoir_Data.x_ds2_index = [];
    Reservoir_Data.y_ds2_index = [];
end

%% ------------ Inflow Cells  ------------ %
% Here we read where the stream gauges are located
if flags.flag_inflow == 1
    % Wshed_Properties.inflow_mask = zeros(size(DEM_raster.Z));
    Wshed_Properties.inflow_cells = zeros(ny,nx,Inflow_Parameters.n_stream_gauges);
    Wshed_Properties.inflow_mask = zeros(ny,nx);
    for i = 1:Inflow_Parameters.n_stream_gauges
        for z = 1:Wshed_Properties.n_inlets(i)
            x = Inflow_Parameters.easting_inlet_cells(z,i);
            y = Inflow_Parameters.northing_inlet_cells(z,i);
            Wshed_Properties.inflow_cells(y,x,i) = 1; % Coordinates of each stream gauge
            Wshed_Properties.inflow_mask(y,x) = 1;
        end
    end
    Wshed_Properties.inflow_mask = logical(Wshed_Properties.inflow_mask);
end


%% ------------ Stage Hydrograph Cells  ------------ %
% Here we read where the stage gauges are located
if flags.flag_stage_hydrograph == 1
    % Wshed_Properties.stage_mask = zeros(size(DEM_raster.Z));
    Wshed_Properties.stage_mask = zeros(ny,nx,Stage_Parameters.n_stage_gauges);
    for i = 1:Stage_Parameters.n_stage_gauges
        for z = 1:Wshed_Properties.n_inlets_stage(i)
            x = Stage_Parameters.easting_inlet_cells(z,i);
            y = Stage_Parameters.northing_inlet_cells(z,i);
            Wshed_Properties.stage_mask(y,x,i) = 1; % Coordinates of each stream gauge
            Wshed_Properties.stage_mask(y,x) = 1;
        end
    end
    Wshed_Properties.stage_mask = logical(Wshed_Properties.stage_mask);
end

%% ------------ Rainfall Matrices ------------ %
if flags.flag_rainfall == 0 % No rainfall
    Wshed_Properties.rainfall_matrix = flags.flag_rainfall*zeros(size(dem));
elseif flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_input_rainfall_map ~= 1
    % Spatial Rainfall Case
    input_table = readtable('Rainfall_Spatial_Input.xlsx');
    % Observations
    n_obs = sum((table2array(input_table(:,2))>=0)); % Number of observations
    n_max_raingauges =  sum(~isnan((table2array(input_table(3,3:end)))));
    Spatial_Rainfall_Parameters.time_step_spatial = table2array(input_table(7,2)) - table2array(input_table(6,2)); % min
    end_rain = (n_obs-1)*Spatial_Rainfall_Parameters.time_step_spatial;
    %     rainfall_spatial_duration = 0:time_step_spatial:(end_rain); % Rainfall data time in minutes
    Spatial_Rainfall_Parameters.rainfall_spatial_duration = 0:Spatial_Rainfall_Parameters.time_step_spatial:(end_rain); % Rainfall data time in minutes
    Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg = 0:running_control.record_time_spatial_rainfall:(end_rain); % Rainfall data time in minutes
    n_spatial_agg = running_control.record_time_spatial_rainfall/Spatial_Rainfall_Parameters.time_step_spatial;
    rainfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),n_spatial_agg);

    % Rainfall Data
    for i = 1:n_max_raingauges
        Spatial_Rainfall_Parameters.rainfall_raingauges(:,i) = table2array(input_table(6:end,3 + (i-1)));
        Spatial_Rainfall_Parameters.coordinates(i,1) = table2array(input_table(3,3 + (i-1)));
        Spatial_Rainfall_Parameters.coordinates(i,2) = table2array(input_table(4,3 + (i-1)));
    end
    Spatial_Rainfall_Parameters.n_raingauges = sum(Spatial_Rainfall_Parameters.rainfall_raingauges(1,:) >= 0); % Number of raingauges
elseif flags.flag_input_rainfall_map == 1
    % We are running rainfall with input GeoTIFF maps (time-series)
    n_obs = Input_Rainfall.num_obs_maps;

    if n_obs < 1
        error('flag_input_rainfall_map=1 but no rainfall rasters were provided in General_Data block.');
    end
    if n_obs < 2 || isnan(Input_Rainfall.time(1))
        error('Need at least 2 rainfall rasters to compute dt. Check Rainfall_Rasters block.');
    end

    Spatial_Rainfall_Parameters.time_step_spatial = Input_Rainfall.time(2) - Input_Rainfall.time(1); % minutes
    end_rain = (n_obs-1) * Spatial_Rainfall_Parameters.time_step_spatial;

    % These two were swapped/inconsistent in your code. Keep them consistent:
    Spatial_Rainfall_Parameters.rainfall_spatial_duration     = 0:Spatial_Rainfall_Parameters.time_step_spatial:end_rain;
    Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg = 0:running_control.record_time_spatial_rainfall:end_rain;

    n_spatial_agg = ceil(running_control.record_time_spatial_rainfall / Spatial_Rainfall_Parameters.time_step_spatial);
    rainfall_spatial_aggregation = zeros(size(dem,1), size(dem,2), n_spatial_agg);

elseif flags.flag_rainfall == 1 && flags.flag_satellite_rainfall == 1
    n_spatial_agg = ceil(running_control.record_time_maps/running_control.record_time_spatial_rainfall);
    if mod(running_control.record_time_maps,running_control.record_time_spatial_rainfall) ~= 0
        error('Please, enter a multiple of the recording time for the input rainfall time-step')
    end
    rainfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),n_spatial_agg);
    % rianfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),saver_memory_maps);

elseif flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~= 1
    % Lumped Rainfall Case
    Wshed_Properties.rainfall_matrix = flags.flag_rainfall*ones(size(dem));
end

%% GW Heads
if flags.flag_baseflow == 1
    Maps.Hydro.GWdepth_save = zeros(size(dem,1),size(dem,2),saver_memory_maps);
end

%% Snowpack
if flags.flag_snow_modeling == 1
    Maps.Hydro.Snowpack = zeros(size(dem,1),size(dem,2),saver_memory_maps);
end

%% Abstraction
if flags.flag_abstraction == 1
    Maps.Hydro.Abstraction = zeros(size(dem,1),size(dem,2),saver_memory_maps);
end

%% ========================================================================
% ETP PREPROCESSING
% Clear separation between:
%   (A) Internal ETP from meteorological forcing
%   (B) Input evaporation/transpiration raster maps
% ========================================================================

% Safe initialization

if flags.flag_ETP ~= 1

    % ------------------------------------------------------------
    % ETP OFF
    % ------------------------------------------------------------
    % Keep empty structs so later code does not crash.
    ETP_Parameters = struct();
    Spatial_ETP_Parameters = struct();

elseif flags.flag_input_ETP_map == 0

    % ------------------------------------------------------------
    % MODE A: INTERNAL ETP FROM METEOROLOGICAL DATA
    % ------------------------------------------------------------
    % This mode uses:
    %   - Tmax, Tmin, Tavg
    %   - u2, RH, G
    %   - DEM
    %   - latitude / longitude
    %   - Krs
    %   - albedo or albedo raster
    %
    % Krs = empirical coefficient used for radiation estimation
    %       (e.g., Hargreaves-Samani type shortwave radiation estimate)
    % ------------------------------------------------------------

    % -------------------------
    % DEM and reference
    % -------------------------
    ETP_Parameters.DEM_etp = single(DEM_raster.Z);   % final aligned DEM

    % Keep spatial reference directly from the in-memory DEM raster.
    % This avoids writing DEM_ETP.tif only to read metadata back.
    ETP_Parameters.SpatialRef = DEM_raster.georef.SpatialRef;

    % Compatibility metadata so downstream legacy code can still use
    % ETP_Parameters.info without requiring a GeoTIFF write/read round-trip.
    ETP_Parameters.info = struct();
    ETP_Parameters.info.Height = size(ETP_Parameters.DEM_etp, 1);
    ETP_Parameters.info.Width  = size(ETP_Parameters.DEM_etp, 2);
    ETP_Parameters.info.SpatialRef = ETP_Parameters.SpatialRef;

    % -------------------------
    % Domain mask
    % -------------------------
    ETP_Parameters.idx_cells = ~idx_nan;
    ETP_Parameters.DEM_etp(idx_nan) = nan;

    % -------------------------
    % ETP empirical parameters
    % -------------------------
    ETP_Parameters.Krs = Krs_ETP * ones(ny, nx);
    ETP_Parameters.Krs(idx_nan) = nan;

    if flags.flag_spatial_albedo == 1 && ~isempty(Albedo_raster)
        ETP_Parameters.alfa_albedo_input = double(Albedo_raster.Z);
    else
        ETP_Parameters.alfa_albedo_input = albedo * ones(ny, nx);
    end
    ETP_Parameters.alfa_albedo_input(idx_nan) = nan;

    % -------------------------
    % Read meteorological forcing table
    % -------------------------
    % Excel mode:
    %   keep legacy default file
    %
    % Bypass mode:
    %   use InputData_Bypass.general.etp_input_spreadsheet
    %   (forwarded by input_data_script into ETP_input_spreadsheet)
    % -------------------------
    if use_inputdata_bypass == 1

        if ~exist('ETP_input_spreadsheet','var') || isempty(ETP_input_spreadsheet)
            try
                input_table = readtable('Input_Data_Sheets/ETP_input_data.xlsx');
            catch
                input_table = readtable('Input_Data_Sheets\ETP_input_data.xlsx');
            end
        else
            if exist(ETP_input_spreadsheet,'file') ~= 2
                error(['Bypass mode requested internal ETP forcing, but the file was not found:\n  %s'], ...
                      ETP_input_spreadsheet);
            end

            input_table = readtable(ETP_input_spreadsheet);
        end

    else

        try
            input_table = readtable('Input_Data_Sheets/ETP_input_data.xlsx');
        catch
            input_table = readtable('Input_Data_Sheets\ETP_input_data.xlsx');
        end

    end

    % -------------------------
    % Time information
    % Clip spreadsheet ETP forcing to model simulation window
    % -------------------------
    etp_data_start_row = 3;
    etp_data_rows_all  = (etp_data_start_row:size(input_table,1))';
    
    raw_time_ETP = datetime(table2array(input_table(etp_data_rows_all,2)));
    raw_time_ETP = raw_time_ETP(:);
    
    raw_n_obs_ETP = numel(raw_time_ETP);
    
    if raw_n_obs_ETP < 2
        error('The ETP forcing table must contain at least 2 time rows.');
    end
    
    if any(isnat(raw_time_ETP))
        error('The ETP forcing table contains invalid datetime values in column 2.');
    end
    
    if any(diff(raw_time_ETP) <= duration(0,0,0))
        error('The ETP forcing timestamps must be strictly increasing.');
    end
    
    model_start_datetime = date_begin;
    model_end_datetime   = date_begin + minutes(running_control.routing_time);
    
    % ---------------------------------------------------------------------
    % Keep one anchor row at or before date_begin.
    % This guarantees that the model has a valid ETP value at model time 0.
    % ---------------------------------------------------------------------
    idx_etp_start = find(raw_time_ETP <= model_start_datetime, 1, 'last');
    
    if isempty(idx_etp_start)
        error(['ETP spreadsheet starts after date_begin. ', ...
               'First ETP time = %s, date_begin = %s. ', ...
               'The ETP time-series must start at or before date_begin.'], ...
               char(raw_time_ETP(1)), char(model_start_datetime));
    end
    
    % ---------------------------------------------------------------------
    % Keep one sentinel row at or after the model end.
    % This prevents the internal ETP updater from reaching the final ETP row
    % too early and treating the forcing as ended.
    % ---------------------------------------------------------------------
    idx_etp_end = find(raw_time_ETP >= model_end_datetime, 1, 'first');
    
    if isempty(idx_etp_end)
        error(['ETP spreadsheet ends before the model simulation ends. ', ...
               'Last ETP time = %s, model end = %s. ', ...
               'The ETP time-series must cover the complete simulation.'], ...
               char(raw_time_ETP(end)), char(model_end_datetime));
    end
    
    % Rows retained from the spreadsheet.
    % This includes:
    %   - the last row at/before date_begin
    %   - all rows inside the model window
    %   - the first row at/after the model end
    idx_etp_keep = idx_etp_start:idx_etp_end;
    etp_rows_to_read = etp_data_rows_all(idx_etp_keep);
    
    ETP_Parameters.time_ETP = raw_time_ETP(idx_etp_keep);
    ETP_Parameters.n_obs_ETP = numel(ETP_Parameters.time_ETP);
    
    if ETP_Parameters.n_obs_ETP < 2
        error('After clipping to the model window, fewer than 2 ETP records remain.');
    end
    
    ETP_Parameters.time_ETP_begin = ETP_Parameters.time_ETP(1);
    
    % Use the actual clipped timestamps relative to date_begin.
    % This is safer than rebuilding an artificial 0:dt:end vector.
    ETP_Parameters.climatologic_spatial_duration = ...
        minutes(ETP_Parameters.time_ETP - model_start_datetime);
    
    % The first retained row may be before date_begin.
    % Assign it to model time 0 so the model has an initial ETP state.
    ETP_Parameters.climatologic_spatial_duration(1) = 0;
    
    % Compatibility fields used elsewhere / for diagnostics
    if ETP_Parameters.n_obs_ETP >= 2
        ETP_Parameters.time_step_etp = median(diff(ETP_Parameters.climatologic_spatial_duration), 'omitnan');
    else
        ETP_Parameters.time_step_etp = NaN;
    end
    
    ETP_Parameters.end_etp = ETP_Parameters.climatologic_spatial_duration(end);
    ETP_Parameters.delta_ETP_date = minutes(ETP_Parameters.time_ETP_begin - date_begin);
    
    % ---------------------------------------------------------------------
    % Console status
    % ---------------------------------------------------------------------
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('ETP spreadsheet time-series clipping\n');
    fprintf('============================================================\n');
    fprintf('Model start date_begin : %s\n', char(model_start_datetime));
    fprintf('Model end datetime     : %s\n', char(model_end_datetime));
    fprintf('Raw ETP first time     : %s\n', char(raw_time_ETP(1)));
    fprintf('Raw ETP last time      : %s\n', char(raw_time_ETP(end)));
    fprintf('Raw ETP records        : %d\n', raw_n_obs_ETP);
    fprintf('Retained ETP records   : %d\n', ETP_Parameters.n_obs_ETP);
    fprintf('Removed before window  : %d\n', idx_etp_start - 1);
    fprintf('Removed after window   : %d\n', raw_n_obs_ETP - idx_etp_end);
    fprintf('Retained first time    : %s -> model time %.3f min\n', ...
            char(ETP_Parameters.time_ETP(1)), ETP_Parameters.climatologic_spatial_duration(1));
    fprintf('Retained last time     : %s -> model time %.3f min\n', ...
            char(ETP_Parameters.time_ETP(end)), ETP_Parameters.climatologic_spatial_duration(end));
    
    if ETP_Parameters.time_ETP(1) < model_start_datetime
        fprintf(['Status note           : first retained ETP row is before date_begin ', ...
                 'and was assigned model time 0 for continuity.\n']);
    end
    
    if ETP_Parameters.time_ETP(end) > model_end_datetime
        fprintf(['Status note           : one ETP row after the model end was retained ', ...
                 'as a sentinel row.\n']);
    end
    
    fprintf('ETP clipping status    : OK\n');
    fprintf('============================================================\n\n');

    % -------------------------
    % Allocate storage
    % -------------------------
    Maps.Hydro.ETP_save = zeros(size(dem,1), size(dem,2), saver_memory_maps);
    Maps.Hydro.ETR_save = zeros(size(dem,1), size(dem,2), saver_memory_maps);

% -------------------------
% Read station data
% Each station block has 6 columns:
% Tmax | Tmin | Tavg | u2 | RH | G
%
% Important:
%   Use etp_rows_to_read so station data is clipped exactly the same way
%   as ETP_Parameters.time_ETP.
% -------------------------
ETP_Parameters.n_max_etp_stations = 50000;
ETP_Parameters.n_stations = 0;

% Clear these fields before filling them, so old workspace values cannot
% survive if fewer rows/stations are retained after clipping.
ETP_Parameters.maxtemp_stations = [];
ETP_Parameters.mintemp_stations = [];
ETP_Parameters.avgtemp_stations = [];
ETP_Parameters.u2_stations      = [];
ETP_Parameters.ur_stations      = [];
ETP_Parameters.G_stations       = [];
ETP_Parameters.coordinates_stations = [];

for i = 1:ETP_Parameters.n_max_etp_stations
    try
        ETP_Parameters.maxtemp_stations(:,i) = table2array(input_table(etp_rows_to_read,6*(i-1) + 3));
        ETP_Parameters.mintemp_stations(:,i) = table2array(input_table(etp_rows_to_read,6*(i-1) + 4));
        ETP_Parameters.avgtemp_stations(:,i) = table2array(input_table(etp_rows_to_read,6*(i-1) + 5));
        ETP_Parameters.u2_stations(:,i)      = table2array(input_table(etp_rows_to_read,6*(i-1) + 6));
        ETP_Parameters.ur_stations(:,i)      = table2array(input_table(etp_rows_to_read,6*(i-1) + 7));
        ETP_Parameters.G_stations(:,i)       = table2array(input_table(etp_rows_to_read,6*(i-1) + 8));

        % Coordinates as stored in your original layout
        ETP_Parameters.coordinates_stations(i,1) = table2array(input_table(1,6*(i-1) + 6));
        ETP_Parameters.coordinates_stations(i,2) = table2array(input_table(1,6*(i-1) + 8));

        ETP_Parameters.n_stations = i;
    catch
        break
    end
end


% Confirm clipped station arrays match clipped time vector
if size(ETP_Parameters.maxtemp_stations,1) ~= ETP_Parameters.n_obs_ETP
    error('Internal ETP clipping failed: station data rows do not match clipped ETP time rows.');
end

fprintf('ETP station data rows  : %d\n', size(ETP_Parameters.maxtemp_stations,1));
fprintf('ETP stations detected  : %d\n', ETP_Parameters.n_stations);
fprintf('ETP station clipping   : OK\n\n');

    if ETP_Parameters.n_stations < 1
        error('No valid ETP meteorological stations were found in the ETP forcing table.');
    end

    % -------------------------
    % Latitude / Longitude from final DEM reference
    % -------------------------
    [height_etp, width_etp] = size(ETP_Parameters.DEM_etp);
    [cols_etp, rows_etp] = meshgrid(1:width_etp, 1:height_etp);

    xWorldLimits = ETP_Parameters.SpatialRef.XWorldLimits;
    yWorldLimits = ETP_Parameters.SpatialRef.YWorldLimits;

    R = imref2d([height_etp, width_etp], xWorldLimits, yWorldLimits);

    [ETP_Parameters.x_etp, ETP_Parameters.y_etp] = intrinsicToWorld(R, cols_etp, rows_etp);

    if isempty(ETP_Parameters.SpatialRef.ProjectedCRS)
        error(['DEM_raster.georef.SpatialRef.ProjectedCRS is empty. ' ...
               'A valid projected CRS is required to compute latitude/longitude for ETP preprocessing.']);
    end

    [ETP_Parameters.lat, ETP_Parameters.lon] = projinv( ...
        ETP_Parameters.SpatialRef.ProjectedCRS, ...
        ETP_Parameters.x_etp, ...
        ETP_Parameters.y_etp);

    ETP_Parameters.lat(idx_nan) = nan;
    ETP_Parameters.lon(idx_nan) = nan;

    Wshed_Properties.pixel_latitude = ETP_Parameters.lat;

elseif flags.flag_input_ETP_map == 1

    % ------------------------------------------------------------
    % MODE B: INPUT EVAPORATION / TRANSPIRATION RASTER MAPS
    % ------------------------------------------------------------
    % Expected input from input_data_script:
    %   Input_Evaporation.time           [minutes]
    %   Input_Evaporation.labels_Directory
    %   Input_Transpiration.time         [minutes]
    %   Input_Transpiration.labels_Directory
    %
    % Excel block says maps are in mm/day.
    % Conversion to model-step depth should be done later as:
    %   mm_per_step = mm_per_day * (time_step_model / 1440)
    % ------------------------------------------------------------

    % -------------------------
    % Sanity checks
    % -------------------------
    if ~exist('Input_Evaporation','var') || ~exist('Input_Transpiration','var')
        error('flag_input_ETP_map=1 but Input_Evaporation and/or Input_Transpiration were not created in input_data_script.');
    end

    if Input_Evaporation.num_obs_maps < 1 || Input_Transpiration.num_obs_maps < 1
        error('flag_input_ETP_map=1 but no evaporation/transpiration rasters were provided.');
    end

    if Input_Evaporation.num_obs_maps ~= Input_Transpiration.num_obs_maps
        error('Evaporation and transpiration raster series must have the same number of maps.');
    end

    if numel(Input_Evaporation.time) < 2 || numel(Input_Transpiration.time) < 2
        error('Need at least 2 evaporation/transpiration maps to compute ETP time step.');
    end

    % -------------------------
    % Time information
    % -------------------------
    Spatial_ETP_Parameters.num_obs_maps = Input_Evaporation.num_obs_maps;
    Spatial_ETP_Parameters.time_input   = Input_Evaporation.time(:)';   % minutes from model start
    Spatial_ETP_Parameters.time_step_spatial = Input_Evaporation.time(2) - Input_Evaporation.time(1);
    Spatial_ETP_Parameters.end_etp = Spatial_ETP_Parameters.time_input(end);

    % Save times / aggregation times
    Spatial_ETP_Parameters.time_save = 0:running_control.record_time_spatial_etp:Spatial_ETP_Parameters.end_etp;

    n_spatial_agg = ceil(running_control.record_time_spatial_etp / Spatial_ETP_Parameters.time_step_spatial);
    etp_spatial_aggregation = zeros(size(dem,1), size(dem,2), n_spatial_agg);

    % Keep compatibility with old naming where needed
    ETP_Parameters.climatologic_spatial_duration = Spatial_ETP_Parameters.time_input;

    % -------------------------
    % Preallocate storage
    % -------------------------
    Maps.Hydro.spatial_evaporation_maps   = zeros(size(dem,1), size(dem,2), saver_memory_maps);
    Maps.Hydro.spatial_transpiration_maps = zeros(size(dem,1), size(dem,2), saver_memory_maps);
end

%% ------------ Recording time of outputs (i.e., flows, concentrations ...) ------------
% Calculations
if flags.flag_real_time_satellite_rainfall == 1
    running_control.routing_time = 360; % to save every six hours, 4 register per day.
    running_control.record_time_maps_2 = 60; % to save the last 6 maps for each register.
    running_control.record_time_hydrographs = 60; % to save the last 6 records for each register.
    running_control.steps = running_control.routing_time/time_step_model; % number of calculation steps
    running_control.number_of_records = floor(running_control.steps*time_step_model/running_control.record_time_maps_2); % number of stored data (size of the vector)
    % Checking if the recording time is working
    if running_control.number_of_records == 0
        error('The recording time is larger than the routing time, please change it')
    end
    running_control.time_records = (0:running_control.record_time_maps_2:running_control.routing_time); % time in minutes
    running_control.time_record_hydrograph = (0:running_control.record_time_hydrographs:running_control.routing_time); % time in minutes
    % running_control.time_change_records = (0:running_control.time_step_change/60:running_control.routing_time); % time in minutes
    % vector to store data
    running_control.time_store = running_control.time_records./time_step_model; % number of steps necessary to reach the record vector
    running_control.time_store(1) = 1; % the zero is the firt time step
    running_control.time_step_save = zeros(running_control.steps,1);
    Courant_Parameters.alfa_save = zeros(running_control.steps,1);
else
    running_control.steps = running_control.routing_time/time_step_model; % number of calculation steps
    running_control.number_of_records = floor(running_control.steps*time_step_model/running_control.record_time_maps); % number of stored data (size of the vector)
    % Checking if the recording time is working
    if running_control.number_of_records == 0
        error('The recording time is larger than the routing time, please change it')
    end
    running_control.time_records = (0:running_control.record_time_maps:running_control.routing_time); % time in minutes
    running_control.time_record_hydrograph = (0:running_control.record_time_hydrographs:running_control.routing_time); % time in minutes
    % running_control.time_change_records = (0:running_control.time_step_change/60:running_control.routing_time); % time in minutes
    % vector to store data
    running_control.time_store = running_control.time_records./time_step_model; % number of steps necessary to reach the record vector
    running_control.time_store(1) = 1; % the zero is the firt time step
    running_control.time_step_save = zeros(running_control.steps,1);
    Courant_Parameters.alfa_save = zeros(running_control.steps,1);
end

%% Inflow and Precipitation data
Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_hydrograph_rate*flags.flag_inflow;
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~=1 % Only for concentrated rainfall
    Rainfall_Parameters.intensity_rainfall = Rainfall_Parameters.intensity_rainfall*flags.flag_rainfall; % If flags.flag_rainfall is zero, no rainfall is considered
end
%% DEM Treatment and Filtering Algorithms
% Fillsinks
% max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
% DEM_filled = fillsinks(DEM,max_depth);
if flags.flag_fill_DEM == 1
    DEM_filled = fillsinks_for_routing(DEM_raster, flags); % Filled DEM
    DIFFDEM = DEM_filled - DEM;
    DIFFDEM.Z(DIFFDEM.Z==0) = nan;
    DEM_raster = DEM_filled;
    % imageschs(DEM_raster,DIFFDEM.Z);
end

% DTM Filter
if flags.flag_DTM == 1
    slope_threshold = GIS_data.slope_DTM; % percentage
    DEM_filtered = DTM_Filter(DEM_raster.Z,Wshed_Properties.Resolution,GIS_data.slope_DTM);
    DEM_raster.Z = DEM_filtered;
end

% Smooth DEM
if flags.flag_smooth_cells == 1
    Vq = imgaussfilt(DEM_raster.Z,'FilterSize',3);
    dem_diff_smooth = Vq;
    dem = Vq; % New dem
    dem(isnan(dem) & ~isnan(DEM_raster.Z)) = DEM_raster.Z(isnan(dem) & ~isnan(DEM_raster.Z));
    DEM_raster.Z = dem;
    close all
end

% Fill Again to Make Sure
% max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
% DEM_filled = fillsinks(DEM,max_depth);
if flags.flag_fill_DEM == 1
    DEM_filled = fillsinks_for_routing(DEM_raster, flags); % Filled DEM
    DIFFDEM = DEM_filled - DEM;
    DIFFDEM.Z(DIFFDEM.Z==0) = nan;
    DEM_raster = DEM_filled;
    % imageschs(DEM_raster,DIFFDEM.Z);
end

% %% D4 Flow Accumulation
% if flags.flag_D8 ~= 1
%     [flowAcc, flowDir] = D4_FlowAccumulation(DEM_raster.Z);
% end
%
% %%
if flags.flag_D8 ~= 1
    [flow_accum, flow_dir] = flow_accumulation_4D(DEM_raster.Z);
end

%% Decrese Elevations in Creeks

% Flow directio 8D raster
FD = FLOWobj(DEM_raster); % Flow direction
% Flow accumulation
As  = flowacc(FD); % Flow Accumulation
% Area of each cell in km2
Wshed_Properties.fac_area = As.Z*(Wshed_Properties.cell_area/1000/1000); % km2
% Mask with flow accumulation larger than the threshold
idx_facc = Wshed_Properties.fac_area >= GIS_data.min_area;
% Drainage Basins of each gauge
% for ii = 1:gauges.num_obs_gauges
%    L = drainagebasins(FD,gauges.northing_obs_gauges(ii),gauges.easting_obs_gauges(ii));
% end

% River Width (only using the first entry)
B = GIS_data.beta_1(1)*Wshed_Properties.fac_area.^(GIS_data.beta_2(1));
B(~idx_facc) = 0;

% River Heigth (only using the first entry)
H = GIS_data.alfa_1(1)*Wshed_Properties.fac_area.^(GIS_data.alfa_2(1));
H(~idx_facc) = 0;

% River inbank area
Flow_Area = B.*H; % m2

% River width raster
Wshed_Properties.River_Width = B;
Wshed_Properties.River_Width(isnan(Wshed_Properties.River_Width)) = 0;
% Wshed_Properties.River_Width(~idx_facc) = 0;

% River height raster
Wshed_Properties.River_Depth = H;
Wshed_Properties.River_Depth(isnan(Wshed_Properties.River_Depth)) = 0;

if flags.flag_D8 ~= 1
    [idx_facc, Wshed_Properties.River_Width, Wshed_Properties.River_Depth] = enforce_4D_flow(idx_facc, Wshed_Properties.River_Width, Wshed_Properties.River_Depth);
    Wshed_Properties.River_Width(isnan(Wshed_Properties.River_Width)) = 0;
    Wshed_Properties.River_Depth(isnan(Wshed_Properties.River_Depth)) = 0;
end

% Wshed_Properties.River_Depth(~idx_facc) = 0;

if flags.flag_subgrid && flags.flag_river_rasters
    Wshed_Properties.River_Width = double(widths_raster.Z);
    Wshed_Properties.River_Width(isnan(Wshed_Properties.River_Width)) = 0;
    Wshed_Properties.River_Depth = double(depths_raster.Z);
    Wshed_Properties.River_Depth(isnan(Wshed_Properties.River_Depth)) = 0;
end

% Reducing Elevation in creeks (if flag reduce DEM is used and not in the subgrid
% is deactivated
if flags.flag_reduce_DEM == 1 && flags.flag_subgrid == 0
    % H_abg = Flow_Area/Wshed_Properties.Resolution; % Old Options
    H_abg = ((Wshed_Properties.River_Width/Wshed_Properties.Resolution) .* Wshed_Properties.River_Depth .^(5/3)) .^ (3/5); % Based on momentum with the same n and S0
    if max(max(H_abg)) > 10^2
        error('Be careful. The decreasing in the DEM to consider the water depths are larger than 100 m. You may change the parameters.')
    end
    DEM_raster.Z = DEM_raster.Z - H_abg; % [m]
end

% New Data
dem = DEM_raster.Z;

%% DEM Smoothening
if flags.flag_smoothening == 1 % Valid only for streams determined with the flow accumulation
    [DEM_raster,DEM,S] = DEM_smoothening(DEM_raster,GIS_data.min_area,flags.flag_trunk,GIS_data.tau,GIS_data.K_value);
    if flags.flag_fill_DEM == 1
        DEM_filled = fillsinks_for_routing(DEM_raster, flags);
        DIFFDEM = DEM_filled - DEM_raster.Z;
        DIFFDEM.Z(DIFFDEM.Z==0) = nan;
        DEM_raster = DEM_filled;
        DEM = DEM_raster.Z;
    end
end

%% Converting to 4D (if required)
% Converting to a 4D flow path if required
if flags.flag_D8 ~= 1
    [flow_mask_4D, B, H] = Transform8Dto4DFlow(double(idx_facc), B, H);

    % River width raster
    Wshed_Properties.River_Width = B;
    Wshed_Properties.River_Width(isnan(Wshed_Properties.River_Width)) = 0;
    % Wshed_Properties.River_Width(~idx_facc) = 0;

    % River height raster
    Wshed_Properties.River_Depth = H;
    Wshed_Properties.River_Depth(isnan(Wshed_Properties.River_Depth)) = 0;
    % Wshed_Properties.River_Depth(~idx_facc) = 0;
end

%% Imposing Minimum Slope - If Required
if flags.flag_diffusive ~= 1 && flags.flag_inertial == 1 && flags.flag_kinematic == 1
    % Impose Mininum Slope
    if flags.flag_smoothening ~=1 % We don't have a S, so we need to calculate it
        FD = FLOWobj(DEM_raster);
        area_km2 = GIS_data.min_area; % km2
        area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
        S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
    end
    DEM_new = imposemin(S,DEM_raster,GIS_data.sl);
    DEM_DIFF = DEM_new - DEM_raster;
    DEM_raster = DEM_new;
    %     imagesc(DEM_DIFF); colorbar; % if you want to plot
end
%% Slope Map
% SLP = arcslope(DEM_raster);
% pause(0.25)
% imagesc(SLP); colorbar
% pause(0.25)

%% Observed Gauges - Catchment Area
if flags.flag_obs_gauges == 1
    % Flow Accumulatin Calculation
    FD = FLOWobj(DEM_raster); % Flow direction
    As  = flowacc(FD); % Flow Accumulation
    Wshed_Properties.fac_area = As.Z*(Wshed_Properties.cell_area/1000/1000); % km2

    % Catchment area of each gauge
    for i = 1:length(gauges.easting_obs_gauges)
        gauges.catchment_area(i,1) = Wshed_Properties.fac_area(gauges.northing_obs_gauges(i,1),gauges.easting_obs_gauges(i,1)); % km2
    end
end

%% DEM with Streams and Observed Points

FD = FLOWobj(DEM_raster);
area_km2 = GIS_data.min_area; % km²
area_cells = area_km2 / ((DEM_raster.cellsize / 1000)^2); % Convert to DEM cells
S = STREAMobj(FD, 'minarea', area_cells); % Generate stream network

% Check if stream object is empty
if isempty(S.x)
    warning('Flow accumulation threshold too high. Lowering threshold to half of domain area.');

    % Lower threshold
    domain_area_cells = numel(DEM_raster.Z);
    area_cells = domain_area_cells / 2;
    S = STREAMobj(FD, 'minarea', area_cells);

    if isempty(S.x)
        error('Stream network still empty after reducing threshold. Check DEM quality or threshold.');
    end
end

MS = STREAMobj2mapstruct(S);

% Ensure output folder exists
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

FileName_String = 'streamnetwork.shp';
FileName = fullfile(folderName, FileName_String);

try
    shapewrite(MS, FileName);
catch ME
    error('Unable to write shapefile: %s\nPossible causes: invalid path, file in use, or insufficient permissions.', ME.message);
end

%% Generate Synthetic Gauges Based on Flow Accumulation Percentiles

if flags.flag_obs_gauges ~= 1 || ...
   (use_inputdata_bypass ~= 1 && ~isfile(model_folder))

    flags.flag_obs_gauges = 1;
    warning('No valid observation file found. Using synthetic stream gauges based on flow accumulation percentiles.');

    % Flow accumulation map (in number of cells)
    FA = flowacc(FD);

    % Extract accumulation values only at stream cells
    FA_stream = FA.Z(S.IXgrid);  % flow accumulation at stream pixels

    % Define desired percentiles (descending: 95%, 90%, ..., 5%)
    target_percentiles = [95:-10:10 5];
    numPoints = numel(target_percentiles);
    total_cells = numel(find(~isnan(FA.Z)));
    sorted_FA = sort(FA_stream, 'descend');

    % Get threshold values at percentiles
    thresholds = prctile(FA_stream, target_percentiles);

    % For each percentile, find closest match in stream cells
    gauges_idx = zeros(1, numPoints);
    for k = 1:numPoints
        [~, i_closest] = min(abs(FA_stream - thresholds(k)));
        gauges_idx(k) = i_closest;
    end

    % Coordinates of selected gauges
    gauges.x_coord_gauges = S.x(gauges_idx);
    gauges.y_coord_gauges = S.y(gauges_idx);
    gauges.num_obs_gauges = numPoints;

    % Labels
    gauges.labels_observed_string = arrayfun(@(p) sprintf('Gauge - %d%%', p), target_percentiles, 'UniformOutput', false);

    % Dummy attributes
    GIS_data.alfa_1 = zeros(numPoints,1);
    GIS_data.alfa_2 = zeros(numPoints,1);
    GIS_data.beta_1 = zeros(numPoints,1);
    GIS_data.beta_2 = zeros(numPoints,1);
    River_Manning = 0.03 * ones(numPoints,1);
    Lateral_Groundwater_Flux = zeros(numPoints,1);

    % Convert absolute to pixel coordinates
    gauges.easting_obs_gauges_absolute = gauges.x_coord_gauges;
    gauges.northing_obs_gauges_absolute = gauges.y_coord_gauges;
    gauges.easting_obs_gauges = round((-GIS_data.xulcorner + gauges.x_coord_gauges) / Wshed_Properties.Resolution);
    gauges.northing_obs_gauges = round((GIS_data.yulcorner - gauges.y_coord_gauges) / Wshed_Properties.Resolution);
end

close all

%% Drainage Basins
try
    D = drainagebasins(FD);
    imageschs(DEM_raster,shufflelabel(D))
    area_km2 = GIS_data.min_area;
    area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
    imagesc(DEM_raster);
    hold on
    S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
    if ~isempty(S)
        plot(S,'linewidth',2,'color','red');
        plot(trunk(S),'linewidth',2,'color','red');
        D = drainagebasins(FD,S);
        imageschs(DEM_raster,shufflelabel(D))
    end

    hold on

    plot(S,'linewidth',2,'color','red');% area_km2 = 500; % km2
    area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
    imagesc(DEM_raster);
    hold on
    S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
    plot(S,'linewidth',2,'color','red');

    % D = drainagebasins(FD,S);
    % D = drainagebasins(FD,gauges.easting_obs_gauges_absolute(2),gauges.northing_obs_gauges(2))

    imageschs(DEM_raster,shufflelabel(D))
    hold on

    plot(S,'linewidth',2,'color','red');

catch
    warning('It was impossible to generate streams with the threshold assumed.')
end


%% Outlet Locations
% ------------ Outlet ------------ %
outlet_index = zeros(zzz);
elevation_nan = dem;
elevation_nan(elevation_nan < min_dem_value) = nan; % Replacing no-info to nan
dem = elevation_nan;
perimeter = zeros(size(dem));
% Identifying Watershed Perimeter
for i = 1:size(dem,1)
    for j = 1:size(dem,2)
        if i == 1 || i == size(dem,1) || j == 1 || j == size(dem,2) % Boundaries of the domain
            if dem(i,j) > min_dem_value % Boundary has values
                perimeter(i,j) = 1;
            else
                perimeter(i,j) = 0;
            end
        end
        if i ~= 1 && i ~= size(dem,1) && j ~= 1 && j ~= size(dem,2)
            if isnan(dem(i,j-1)) && isnan(dem(i,j+1)) && isnan(dem(i-1,j)) && isnan(dem(i+1,j))
                perimeter(i,j) = 0;
            elseif  dem(i,j) > min_dem_value && (isnan(dem(i,j-1)) || isnan(dem(i,j+1)) || isnan(dem(i-1,j)) || isnan(dem(i+1,j)))
                perimeter(i,j) = 1;
            else
                perimeter(i,j) = 0;
            end
        end
    end
end
[row_boundary, col_boundary] = find(perimeter > 0); % Boundary Cells
max_length = (max(col_boundary) - min(col_boundary))*Wshed_Properties.Resolution - (max(row_boundary) - min(row_boundary))*Wshed_Properties.Resolution; % Length
max_width  = (max(row_boundary) - min(row_boundary)) * Wshed_Properties.Resolution;
max_length = (max(col_boundary) - min(col_boundary)) * Wshed_Properties.Resolution;

% Number of Outlets
idx_outlet = elevation_nan == min(min(elevation_nan)); % Finding cells with lowest elevations
[Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(idx_outlet == 1); % finding rows and cols of these cells
outlet_index(idx_outlet) = 1;

% Perimeter
Wshed_Properties.perimeter = perimeter;
Wshed_Properties.domain = dem > min_dem_value; % Cells that are part of the domain


%% Outlet Calculations - 1  -- REVISED

% Number of requested outlet cells
if ~exist('n_outlets_data','var') || isempty(n_outlets_data)
    n_outlets_data = 1;
end
n_outlets = max(1, round(n_outlets_data));

% Only perimeter cells can be outlets
valid_perimeter = perimeter > 0 & ~isnan(dem);

if ~any(valid_perimeter(:))
    error('No valid perimeter cells found for outlet placement.');
end

outlet_from_csv = false;
outlet_csv_path = '';
if exist('InputPaths','var') && isstruct(InputPaths) && ...
        isfield(InputPaths,'Outlet_Cells_CSV') && ~isempty(InputPaths.Outlet_Cells_CSV) && ...
        exist(InputPaths.Outlet_Cells_CSV,'file') == 2
    outlet_csv_path = InputPaths.Outlet_Cells_CSV;
    try
        outlet_table = readtable(outlet_csv_path);
        outlet_vars = outlet_table.Properties.VariableNames;
        x_candidates = {'x_m','x','easting_m','easting','Easting','X'};
        y_candidates = {'y_m','y','northing_m','northing','Northing','Y'};
        x_name = '';
        y_name = '';
        for ii = 1:numel(x_candidates)
            pos = find(strcmpi(outlet_vars, x_candidates{ii}), 1, 'first');
            if ~isempty(pos)
                x_name = outlet_vars{pos};
                break
            end
        end
        for ii = 1:numel(y_candidates)
            pos = find(strcmpi(outlet_vars, y_candidates{ii}), 1, 'first');
            if ~isempty(pos)
                y_name = outlet_vars{pos};
                break
            end
        end
        if isempty(x_name) || isempty(y_name)
            error('Outlet CSV must contain x_m/y_m or equivalent coordinate columns.');
        end
        x_abs = double(outlet_table.(x_name));
        y_abs = double(outlet_table.(y_name));
        col_csv = round((-GIS_data.xulcorner + x_abs) ./ Wshed_Properties.Resolution);
        row_csv = round((GIS_data.yulcorner - y_abs) ./ Wshed_Properties.Resolution);
        keep = isfinite(row_csv) & isfinite(col_csv) & ...
               row_csv >= 1 & row_csv <= size(dem,1) & ...
               col_csv >= 1 & col_csv <= size(dem,2);
        row_csv = row_csv(keep);
        col_csv = col_csv(keep);
        if ~isempty(row_csv)
            csv_linear = sub2ind(size(dem), row_csv, col_csv);
            csv_linear = unique(csv_linear(:));
            csv_linear = csv_linear(~isnan(dem(csv_linear)));
            if ~isempty(csv_linear)
                outlet_index = false(size(dem));
                outlet_index(csv_linear) = true;
                outlet_from_csv = true;
                fprintf('Using %d outlet cells from CSV: %s\n', nnz(outlet_index), outlet_csv_path);
            end
        end
        if ~outlet_from_csv
            warning('Outlet CSV was found but no valid cells matched the DEM grid; falling back to automatic outlet placement.');
        end
    catch ME
        warning('Could not use outlet CSV "%s": %s. Falling back to automatic outlet placement.', outlet_csv_path, ME.message);
    end
end

if ~outlet_from_csv
    % Start from the lowest DEM cell ON THE PERIMETER
    perim_vals = dem(valid_perimeter);
    min_perim_el = min(perim_vals);
    
    candidate_seed = valid_perimeter & dem == min_perim_el;
    [row_seed, col_seed] = find(candidate_seed, 1, 'first');
    
    % Reset outlet mask completely
    % This is important: otherwise interior minimum cells remain outlets.
    outlet_index = false(size(dem));
    outlet_index(row_seed, col_seed) = true;
    
    % Expand outlets along connected perimeter cells
    neighbor_kernel = ones(3,3);
    
    while nnz(outlet_index) < n_outlets
    
        neighbor_mask = conv2(double(outlet_index), neighbor_kernel, 'same') > 0;
    
        candidate_mask = neighbor_mask & valid_perimeter & ~outlet_index;
    
        if ~any(candidate_mask(:))
            warning('Requested %d outlets, but only %d connected perimeter outlets could be placed.', ...
                    n_outlets, nnz(outlet_index));
            break
        end
    
        % Prefer the lowest neighboring perimeter cell
        candidate_vals = dem;
        candidate_vals(~candidate_mask) = nan;
    
        min_candidate_el = min(candidate_vals(candidate_mask));
    
        next_mask = candidate_mask & candidate_vals == min_candidate_el;
        [r_next, c_next] = find(next_mask, 1, 'first');
    
        outlet_index(r_next, c_next) = true;
    end
end

% Final outlet variables
idx_outlet = logical(outlet_index);

[Wshed_Properties.row_outlet, Wshed_Properties.col_outlet] = find(idx_outlet);

Wshed_Properties.el_outlet = dem;
Wshed_Properties.el_outlet(~idx_outlet) = nan;

outlet_vals = Wshed_Properties.el_outlet(idx_outlet);
outlet_vals = outlet_vals(~isnan(outlet_vals));

Wshed_Properties.stage_min = min(outlet_vals);

[Wshed_Properties.row_min, Wshed_Properties.col_min] = ...
    find(idx_outlet & dem == Wshed_Properties.stage_min);

% % Save Map
FileName = 'Outlet_Mask.tif';
FileName = fullfile(Paths.Shapes,FileName);
% Exporting Outlet_Index as a Mask
raster_to_export = DEM_raster; % Just to get the properties
raster_to_export.Z = outlet_index; % Adding Outlets
try
    raster_to_export.Z(isnan(DEM_raster.Z)) = nan;
end
try
    GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
catch me
    warning('Outlet mask not exported.')
end

% Inlet Mask
if flags.flag_inflow == 1
    FileName = 'Inlet_Mask.tif';
    FileName = fullfile(Paths.Shapes,FileName);
    % Exporting Outlet_Index as a Mask
    raster_to_export = DEM_raster; % Just to get the properties
    temp = 0*raster_to_export.Z;
    for i = 1:Inflow_Parameters.n_stream_gauges
        temp = Wshed_Properties.inflow_cells(:,:,i) + temp;
    end
    temp = double(temp>0);
    raster_to_export.Z = temp; % Adding Outlets
    raster_to_export.Z(isnan(DEM_raster.Z)) = nan;
    try
        GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
    catch me
        warning('Inlet mask not exported.')
    end
end

%% Second Way - Outlet Calculation
% Here we find the k-lowest points in the perimeter of the catchment
% min_el = min(min(dem)); % Minimum elevation
% min_el_out = min_el;
% for i = 1:n_outlets_data
%     [row_new,col_new] = find(dem > min_el_out,1,'first');
%     outlet_index(row_new,col_new) = 1;
%     min_el_out = dem(row_new,col_new);
% end
% [Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1);  % Refreshing outlet
%
%
% % % Filling Inside Gaps
% [row_outlets,col_outlets] = find(outlet_index == 1);
% for i = 1:sum(sum(outlet_index))
%     row_out = row_outlets(i);
%     col_out = col_outlets(i);
%     % 0, 45, 90, 135, 180, 225, 270, 315, 360
%     r = row_out; c = col_out + 1;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 0
%         outlet_index(r,c) = 1;
%     end
%     r = row_out-1; c = col_out+1;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 45
%         outlet_index(r,c) = 1;
%     end
%     r = row_out-1; c = col_out;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 90
%         outlet_index(r,c) = 1;
%     end
%     r = row_out-1; c = col_out-1;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 135
%         outlet_index(r,c) = 1;
%     end
%     r = row_out; c = col_out-1;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 180
%         outlet_index(r,c) = 1;
%     end
%     r = row_out+1; c = col_out-1;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 225
%         outlet_index(r,c) = 1;
%     end
%     r = row_out+1; c = col_out;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 270
%         outlet_index(r,c) = 1;
%     end
%     r = row_out+1; c = col_out+1;
%     if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 315
%         outlet_index(r,c) = 1;
%     end
% end
%
%
% % outlet_index = perimeter; % DELETE HERE
% [Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1); % finding rows and cols of these cells
% Wshed_Properties.el_outlet = dem;
% Wshed_Properties.el_outlet(idx_outlet < 1) = nan;
% Wshed_Properties.stage_min = min(min(Wshed_Properties.el_outlet));
% [Wshed_Properties.row_min, Wshed_Properties.col_min] = find(dem == Wshed_Properties.stage_min);
%
% % Creating Modeling Results Folder
% folderName = 'Modeling_Results';
%
% % Check if the folder already exists
% if ~exist(folderName, 'dir')
%     % If it doesn't exist, create the folder
%     mkdir(folderName);
%     disp('Folder "Modeling_Results" created successfully!');
% else
%     disp('Data sucessfully exported in Modeling_Results Folder');
% end
%
% % Save Map
% FileName = 'Outlet_Mask.tif';
% FileName = fullfile(folderName,FileName);
% % Exporting Outlet_Index as a Mask
% raster_to_export = DEM_raster; % Just to get the properties
% raster_to_export.Z = outlet_index; % Adding Outlets
% raster_to_export.Z(isnan(DEM_raster.Z)) = nan;
% GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
%
% % Save Map
% FileName = 'Perimeter_Mask.tif';
% FileName = fullfile(folderName,FileName);
% % Exporting Outlet_Index as a Mask
% raster_to_export = DEM_raster; % Just to get the properties
% raster_to_export.Z = perimeter; % Adding Outlets
% raster_to_export.Z(isnan(DEM_raster.Z)) = nan;
% GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map

%% Grid Domain
%%%%%% ORIGINAL GRID %%%%%%
xmin = 1; % initial position x in the grid (collums)
ymin = 1; % lines (up to down)
xmax = zzz(2);
ymax = zzz(1);
% ------------ Cutting Cells ------------
if flags.flag_inflow == 1
    Wshed_Properties.inflow_cells = Wshed_Properties.inflow_cells(ymin:ymax,xmin:xmax,:);
end
if flags.flag_stage_hydrograph == 1
    Wshed_Properties.stage_mask = Wshed_Properties.stage_mask(ymin:ymax,xmin:xmax,:);
end
outlet_index = outlet_index(ymin:ymax,xmin:xmax);
elevation = dem(ymin:ymax,xmin:xmax); % using only specified grid
spatial_domain = zeros(size(dem));
dem = dem(ymin:ymax,xmin:xmax);
lulc_matrix = imp(ymin:ymax,xmin:xmax); % using only specified grid
soil_matrix  = soil(ymin:ymax,xmin:xmax); % using only specified grid;
Wshed_Properties.rainfall_matrix = double(elevation >= min_dem_value) ; % Elevation could be negative
Wshed_Properties.rainfall_matrix(Wshed_Properties.rainfall_matrix == 0) = nan;

% ------------ GRID ------------- %
% --- Only for Plane Watershed. Deleted it afterwards ---
x_grid = 0 + Wshed_Properties.Resolution*[1:1:size(dem,2)];
y_grid = 0 + Wshed_Properties.Resolution*[1:1:size(dem,1)];

% Lower Left Corner
xllcorner = GIS_data.xulcorner; % Same
yllcorner = GIS_data.yulcorner - Wshed_Properties.Resolution*(size(dem,1)); % Subtract Y distance

%% Determining soil properties for each cell, according to the imperviousness and Cutting all full-domain matrices
% Testing if Imp and DEM have the same -9999 values
idx_elevation = isnan(elevation); % find where elevation is nan (we already treated elevation to nans)
idx_LULC = lulc_matrix < 0; % find where lulc is -9999 or negative
idx_SOIL = soil_matrix < 0; % find where soil is -9999 or negative
% ----------------- Checking if Data are Equivalent

%% Pre-allocating Arrays and Checking no-data Values
% ------------ New-way ------------:  Most of the time these data have different resolutions and
% no-data values. In this case, if at least one of the bothr rasters (i.e.,
% DEM and IMP have a missing data, we assume both don't have)
idx = logical(idx_elevation + idx_LULC + idx_SOIL); % Logical equation here
idx = idx > 0; % Assuming only values of 0 and 1
idx_nan = idx; %  saving nan_matrix

% Masking Inputs
DEM_raster.Z(idx_nan) = nan;
SOIL_raster.Z(idx_nan) = nan;
LULC_raster.Z(idx_nan) = nan;

soil_matrix(idx_nan) = nan;
lulc_matrix(idx_nan) = nan;

idx_not_nan = idx_nan < 1;
if flags.flag_waterquality == 1
    idx_nan_5 = zeros(zzz(1),zzz(2),5);

    for i = 1:5 % Left, right, up, down, outlet
        idx_nan_5(:,:,i) = idx_nan; % This is used in CA model
    end
    idx_nan_5 = logical(idx_nan_5); % Converting to logical
end
elevation(idx) = nan; % change those values for nan
lulc_matrix(idx) = nan; % change those values for nan
soil_matrix(idx) = nan; % change those values for nan

spatial_domain(idx) = nan; % change those values for inf

%% Struct Variables
% ----------- Struct Cells ------------- %
outflow_rate = struct('qout_left_t',spatial_domain,'qout_right_t',spatial_domain,'qout_up_t',spatial_domain,'qout_down_t',spatial_domain); % Struct array
% depths = struct('d_0',spatial_domain,'d_avg',spatial_domain,'d_left_cell',spatial_domain,'d_right_cell',spatial_domain,'d_up_cell',spatial_domain,'d_down_cell',spatial_domain,'d_tot',spatial_domain,'depth_wse',depth_wse); % struct
depths = struct('d_0',spatial_domain,'d_avg',spatial_domain,'d_tot',spatial_domain,'depth_wse',depth_wse); % struct
velocities = struct('vel_left',spatial_domain,'vel_right',spatial_domain,'vel_up',spatial_domain,'vel_down',spatial_domain,'velocity_raster',spatial_domain,'vmax_final',spatial_domain); % struct

Soil_Properties = struct('ksat',spatial_domain,'n_vg',spatial_domain,'alpha_vg',spatial_domain,'theta_sat',spatial_domain,'theta_r',spatial_domain,'theta_i',spatial_domain,'Sy',spatial_domain,'ksat_gw',spatial_domain,'I_0',spatial_domain,'Soil_Depth',spatial_domain); % struct
LULC_Properties = struct('roughness',spatial_domain,'h_0',spatial_domain,'C_1',spatial_domain,'C_2',spatial_domain,'C_3',spatial_domain,'C_4',spatial_domain,'B_0',spatial_domain,'ADD',0,'root_depth_m',spatial_domain); % struct
Soil_Properties.Ks_multiplier_near_surface = ones(size(spatial_domain));
Soil_Properties.Ks_multiplier_root_zone = ones(size(spatial_domain));
Soil_Properties.Ks_multiplier_transmission = ones(size(spatial_domain));

% Elevation_Properties = struct('elevation_cell',spatial_domain,'elevation_left_t',spatial_domain,'elevation_right_t',spatial_domain,'elevation_up_t',spatial_domain,'elevation_down_t',spatial_domain); % struct

% ------------- Preallocating arrays to avoid large computational efforts
y_1 = length(ymin:ymax);
x_1 = length(xmin:xmax);
if flags.flag_waterquality == 1 % Build-up and Wash-off matrices
    LULC_Properties.C_1 = spatial_domain;
    LULC_Properties.C_2 = spatial_domain;
    LULC_Properties.C_3 = spatial_domain;
    LULC_Properties.C_4 = spatial_domain;
    WQ_States.B_0 = spatial_domain;
end
%% ----------------- Fill Variables ----------------- %%
% Correcting LULC_index
LULC_Properties.n_lulc = n_lulc;
LULC_Properties.ADD = ADD;
LULC_Properties.min_Bt = min_Bt;
LULC_Properties.Bmin = Bmin;
LULC_Properties.Bmax = Bmax;
LULC_Properties.Pol_min = Pol_min;

lulc_matrix(idx_nan) = -9999;
LULC_Properties.idx_lulc = zeros(size(DEM,1),size(dem,2),LULC_Properties.n_lulc);
imp_index_matrix = [];
% Assuming that we know where is the impervious areas
for i = 1:LULC_Properties.n_lulc
    index = lulc_matrix == LULC_index(i,1);
    LULC_Properties.idx_lulc(:,:,i) = index;
    if LULC_index(i,1) == imp_index
        imp_index_matrix = i;
    end
end
LULC_Properties.idx_imp = LULC_Properties.idx_lulc(:,:,imp_index_matrix);
% idx = cat(3,idx_1,idx_2,idx_3,idx_4,idx_5,idx_6); % Concatenating all of them
impervious_cells = sum(sum(LULC_Properties.idx_imp));
pervious_cells = sum(sum(sum(LULC_Properties.idx_lulc))) - impervious_cells;
% Converting to Logical Values
LULC_Properties.idx_lulc = logical(LULC_Properties.idx_lulc);
LULC_Properties.idx_imp = logical(LULC_Properties.idx_imp);
LULC_Properties.ADD = ADD;
for i = 1:LULC_Properties.n_lulc % Types of LULC
    % Only Roughness and h_0
    index = LULC_index(i,1);
    LULC_Properties.roughness(LULC_Properties.idx_lulc(:,:,i)) = lulc_parameters(i,1); % assigning values for roughness at impervious areas
    LULC_Properties.h_0(LULC_Properties.idx_lulc(:,:,i)) = lulc_parameters(i,2); % Initial Abstraction
    if size(lulc_parameters,2) >= 9 && isfinite(lulc_parameters(i,9))
        root_depth_i = lulc_parameters(i,9);
    elseif LULC_index(i,1) == imp_index
        root_depth_i = 0;
    else
        root_depth_i = 1.0;
    end
    LULC_Properties.root_depth_m(LULC_Properties.idx_lulc(:,:,i)) = max(root_depth_i, 0);

    % ------------ Warm-up Data ------------:
    if flags.flag_warmup == 1
        if i == 1
            % Identifying positive values of the Warmup
            try
                Warmup_Raster = GRIDobj(Warmup_Depth_path);
                if sum(size(Warmup_Raster.Z)) ~= sum(size(DEM))
                    Warmup_Raster = resample(Warmup_Raster,DEM_raster);
                end
                Warmup_Depths = double(Warmup_Raster.Z);
                idx_warmup = Warmup_Depths >= 0;
                depths.d_0 =  idx_warmup.*Warmup_Depths*1000; % Depths in m converting to mm
                % Treat inf values
                depths.d_0(idx_nan) = nan; % carefull here
                depths.d_0(idx_nan == 0 & isnan(depths.d_0)) = 0;
                depths.d_t = depths.d_0;
            catch
                warning('Warm-up raster not found. Assuming no initial depths.')
                depths.d_t = depths.d_0;
            end
        end
    else
        depths.d_0(LULC_Properties.idx_lulc(:,:,i)) = lulc_parameters(i,3);
    end
    if flags.flag_waterquality == 1 && flags.flag_WQ_Rasters ~= 1 % Only if water quality is being modeled
        LULC_Properties.C_1(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,4);
        LULC_Properties.C_2(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,5);
        LULC_Properties.C_3(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,6);
        LULC_Properties.C_4(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,7);

    elseif flags.flag_waterquality == 1 && flags.flag_WQ_Rasters == 1 % We are using rasters
        temp = GRIDobj(B1_path); temp(isnan(DEM_raster.Z)) = nan;
        LULC_Properties.C_1 = temp.Z;
        temp = GRIDobj(B2_path); temp(isnan(DEM_raster.Z)) = nan;
        LULC_Properties.C_2 = temp.Z;
        temp = GRIDobj(W1_path); temp(isnan(DEM_raster.Z)) = nan;
        LULC_Properties.C_3 = temp.Z;
        temp = GRIDobj(W2_path); temp(isnan(DEM_raster.Z)) = nan;
        LULC_Properties.C_4 = temp.Z;
    end
    if i == LULC_Properties.n_lulc % Last land use
        % Do we add another mass in the initial buildup or not?
        flags.flag_mass_sum = 0;
        if flags.flag_initial_buildup == 1 && flags.flag_mass_sum == 1 && flags.flag_waterquality == 1
            Initial_Buildup_Raster = GRIDobj(Initial_Buildup_path);
            WQ_States.B_0 = Initial_Buildup_Raster.Z; % kg/cell
            % Treat inf values
            WQ_States.B_0(idx_nan) = inf; % carefull here
            % Another filter
            WQ_States.B_0(WQ_States.B_0 < 0) = 0; % Attention here
            B_0_add = C_1.*(1 - exp(1).^(-C_2*LULC_Properties.ADD )); % kg/ha
            WQ_States.B_t = B_0 + B_0_add*Wshed_Properties.cell_area/10^4; % kg
            Maps.WQ_States.initial_buildup_map = B_t; % kg
        elseif flags.flag_initial_buildup == 1 && flags.flag_mass_sum ~= 1 && flags.flag_waterquality == 1
            Initial_Buildup_Raster = GRIDobj(Initial_Buildup_path);
            WQ_States.B_0 = Initial_Buildup_Raster.Z; % kg/cell
            % Treat inf values
            WQ_States.B_0(idx_nan) = inf; % carefull here
            % Another filter
            WQ_States.B_0(WQ_States.B_0 < 0) = 0; % Attention here
            WQ_States.B_t = B_0;
            Maps.WQ_States.initial_buildup_map = B_t; % kg
        elseif flags.flag_waterquality == 1 % Let's calculate iT
            % Calculate Build-up using C1 and C2
            WQ_States.B_0 = LULC_Properties.C_1.*(1 - exp(1).^(-LULC_Properties.C_2*LULC_Properties.ADD )); % kg/ha
            WQ_States.B_t = WQ_States.B_0*Wshed_Properties.cell_area/10^4.*idx_not_nan; % kg
            Maps.WQ_States.initial_buildup_map = WQ_States.B_t; % kg
        end
    end
end

LULC_Properties.LULC = LULC;

% Soil Indexes
Soil_Properties.n_soil = n_soil;
idx_soil = zeros(size(dem,1),size(dem,2),Soil_Properties.n_soil);
for i = 1:Soil_Properties.n_soil
    index = soil_matrix == SOIL_index(i,1);
    idx_soil(:,:,i) = index;
end

idx_soil = logical(idx_soil);

% ========================================================================
% Soil Parameter Assigning (Updated for new Excel format)
% ========================================================================

Soil_Properties.Soil = SOIL;

for i = 1:Soil_Properties.n_soil

    Soil_Properties.ksat(idx_soil(:,:,i))        = soil_parameters(i,1); % mm/h
    Soil_Properties.n_vg(idx_soil(:,:,i))       = soil_parameters(i,2);
    Soil_Properties.alpha_vg(idx_soil(:,:,i))   = soil_parameters(i,3); % 1/m
    Soil_Properties.theta_sat(idx_soil(:,:,i))  = soil_parameters(i,4);
    Soil_Properties.theta_r(idx_soil(:,:,i))    = soil_parameters(i,5);
    Soil_Properties.theta_i(idx_soil(:,:,i))    = soil_parameters(i,6);
    Soil_Properties.Sy(idx_soil(:,:,i))         = soil_parameters(i,7);
    Soil_Properties.ksat_gw(idx_soil(:,:,i))    = soil_parameters(i,8);
    if size(soil_parameters,2) >= 10 && isfinite(soil_parameters(i,10))
        Soil_Properties.Ks_multiplier_near_surface(idx_soil(:,:,i)) = soil_parameters(i,10);
    else
        Soil_Properties.Ks_multiplier_near_surface(idx_soil(:,:,i)) = 1;
    end
    if size(soil_parameters,2) >= 11 && isfinite(soil_parameters(i,11))
        Soil_Properties.Ks_multiplier_root_zone(idx_soil(:,:,i)) = soil_parameters(i,11);
    else
        Soil_Properties.Ks_multiplier_root_zone(idx_soil(:,:,i)) = 1;
    end
    if size(soil_parameters,2) >= 12 && isfinite(soil_parameters(i,12))
        Soil_Properties.Ks_multiplier_transmission(idx_soil(:,:,i)) = soil_parameters(i,12);
    else
        Soil_Properties.Ks_multiplier_transmission(idx_soil(:,:,i)) = 1;
    end
    Soil_Properties.m_vg = 1 - 1 ./ Soil_Properties.n_vg;

end

% Soil depth (DTB) assignment
if isempty(DTB_raster)
    for i = 1:Soil_Properties.n_soil
        Soil_Properties.Soil_Depth(idx_soil(:,:,i)) = soil_parameters(i,9); % [m] soil depth per soil class
    end
else
    Soil_Properties.Soil_Depth = DTB_raster.Z;
end

LULC_Properties.River_K_coeff = River_K_coeff;

% Mask in Impervious Areas. Cells that are built areas have no infiltration
Soil_Properties.ksat(LULC_Properties.idx_imp) = 0; % Impervious areas
Soil_Properties.soil_matrix = soil_matrix;
Soil_Properties.idx_soil = idx_soil;
% Soil initial storage from Initial_Soil_Moisture_path is handled later
% in the unified "Initial Vadose Zone Volume" block, after zwt_init is known.

% Initial Moisture
Soil_Properties.I_0(LULC_Properties.idx_imp) = 0; % Impervious areas

% Replenishing Coefficient
Soil_Properties.kr = (1/75*(sqrt(Soil_Properties.ksat/25.4))); % Replenishing rate (1/hr) (Check Rossman, pg. 110)
Soil_Properties.Tr = 4.5./sqrt(Soil_Properties.ksat/25.4); %  Recovery time (hr) Check rossman pg. 110
Soil_Properties.Lu = 4.*sqrt(Soil_Properties.ksat/25.4); % Inches - Uppermost layer of the soil
Soil_Properties.Lu = Soil_Properties.Lu*2.54/100; % Meters
Soil_Properties.k_out = (Soil_Properties.theta_sat - Soil_Properties.theta_r).*Soil_Properties.kr.*Soil_Properties.Lu*1000; % Rate of replenishing exfiltration from the saturated zone during recoverying times (mm/hr)
% clear idx % clear this matrices to avoid huge storage data

%%%% New xmax and ymax
xmin = 1; %initial position x in the grid (collums)
ymin = 1; % lines (up to down)
zzz = size(elevation);
xmax = zzz(2);
ymax = zzz(1);

%% Neglecting infiltration in rivers
if flags.flag_neglect_infiltration_river == 1
    LULC_Properties.h_0(idx_facc) = 0; % No initial abstraction in rivers
    Soil_Properties.ksat(idx_facc) = 0; % No infiltration in rivers
end

 %% Native inflow forcing storage (NO discretization to model time-step)
if flags.flag_inflow == 1
    % Keep inflow in its native forcing grid.
    % Inflow_Parameters.inflow_hydrograph_rate is already in mm/h.
    % Inflow_Parameters.time_inflow is in minutes.
    % Inflow_Parameters.time_step_inflow is the native forcing interval in minutes.

    % Store end-stamped forcing times for update_spatial_BC.
    % Assumption:
    %   value at row i is valid from previous forcing time up to time_inflow(i)
    BC_States.time_inflow_forcing = Inflow_Parameters.time_inflow(:);

    % Native rate series [nTimes x nGauges] -> transpose for convenience later if needed
    BC_States.inflow_rate_native = Inflow_Parameters.inflow_hydrograph_rate;
end

%% Native lumped rainfall forcing storage (NO discretization to model time-step)
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~= 1
    % Keep rainfall in its native forcing grid.
    % Rainfall_Parameters.intensity_rainfall is already in mm/h.
    % Rainfall_Parameters.time_rainfall is in minutes.
    BC_States.time_rainfall_forcing = Rainfall_Parameters.time_rainfall(:);
    BC_States.rainfall_rate_native  = Rainfall_Parameters.intensity_rainfall(:);
end
%% Determination of grid parameters and outlet coordinates (Whole Domain)

% Calculates the number of non-inf cells to determine the watershed area
matrix_nan = isnan(elevation);
number_nan = sum(sum(matrix_nan));
Wshed_Properties.drainage_area = sum(sum((~isnan(DEM_raster.Z))))*Wshed_Properties.Resolution^2; % m2
Edges = bwperim(dem>=0,8);
Wshed_Properties.watershed_perimeter = sum(sum(Edges))*Wshed_Properties.Resolution/1000;  % km
Wshed_Properties.impervious_area = impervious_cells * Wshed_Properties.cell_area;
Wshed_Properties.pervious_area   = pervious_cells   * Wshed_Properties.cell_area;
Wshed_Properties.impervious_rate = Wshed_Properties.impervious_area/Wshed_Properties.drainage_area;
Wshed_Properties.pervious_rate = Wshed_Properties.pervious_area / Wshed_Properties.drainage_area;

% Geometrical Properties
if ~exist('S','var') % S does not exists
    FD = FLOWobj(DEM_raster);
    area_km2 = GIS_data.min_area; % km2
    area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
    S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
end
Wshed_Properties.avg_river_length = S.distance(end); % Lenth of the river (m)
Wshed_Properties.compactness_coefficient = 0.28*Wshed_Properties.watershed_perimeter*1000/sqrt(Wshed_Properties.drainage_area);
Wshed_Properties.form_factor = Wshed_Properties.drainage_area/(Wshed_Properties.avg_river_length^2); % A / L^2
Wshed_Properties.circularity_index = 12.57*Wshed_Properties.drainage_area/((Wshed_Properties.watershed_perimeter*1000)^2); % more close to 1, closer to a circle
% Pollutant Mass
if flags.flag_waterquality == 1
    initial_mass = sum(sum(WQ_States.B_t(~isinf(WQ_States.B_t)))); % kg of pollutant
end
%% Inflow forcing metadata for native-time integration
if flags.flag_inflow == 1
    % Cursor starts at first forcing interval
    BC_States.cursor_inflow = 1;

    % Optional convenience field
    BC_States.n_inflow_times = size(BC_States.inflow_rate_native, 1);
end
%% Lumped rainfall forcing metadata for native-time integration
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~= 1
    BC_States.cursor_lumped_rain = 1;
    BC_States.n_rain_times = numel(BC_States.rainfall_rate_native);

    % keep this if used elsewhere
    outlet_runoff_volume = zeros(size(BC_States.rainfall_rate_native));
end
%% Pre allocating more arrays
Total_Inflow = 0;
if flags.flag_warmup ~= 1
    depths.d_t = spatial_domain;
end
outflow_rate.outflow_cms_t = spatial_domain;
% ----------------- Space dependent arrays -----------------

depths.d_p = spatial_domain;
outflow_rate.qin_t = spatial_domain;
% inundated_cells = spatial_domain;
% ----------------- Time dependent arrays -----------------
time_size = length(running_control.time_store);
time_size_2 = saver_memory_maps; % for all matrices in order to save memory space
WQ_States.EMC_outlet = zeros(time_size,1);
Maps.Hydro.d = zeros(ny,nx,time_size_2);
if flags.flag_human_instability ==1
    Maps.Hydro.risk = zeros(ny,nx,time_size_2);
elseif flags.flag_human_instability == 2
    Maps.Hydro.risk_bti = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_fti = zeros(ny,nx,time_size_2);
elseif flags.flag_human_instability == 3
    Maps.Hydro.risk_cm = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_tm = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_am = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_om = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_cf = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_tf = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_af = zeros(ny,nx,time_size_2);
    Maps.Hydro.risk_of = zeros(ny,nx,time_size_2);
end
Maps.Hydro.I_t = zeros(ny,nx,time_size_2);
Maps.Hydro.C = Maps.Hydro.I_t;
Maps.Hydro.f = Maps.Hydro.I_t;
outet_hydrograph = zeros(time_size,1);
time_hydrograph = zeros(time_size,1);
if flags.flag_waterquality == 1
    Maps.WQ_States.Pol_Conc_Map = zeros(ny,nx,time_size_2);
    Maps.WQ_States.Pol_Mass_Map = zeros(ny,nx,time_size_2);
    Maps.WQ_States.Pol_Load_Map = zeros(ny,nx,time_size_2);
    Maps.WQ_States.outet_pollutograph = zeros(time_size,1);
end
%% Clearing a few variables
clear accum_precipitation  precipitation imp time_size idx_soil  Land_Cover_Data Elevation_DATA idx_ idx_1 matrix_nan idx accum_inflow col col_check d_0_imp d_0_per h_0_imp h_0_per I_0_per I_0_imp 

%% ---------------------- Main Routing (GA + 4D/8D CA + BW models) ------------------%%
tic % Start counting time

%%%% Minimum Soil Moisture
eps = 1e-3;
min_soil_moisture = eps * double(Soil_Properties.I_0 > 0);
min_soil_moisture(Soil_Properties.I_0 == 0) = eps;
min_soil_moisture(idx_nan) = nan;

% ----------------- Initialize Variables -----------------
time_calculation_routing = 0;
k = 1; % Start counter of the while loop

recording_parameters.actual_record_state   = 1;
recording_parameters.last_record_maps      = 1;
recording_parameters.last_record_hydrograph = 1;

if flags.flag_satellite_rainfall == 1
    recording_parameters.last_record_maps_rainfall = 1;
end

if flags.flag_inflow == 1
    BC_States.delta_inflow_agg = zeros(Inflow_Parameters.n_stream_gauges,1);
end

% ----------------- Initial simulation clock -----------------
t = 0;   % minutes at initialization

% ----------------- Initialize rainfall forcing -----------------
% Rainfall depth for each step is now computed inside update_spatial_BC
% from the native forcing series over [t_previous, t].
BC_States.delta_p_agg = 0;

% ----------------- Initialize ETP forcing -----------------
BC_States.delta_e_agg  = zeros(size(dem));   % evaporation depth over model step [mm]
BC_States.delta_tr_agg = zeros(size(dem));   % transpiration depth over model step [mm]

% ----------------- Grid for spatial rainfall interpolation -----------------
Spatial_Rainfall_Parameters.x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution * (1:size(DEM_raster.Z,2))'; % Pixel eastings
Spatial_Rainfall_Parameters.y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution * (1:size(DEM_raster.Z,1))'; % Pixel northings

% ======================================================================
% INITIAL SPATIAL RAINFALL FROM GAUGES
% ======================================================================
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && ...
        flags.flag_input_rainfall_map ~= 1 && ...
        flags.flag_real_time_satellite_rainfall ~= 1 && ...
        flags.flag_satellite_rainfall ~= 1

    nsteps_spatial_rainfall = floor(running_control.routing_time / running_control.record_time_spatial_rainfall);

    if nsteps_spatial_rainfall < 1
        error('Please enter a feasible time-step for the aggregation of spatial rainfall.');
    end

    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1), size(dem,2), saver_memory_maps);

    % At initialization use first gauge rainfall snapshot
    z = 1;

    Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges, 1);
    Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges, 2);

    rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(z, 1:Spatial_Rainfall_Parameters.n_raingauges)';

    idx_rainfall = isnan(rainfall);
    rainfall(idx_rainfall) = [];
    Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = [];
    Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = [];

    spatial_rainfall = Rainfall_Interpolator( ...
        Spatial_Rainfall_Parameters.x_coordinate, ...
        Spatial_Rainfall_Parameters.y_coordinate, ...
        rainfall, ...
        Spatial_Rainfall_Parameters.x_grid, ...
        Spatial_Rainfall_Parameters.y_grid);

    spatial_rainfall(idx_nan) = nan;

    % rainfall is in mm/h -> convert to mm over model step
    BC_States.delta_p_agg = spatial_rainfall * (time_step_model / 60);

    Maps.Hydro.spatial_rainfall_maps(:,:,1) = spatial_rainfall;
    BC_States.average_spatial_rainfall(1,1) = mean(mean(spatial_rainfall, 'omitnan'));
end

% ----------------- Initial numerical time step -----------------
time_step = running_control.min_time_step / 60;   % hours? kept as original logic
t = time_step_model;                              % minutes

% ======================================================================
% INITIAL RAINFALL TILE (INPUT MAPS)
% ======================================================================
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && ...
        flags.flag_input_rainfall_map == 1 && ...
        flags.flag_satellite_rainfall == 0 && ...
        flags.flag_real_time_satellite_rainfall == 0

    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1), size(dem,2), saver_memory_maps);

    f0 = string(Input_Rainfall.labels_Directory{1});

    try
        % Geographic raster
        [A, Rr] = readgeoraster(f0, 'CoordinateSystemType', 'geographic');
        Rr.GeographicCRS = geocrs(4326);
        A = raster_cutter(DEM_raster, Rr, A, 1);
    catch
        % Projected raster
        [A, Rr] = readgeoraster(f0);
        A = raster_cutter(DEM_raster, Rr, A, 0);
    end

    input_rainfall = double(A.Z);
    input_rainfall(idx_nan) = nan;

    Maps.Hydro.spatial_rainfall_maps(:,:,1) = input_rainfall;

    % Input rainfall maps are assumed to be in mm/h
    BC_States.delta_p_agg = input_rainfall * (time_step_model / 60);   % mm per model step

    BC_States.average_spatial_rainfall(1,1) = mean(mean(input_rainfall, 'omitnan'));
end

% ======================================================================
% INITIAL ETP TILES (INPUT EVAPORATION / TRANSPIRATION MAPS)
% ======================================================================
if flags.flag_ETP == 1 && flags.flag_input_ETP_map == 1

    Maps.Hydro.spatial_evaporation_maps   = zeros(size(dem,1), size(dem,2), saver_memory_maps);
    Maps.Hydro.spatial_transpiration_maps = zeros(size(dem,1), size(dem,2), saver_memory_maps);

    fe0 = string(Input_Evaporation.labels_Directory{1});
    ft0 = string(Input_Transpiration.labels_Directory{1});

    try
        % Geographic rasters
        [input_evaporation, rR]   = readgeoraster(fe0, 'CoordinateSystemType', 'geographic');
        [input_transpiration, ~]  = readgeoraster(ft0, 'CoordinateSystemType', 'geographic');

        rR.GeographicCRS = geocrs(4326);

        input_evaporation   = raster_cutter(DEM_raster, rR, input_evaporation,   1);
        input_transpiration = raster_cutter(DEM_raster, rR, input_transpiration, 1);

    catch
        % Projected rasters
        [input_evaporation, rR]   = readgeoraster(fe0);
        [input_transpiration, ~]  = readgeoraster(ft0);

        input_evaporation   = raster_cutter(DEM_raster, rR, input_evaporation,   0);
        input_transpiration = raster_cutter(DEM_raster, rR, input_transpiration, 0);
    end

    % Extract matrices
    input_evaporation   = double(input_evaporation.Z);
    input_transpiration = double(input_transpiration.Z);

    % Mask outside valid domain
    input_evaporation(idx_nan)   = nan;
    input_transpiration(idx_nan) = nan;

    % Save first maps
    Maps.Hydro.spatial_evaporation_maps(:,:,1)   = input_evaporation;
    Maps.Hydro.spatial_transpiration_maps(:,:,1) = input_transpiration;

    % IMPORTANT:
    % Input ET component maps are assumed to be in mm/day
    % Convert to mm over one model time step:
    %   mm_per_step = mm_per_day * (time_step_model / 1440)
    BC_States.delta_e_agg  = input_evaporation   * (time_step_model / 1440);
    BC_States.delta_tr_agg = input_transpiration * (time_step_model / 1440);

    BC_States.average_spatial_evaporation(1,1)   = mean(mean(input_evaporation, 'omitnan'));
    BC_States.average_spatial_transpiration(1,1) = mean(mean(input_transpiration, 'omitnan'));
end

% ======================================================================
% INITIAL SATELLITE RAINFALL
% ======================================================================
register = [];

if flags.flag_satellite_rainfall == 1
    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1), size(dem,2), saver_memory_maps);

    register = 0;
    product = 'PDIRNow1hourly';

    [rainfall_raster, register, register_data, ~] = Satellite_rainfall_processing( ...
        [], [], register, product, date_begin, date_end, ...
        flags.flag_satellite_rainfall, flags.flag_real_time_satellite_rainfall, DEM_raster);

    input_rainfall = double(rainfall_raster.Z);
    input_rainfall(idx_nan) = nan;

    Maps.Hydro.spatial_rainfall_maps(:,:,1) = input_rainfall;

    % Satellite rainfall assumed in mm/h
    BC_States.delta_p_agg = input_rainfall * (time_step_model / 60);

    Spatial_Rainfall_Parameters.rainfall_spatial_duration = Input_Rainfall.time;
    BC_States.average_spatial_rainfall(1,1) = mean(mean(input_rainfall, 'omitnan'));
end

% ======================================================================
% INITIAL REAL-TIME SATELLITE RAINFALL
% ======================================================================
if flags.flag_real_time_satellite_rainfall == 1
    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1), size(dem,2), saver_memory_maps);

    register = 0;
    product = 'PDIRNow1hourly';

    [rainfall_raster, register] = Satellite_rainfall_processing( ...
        register, product, date_begin, date_end, ...
        flags.flag_satellite_rainfall, flags.flag_real_time_satellite_rainfall, DEM_raster);

    input_rainfall = double(rainfall_raster.Z);
    input_rainfall(idx_nan) = nan;

    Maps.Hydro.spatial_rainfall_maps(:,:,1) = input_rainfall;

    % Real-time rainfall assumed in mm/h
    BC_States.delta_p_agg = input_rainfall * (time_step_model / 60);

    Spatial_Rainfall_Parameters.rainfall_spatial_duration = Input_Rainfall.time;
    BC_States.average_spatial_rainfall(1,1) = mean(mean(input_rainfall, 'omitnan'));
end

% ======================================================================
% GENERAL INITIALIZATION AFTER FORCING
% ======================================================================
running_control.time_step_model   = time_step_model;
running_control.time_save_previous = 0; % minutes

Previous_Volume = 0;
t_previous      = 0;
step_error      = zeros(1, running_control.steps);

% CA_States.I_cell = zeros(ny,nx,5);

WQ_States.P_conc      = 0;
tmin_wq               = running_control.max_time_step;
dmax_final            = zeros(size(elevation));
WQ_States.mass_lost   = 0;
WQ_States.mass_outlet = 0;
Out_Conc              = 0;
WQ_States.vol_outlet  = 0;
outlet_runoff_volume  = 0;
WQ_States.Tot_Washed  = spatial_domain;

BC_States.outflow_volume = 0;
BC_States.inflow_volume  = 0;

CA_States.I_tot_end_cell = zeros(size(spatial_domain));

Hydro_States.ETP = zeros(size(DEM));
Hydro_States.ETR = zeros(size(DEM));
Hydro_States.S   = zeros(size(DEM));

% ======================================================================
% CELL AREA / SUBGRID EFFECTIVE AREA
% ======================================================================
if flags.flag_subgrid == 1
    if flags.flag_overbanks == 1
        warning('HydroPol2D:LegacySubgridOverbanks', ...
            ['flag_subgrid=1 with flag_overbanks=1 uses the legacy River_Width/River_Depth overbank pathway. ', ...
             'Phase 1 validation covers lookup-table subgrid hydraulics with flag_overbanks=0.']);
    end

    C_a = Wshed_Properties.cell_area * ones(ny, nx);
    C_a(isnan(dem)) = nan;

    if flags.flag_overbanks == 1
        % Legacy overbank pathway: effective area comes from explicit river
        % width/depth rasters. Lookup-table subgrid storage is handled by
        % SubgridTables and keeps C_a at the coarse-cell area.
        index = (depths.d_0 / 1000 < Wshed_Properties.River_Depth) & (Wshed_Properties.River_Width > 0);
        C_a(index) = Wshed_Properties.River_Width(index) * Wshed_Properties.Resolution;
    end

    % Make sure inflow cells use full raster cell area
    if flags.flag_inflow == 1
        C_a(Wshed_Properties.inflow_mask) = Wshed_Properties.cell_area;
        if flags.flag_overbanks == 1
            Wshed_Properties.River_Width(Wshed_Properties.inflow_mask) = 0;
            Wshed_Properties.River_Depth(Wshed_Properties.inflow_mask) = 0;
        end
    end
else
    C_a = Wshed_Properties.cell_area * ones(ny, nx);
    C_a(isnan(dem)) = nan;
end

% ======================================================================
% INFLOW BOUNDARY CONDITION
% ======================================================================
if flags.flag_inflow == 1
    BC_States.inflow = spatial_domain; % spatially distributed inflow matrix
    for i = 1:Inflow_Parameters.n_stream_gauges
        BC_States.inflow = BC_States.inflow + ...
            BC_States.delta_inflow_agg(i,1) * Wshed_Properties.inflow_cells(:,:,i); % mm
    end
else
    BC_States.inflow = 0 * spatial_domain;
end

%%%% ELEVATIONS %%%
%%%% ASSIGNING VALUES %%%

Elevation_Properties.elevation_cell = elevation;

if flags.flag_D8 == 1
    zero_matrix = NaN(size(Elevation_Properties.elevation_cell));

    Elevation_Properties.elevation_left_t  = [zeros(ny,1), elevation(:,1:(nx-1))];
    Elevation_Properties.elevation_right_t = [elevation(:,2:nx), zeros(ny,1)];
    Elevation_Properties.elevation_up_t    = [zeros(1,nx); elevation(1:(ny-1),:)];
    Elevation_Properties.elevation_down_t  = [elevation(2:ny,:); zeros(1,nx)];

    Elevation_Properties.elevation_NE_t = zero_matrix;
    Elevation_Properties.elevation_SE_t = zero_matrix;
    Elevation_Properties.elevation_SW_t = zero_matrix;
    Elevation_Properties.elevation_NW_t = zero_matrix;

    depths.d_NE_t = zero_matrix;
    depths.d_SE_t = zero_matrix;
    depths.d_SW_t = zero_matrix;
    depths.d_NW_t = zero_matrix;

    Elevation_Properties.elevation_NE_t(2:ny,1:(nx-1)) = Elevation_Properties.elevation_cell(1:(ny-1),2:nx);
    Elevation_Properties.elevation_SE_t(1:(ny-1),1:(nx-1)) = Elevation_Properties.elevation_cell(2:ny,2:nx);
    Elevation_Properties.elevation_SW_t(1:(ny-1),2:nx) = Elevation_Properties.elevation_cell(2:ny,1:(nx-1));
    Elevation_Properties.elevation_NW_t(2:ny,2:nx) = Elevation_Properties.elevation_cell(1:(ny-1),1:(nx-1));
end


%% Checking if all cells of the domain have data
Wshed_Properties.rainfall_matrix(idx_nan) = nan;
zzz = Wshed_Properties.rainfall_matrix; zzz = zzz > 0;
idx_cells = logical(zzz);
if min(min(LULC_Properties.roughness(idx_cells))) == 0
    warning('Cells with not associated LULC parameters')
    warning('Assuming n = 0.03 for these areas. Also assuming h0 = 0 for them.')
    idx_not_assigned = LULC_Properties.roughness == 0 & idx_cells == 1;
    LULC_Properties.roughness(idx_not_assigned) = nanmean(nanmean(LULC_Properties.roughness));
    LULC_Properties.h_0(idx_not_assigned) = 0;
    pause(1)
end

if min(min(Soil_Properties.ksat(and(idx_cells, ~LULC_Properties.idx_imp)))) == 0
    warning('Cells with not associated SOIL parameters or K = 0 that are not impervious areas')
    warning('Assuming K = 0 for these areas.')
    idx_not_assigned = Soil_Properties.ksat == 0 & idx_cells == 1 & ~LULC_Properties.idx_imp;
    Soil_Properties.ksat(idx_not_assigned) = nanmean(nanmean(Soil_Properties.ksat(Soil_Properties.ksat>0)));
    Soil_Properties.theta_sat(idx_not_assigned) = nanmean(nanmean(Soil_Properties.theta_sat(Soil_Properties.theta_sat>0)));
    Soil_Properties.theta_r(idx_not_assigned) = nanmean(nanmean(Soil_Properties.theta_r(Soil_Properties.theta_r>0)));
    Soil_Properties.ksat_gw(idx_not_assigned) = nanmean(nanmean(Soil_Properties.ksat_gw(Soil_Properties.ksat_gw>0)));  % Attention here
    Soil_Properties.Sy(idx_not_assigned) = nanmean(nanmean(Soil_Properties.Sy(Soil_Properties.Sy>0)));  % Attention here
    pause(1)
end

% Making sure all GW cells have properties
if nansum(nansum((and(Soil_Properties.Sy == 0, ~isnan(DEM_raster.Z))))) > 0
    idx_gw = (and(Soil_Properties.Sy == 0, ~isnan(DEM_raster.Z)));
    Soil_Properties.ksat_gw(idx_gw) = nanmean(nanmean(Soil_Properties.ksat_gw(Soil_Properties.ksat_gw>0)));  % Attention here
    Soil_Properties.Sy(idx_gw) = nanmean(nanmean(Soil_Properties.Sy(Soil_Properties.Sy>0)));  % Attention here
end


%% CA-8D Matrices
% Slope
if flags.flag_GPU
    dim1 = size(Elevation_Properties.elevation_cell,1);
    dim2 = size(Elevation_Properties.elevation_cell,2);
    if flags.flag_D8 == 1
        dim3 = 9;
        dim3_ = 5;
    else
        dim3 = 5;
        dim3_ = 3; % Now bates has 3 dimensions
    end
    wse_slope_zeros = gpuArray(zeros(dim1,dim2,dim3));
    Distance_Matrix = gpuArray(zeros(dim1,dim2));
    outflow_bates = gpuArray(zeros(dim1,dim2,dim3_));
    Qc = gpuArray(zeros(dim1,dim2,dim3_ - 1));
    Qf = gpuArray(zeros(dim1,dim2,dim3_ - 1));
    Qci = gpuArray(zeros(dim1,dim2,dim3_ - 1));
    Qfi = gpuArray(zeros(dim1,dim2,dim3_ - 1));
else
    dim1 = size(Elevation_Properties.elevation_cell,1);
    dim2 = size(Elevation_Properties.elevation_cell,2);
    if flags.flag_D8 == 1
        dim3 = 9;
        dim3_ = 5;
    else
        dim3 = 5;
        dim3_ = 3;
    end
    wse_slope_zeros = (zeros(dim1,dim2,dim3));
    Distance_Matrix = (zeros(dim1,dim2));
    outflow_bates = (zeros(dim1,dim2,dim3_));
    Qc = (zeros(dim1,dim2,dim3_ - 1));
    Qf = (zeros(dim1,dim2,dim3_ - 1));
    Qci = (zeros(dim1,dim2,dim3_ - 1));
    Qfi = (zeros(dim1,dim2,dim3_ - 1));
end


%% Domain Values for Hydrological Inputs and Kernel Filter
if ~isempty(LAI_raster)
    try
        LAI_raster.Z(idx_nan) = nan;
        kernell = applyKernelFilter(LAI_raster.Z, 3, 'gaussian');
        kernell(isnan(kernell) & ~isnan(LAI_raster.Z)) = LAI_raster.Z(isnan(kernell) & ~isnan(LAI_raster.Z));
        LAI_raster.Z = kernell;
    end
end
if ~isempty(Albedo_raster)
    try
        Albedo_raster.Z(idx_nan) = nan;
        kernell = applyKernelFilter(Albedo_raster.Z, 3, 'gaussian');
        kernell(isnan(kernell) & ~isnan(Albedo_raster.Z)) = Albedo_raster.Z(isnan(kernell) & ~isnan(Albedo_raster.Z));
        Albedo_raster.Z = kernell;
    end
end

if ~isempty(DTB_raster)
    try
        DTB_raster.Z(idx_nan) = nan;
        kernell = applyKernelFilter(DTB_raster.Z, 3, 'gaussian');
        kernell(isnan(kernell) & ~isnan(DTB_raster.Z)) = DTB_raster.Z(isnan(kernell) & ~isnan(DTB_raster.Z));
        DTB_raster.Z = kernell;
    end
end


% Soil Thickness Input
try
    Soil_Properties.Soil_Depth = double(DTB_raster.Z);
    
    % Fill only missing values inside the domain.
    Soil_Properties.Soil_Depth(~idx_nan & isnan(Soil_Properties.Soil_Depth)) = 1;
    
    % Keep a small minimum depth for numerical safety.
    % This preserves Case 2 because Case 2 has DTB = 0.30 m.
    Soil_Properties.Soil_Depth(Soil_Properties.Soil_Depth < 0.1) = 0.1;

    % Assuming the average of soil depth
    avg_soil_depth = nanmean(nanmean(Soil_Properties.Soil_Depth));

    % Activate this to use an average soil depth
    % Soil_Properties.Soil_Depth(Soil_Properties.Soil_Depth>=0) = avg_soil_depth;

    % Assuming a very shallow aquifer in rivers
    if flags.flag_subgrid == 1
        idx_GW_river = Wshed_Properties.River_Width > 0;
    else
        idx_GW_river = idx_rivers;
    end

    % % Activate this to ensure very shallow depth in rivers
    % Soil_Properties.Soil_Depth(idx_GW_river) = 0.05; % m
end

if isempty(DTB_path)
    Soil_Properties.Soil_Depth = elevation ./ elevation; % 1 meter soil
end

C = 0; k = 1; Rainfall_Parameters.index_aggregation = 1;

%% Imposing that we are not doing automatic calibration
flags.flag_automatic_calibration = 0;
% % If you want to save the pre-processing data, use the code below:
% save('HydroPol2D_preprocessing_input_data.mat');

%% Deleting Water Quality States
if flags.flag_waterquality ~= 1
    clear WQ_States % We delete WQ_States if we are not modeling it
    WQ_States = []; % Just adding an empty array to use in other functions
end

%% Calculating River Matrix to Estimate Groundwater
if flags.flag_groundwater_modeling == 1
    FD = FLOWobj(DEM_raster); % Flow direction
    As  = flowacc(FD); % Flow Accumulation
    Wshed_Properties.fac_area = As.Z*(Wshed_Properties.cell_area/1000/1000); % km2
else
    Lateral_Groundwater_Flux = 0; % m3/s/km of river
end
idx_rivers = Wshed_Properties.fac_area >= GIS_data.min_area;  % Logical Matrix with 1 being pixels with rivers
Wshed_Properties.idx_rivers = idx_rivers;

if flags.flag_baseflow == 1
    % Activate this to ensure very shallow depth in rivers
    % Soil_Properties.Soil_Depth(Wshed_Properties.idx_rivers) = 0.5; % m
end

%% Groundwater States

% Bottom of soil column / bedrock elevation [m]
BC_States.z0 = elevation - Soil_Properties.Soil_Depth;
BC_States.z0(idx_nan) = nan;

% Fallback initial groundwater depth above bedrock [m]
h_depth = 0.01;

% Initial groundwater table elevation [m]
% Priority:
%   1) Use GW_table_raster if available after preprocessing/resampling.
%   2) If not available, fall back to z0 + 0.01 m.
if exist('GW_table_raster','var') && ~isempty(GW_table_raster)

    BC_States.h_0 = double(GW_table_raster.Z);

    % Constraint at DEM elevation
    min_zwt = 0.05;                                % Minimum vadose depth [m]
    BC_States.h_0 = min(BC_States.h_0, elevation - min_zwt);

    % Safety check: by this point GW_table_raster should already be aligned
    % to the DEM/model domain.
    if ~isequal(size(BC_States.h_0), size(BC_States.z0))
        error(['GW_table_raster exists, but its size does not match the model domain. ', ...
               'Check raster loading/resampling/alignment before Groundwater States.']);
    end

    % Apply model-domain mask
    BC_States.h_0(idx_nan) = nan;

    % If the raster has missing values inside the valid domain, fill only
    % those cells with the fallback groundwater depth above bedrock.
    idx_missing_gw = ~idx_nan & (isnan(BC_States.h_0) | isinf(BC_States.h_0));

    if any(idx_missing_gw(:))
        BC_States.h_0(idx_missing_gw) = BC_States.z0(idx_missing_gw) + h_depth;

        warning(['GW_table_raster was loaded, but it contains missing/invalid values inside the model domain. ', ...
                 'Those cells were filled using z0 + ', num2str(h_depth), ' m.']);
    else
        disp('Initial groundwater table loaded from GW_table_raster.');
    end

else

    % No GW raster was provided or loaded.
    % Fall back to constant initial groundwater depth above bedrock.
    BC_States.h_0 = BC_States.z0 + h_depth;
    BC_States.h_0(idx_nan) = nan;

    warning(['GW_table_raster not found. Falling back to constant initial groundwater depth of ', ...
             num2str(h_depth), ' m above bedrock.']);

end

BC_States.h_t = BC_States.h_0;

% Example for visual check
% surfmap(BC_States.h_t - elevation);

%% Initial Vadose Zone Volume
% -------------------------------------------------------------------------
% Initial unsaturated thickness from groundwater table
% -------------------------------------------------------------------------
GW_Depth_init = BC_States.h_t - (elevation - Soil_Properties.Soil_Depth);

zwt_init = Soil_Properties.Soil_Depth - GW_Depth_init;
zwt_init = max(zwt_init, 0);
zwt_init(idx_nan) = nan;

% -------------------------------------------------------------------------
% Unified initial soil water storage [mm]
%
% Case A:
%   No warmup soil-storage raster:
%       I_0 = (theta_i - theta_r) * zwt_init * 1000
%
% Case B:
%   flag_warmup = 1 and Initial_Soil_Moisture_path exists:
%       Initial_Soil_Moisture_path is interpreted as I0 [mm]
%       I_0 is read directly from raster
%       theta_i is back-calculated for consistency:
%           theta_i = theta_r + I_0 / (zwt_init * 1000)
%
% IMPORTANT:
%   I_0 is storage ABOVE residual water content [mm].
% -------------------------------------------------------------------------

use_initial_soil_moisture_raster = false;

if flags.flag_warmup == 1 && ...
        exist('Initial_Soil_Moisture_path','var') && ...
        strlength(string(Initial_Soil_Moisture_path)) > 0 && ...
        isfile(char(Initial_Soil_Moisture_path))

    use_initial_soil_moisture_raster = true;
end

if use_initial_soil_moisture_raster

    % ------------------------------------------------------------
    % Warmup raster gives initial soil bucket storage I0 [mm]
    % NOTE:
    % Despite the legacy name Initial_Soil_Moisture_path, the raster values
    % are interpreted here as soil water storage above residual [mm].
    % ------------------------------------------------------------
    try
        Initial_Soil_Moisture_raster = GRIDobj(char(Initial_Soil_Moisture_path));
    
        % Align to DEM grid
        Initial_Soil_Moisture_raster = resample( ...
            Initial_Soil_Moisture_raster, DEM_raster, 'bilinear');
    
        I0_raster_mm = double(Initial_Soil_Moisture_raster.Z);
    
        % Clean invalid cells
        I0_raster_mm(idx_nan) = nan;
        I0_raster_mm(I0_raster_mm < 0) = 0;
    
        % Maximum active storage allowed by theta_sat and theta_r
        I0_max_mm = max(Soil_Properties.theta_sat - Soil_Properties.theta_r, 0) .* ...
            zwt_init .* 1000;
    
        % Do not allow raster storage to exceed physically available pore space
        I0_raster_mm = min(I0_raster_mm, I0_max_mm);
    
        % Assign storage directly
        Soil_Properties.I_0 = I0_raster_mm;
    
        % Back-calculate equivalent theta_i for consistency
        theta_i_from_I0 = Soil_Properties.theta_r + ...
            Soil_Properties.I_0 ./ max(zwt_init .* 1000, eps);
    
        % Clip theta_i physically
        theta_i_from_I0 = max(theta_i_from_I0, Soil_Properties.theta_r);
        theta_i_from_I0 = min(theta_i_from_I0, Soil_Properties.theta_sat);
    
        theta_i_from_I0(idx_nan) = nan;
    
        Soil_Properties.theta_i = theta_i_from_I0;
    
        Soil_Properties.initial_soil_storage_source = 'Initial_Soil_Moisture_path_I0_mm';
    catch
        warning('No initial soil moisture raster found. Assuming initial soil moisture from tables.')
    end

else

    % ------------------------------------------------------------
    % Class-based theta_i from SOIL_table.
    % This may be entered directly, computed field capacity, or manual
    % field capacity depending on input_data_bypass_script settings.
    % ------------------------------------------------------------
    Soil_Properties.I_0 = ...
        max(Soil_Properties.theta_i - Soil_Properties.theta_r, 0) .* ...
        zwt_init .* 1000;

    Soil_Properties.initial_soil_storage_source = 'SOIL_table_theta_i';

end

% Domain mask
Soil_Properties.I_0(idx_nan) = nan;

% Impervious cells have no soil infiltration storage
Soil_Properties.I_0(LULC_Properties.idx_imp) = 0;

% Apply minimum soil moisture only to valid, non-impervious cells
valid_soil_storage = ~idx_nan & ~LULC_Properties.idx_imp & isfinite(Soil_Properties.I_0);

if isscalar(min_soil_moisture)
    Soil_Properties.I_0(valid_soil_storage) = ...
        max(Soil_Properties.I_0(valid_soil_storage), min_soil_moisture);
else
    Soil_Properties.I_0(valid_soil_storage) = ...
        max(Soil_Properties.I_0(valid_soil_storage), min_soil_moisture(valid_soil_storage));
end

% Final states
Soil_Properties.I_p = Soil_Properties.I_0;
Soil_Properties.I_t = Soil_Properties.I_0;

layer_options = struct('near_surface_depth_m', 0.10, 'min_layer_thickness_m', 0.005);
[Soil_Properties, LULC_Properties] = derive_layered_soil_profile( ...
    Soil_Properties, LULC_Properties, BC_States, elevation, idx_nan, layer_options);
Soil_Properties = initialize_layered_soil_storage(Soil_Properties, idx_nan);

GW_Depth_check = BC_States.h_t - (elevation - Soil_Properties.Soil_Depth);
zwt_check = Soil_Properties.Soil_Depth - GW_Depth_check;

fprintf('\n===== SUBSURFACE INITIALIZATION CHECK =====\n');
fprintf('Soil I0 source  = %s\n', Soil_Properties.initial_soil_storage_source);
fprintf('Mean Soil_Depth = %.4f m\n', mean(Soil_Properties.Soil_Depth(:), 'omitnan'));
fprintf('Mean h_t        = %.4f m\n', mean(BC_States.h_t(:), 'omitnan'));
fprintf('Mean elevation  = %.4f m\n', mean(elevation(:), 'omitnan'));
fprintf('Mean theta_i    = %.4f\n', mean(Soil_Properties.theta_i(:), 'omitnan'));
fprintf('Mean I_t        = %.4f mm\n', mean(Soil_Properties.I_t(:), 'omitnan'));
fprintf('==========================================\n');
fprintf('==========================================\n');

%% Characterizing River Roughness
% Inbank Manning (only using the first entry)
if flags.flag_obs_gauges == 1
    Wshed_Properties.Inbank_Manning = LULC_Properties.roughness; % Only first entry used
    try
        Wshed_Properties.Inbank_Manning(idx_rivers) = River_Manning(1); % Only first entry used
    catch
        Wshed_Properties.Inbank_Manning(idx_rivers) = 0.035;
    end
    % Overbank Manning (assuming the same of the LULC) )
    Wshed_Properties.Overbank_Manning = LULC_Properties.roughness; % Can be altered
else
    Wshed_Properties.Inbank_Manning = LULC_Properties.roughness; % Only first entry used
    Wshed_Properties.Inbank_Manning(idx_rivers) = LULC_Parameters.River_Manning; % Only first entry used
    % Overbank Manning (assuming the same of the LULC) )
    Wshed_Properties.Overbank_Manning = LULC_Properties.roughness; % Can be altered
end

LULC_Properties.River_K_coeff = River_K_coeff; % River recharge coefficient

% Enforcing chosen roughness coefficient for river cells
if flags.flag_reduce_DEM == 1 && flags.flag_subgrid == 0
    LULC_Properties.roughness(idx_facc) = LULC_Parameters.River_Manning;
end

%% Subgrid Functions
SubgridTables = [];
if flags.flag_subgrid == 1 && flags.flag_resample == 1
    % DEM Treatment and Filtering Algorithms
    % Fillsinks
    % max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
    % DEM_filled = fillsinks(DEM,max_depth);
    if flags.flag_fill_DEM == 1
        DEM_filled = fillsinks_for_routing(DEM_raster, flags); % Filled DEM
        DIFFDEM = DEM_filled - DEM_raster.Z;
        DIFFDEM.Z(DIFFDEM.Z==0) = nan;
        DEM_raster = DEM_filled;
        % imageschs(DEM_raster,DIFFDEM.Z);
    end

    % DTM Filter
    flags.flag_DTM = 1;
    if flags.flag_DTM == 1
        slope_threshold = GIS_data.slope_DTM; % percentage
        DEM_filtered = DTM_Filter(DEM_raster.Z,Wshed_Properties.Resolution,GIS_data.slope_DTM);
        DEM_raster.Z = DEM_filtered;
    end

    % Smooth DEM
    if flags.flag_smooth_cells == 1
        Vq = imgaussfilt(DEM_raster.Z,'FilterSize',3);
        dem_diff_smooth = Vq;
        dem = Vq; % New dem
        DEM_raster.Z = dem;
        close all
    end

    % Fill Again to Make Sure
    % max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
    % DEM_filled = fillsinks(DEM,max_depth);
    if flags.flag_fill_DEM == 1
        DEM_filled = fillsinks_for_routing(DEM_raster, flags); % Filled DEM
        DIFFDEM = DEM_filled - DEM_raster.Z;
        DIFFDEM.Z(DIFFDEM.Z==0) = nan;
        DEM_raster = DEM_filled;
        % imageschs(DEM_raster,DIFFDEM.Z);
    end
    % Saving high resolution DEM
    DEM_raster_high_resolution = GRIDobj(Subgrid_DEM_path);
    %%
    if flags.flag_overbanks ~= 1
        % Polynomial Subgrid
        % Polynomial order
        % poly_order = 2;
        % SubgridTables = [];
        % [Subgrid_Properties.A_spline, Subgrid_Properties.V_spline, Subgrid_Properties.Rh_east_spline, Subgrid_Properties.Rh_north_spline, Subgrid_Properties.Subgrid_Properties.W_east_spline, Subgrid_Properties.W_north_spline, Subgrid_Properties.A_east_spline, Subgrid_Properties.A_north_spline , Subgrid_Properties.Poly_NSE, Subgrid_Properties.invert_el] = Subgrid_Properties_Function(DEM_raster_high_resolution, DEM_raster, GIS_data.resolution_resample, poly_order);
        % Lookup Subgrid
        Roughness_raster = DEM_raster; % clone structure
        Roughness_raster.Z = LULC_Properties.roughness;
        subgrid_nr_levels = 20;
        subgrid_dz_m = 0.02;
        subgrid_max_depth_m = 10.0;
        if exist('InputData_Bypass', 'var') && isfield(InputData_Bypass, 'subgrid')
            if isfield(InputData_Bypass.subgrid, 'subgrid_nr_levels')
                subgrid_nr_levels = InputData_Bypass.subgrid.subgrid_nr_levels;
            end
            if isfield(InputData_Bypass.subgrid, 'subgrid_dz_m')
                subgrid_dz_m = InputData_Bypass.subgrid.subgrid_dz_m;
            end
            if isfield(InputData_Bypass.subgrid, 'subgrid_max_depth_m')
                subgrid_max_depth_m = InputData_Bypass.subgrid.subgrid_max_depth_m;
            end
        end
        [SubgridTables, Subgrid_Properties.invert_el] = Subgrid_Properties_Lookup( ...
            DEM_raster_high_resolution, ...
            Roughness_raster, ...
            DEM_raster, ...
            Wshed_Properties.Resolution, ...
            'nr_levels', subgrid_nr_levels, ...
            'dz', subgrid_dz_m, ...
            'maxDepth', subgrid_max_depth_m);
        hp2d_validate_subgrid_tables(SubgridTables);
        % Calculate and print size of SubgridTables in MB
        info = whos('SubgridTables');                   % Get memory info for the variable
        size_MB = info.bytes / 1024 / 1024;             % Convert from bytes to megabytes
        fprintf('SubgridTables size: %.2f MB\n', size_MB);

    else
        Subgrid_Properties = [];
    end
elseif flags.flag_subgrid == 1
    if flags.flag_overbanks ~= 1
        % Saving high resolution DEM
        DEM_raster_high_resolution = GRIDobj(Subgrid_DEM_path);
        % Polynomial Subgrid
        % Polynomial order
        % poly_order = 2;
        % SubgridTables = [];
        % [Subgrid_Properties.A_spline, Subgrid_Properties.V_spline, Subgrid_Properties.Rh_east_spline, Subgrid_Properties.Rh_north_spline, Subgrid_Properties.Subgrid_Properties.W_east_spline, Subgrid_Properties.W_north_spline, Subgrid_Properties.A_east_spline, Subgrid_Properties.A_north_spline , Subgrid_Properties.Poly_NSE, Subgrid_Properties.invert_el] = Subgrid_Properties_Function(DEM_raster_high_resolution, DEM_raster, GIS_data.resolution_resample, poly_order);
        % Lookup Subgrid
        Roughness_raster = DEM_raster; % clone structure
        Roughness_raster.Z = LULC_Properties.roughness;  % your Manning map [s/m^(1/3)]
        subgrid_nr_levels = 20;
        subgrid_dz_m = 0.02;
        subgrid_max_depth_m = 10.0;
        if exist('InputData_Bypass', 'var') && isfield(InputData_Bypass, 'subgrid')
            if isfield(InputData_Bypass.subgrid, 'subgrid_nr_levels')
                subgrid_nr_levels = InputData_Bypass.subgrid.subgrid_nr_levels;
            end
            if isfield(InputData_Bypass.subgrid, 'subgrid_dz_m')
                subgrid_dz_m = InputData_Bypass.subgrid.subgrid_dz_m;
            end
            if isfield(InputData_Bypass.subgrid, 'subgrid_max_depth_m')
                subgrid_max_depth_m = InputData_Bypass.subgrid.subgrid_max_depth_m;
            end
        end
        [SubgridTables, Subgrid_Properties.invert_el] = Subgrid_Properties_Lookup( ...
            DEM_raster_high_resolution, ...   % fine DEM
            Roughness_raster, ...             % fine Manning raster
            DEM_raster, ...                   % coarse reference grid
            Wshed_Properties.Resolution, ...  % coarse resolution [m]
            'nr_levels', subgrid_nr_levels, ...
            'dz', subgrid_dz_m, ...
            'maxDepth', subgrid_max_depth_m);
        hp2d_validate_subgrid_tables(SubgridTables);
        % Calculate and print size of SubgridTables in MB
        info = whos('SubgridTables');                   % Get memory info for the variable
        size_MB = info.bytes / 1024 / 1024;             % Convert from bytes to megabytes
        fprintf('SubgridTables size: %.2f MB\n', size_MB);
    end
end

%% Dam Break Parameters
% If dam breaking is simulated
% if flags.flag_dam_break == 1
%     input_table_breaker = readtable('Damns_breaks_points.xlsx');
%     breakers.easting_absolute = table2array(input_table_breaker(:,2)); % Easting Coordinates
%     breakers.northing_absolute = table2array(input_table_breaker(:,3)); % Northing Coordinates
%     breakers.ID = table2array(input_table_breaker(:,1)); % ID
%     % --- Converting coordinates to local coordinates in pixels
%     breakers.easting= round((-GIS_data.xulcorner + breakers.easting_absolute)/Wshed_Properties.Resolution);
%     breakers.northing = round((GIS_data.yulcorner - breakers.northing_absolute)/Wshed_Properties.Resolution);
%     flag_break_1 = 1; flag_break_2 = 1;
% end

% We need to finish this part. Also, we need to fix it in the main while

%% Outlet Perimeter B.C.
% If you want to assign outlet normal flow boundary conditions to all
% domain perimeter, activate the following code
% outlet_index = Wshed_Properties.perimeter;
if flags.flag_boundary == 1
    outlet_index = Wshed_Properties.perimeter; % Perimeter of the domain is the outlet
    % Taking away cells that are not in the perimeter but are inflow cells
    if flags.flag_inflow == 1
        outlet_index(logical(Wshed_Properties.inflow_mask)) = 0;
    end

    idx_outlet = logical(outlet_index);
    [Wshed_Properties.row_outlet, Wshed_Properties.col_outlet] = find(idx_outlet);

    Wshed_Properties.el_outlet = dem;
    Wshed_Properties.el_outlet(~idx_outlet) = nan;
    if any(idx_outlet(:))
        Wshed_Properties.stage_min = min(Wshed_Properties.el_outlet(idx_outlet), [], 'omitnan');
        [Wshed_Properties.row_outlet_stage_min, Wshed_Properties.col_outlet_stage_min] = ...
            find(idx_outlet & dem == Wshed_Properties.stage_min);
    end

    if exist('Paths','var') && isfield(Paths,'Shapes') && exist('DEM_raster','var')
        FileName = fullfile(Paths.Shapes, 'Outlet_Mask.tif');
        raster_to_export = DEM_raster;
        raster_to_export.Z = double(outlet_index);
        try
            raster_to_export.Z(isnan(DEM_raster.Z)) = nan;
        catch
        end
        try
            GRIDobj2geotiff(raster_to_export, FileName);
        catch
            warning('Final perimeter outlet mask not exported.');
        end
    end
end

%% Converting Arrays to Single Precision
if flags.flag_single == 1

    % ------------------------------------------------------------
    % Preserve CPU-only / nonnumeric fields before casting
    % ------------------------------------------------------------
    if ~exist('extra_parameters','var') || ~isstruct(extra_parameters)
        extra_parameters = struct();
    end

    % -------------------------------
    % Gauges labels MUST stay on CPU
    % -------------------------------
    if isfield(flags,'flag_obs_gauges') && flags.flag_obs_gauges == 1
        if ~isfield(extra_parameters,'gauges') || ~isstruct(extra_parameters.gauges)
            extra_parameters.gauges = struct();
        end

        if exist('gauges','var') && isstruct(gauges)
            if isfield(gauges,'labels_observed_string')
                extra_parameters.gauges.labels_observed_string = gauges.labels_observed_string;
                gauges.labels_observed_string = [];
            else
                extra_parameters.gauges.labels_observed_string = {};
            end
        end
    end

    % -------------------------------
    % ETP datetime fields stay on CPU
    % -------------------------------
    if isfield(flags,'flag_ETP') && flags.flag_ETP == 1
        if ~isfield(extra_parameters,'ETP') || ~isstruct(extra_parameters_ETP)
            extra_parameters_ETP = struct();
        end

        if exist('ETP_Parameters','var') && isstruct(ETP_Parameters)
            if isfield(ETP_Parameters,'time_ETP_begin')
                extra_parameters_ETP.time_ETP_begin = ETP_Parameters.time_ETP_begin;
                ETP_Parameters.time_ETP_begin = [];
            end
            if isfield(ETP_Parameters,'time_ETP')
                extra_parameters_ETP.time_ETP = ETP_Parameters.time_ETP;
                ETP_Parameters.time_ETP = [];
            end
        end
    end

    % -------------------------------
    % Human instability text fields stay on CPU
    % -------------------------------
    if isfield(flags,'flag_human_instability') && flags.flag_human_instability > 0
        if exist('Human_Instability','var') && isstruct(Human_Instability)
            if ~isfield(extra_parameters,'Human_Instability') || ~isstruct(extra_parameters.Human_Instability)
                extra_parameters.Human_Instability = struct();
            end

            if isfield(Human_Instability,'list')
                extra_parameters.Human_Instability.list = Human_Instability.list;
                Human_Instability.list = [];
            end
            if isfield(Human_Instability,'names')
                extra_parameters.Human_Instability.names = Human_Instability.names;
                Human_Instability.names = [];
            end
        end
    end

    % ------------------------------------------------------------
    % Convert struct arrays to single
    % ------------------------------------------------------------
    BC_States = structfun(@single, BC_States, 'UniformOutput', false);
    CA_States = structfun(@single, CA_States, 'UniformOutput', false);
    Courant_Parameters = structfun(@single, Courant_Parameters, 'UniformOutput', false);
    depths = structfun(@single, depths, 'UniformOutput', false);
    Elevation_Properties = structfun(@single, Elevation_Properties, 'UniformOutput', false);

    if exist('Subgrid_Properties','var') && isstruct(Subgrid_Properties) && ~isempty(Subgrid_Properties)
        Subgrid_Properties = structfun(@single, Subgrid_Properties, 'UniformOutput', false);
    end

    flags = structfun(@single, flags, 'UniformOutput', false);
    GIS_data = structfun(@single, GIS_data, 'UniformOutput', false);
    Hydro_States = structfun(@single, Hydro_States, 'UniformOutput', false);
    Inflow_Parameters = structfun(@single, Inflow_Parameters, 'UniformOutput', false);
    LULC_Properties = structfun(@single, LULC_Properties, 'UniformOutput', false);
    Rainfall_Parameters = structfun(@single, Rainfall_Parameters, 'UniformOutput', false);
    recording_parameters = structfun(@single, recording_parameters, 'UniformOutput', false);
    running_control = structfun(@single, running_control, 'UniformOutput', false);
    Soil_Properties = convert_struct_fields_recursive(Soil_Properties, @single);
    Wshed_Properties = structfun(@single, Wshed_Properties, 'UniformOutput', false);

    if exist('gauges','var') && isstruct(gauges)
        gauges = structfun(@single, gauges, 'UniformOutput', false);
    end

    % if exist('Maps','var') && isstruct(Maps) && ~isempty(Maps)
    %     Maps = structfun(@single, Maps, 'UniformOutput', false);
    % end

    if isfield(flags,'flag_waterquality') && flags.flag_waterquality == 1
        WQ_States = structfun(@single, WQ_States, 'UniformOutput', false);
    end

    if isfield(flags,'flag_reservoir') && flags.flag_reservoir == 1
        Reservoir_Data = structfun(@single, Reservoir_Data, 'UniformOutput', false);
    end

    if isfield(flags,'flag_spatial_rainfall') && flags.flag_spatial_rainfall == 1
        Spatial_Rainfall_Parameters = structfun(@single, Spatial_Rainfall_Parameters, 'UniformOutput', false);
    end

    if isfield(flags,'flag_snow_modeling') && flags.flag_snow_modeling == 1
        Snow_Properties = structfun(@single, Snow_Properties, 'UniformOutput', false);
    end

%     if isfield(flags,'flag_ETP') && flags.flag_ETP == 1
%         ETP_Parameters = structfun(@single, ETP_Parameters, 'UniformOutput', false);
%     end

    if isfield(flags,'flag_ETP') && flags.flag_ETP == 1
        fn = fieldnames(ETP_Parameters);
        for k = 1:numel(fn)
            v = ETP_Parameters.(fn{k});

            if isnumeric(v) || islogical(v)
                ETP_Parameters.(fn{k}) = gpuArray(v);
            else
                fprintf('ETP_Parameters.%s left on CPU (class: %s)\n', fn{k}, class(v));
            end
        end
    end

    if isfield(flags,'flag_ETP') && flags.flag_ETP == 1 && ...
       isfield(flags,'flag_input_ETP_map') && flags.flag_input_ETP_map == 1
        if exist('Spatial_ETP_Parameters','var') && isstruct(Spatial_ETP_Parameters) && ~isempty(Spatial_ETP_Parameters)
            fn = fieldnames(Spatial_ETP_Parameters);
            for k = 1:numel(fn)
                v = Spatial_ETP_Parameters.(fn{k});

                if isnumeric(v) || islogical(v)
                    Spatial_ETP_Parameters.(fn{k}) = gpuArray(v);
                else
                    fprintf('Spatial_ETP_Parameters.%s left on CPU (class: %s)\n', fn{k}, class(v));
                end
            end
        end
    end

    if isfield(flags,'flag_human_instability') && flags.flag_human_instability > 0
        Human_Instability.slope = arcslope(DEM_raster,'degree');
        Human_Instability.slope = Human_Instability.slope.Z;
        Human_Instability = structfun(@single, Human_Instability, 'UniformOutput', false);
    end

    % ------------------------------------------------------------
    % Restore logical fields that should never stay single
    % ------------------------------------------------------------
    if isfield(LULC_Properties,'idx_lulc')
        LULC_Properties.idx_lulc = logical(LULC_Properties.idx_lulc);
    end
    if isfield(LULC_Properties,'idx_imp')
        LULC_Properties.idx_imp = logical(LULC_Properties.idx_imp);
    end
    if isfield(Soil_Properties,'idx_soil')
        Soil_Properties.idx_soil = logical(Soil_Properties.idx_soil);
    end
    if exist('idx_nan','var')
        idx_nan = logical(idx_nan);
    end
    if exist('idx_nan_5','var')
        if isfield(flags,'flag_waterquality') && flags.flag_waterquality == 1
            idx_nan_5 = logical(idx_nan_5);
        else
            idx_nan_5 = [];
        end
    end
    if exist('idx_outlet','var')
        idx_outlet = logical(idx_outlet);
    end

    % ------------------------------------------------------------
    % Convert standalone arrays to single
    % ------------------------------------------------------------
    if flags.flag_subgrid == 1
    elevation = single(elevation);

    fn = fieldnames(SubgridTables);
    for k = 1:numel(fn)
        v = SubgridTables.(fn{k});
        if isnumeric(v)
            SubgridTables.(fn{k}) = single(v);
        end
    end
    end
    if exist('C_a','var') && ~isempty(C_a)
        C_a = single(C_a);
    end
    if exist('wse_slope_zeros','var') && ~isempty(wse_slope_zeros)
        wse_slope_zeros = single(wse_slope_zeros);
    end
    if exist('Distance_Matrix','var') && ~isempty(Distance_Matrix)
        Distance_Matrix = single(Distance_Matrix);
    end
    if exist('outflow_bates','var') && ~isempty(outflow_bates)
        outflow_bates = single(outflow_bates);
    end
    if exist('Qc','var') && ~isempty(Qc)
        Qc = single(Qc);
    end
    if exist('Qf','var') && ~isempty(Qf)
        Qf = single(Qf);
    end
    if exist('Qci','var') && ~isempty(Qci)
        Qci = single(Qci);
    end
    if exist('Qfi','var') && ~isempty(Qfi)
        Qfi = single(Qfi);
    end
    if exist('rainfall_spatial_aggregation','var') && ~isempty(rainfall_spatial_aggregation)
        rainfall_spatial_aggregation = single(rainfall_spatial_aggregation);
    end
    if exist('etp_spatial_aggregation','var') && ~isempty(etp_spatial_aggregation)
        etp_spatial_aggregation = single(etp_spatial_aggregation);
    end
    if exist('input_evaporation','var') && ~isempty(input_evaporation)
        input_evaporation = single(input_evaporation);
    end
    if exist('input_transpiration','var') && ~isempty(input_transpiration)
        input_transpiration = single(input_transpiration);
    end
    if exist('register','var') && isnumeric(register) && ~isempty(register)
        register = single(register);
    end

    k = single(k);
    nx = single(nx);
    ny = single(ny);

    if exist('Out_Conc','var') && ~isempty(Out_Conc)
        Out_Conc = single(Out_Conc);
    end
    if exist('outlet_index','var') && ~isempty(outlet_index)
        outlet_index = single(outlet_index);
    end
    if exist('outlet_runoff_volume','var') && ~isempty(outlet_runoff_volume)
        outlet_runoff_volume = single(outlet_runoff_volume);
    end
    if exist('outlet_type','var') && ~isempty(outlet_type)
        outlet_type = single(outlet_type);
    end
    if exist('slope_outlet','var') && ~isempty(slope_outlet)
        slope_outlet = single(slope_outlet);
    end
    if exist('spatial_domain','var') && ~isempty(spatial_domain)
        spatial_domain = single(spatial_domain);
    end
    if exist('t','var') && ~isempty(t)
        t = single(t);
    end
    if exist('t_previous','var') && ~isempty(t_previous)
        t_previous = single(t_previous);
    end
    if exist('time_calculation_routing','var') && ~isempty(time_calculation_routing)
        time_calculation_routing = single(time_calculation_routing);
    end
    if exist('time_step','var') && ~isempty(time_step)
        time_step = single(time_step);
    end
    if exist('time_step_model','var') && ~isempty(time_step_model)
        time_step_model = single(time_step_model);
    end
    if exist('tmin_wq','var') && ~isempty(tmin_wq)
        tmin_wq = single(tmin_wq);
    end
    if exist('C','var') && ~isempty(C)
        C = single(C);
    end
    if exist('min_soil_moisture','var') && ~isempty(min_soil_moisture)
        min_soil_moisture = single(min_soil_moisture);
    end

    % Force k back to scalar single 1
    k = single(1);

end

%% Plotting Input Rasters
try
    if strcmpi(strtrim(getenv('HYDROPOL2D_SKIP_INITIAL_PLOTS')), '1')
        warning('HydroPol2D:InitialMapPlotSkipped', ...
            'Initial input-map plotting skipped by HYDROPOL2D_SKIP_INITIAL_PLOTS=1.');
    else
        Plot_Initial_Maps; % Script to plot initial maps
    end
catch plot_initial_maps_error
    warning('HydroPol2D:InitialMapPlotSkipped', ...
        ['Initial input-map plotting failed and was skipped. ', ...
        'Hydraulic preprocessing will continue. Original error: %s'], ...
        plot_initial_maps_error.message);
end

%% Manifest
% WRITE MODEL PARAMETER MANIFEST
% =========================
try
    run_manifest = struct();

    % -------------------------
    % Run identity
    % -------------------------
    run_manifest.run_name = "hydropol2d_preprocessing";

    if use_inputpaths_bypass == 1
        run_manifest.mode = "bypass";
    else
        run_manifest.mode = "excel";
    end

    % -------------------------
    % Real-world timing
    % -------------------------
    run_manifest.time = struct();
    
    if ~exist('run_start_datetime','var') || isempty(run_start_datetime)
        run_start_datetime = datetime('now');
    end
    
    run_manifest.time.run_start = char(datetime(run_start_datetime, ...
        'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
    run_manifest.time.run_start_datetime = run_start_datetime;

    if exist('date_begin','var')
        try
            run_manifest.time.date_begin = char(datetime(date_begin, ...
                'Format','yyyy-MM-dd HH:mm:ss'));
        catch
            run_manifest.time.date_begin = date_begin;
        end
    else
        run_manifest.time.date_begin = "";
    end

    if exist('date_end','var')
        try
            run_manifest.time.date_end = char(datetime(date_end, ...
                'Format','yyyy-MM-dd HH:mm:ss'));
        catch
            run_manifest.time.date_end = date_end;
        end
    else
        run_manifest.time.date_end = "";
    end

    try
        db = datetime(date_begin);
        de = datetime(date_end);
        run_manifest.time.simulation_duration_hours = hours(de - db);
    catch
        run_manifest.time.simulation_duration_hours = [];
    end

    % -------------------------
    % Output / environment
    % -------------------------
    run_manifest.results_dir = results_dir;
    run_manifest.temp_dir = temp_dir;
    run_manifest.folderName = folderName;
    run_manifest.folderName_2 = folderName_2;
    run_manifest.saver_memory_maps = saver_memory_maps;
    run_manifest.use_inputpaths_bypass = use_inputpaths_bypass;

    if exist('model_folder','var')
        run_manifest.model_folder = model_folder;
    else
        run_manifest.model_folder = "";
    end

    % -------------------------
    % Core tool paths
    % -------------------------
    run_manifest.topo_path = topo_path;
    run_manifest.hydropol2d_tools = hydropol2d_tools;

    % -------------------------
    % Resolved input paths
    % -------------------------
    run_manifest.paths = struct();

    % Required rasters
    run_manifest.paths.DEM_path  = DEM_path;
    run_manifest.paths.LULC_path = LULC_path;
    run_manifest.paths.SOIL_path = SOIL_path;

    % Optional rasters
    run_manifest.paths.Warmup_Depth_path          = Warmup_Depth_path;
    run_manifest.paths.Initial_Buildup_path       = Initial_Buildup_path;
    run_manifest.paths.Initial_Soil_Moisture_path = Initial_Soil_Moisture_path;
    run_manifest.paths.Albedo_path                = Albedo_path;
    run_manifest.paths.LAI_path                   = LAI_path;
    run_manifest.paths.DTB_path                   = DTB_path;
    run_manifest.paths.B1_path                    = B1_path;
    run_manifest.paths.B2_path                    = B2_path;
    run_manifest.paths.W1_path                    = W1_path;
    run_manifest.paths.W2_path                    = W2_path;
    run_manifest.paths.Subgrid_DEM_path           = Subgrid_DEM_path;
    run_manifest.paths.RiverWidths_path           = RiverWidths_path;
    run_manifest.paths.RiverDepths_path           = RiverDepths_path;

    % -------------------------
    % Flags
    % -------------------------
    if exist('flags','var') && isstruct(flags)
        run_manifest.flags = flags;
    else
        run_manifest.flags = struct();
    end

    % -------------------------
    % Input source snapshots
    % -------------------------
    if exist('InputPaths','var') && isstruct(InputPaths)
        run_manifest.InputPaths = InputPaths;
    end

    if exist('InputData_Bypass','var') && isstruct(InputData_Bypass)
        run_manifest.InputData_Bypass_fields = fieldnames(InputData_Bypass);
    end

    % -------------------------
    % Raster summary
    % -------------------------
    run_manifest.rasters = struct();
    run_manifest.rasters.DEM_size  = size(DEM_raster.Z);
    run_manifest.rasters.LULC_size = size(LULC_raster.Z);
    run_manifest.rasters.SOIL_size = size(SOIL_raster.Z);

    run_manifest.rasters.has_DTB    = ~isempty(DTB_raster);
    run_manifest.rasters.has_LAI    = ~isempty(LAI_raster);
    run_manifest.rasters.has_Albedo = ~isempty(Albedo_raster);
    run_manifest.rasters.has_widths = ~isempty(widths_raster);
    run_manifest.rasters.has_depths = ~isempty(depths_raster);

    try
        run_manifest.rasters.cellsize = DEM_raster.cellsize;
    catch
        run_manifest.rasters.cellsize = [];
    end

    % -------------------------
    % GIS / model metadata if available
    % -------------------------
    if exist('GIS_data','var')
        run_manifest.GIS_data = GIS_data;
    end

    % -------------------------
    % MATLAB metadata
    % -------------------------
    run_manifest.matlab = struct();
    run_manifest.matlab.version = version;
    run_manifest.matlab.release = version('-release');
    run_manifest.matlab.computer = computer;

    % -------------------------
    % Optional git hash
    % -------------------------
    try
        [status, git_hash] = system('git rev-parse HEAD');
        if status == 0
            run_manifest.git_commit = strtrim(git_hash);
        else
            run_manifest.git_commit = '';
        end
    catch
        run_manifest.git_commit = '';
    end

    % -------------------------
    % Write files
    % -------------------------
    manifest_paths_pre = write_run_manifest(results_dir, ...
        run_manifest, "hydropol2d_preprocessing");

catch ME
    warning('Could not write preprocessing manifest: %s')
end

%% Clearing Variables

% clearvars  -except register saver_memory_maps idx_rivers rainfall_spatial_aggregation model_folder Input_Rainfall Reservoir_Data wse_slope_zeros Distance_Matrix depths Maps Spatial_Rainfall_Parameters GIS_data Inflow_Parameters ETP_Parameters Rainfall_Parameters CA_States BC_States Wshed_Properties Wshed_Properties Human_Instability gauges Hydro_States recording_parameters Courant_Parameters running_control Elevation_Properties inflow_volume idx_outlet outflow_volume outlet_runoff_volume I_t num_obs_gauges drainage_area northing_obs_gauges easting_obs_gauges depths time_record_hydrograph last_record_hydrograph initial_mass delta_p WQ_States routing_time flags LULC_Properties Soil_Properties topo_path idx_lulc idx_imp idx_soil d steps 	alfa_albedo_input 	alfa_max 	alfa_min 	alfa_save 	avgtemp_stations 	B_t   	C  	Cd 	cell_area 	climatologic_spatial_duration 	col_outlet 	coordinate_x 	coordinate_y 	coordinates_stations d_t  d_p 	date_begin  date_end	delta_p_agg  	DEM_etp 	DEM_raster 	depth_tolerance 	elevation    	ETP 	ETP_save 	factor_cells		flow_tolerance	flows_cells	G_stations	gravity	I_tot_end_cell	idx_nan	idx_nan_5	inflow	inflow_cells	k	k_out	Krs	ksat_fulldomain	last_record_maps	lat	mass_lost	mass_outlet	running_control.max_time_step	maxtemp_stations	min_time_step	mintemp_stations	mu	Inflow_Parameters.n_stream_gauges	nx	ny	Out_Conc	outlet_index	outlet_index_fulldomain	outlet_type	P_conc	psi_fulldomain	rainfall_matrix	rainfall_matrix_full_domain	Resolution	ro_water	roughness	roughness_fulldomain	row_outlet	slope_alfa	slope_outlet	spatial_domain	t	t_previous	theta_r_fulldomain	theta_sat	theta_sat_fulldomain	time_calculation_routing	time_change_matrices	time_change_records	time_deltap	time_ETP	time_records	time_save_previous	time_step	time_step_change	time_step_increments	time_step_model	time_step_save	tmin_wq	Tot_Washed	Tr	u2_stations	ur_stations	v_threshold	vel_down	vel_left	vel_right	vel_up	vol_outlet	weight_person	width1_person	width2_person
clearvars  -except run_start_datetime hydropol2d_tools temp_dir folderName_2 folderName run_start_str results_dir use_inputpaths_bypass use_inputdata_bypass SubgridTables extra_parameters_ETP export_root_dir enable_logging input_excel_file input_sheets_folder add_input_sheets_to_path clean_output_folder run_postprocessing Paths InputPaths ValidationCase Cases RunSummary SummaryRow icase case_id case_dir case_root static_dir full_root summary_dir config_dir output_root input_data_bypass_script_path DEM_raster_high_resolution LAI_raster NDVI_raster DTB_raster Albedo_raster C_a input_evaporation input_transpiration Stage_Parameters Qc Qf Qci Qfi register saver_memory_maps extra_parameters Lateral_Groundwater_Flux idx_rivers min_soil_moisture rainfall_spatial_aggregation  model_folder Input_Rainfall Input_Evaporation Input_Transpiration Reservoir_Data wse_slope_zeros Distance_Matrix depths Maps Spatial_Rainfall_Parameters Spatial_ETP_Parameters GIS_data Inflow_Parameters Snow_Properties ETP_Parameters Rainfall_Parameters CA_States BC_States Wshed_Properties Wshed_Properties Human_Instability gauges Hydro_States recording_parameters Courant_Parameters running_control Subgrid_Properties Elevation_Properties inflow_volume idx_outlet outflow_volume outlet_runoff_volume I_t num_obs_gauges drainage_area northing_obs_gauges easting_obs_gauges depths time_record_hydrograph last_record_hydrograph initial_mass delta_p WQ_States routing_time flags LULC_Properties Soil_Properties topo_path idx_lulc idx_imp idx_soil d steps 	alfa_albedo_input 	alfa_max 	alfa_min 	alfa_save 	avgtemp_stations 	B_t   	C  	Cd 	cell_area 	climatologic_spatial_duration 	col_outlet 	coordinate_x 	coordinate_y 	coordinates_stations d_t  d_p 	date_begin  date_end	delta_p_agg  delta_e_agg delta_tr_agg	DEM_etp 	DEM_raster 	depth_tolerance 	elevation    	ETP 	ETP_save 	factor_cells		flow_tolerance	flows_cells	G_stations	gravity	I_tot_end_cell	idx_nan	idx_nan_5	inflow	inflow_cells	k	k_out	Krs	ksat_fulldomain	last_record_maps	lat	mass_lost	mass_outlet	running_control.max_time_step	maxtemp_stations	min_time_step	mintemp_stations	mu	Inflow_Parameters.n_stream_gauges	nx	ny	Out_Conc	outlet_index	outlet_index_fulldomain	outlet_type	P_conc	psi_fulldomain	rainfall_matrix	rainfall_matrix_full_domain	Resolution	ro_water	roughness	roughness_fulldomain	row_outlet	slope_alfa	slope_outlet	spatial_domain	t	t_previous	theta_r_fulldomain	theta_sat	theta_sat_fulldomain	time_calculation_routing	time_change_matrices	time_change_records	time_deltap	time_ETP	time_records	time_save_previous	time_step	time_step_change	time_step_increments	time_step_model	time_step_save	tmin_wq	Tot_Washed	Tr	u2_stations	ur_stations	v_threshold	vel_down	vel_left	vel_right	vel_up	vol_outlet	weight_person	width1_person	width2_person outflow_bates


%% Converting Arrays to GPU Arrays, if required

% DEM info
if flags.flag_ETP == 1
    ETP_info = ETP_Parameters.SpatialRef;
    ETP_Parameters.info = [];
end

if flags.flag_GPU == 1

    % ------------------------------------------------------------
    % Make sure CPU stash exists for nonnumeric/text fields
    % ------------------------------------------------------------
    if ~exist('extra_parameters','var') || ~isstruct(extra_parameters)
        extra_parameters = struct();
    end

    % ------------------------------------------------------------
    % Preserve CPU-only / nonnumeric fields before gpuArray casting
    % ------------------------------------------------------------

    % Gauges labels
    if isfield(flags,'flag_obs_gauges') && flags.flag_obs_gauges == 1
        if ~isfield(extra_parameters,'gauges') || ~isstruct(extra_parameters.gauges)
            extra_parameters.gauges = struct();
        end

        if exist('gauges','var') && isstruct(gauges)
            if isfield(gauges,'labels_observed_string')
                extra_parameters.gauges.labels_observed_string = gauges.labels_observed_string;
                gauges.labels_observed_string = [];
            else
                extra_parameters.gauges.labels_observed_string = {};
            end
        end
    end

    % ETP datetime fields
    if isfield(flags,'flag_ETP') && flags.flag_ETP == 1

        if exist('ETP_Parameters','var') && isstruct(ETP_Parameters)
            if isfield(ETP_Parameters,'time_ETP_begin') && ~isempty(ETP_Parameters.time_ETP_begin)
                extra_parameters_ETP.time_ETP_begin = ETP_Parameters.time_ETP_begin;
                ETP_Parameters.time_ETP_begin = [];
            end
            if isfield(ETP_Parameters,'time_ETP') && ~isempty(ETP_Parameters.time_ETP)
                extra_parameters_ETP.time_ETP = ETP_Parameters.time_ETP;
                ETP_Parameters.time_ETP = [];
            end
        end
    end

    % Human instability text fields
    if isfield(flags,'flag_human_instability') && flags.flag_human_instability > 0
        if exist('Human_Instability','var') && isstruct(Human_Instability)
            if ~isfield(extra_parameters,'Human_Instability') || ~isstruct(extra_parameters.Human_Instability)
                extra_parameters.Human_Instability = struct();
            end

            if isfield(Human_Instability,'list')
                extra_parameters.Human_Instability.list = Human_Instability.list;
                Human_Instability.list = [];
            end
            if isfield(Human_Instability,'names')
                extra_parameters.Human_Instability.names = Human_Instability.names;
                Human_Instability.names = [];
            end
        end
    end

    % ------------------------------------------------------------
    % Convert struct arrays to GPU
    % ------------------------------------------------------------
    BC_States = structfun(@gpuArray, BC_States, 'UniformOutput', false);
    CA_States = structfun(@gpuArray, CA_States, 'UniformOutput', false);
    Courant_Parameters = structfun(@gpuArray, Courant_Parameters, 'UniformOutput', false);
    depths = structfun(@gpuArray, depths, 'UniformOutput', false);
    Elevation_Properties = structfun(@gpuArray, Elevation_Properties, 'UniformOutput', false);

    if exist('Subgrid_Properties','var') && isstruct(Subgrid_Properties) && ~isempty(Subgrid_Properties)
        Subgrid_Properties = structfun(@gpuArray, Subgrid_Properties, 'UniformOutput', false);
    end

    flags = structfun(@gpuArray, flags, 'UniformOutput', false);
    GIS_data = structfun(@gpuArray, GIS_data, 'UniformOutput', false);
    Hydro_States = structfun(@gpuArray, Hydro_States, 'UniformOutput', false);
    Inflow_Parameters = structfun(@gpuArray, Inflow_Parameters, 'UniformOutput', false);
    LULC_Properties = structfun(@gpuArray, LULC_Properties, 'UniformOutput', false);
    Rainfall_Parameters = structfun(@gpuArray, Rainfall_Parameters, 'UniformOutput', false);
    recording_parameters = structfun(@gpuArray, recording_parameters, 'UniformOutput', false);
    running_control = structfun(@gpuArray, running_control, 'UniformOutput', false);
    Soil_Properties = convert_struct_fields_recursive(Soil_Properties, @gpuArray);
    Wshed_Properties = structfun(@gpuArray, Wshed_Properties, 'UniformOutput', false);

    if exist('gauges','var') && isstruct(gauges)
        gauges = structfun(@gpuArray, gauges, 'UniformOutput', false);
    end

%     if exist('Maps','var') && isstruct(Maps) && ~isempty(Maps)
%         Maps = structfun(@gpuArray, Maps, 'UniformOutput', false);
%     end

    if isfield(flags,'flag_waterquality') && flags.flag_waterquality == 1
        WQ_States = structfun(@gpuArray, WQ_States, 'UniformOutput', false);
    else
        WQ_States = [];
    end

    if isfield(flags,'flag_reservoir') && flags.flag_reservoir == 1
        Reservoir_Data = structfun(@gpuArray, Reservoir_Data, 'UniformOutput', false);
    end

    if isfield(flags,'flag_spatial_rainfall') && flags.flag_spatial_rainfall == 1
        Spatial_Rainfall_Parameters = structfun(@gpuArray, Spatial_Rainfall_Parameters, 'UniformOutput', false);
    end

    if isfield(flags,'flag_snow_modeling') && flags.flag_snow_modeling == 1
        Snow_Properties = structfun(@gpuArray, Snow_Properties, 'UniformOutput', false);
    end

    fn = fieldnames(ETP_Parameters);
    for k = 1:numel(fn)
        v = ETP_Parameters.(fn{k});
    
        if isnumeric(v) || islogical(v)
            ETP_Parameters.(fn{k}) = gpuArray(v);
        else
            fprintf('ETP_Parameters.%s left on CPU (class: %s)\n', fn{k}, class(v));
        end
    end

    if isfield(flags,'flag_ETP') && flags.flag_ETP == 1 && ...
       isfield(flags,'flag_input_ETP_map') && flags.flag_input_ETP_map == 1
        if exist('Spatial_ETP_Parameters','var') && isstruct(Spatial_ETP_Parameters) && ~isempty(Spatial_ETP_Parameters)
            Spatial_ETP_Parameters = structfun(@gpuArray, Spatial_ETP_Parameters, 'UniformOutput', false);
        end
    end

    if isfield(flags,'flag_human_instability') && flags.flag_human_instability > 0
        Human_Instability.slope = arcslope(DEM_raster,'degree');
        Human_Instability.slope = Human_Instability.slope.Z;
        Human_Instability = structfun(@gpuArray, Human_Instability, 'UniformOutput', false);
    end

    % ------------------------------------------------------------
    % Restore logical fields after structfun(@gpuArray,...)
    % ------------------------------------------------------------
    if isfield(LULC_Properties,'idx_lulc')
        LULC_Properties.idx_lulc = logical(LULC_Properties.idx_lulc);
    end
    if isfield(LULC_Properties,'idx_imp')
        LULC_Properties.idx_imp = logical(LULC_Properties.idx_imp);
    end
    if isfield(Soil_Properties,'idx_soil')
        Soil_Properties.idx_soil = logical(Soil_Properties.idx_soil);
    end

    % ------------------------------------------------------------
    % Convert standalone arrays to GPU
    % ------------------------------------------------------------
    if exist('SubgridTables','var') && ~isempty(SubgridTables)
        if isnumeric(SubgridTables) || islogical(SubgridTables)
            SubgridTables = gpuArray(SubgridTables);
        elseif isstruct(SubgridTables)
            fn = fieldnames(SubgridTables);
            for k = 1:numel(fn)
                v = SubgridTables.(fn{k});
                if isnumeric(v) || islogical(v)
                    SubgridTables.(fn{k}) = gpuArray(v);
                end
            end
        end
    end
    if exist('C_a','var') && ~isempty(C_a)
        C_a = gpuArray(C_a);
    end
    if exist('wse_slope_zeros','var') && ~isempty(wse_slope_zeros)
        wse_slope_zeros = gpuArray(wse_slope_zeros);
    end
    if exist('Distance_Matrix','var') && ~isempty(Distance_Matrix)
        Distance_Matrix = gpuArray(Distance_Matrix);
    end
    if exist('outflow_bates','var') && ~isempty(outflow_bates)
        outflow_bates = gpuArray(outflow_bates);
    end
    if exist('Qc','var') && ~isempty(Qc)
        Qc = gpuArray(Qc);
    end
    if exist('Qf','var') && ~isempty(Qf)
        Qf = gpuArray(Qf);
    end
    if exist('Qci','var') && ~isempty(Qci)
        Qci = gpuArray(Qci);
    end
    if exist('Qfi','var') && ~isempty(Qfi)
        Qfi = gpuArray(Qfi);
    end
    if exist('rainfall_spatial_aggregation','var') && ~isempty(rainfall_spatial_aggregation)
        rainfall_spatial_aggregation = gpuArray(rainfall_spatial_aggregation);
    end
    if exist('etp_spatial_aggregation','var') && ~isempty(etp_spatial_aggregation)
        etp_spatial_aggregation = gpuArray(etp_spatial_aggregation);
    end
    if exist('input_evaporation','var') && ~isempty(input_evaporation)
        input_evaporation = gpuArray(input_evaporation);
    end
    if exist('input_transpiration','var') && ~isempty(input_transpiration)
        input_transpiration = gpuArray(input_transpiration);
    end
    if exist('register','var') && isnumeric(register) && ~isempty(register)
        register = gpuArray(register);
    end

    elevation = gpuArray(elevation);
    idx_nan = gpuArray(idx_nan);

    if isfield(flags,'flag_waterquality') && flags.flag_waterquality == 1
        idx_nan_5 = gpuArray(idx_nan_5);
    end

    idx_outlet = gpuArray(idx_outlet);
    k = gpuArray(k);
    nx = gpuArray(nx);
    ny = gpuArray(ny);

    if exist('Out_Conc','var') && ~isempty(Out_Conc)
        Out_Conc = gpuArray(Out_Conc);
    end
    if exist('outlet_index','var') && ~isempty(outlet_index)
        outlet_index = gpuArray(outlet_index);
    end
    if exist('outlet_runoff_volume','var') && ~isempty(outlet_runoff_volume)
        outlet_runoff_volume = gpuArray(outlet_runoff_volume);
    end
    if exist('outlet_type','var') && ~isempty(outlet_type)
        outlet_type = gpuArray(outlet_type);
    end
    if exist('slope_outlet','var') && ~isempty(slope_outlet)
        slope_outlet = gpuArray(slope_outlet);
    end
    if exist('spatial_domain','var') && ~isempty(spatial_domain)
        spatial_domain = gpuArray(spatial_domain);
    end
    if exist('t','var') && ~isempty(t)
        t = gpuArray(t);
    end
    if exist('t_previous','var') && ~isempty(t_previous)
        t_previous = gpuArray(t_previous);
    end
    if exist('time_calculation_routing','var') && ~isempty(time_calculation_routing)
        time_calculation_routing = gpuArray(time_calculation_routing);
    end
    if exist('time_step','var') && ~isempty(time_step)
        time_step = gpuArray(time_step);
    end
    if exist('time_step_model','var') && ~isempty(time_step_model)
        time_step_model = gpuArray(time_step_model);
    end
    if exist('tmin_wq','var') && ~isempty(tmin_wq)
        tmin_wq = gpuArray(tmin_wq);
    end
    if exist('C','var') && ~isempty(C)
        C = gpuArray(C);
    end
    if exist('min_soil_moisture','var') && ~isempty(min_soil_moisture)
        min_soil_moisture = gpuArray(min_soil_moisture);
    end

    % Force k back to scalar gpuArray(1)
    k = gpuArray(1);

    % ------------------------------------------------------------
    % Restore CPU-only text fields after numeric GPU casting
    % ------------------------------------------------------------
    if isfield(flags,'flag_obs_gauges') && gather(flags.flag_obs_gauges) == 1
        gNames = fieldnames(gauges);
        for ii = 1:numel(gNames)
            v = gauges.(gNames{ii});
            if isa(v,'gpuArray')
                gauges.(gNames{ii}) = gather(v);
            end
        end

        if exist('extra_parameters','var') && isfield(extra_parameters,'gauges') && ...
                isfield(extra_parameters.gauges,'labels_observed_string')
            gauges.labels_observed_string = extra_parameters.gauges.labels_observed_string;
        else
            gauges.labels_observed_string = {};
        end
    end

    if isfield(flags,'flag_human_instability') && gather(flags.flag_human_instability) > 0
        if exist('extra_parameters','var') && isfield(extra_parameters,'Human_Instability')
            if isfield(extra_parameters.Human_Instability,'list')
                Human_Instability_text.list = extra_parameters.Human_Instability.list;
            end
            if isfield(extra_parameters.Human_Instability,'names')
                Human_Instability_text.names = extra_parameters.Human_Instability.names;
            end
        end
    end
end

clear spatial_domain

%%

function T = read_block_table(GD, blockName, requiredHeaders)

    % Find the anchor cell
    [r0,c0] = find(strcmp(string(GD), blockName), 1, 'first');
    if isempty(r0)
        error('General_Data block "%s" was not found.', blockName);
    end

    % Candidate header row 1: same row as anchor
    hdr1 = string(GD(r0, c0:min(size(GD,2), c0 + numel(requiredHeaders)-1)));

    % Candidate header row 2: next row after anchor
    if r0 + 1 <= size(GD,1)
        hdr2 = string(GD(r0+1, c0:min(size(GD,2), c0 + numel(requiredHeaders)-1)));
    else
        hdr2 = strings(1,numel(requiredHeaders));
    end

    % Decide which row is the true header row
    if all(strcmp(strtrim(hdr1), requiredHeaders))
        hr = r0;
    elseif all(strcmp(strtrim(hdr2), requiredHeaders))
        hr = r0 + 1;
    else
        error('General_Data block missing required headers for "%s". Expected: %s', ...
              blockName, strjoin(requiredHeaders, ', '));
    end

    % Read data rows until blank first column
    data = {};
    rr = hr + 1;
    while rr <= size(GD,1)
        firstCell = GD{rr,c0};
        if (isstring(firstCell) || ischar(firstCell)) && strlength(strtrim(string(firstCell))) == 0
            break
        elseif isempty(firstCell) || (isnumeric(firstCell) && isnan(firstCell))
            break
        end

        data(end+1,1:numel(requiredHeaders)) = GD(rr, c0:c0+numel(requiredHeaders)-1); %#ok<AGROW>
        rr = rr + 1;
    end

    T = cell2table(data, 'VariableNames', matlab.lang.makeValidName(requiredHeaders));
    T.Properties.VariableNames = requiredHeaders;
end

function [rr,cc] = find_cell(GD, key)
S = strings(size(GD));
for r = 1:size(GD,1)
    for c = 1:size(GD,2)
        if ischar(GD{r,c}) || isstring(GD{r,c})
            S(r,c) = string(GD{r,c});
        end
    end
end
[rr,cc] = find(strcmpi(strtrim(S), key), 1, 'first');
if isempty(rr)
    error("General_Data: key/header '%s' not found.", key);
end
end



function v = xlget(GD, key)
%XLGET  Find label in a readcell() grid and return value to the right.
% Case-insensitive. Errors if key not found.

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
    error("General_Data: key '%s' found in last column; no value to the right.", key);
end

v = GD{rr, cc+1};
end

function p = xlpath(GD, key, required)
%XLPATH  Read a path from the General_Data cell array by key.
% If missing and not required -> returns [].
%
% Robust to values coming from readcell/readtable:
%   - cell scalars, char, string scalars/arrays, missing, etc.
% Prevents && / || errors by ensuring scalar logical tests.

if nargin < 3
    required = false;
end

% --- Your existing lookup (keep as-is) ---
v = xlget(GD, key);

% ------------------------------------------------------------
% 1) Normalize v so all tests are scalar-safe
% ------------------------------------------------------------
% If xlget returns a cell (common with readcell), extract first element.
if iscell(v)
    if isempty(v)
        v = [];
    else
        v = v{1};
    end
end

% If xlget returns a string array, force scalar.
if isstring(v) && ~isscalar(v)
    v = v(1);
end

% If xlget returns a char array, keep it as-is (char is fine).

% ------------------------------------------------------------
% 2) Safe empty / missing checks (no non-scalar logicals)
% ------------------------------------------------------------
isEmpty = isempty(v) || ...
          (ischar(v)   && numel(strtrim(v)) == 0) || ...
          (isstring(v) && isscalar(v) && strlength(strtrim(v)) == 0);

% ismissing can return non-scalar if s is not scalar; we guard with isscalar
isMiss  = (isstring(v) && isscalar(v) && ismissing(v));

if isEmpty || isMiss
    if required
        error('Missing required path for "%s" in General_Data.', key);
    else
        p = [];
        return
    end
end

% ------------------------------------------------------------
% 3) Convert to char path (robust)
% ------------------------------------------------------------
% If it's already char, just trim. If it's numeric, try to stringify.
if ischar(v)
    s = string(v);
elseif isstring(v)
    s = v; % already scalar (forced above)
else
    % numeric / other types: attempt conversion
    try
        s = string(v);
    catch
        s = "";
    end
end

% Scalar-safe blank/missing check on s
if ~(isstring(s) && isscalar(s))
    s = string(s);
    if ~isscalar(s), s = s(1); end
end

if ismissing(s) || strlength(strtrim(s)) == 0
    if required
        error('Missing required path for "%s" in General_Data.', key);
    else
        p = [];
        return
    end
end

p = char(strtrim(s));   % return as char (GRIDobj likes char better)

% ------------------------------------------------------------
% 4) Validate only if required
% ------------------------------------------------------------
if required && ~(isfile(p) || isfolder(p))
    error('Path for "%s" does not exist: %s', key, p);
end

end

function s = toStringScalar(v)
% Convert whatever readcell gave you into a scalar string safely.

if iscell(v)
    if isempty(v)
        s = string(missing);
        return
    end
    v = v{1};  % take first cell
end

if isstring(v)
    if isempty(v)
        s = string(missing);
    else
        s = v(1); % force scalar
    end
    return
end

if ischar(v)
    s = string(v);
    return
end

if isnumeric(v)
    if isempty(v) || all(isnan(v(:)))
        s = string(missing);
    else
        s = string(v(1));
    end
    return
end

% fallback
try
    s = string(v);
    if numel(s) > 1, s = s(1); end
catch
    s = string(missing);
end
end

function R = read_raster_or_empty(pathStr)

R = [];

try
    if nargin == 0
        return
    end

    % Handle cell/string/char
    if iscell(pathStr)
        if isempty(pathStr) || isempty(pathStr{1})
            return
        end
        pathStr = pathStr{1};
    end

    pathStr = string(pathStr);

    if ismissing(pathStr) || strlength(strtrim(pathStr)) == 0
        return
    end

    if ~isfile(pathStr)
        return
    end

    % ===== YOUR READER HERE =====
    % You are using TopoToolbox:
    R = GRIDobj(char(pathStr));

catch
    R = [];
end

end

function flags = read_flags_table(FlagsGrid)
%READ_FLAGS_SHEET Builds flags struct by scanning for any cell 'flag_*'
% and reading the value in the cell to the right.
%
% Supports:
% - Flags sheet is a readcell() grid
% - flag names anywhere in the grid (case-insensitive)
% - value is assumed to be in the cell immediately to the right
% - values can be numeric, logical, char, string, empty, missing

    flags = struct();

    nR = size(FlagsGrid,1);
    nC = size(FlagsGrid,2);

    for r = 1:nR
        for c = 1:nC

            % ---- read potential key cell safely ----
            keyRaw = FlagsGrid{r,c};

            % unwrap cell-in-cell
            if iscell(keyRaw) && numel(keyRaw) >= 1
                keyRaw = keyRaw{1};
            end

            % convert to string for checking (safe)
            try
                keyStr = string(keyRaw);
            catch
                continue
            end

            if numel(keyStr) == 0
                continue
            end
            keyStr = keyStr(1);                 % scalar
            keyStr = strtrim(keyStr);

            if strlength(keyStr) == 0
                continue
            end

            % only accept keys like "flag_*"
            if ~startsWith(lower(keyStr), "flag_")
                continue
            end

            % ---- value is the cell to the right ----
            if c == nC
                % flag in last column: no value cell
                valRaw = [];
            else
                valRaw = FlagsGrid{r,c+1};
            end

            % unwrap cell-in-cell
            if iscell(valRaw) && numel(valRaw) >= 1
                valRaw = valRaw{1};
            end

            % ---- normalize value to numeric (0/1) or NaN ----
            val = NaN;

            if isempty(valRaw)
                val = NaN;

            elseif ismissing(string(valRaw))  % handles <missing>, NaT, etc
                val = NaN;

            elseif islogical(valRaw)
                val = double(valRaw);

            elseif isnumeric(valRaw)
                if isscalar(valRaw)
                    val = double(valRaw);
                else
                    % if someone pasted a vector by mistake, take first
                    val = double(valRaw(1));
                end

            else
                % char/string like '1', '0', 'true', 'false', 'yes', 'no'
                s = strtrim(string(valRaw));
                if strlength(s) == 0 || ismissing(s)
                    val = NaN;
                else
                    sLow = lower(s);

                    if any(sLow == ["true","t","yes","y","on"])
                        val = 1;
                    elseif any(sLow == ["false","f","no","n","off"])
                        val = 0;
                    else
                        tmp = str2double(sLow);
                        if ~isnan(tmp)
                            val = tmp;
                        else
                            val = NaN;
                        end
                    end
                end
            end

            % store field (make valid)
            fieldName = matlab.lang.makeValidName(char(keyStr));
            flags.(fieldName) = val;

        end
    end

    % Optional default (prevents crash if user forgot it)
    if ~isfield(flags,'flag_warmup')
        flags.flag_warmup = 0;
    end

    % Read GIS_data.resolution_resample
    
end


% ========================================================================
% HELPERS (put at the END of HydroPol2D_preprocessing.m)
% ========================================================================

function out = xlgetstr(GD, key, defaultValue)
%XLGETSTR  Read a value from General_Data by key and return it as a string scalar.
% Robust to key being anywhere in GD (not just col 1), and to values being
% char/string/numeric/missing. If not found or empty -> defaultValue.
%
% Usage:
%   topo_path = xlgetstr(GD,'topo_path',"");
%
    if nargin < 3, defaultValue = ""; end

    v = xlget(GD, key);  % may return [] if not found
    if isempty(v)
        out = string(defaultValue);
        return;
    end

    % Convert to scalar string safely
    out = toStringScalar(v);

    % Treat "missing-like" content as empty
    if strlength(strtrim(out)) == 0 || any(strcmpi(out, ["nan","none","null","na"]))
        out = string(defaultValue);
    end
end

function G = ensure_projected_crs_from_geotiff(G, geotiff_path)
%ENSURE_PROJECTED_CRS_FROM_GEOTIFF  Attach ProjectedCRS to a TopoToolbox GRIDobj if missing.
%
%   G = ensure_projected_crs_from_geotiff(G, geotiff_path)
%
% If G.georef.SpatialRef.ProjectedCRS is empty/missing, read the GeoTIFF
% spatial reference and copy its ProjectedCRS into G.

    if nargin < 2 || isempty(geotiff_path)
        return
    end

    try
        hasCRS = isfield(G, 'georef') && isfield(G.georef, 'SpatialRef') && ...
                 isprop(G.georef.SpatialRef, "ProjectedCRS") && ~isempty(G.georef.SpatialRef.ProjectedCRS);
        if hasCRS
            return
        end
    catch
        % If any weird structure issue, we'll try to set it anyway
    end

    % Read CRS from the file
    try
        [~, R] = readgeoraster(geotiff_path);
    catch ME
        error("Failed to read GeoTIFF '%s' to obtain CRS: %s", geotiff_path, ME.message);
    end

    % Only projected CRS makes sense for HydroPol2D grids
    if ~isprop(R, "ProjectedCRS") || isempty(R.ProjectedCRS)
        error("GeoTIFF '%s' does not have a ProjectedCRS. Reproject the raster to a projected CRS (UTM, etc.).", geotiff_path);
    end

    % Attach into GRIDobj
    try
        G.georef.SpatialRef.ProjectedCRS = R.ProjectedCRS;
    catch
        % In case georef/SR structure is missing, try to preserve existing and set field
        try
            SR = G.georef.SpatialRef;
            SR.ProjectedCRS = R.ProjectedCRS;
            G.georef.SpatialRef = SR;
        catch ME
            error("Could not attach ProjectedCRS to GRIDobj: %s", ME.message);
        end
    end
end

function [DEM, Rdem, DEM_CRS] = ensure_projected_crs_like(DEM, crs_input)
%ENSURE_PROJECTED_CRS_LIKE_GRIDOBJ
% Ensures a TopoToolbox GRIDobj has a Mapping Toolbox SpatialRef with ProjectedCRS.
%
% INPUTS
%   DEM       : TopoToolbox GRIDobj (must have .refmat and .size)
%   crs_input : (optional) EPSG code (e.g., 32643) OR projcrs object OR proj string.
%               If DEM already has a ProjectedCRS attached, this can be [].
%
% OUTPUTS
%   DEM     : updated GRIDobj with DEM.georef.SpatialRef
%   Rdem    : map.rasterref.MapCellsReference built from DEM.refmat
%   DEM_CRS : projcrs object (ProjectedCRS)
%
% NOTES
% - CRS cannot be inferred from refmat alone. If DEM doesn't already contain one,
%   you must provide crs_input (EPSG/proj).
%
% Maria — adapted for your DEM GRIDobj structure.

    arguments
        DEM
        crs_input = []
    end

    %-----------------------------
    % 1) Build a SpatialRef (Rdem) from refmat + size
    %-----------------------------
    if ~isfield(DEM,'refmat') || isempty(DEM.refmat)
        error('GRIDobj DEM has no refmat. Cannot build SpatialRef.');
    end

    sz = DEM.size;  % [nRows nCols]
    if numel(sz) ~= 2
        error('DEM.size must be [nRows nCols].');
    end

    % Convert referencing matrix -> MapCellsReference
    % This is Mapping Toolbox:
    try
        Rdem = refmatToMapRasterReference(DEM.refmat, sz);
    catch ME
        error(['Failed to convert DEM.refmat to SpatialRef. ' ...
               'Do you have Mapping Toolbox? Original error:\n%s'], ME.message);
    end

    %-----------------------------
    % 2) Detect existing ProjectedCRS (if already attached)
    %-----------------------------
    DEM_CRS = [];
    if isfield(DEM,'georef') && isstruct(DEM.georef) && isfield(DEM.georef,'SpatialRef')
        SR = DEM.georef.SpatialRef;
        try
            if isprop(SR,'ProjectedCRS') && ~isempty(SR.ProjectedCRS)
                DEM_CRS = SR.ProjectedCRS;
            end
        catch
            % ignore (some SR might not support ProjectedCRS)
        end
    end

    %-----------------------------
    % 3) If no CRS attached, create one from crs_input
    %-----------------------------
    if isempty(DEM_CRS)
        if isempty(crs_input)
            error(['DEM has no attached ProjectedCRS. ' ...
                   'Provide crs_input (EPSG like 326xx, a projcrs object, or a proj string).']);
        end

        DEM_CRS = local_make_projcrs(crs_input);
    end

    %-----------------------------
    % 4) Attach CRS to Rdem and store in DEM.georef
    %-----------------------------
    % Attach projected CRS to Rdem (MapCellsReference supports this)
    try
        Rdem.ProjectedCRS = DEM_CRS;
    catch ME
        error('Could not assign ProjectedCRS to MapCellsReference: %s', ME.message);
    end

    if ~isfield(DEM,'georef') || ~isstruct(DEM.georef) || isempty(DEM.georef)
        DEM.georef = struct();
    end
    DEM.georef.SpatialRef = Rdem;

end


function CRS = local_make_projcrs(crs_input)
% Create a projcrs object from EPSG / proj string / projcrs itself.

    if isa(crs_input,'projcrs')
        CRS = crs_input;
        return;
    end

    if isnumeric(crs_input) && isscalar(crs_input)
        CRS = projcrs(crs_input);  % EPSG
        return;
    end

    if ischar(crs_input) || (isstring(crs_input) && isscalar(crs_input))
        s = char(crs_input);

        % Try EPSG in string form e.g., "EPSG:32643"
        tok = regexp(s,'EPSG\:(\d+)','tokens','once');
        if ~isempty(tok)
            CRS = projcrs(str2double(tok{1}));
            return;
        end

        % Otherwise assume it is a WKT / PROJ string if your MATLAB supports it
        % Mapping Toolbox supports projcrs from WKT in recent versions:
        try
            CRS = projcrs(s);
            return;
        catch
            error(['crs_input string not recognized as EPSG:#### or valid projcrs input.\n' ...
                   'Provide numeric EPSG (recommended).']);
        end
    end

    error('Unsupported crs_input type. Use numeric EPSG, projcrs, or "EPSG:####".');
end

% ======================================================================
% Robust getters/setters
% ======================================================================
function R = getSpatialRefRobust(obj)
    % If user passes the ref directly
    if isa(obj,"map.rasterref.MapCellsReference") || isa(obj,"map.rasterref.MapPostingsReference") || ...
       isa(obj,"map.rasterref.GeographicCellsReference") || isa(obj,"map.rasterref.GeographicPostingsReference")
        R = obj; return;
    end

    if ~isstruct(obj), R = []; return; end

    candidates = { ...
        {'georef','SpatialRef'}, ...
        {'georef','R'}, ...
        {'R'}, ...
        {'Ref'}, ...
        {'SpatialRef'}, ...
        {'spatialRef'}, ...
        {'rasterRef'}, ...
        {'RasterRef'} ...
    };

    R = [];
    for k = 1:numel(candidates)
        path = candidates{k};
        v = getNested(obj, path);
        if isa(v,"map.rasterref.MapCellsReference") || isa(v,"map.rasterref.MapPostingsReference") || ...
           isa(v,"map.rasterref.GeographicCellsReference") || isa(v,"map.rasterref.GeographicPostingsReference")
            R = v; return;
        end
    end
end

function [Z, R] = getRasterAndRefRobust(r)
    % Accept either a struct (HydroPol2D style) or a numeric+ref pair (not here)
    if ~isstruct(r)
        error("in_raster must be a struct with data + SpatialRef (e.g., .Z + .R/.Ref/.georef.SpatialRef).");
    end

    % Data
    if isfield(r,'Z')
        Z = r.Z;
    elseif isfield(r,'A')
        Z = r.A;
    elseif isfield(r,'Data')
        Z = r.Data;
    else
        error("in_raster must contain raster data in .Z (or .A/.Data).");
    end

    % Ref
    R = getSpatialRefRobust(r);
    if isempty(R)
        error("in_raster does not contain a recognizable SpatialRef (e.g., .R, .Ref, .SpatialRef, .georef.SpatialRef).");
    end
end

function r = setRasterZRobust(r, Znew)
    if isfield(r,'Z')
        r.Z = Znew;
    elseif isfield(r,'A')
        r.A = Znew;
    elseif isfield(r,'Data')
        r.Data = Znew;
    else
        r.Z = Znew; % fallback
    end
end

function r = setSpatialRefLikeDEM(r, DEM_raster, Rdem)
    % Write output ref in the same "style" as DEM_raster uses (if DEM is struct),
    % otherwise default to r.georef.SpatialRef.
    if ~isstruct(DEM_raster)
        % DEM passed as ref object
        if isfield(r,'georef') && isstruct(r.georef)
            r.georef.SpatialRef = Rdem;
        else
            r.georef = struct('SpatialRef', Rdem);
        end
        return;
    end

    % Prefer to mirror DEM’s storage path
    if isfield(DEM_raster,'georef') && isstruct(DEM_raster.georef)
        if isfield(DEM_raster.georef,'SpatialRef')
            if ~isfield(r,'georef') || ~isstruct(r.georef), r.georef = struct(); end
            r.georef.SpatialRef = Rdem; return;
        elseif isfield(DEM_raster.georef,'R')
            if ~isfield(r,'georef') || ~isstruct(r.georef), r.georef = struct(); end
            r.georef.R = Rdem; return;
        end
    end

    if isfield(DEM_raster,'R')
        r.R = Rdem; return;
    elseif isfield(DEM_raster,'Ref')
        r.Ref = Rdem; return;
    elseif isfield(DEM_raster,'SpatialRef')
        r.SpatialRef = Rdem; return;
    elseif isfield(DEM_raster,'spatialRef')
        r.spatialRef = Rdem; return;
    end

    % fallback
    if ~isfield(r,'georef') || ~isstruct(r.georef), r.georef = struct(); end
    r.georef.SpatialRef = Rdem;
end

function v = getNested(s, path)
    v = [];
    try
        tmp = s;
        for i = 1:numel(path)
            key = path{i};
            if ~isstruct(tmp) || ~isfield(tmp, key)
                v = []; return;
            end
            tmp = tmp.(key);
        end
        v = tmp;
    catch
        v = [];
    end
end

% =======================================================================
% Helper: robust world grid creation from raster reference
% =======================================================================
function [X, Y] = worldGridFromRef(R)
    if exist("worldGrid","file") == 2
        [X, Y] = worldGrid(R);
        return;
    end

    nRows = R.RasterSize(1);
    nCols = R.RasterSize(2);

    dx = diff(R.XWorldLimits) / nCols;
    dy = diff(R.YWorldLimits) / nRows;

    xCenters = R.XWorldLimits(1) + dx*(0.5 : 1 : nCols-0.5);
    yTop     = R.YWorldLimits(2) - dy*0.5;
    yCenters = yTop - dy*(0 : 1 : nRows-1);

    [X, Y] = meshgrid(xCenters, yCenters);
end

function [LULC_raster, SOIL_raster, DTB_raster, LAI_raster, Albedo_raster, widths_raster, depths_raster] = ...
    align_all_to_dem(DEM_raster, LULC_raster, SOIL_raster, DTB_raster, LAI_raster, Albedo_raster, widths_raster, depths_raster)
%ALIGN_ALL_TO_DEM  Resample all rasters to DEM grid if needed.
% Simple version based on your old preprocessing logic.

    % -----------------------------
    % LULC
    % -----------------------------
    if ~isempty(LULC_raster)
        if any(size(LULC_raster.Z) ~= size(DEM_raster.Z))
            LULC_raster = resample(LULC_raster, DEM_raster, 'nearest');
            LULC_raster.Z = round(LULC_raster.Z);
        end
    end

    % -----------------------------
    % SOIL
    % -----------------------------
    if ~isempty(SOIL_raster)
        if any(size(SOIL_raster.Z) ~= size(DEM_raster.Z))
            SOIL_raster = resample(SOIL_raster, DEM_raster, 'nearest');
            SOIL_raster.Z = round(SOIL_raster.Z);
        end
    end

    % -----------------------------
    % DTB
    % -----------------------------
    if ~isempty(DTB_raster)
        if any(size(DTB_raster.Z) ~= size(DEM_raster.Z))
            DTB_raster = resample(DTB_raster, DEM_raster, 'bilinear');
        end
    end

    % -----------------------------
    % LAI
    % -----------------------------
    if ~isempty(LAI_raster)
        if any(size(LAI_raster.Z) ~= size(DEM_raster.Z))
            LAI_raster = resample(LAI_raster, DEM_raster, 'bilinear');
        end
    end

    % -----------------------------
    % Albedo
    % -----------------------------
    if ~isempty(Albedo_raster)
        if any(size(Albedo_raster.Z) ~= size(DEM_raster.Z))
            Albedo_raster = resample(Albedo_raster, DEM_raster, 'bilinear');
        end
    end

    % -----------------------------
    % River widths
    % -----------------------------
    if ~isempty(widths_raster)
        if any(size(widths_raster.Z) ~= size(DEM_raster.Z))
            widths_raster = resample(widths_raster, DEM_raster, 'nearest');
        end
    end

    % -----------------------------
    % River depths
    % -----------------------------
    if ~isempty(depths_raster)
        if any(size(depths_raster.Z) ~= size(DEM_raster.Z))
            depths_raster = resample(depths_raster, DEM_raster, 'nearest');
        end
    end
end

function x = xlnum(GD, key)
%XLNUM numeric version of xlget
    v = xlget(GD, key);
    if isempty(v) || (isstring(v) && strlength(v)==0)
        x = NaN;
        return
    end
    x = double(v);
end

function S = cast_struct_numeric_to_single(S)
    fn = fieldnames(S);
    for i = 1:numel(fn)
        f = fn{i};
        v = S.(f);

        if isnumeric(v) || islogical(v)
            S.(f) = single(v);
        elseif isstruct(v)
            S.(f) = cast_struct_numeric_to_single(v);
        else
            % leave chars, strings, datetime, cells, objects unchanged
        end
    end
end

function S = convert_struct_fields_recursive(S, converter)
    fn = fieldnames(S);
    for i = 1:numel(fn)
        f = fn{i};
        v = S.(f);

        if isnumeric(v) || islogical(v)
            S.(f) = converter(v);
        elseif isstruct(v)
            S.(f) = convert_struct_fields_recursive(v, converter);
        else
            % leave chars, strings, datetime, cells, objects unchanged
        end
    end
end

function val = get_inputpaths_field(InputPaths, fieldName, defaultValue)

    if nargin < 3
        defaultValue = '';
    end

    if isstruct(InputPaths) && isfield(InputPaths, fieldName)
        v = InputPaths.(fieldName);

        if isstring(v)
            if isscalar(v)
                val = char(v);
            else
                val = char(v(1));
            end
        elseif ischar(v)
            val = v;
        else
            try
                val = char(string(v));
            catch
                val = defaultValue;
            end
        end

        if isempty(val)
            val = defaultValue;
        end
    else
        val = defaultValue;
    end
end
