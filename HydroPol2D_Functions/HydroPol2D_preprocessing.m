
% % Creating Modeling Results Folder
folderName = 'Modeling_Results';
folderName_2 = 'Temporary_Files';
test_dir = dir;
% Check if the folder already exists]
% if exist(folderName, 'dir') == 0
try
    % If it doesn't exist, create the folder
    mkdir(folderName);
    mkdir(folderName_2);
    disp('Folder "Modeling_Results" created successfully!');
catch ME
    disp('Impossible to create the folder for some reason');
end

% variable to specify the size of matrices to store out of memory
saver_memory_maps = 12;

%input_table = readtable('General_Data_HydroPol2D.xlsx');
% Load TopoToolBox Tools
topo_path = string(table2cell(input_table(1,31)));
addpath(genpath(char(topo_path)));

% Load Model Functions
HydroPol2D_tools = char(table2cell(input_table(9,31)));
addpath(genpath(char(HydroPol2D_tools)));

% Input Data Paths
DEM_path = char(table2cell((input_table(3,31))));
LULC_path = char((table2cell(input_table(4,31))));
SOIL_path = char(table2cell(input_table(5,31)));
Warmup_Depth_path = char(table2cell(input_table(6,31)));
Initial_Buildup_path = char(table2cell(input_table(7,31)));
Initial_Soil_Moisture_path = char(table2cell(input_table(8,31)));

% Rasters
fname_LULC = LULC_path; fname_DEM = DEM_path;
fname_SOIL = SOIL_path;



LULC_raster = GRIDobj(fname_LULC); % Land Use and Land Cover Classification
DEM_raster = GRIDobj(fname_DEM); % Digital Elevation Model (m)
SOIL_raster = GRIDobj(fname_SOIL); % Soil Map
min_dem_value = -200; % min value that a dem can have

% Checking CRS
try
if isempty(DEM_raster.georef.SpatialRef.ProjectedCRS)
    % Enter with CRS code for the raster
%     prompt = "What is the CRS code for your raster data? ";
%     code = input(prompt);
%     projected_data = projcrs(code);
%     gtiffinfo.SpatialRef.ProjectedCRS = projected_data;
%     DEM_raster.georef.SpatialRef  = gtiffinfo.SpatialRef;
%     LULC_raster.georef.SpatialRef  = gtiffinfo.SpatialRef;
%     SOIL_raster.georef.SpatialRef  = gtiffinfo.SpatialRef;
    [~,R] = readgeoraster(fname_DEM);
    proj = R.ProjectedCRS;
    DEM_raster.georef.SpatialRef.ProjectedCRS  = proj;
    LULC_raster.georef.SpatialRef.ProjectedCRS  = proj;
    SOIL_raster.georef.SpatialRef.ProjectedCRS  = proj;
end
catch me
    if isempty(DEM_raster.georef)
       % some code?
    end
end


%Changing CRS


%% Checking Extent Problem
crs_save = DEM_raster.georef.SpatialRef.ProjectedCRS;
% croping all nan rows and columns
DEM_raster = crop(DEM_raster);
LULC_raster = crop(LULC_raster);
SOIL_raster = crop(SOIL_raster);
% 
DEM_raster.georef.SpatialRef.ProjectedCRS = crs_save;
LULC_raster.georef.SpatialRef.ProjectedCRS = crs_save;
SOIL_raster.georef.SpatialRef.ProjectedCRS = crs_save;

% if min(min(DEM_raster.Z)) <= 0
%     error('Please make sure that you actually have negative or 0 values in the DEM. Otherwise, treat non-value points as NaN or -9999.')
% end
if sum(size(DEM_raster.Z)) == sum(size(LULC_raster.Z)) && sum(size(DEM_raster.Z)) == sum(size(SOIL_raster.Z))
else
    if sum(size(DEM_raster.Z)) >= sum(size(LULC_raster.Z)) && sum(size(DEM_raster.Z)) >= sum(size(SOIL_raster.Z)) % DEM is larger
        raster_resample = DEM_raster;
        % Resample other two rasters
        % ---- Constraint at LULC Raster
        LULC_raster = resample(LULC_raster,raster_resample,'nearest');
        % LULC_raster.Z = round(LULC_raster.Z,'nearest'); % Only Integers
        % ---- Constraint at SOIL Raster
        SOIL_raster = resample(SOIL_raster,raster_resample,'nearest');
        % SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
    end

    if sum(size(SOIL_raster.Z)) >= sum(size(DEM_raster.Z)) && sum(size(SOIL_raster.Z)) >= sum(size(LULC_raster.Z))  % SOIL is larger
        raster_resample = SOIL_raster;
        % Resample other two rasters
        LULC_raster = resample(LULC_raster,raster_resample,'nearest');
        % ---- Constraint at LULC Raster
        LULC_raster = resample(LULC_raster,raster_resample,'nearest');
        % LULC_raster.Z = round(LULC_raster.Z); % Only Integers
        DEM_raster = resample(DEM_raster,raster_resample,'bilinear');
    end

    if sum(size(LULC_raster.Z)) >= sum(size(DEM_raster.Z)) && sum(size(LULC_raster.Z)) >= sum(size(SOIL_raster.Z))  % SOIL is larger
        raster_resample = LULC_raster;
        % Resample other two rasters
        % ---- Constraint at SOIL Raster
        SOIL_raster = resample(SOIL_raster,raster_resample,'nearest');
        % SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
        DEM_raster = resample(DEM_raster,raster_resample,'bilinear');
    end
end

% Checking if there are nan cells within the study area to avoid numerical
% instability 
if any(~isnan(DEM_raster.Z).*isnan(LULC_raster.Z)) == 1
    imagesc(~isnan(DEM_raster.Z).*isnan(LULC_raster.Z));
    error('Please, check your DEM and LULC rasters. There are cells with no information which will produce numerical instability')
elseif any(~isnan(DEM_raster.Z).*isnan(SOIL_raster.Z)) == 1 
    imagesc(~isnan(DEM_raster.Z).*isnan(SOIL_raster.Z));
    error('Please, check your DEM and SOIL rasters. There are cells with no information which will produce numerical instability')
end

% Raster Extent
GIS_data.xulcorner = DEM_raster.refmat(3,1); % Up Left Corner
GIS_data.yulcorner = DEM_raster.refmat(3,2);
% - Extent is already solved, we can login the input data
Wshed_Properties.Resolution = DEM_raster.cellsize; % m

%% Load General Data Input
input_data_script;  % Load general data, soil, and LULC parameters

%% Resampling Maps
% In case we want to resample the DEM
if flags.flag_resample == 1
    resolution = GIS_data.resolution_resample; % m
    % DEM
    DEM_raster = resample(DEM_raster,resolution,'bilinear');
    % LULC
    LULC_raster = resample(LULC_raster,resolution);
    LULC_raster.Z = round(LULC_raster.Z);
    % SOIL
    SOIL_raster = resample(SOIL_raster,resolution);
    SOIL_raster.Z = round(SOIL_raster.Z);

    % Extent Problem
    if sum(size(DEM_raster.Z)) > sum(size(LULC_raster.Z)) && sum(size(DEM_raster.Z)) >= sum(size(SOIL_raster.Z)) % DEM is larger
        raster_resample = DEM_raster;
        % Resample other two rasters
        % ---- Constraint at LULC Raster
        LULC_raster = resample(LULC_raster,raster_resample);
        LULC_raster.Z = round(LULC_raster.Z); % Only Integers
        % ---- Constraint at SOIL Raster
        SOIL_raster = resample(SOIL_raster,raster_resample);
        SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
    end

    if sum(size(SOIL_raster.Z)) > sum(size(DEM_raster.Z)) && sum(size(SOIL_raster.Z)) >= sum(size(LULC_raster.Z))  % SOIL is larger
        raster_resample = SOIL_raster;
        % Resample other two rasters
        LULC_raster = resample(LULC_raster,raster_resample);
        % ---- Constraint at LULC Raster
        LULC_raster = resample(LULC_raster,raster_resample);
        LULC_raster.Z = round(LULC_raster.Z); % Only Integers
        DEM_raster = resample(DEM_raster,raster_resample);
    end

    if sum(size(LULC_raster.Z)) > sum(size(DEM_raster.Z)) && sum(size(LULC_raster.Z)) >= sum(size(SOIL_raster.Z))  % SOIL is larger
        raster_resample = DEM_raster;
        % Resample other two rasters
        % ---- Constraint at SOIL Raster
        SOIL_raster = resample(SOIL_raster,raster_resample);
        SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
        DEM_raster = resample(DEM_raster,raster_resample);
    end

    % Raster Extent
    GIS_data.xulcorner = DEM_raster.refmat(3,1); % Up Left Corner
    GIS_data.yulcorner = DEM_raster.refmat(3,2);
    % - Extent is already solved, we can login the input data
    Wshed_Properties.Resolution = DEM_raster.cellsize; % m
end

%% Observed Gauges with Resampled DEM
% Observed Gauges
if flags.flag_obs_gauges ==1
    input_table = readtable(model_folder);
    input_table_observed = (input_table(2:end,37:39));
    obs_gauges = table2array(input_table_observed(:,1)); % Indexes
    gauges.num_obs_gauges = sum(~isnan(obs_gauges));
    gauges.easting_obs_gauges_absolute = table2array(input_table_observed(1:gauges.num_obs_gauges,2)); % Easting Coordinates
    gauges.northing_obs_gauges_absolute = table2array(input_table_observed(1:gauges.num_obs_gauges,3)); % Northing Coordinates
    % --- Converting coordinates to local coordinates in pixels
    gauges.easting_obs_gauges = round((-GIS_data.xulcorner + gauges.easting_obs_gauges_absolute)/Wshed_Properties.Resolution);
    gauges.northing_obs_gauges = round((GIS_data.yulcorner - gauges.northing_obs_gauges_absolute)/Wshed_Properties.Resolution);
    
    % Getting the labels
    labels_observed = (input_table(2:(2+gauges.num_obs_gauges-1),40));
    if flags.flag_obs_gauges == 1
        for i = 1:gauges.num_obs_gauges
            gauges.labels_observed_string{i,:} = string(labels_observed{i,:});
        end
    end
end
%% ETP with new resampled DEM
if flags.flag_ETP == 1
    GRIDobj2geotiff(DEM_raster,'DEM_ETP.tif');
    [DEM_etp,R_etp] = readgeoraster('DEM_ETP.tif'); % Getting Raster Information
    ETP_Parameters.DEM_etp = double(DEM_etp);
    ETP_Parameters.info = geotiffinfo('DEM_ETP.tif');
end

%% ----- Transforming Raster into Matrix with Values ----- %

LULC = double(LULC_raster.Z);
DEM = double(DEM_raster.Z);
SOIL = double(SOIL_raster.Z);

% temporal modification, only to consider negative values in DEM due to
% coastal negative values
if flags.flag_input_rainfall_map
    inf_nan_MAPS = isinf(DEM) + isnan(DEM); % Logical array
else
    neg_DEM = DEM < min_dem_value;
    neg_LULC = LULC < 0;
    neg_SOIL = SOIL < 0;
    inf_nan_MAPS = isinf(DEM) + isnan(DEM) + neg_DEM + isnan(LULC) + isnan(SOIL) + neg_LULC + neg_SOIL + isinf(LULC) + isinf(SOIL); % Logical array
end
idx = inf_nan_MAPS > 0;

% Rebuilding Rasters to the Lowest Extent
LULC_raster.Z = LULC;% Land Use and Land Cover Classification
DEM_raster.Z = DEM; % Digital Elevation Model
SOIL_raster.Z = SOIL; % Soil Map

% Replacing Values with Issues
DEM_raster.Z(idx) = nan;
LULC_raster.Z(idx) = nan;
SOIL_raster.Z(idx) = nan;

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
    Wshed_Properties.inflow_cells = zeros(ny,nx,Inflow_Parameters.n_stream_gauges);
    for i = 1:Inflow_Parameters.n_stream_gauges
        for z = 1:Wshed_Properties.n_inlets(i)
            x = easting_inlet_cells(z,i);
            y = northing_inlet_cells(z,i);
            Wshed_Properties.inflow_cells(y,x,i) = 1; % Coordinates of each stream gauge
        end
    end
end

%% ------------ Rainfall Matrices ------------ %
if flags.flag_rainfall == 0 % No rainfall
    Wshed_Properties.rainfall_matrix = flags.flag_rainfall*zeros(size(dem));
elseif flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && flags.flag_real_time_satellite_rainfall ~= 1 && flags.flag_satellite_rainfall ~= 1 && flags.flag_input_rainfall_map ~= 1
    % Spatial Rainfall Case
    input_table = readtable('Rainfall_Spatial_Input.xlsx');
    % Observations
    n_obs = sum((table2array(input_table(:,2))>=0)); % Number of observations
    n_max_raingauges = 50;
    Spatial_Rainfall_Parameters.time_step_spatial = table2array(input_table(7,2)) - table2array(input_table(6,2)); % min
    end_rain = (n_obs-1)*Spatial_Rainfall_Parameters.time_step_spatial;
    %     rainfall_spatial_duration = 0:time_step_spatial:(end_rain); % Rainfall data time in minutes
    Spatial_Rainfall_Parameters.rainfall_spatial_duration = 0:running_control.record_time_spatial_rainfall:(end_rain); % Rainfall data time in minutes
    Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg = 0:Spatial_Rainfall_Parameters.time_step_spatial:(end_rain); % Rainfall data time in minutes
    n_spatial_agg = 1 + running_control.record_time_spatial_rainfall/Spatial_Rainfall_Parameters.time_step_spatial;
    rainfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),n_spatial_agg);
    % rianfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),saver_memory_maps);

    % Rainfall Data
    for i = 1:n_max_raingauges
        Spatial_Rainfall_Parameters.rainfall_raingauges(:,i) = table2array(input_table(6:end,3 + (i-1)));
        Spatial_Rainfall_Parameters.coordinates(i,1) = table2array(input_table(3,3 + (i-1)));
        Spatial_Rainfall_Parameters.coordinates(i,2) = table2array(input_table(4,3 + (i-1)));
    end
    Spatial_Rainfall_Parameters.n_raingauges = sum(Spatial_Rainfall_Parameters.rainfall_raingauges(1,:) >= 0); % Number of raingauges
elseif flags.flag_input_rainfall_map == 1 
    % We are running rainfall with input tiff maps
    % Observations
    n_obs = sum((table2array(input_table_rainfall(:,1))>=0)); % Number of observations
    Spatial_Rainfall_Parameters.time_step_spatial = table2array(input_table_rainfall(2,1)) - table2array(input_table_rainfall(1,1)); % min
    end_rain = (n_obs-1)*Spatial_Rainfall_Parameters.time_step_spatial;
    Spatial_Rainfall_Parameters.rainfall_spatial_duration = 0:running_control.record_time_spatial_rainfall:(end_rain); % Rainfall data time in minutes
    Spatial_Rainfall_Parameters.rainfall_spatial_duration_agg = 0:Spatial_Rainfall_Parameters.time_step_spatial:(end_rain); % Rainfall data time in minutes
    n_spatial_agg = 1 + running_control.record_time_spatial_rainfall/Spatial_Rainfall_Parameters.time_step_spatial;
    rainfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),n_spatial_agg);
    % rianfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),saver_memory_maps);

elseif flags.flag_rainfall == 1 && flags.flag_satellite_rainfall == 1
    n_spatial_agg = 1 + round(running_control.record_time_maps/running_control.record_time_spatial_rainfall,1);
    if mod(running_control.record_time_maps,running_control.record_time_spatial_rainfall) ~= 0 
        error('Please, enter a multiple of the recording time for the input rainfall time-step')
    end
    rainfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),n_spatial_agg);
    % rianfall_spatial_aggregation = zeros(size(dem,1),size(dem,2),saver_memory_maps);

elseif flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~= 1
    % Lumped Rainfall Case
    Wshed_Properties.rainfall_matrix = flags.flag_rainfall*ones(size(dem));
end

%% ------------ ETP Matrices ------------ %%
if flags.flag_ETP == 1
    % We are running ETP
    ETP_Parameters.Krs = 0.16; % Parameter
    ETP_Parameters.alfa_albedo_input = 0.23; % Parameter
    input_table = readtable('ETP_input_data.xlsx');
    % Observations
    n_obs_ETP = min(size(input_table,1) - 2,days(date_end - date_begin));
    ETP_Parameters.n_obs_ETP = size(input_table,1) - 2;
    ETP_Parameters.n_max_etp_stations = 50;
    ETP_Parameters.time_step_etp = minutes(datetime(table2array(input_table(5,2)) ) - datetime(table2array(input_table(4,2)) )); % min
    ETP_Parameters.end_etp = min((ETP_Parameters.n_obs_ETP-1)*ETP_Parameters.time_step_etp,running_control.routing_time);
    ETP_Parameters.time_ETP = datetime(table2array(input_table(3:end,2)));
    ETP_Parameters.time_ETP_begin = ETP_Parameters.time_ETP(1);

    % Check Initial ETP Dates
    ETP_Parameters.delta_ETP_date = minutes(ETP_Parameters.time_ETP_begin - date_begin); % minutes positive or negative
    ETP_Parameters.climatologic_spatial_duration = 0:ETP_Parameters.time_step_etp:(ETP_Parameters.end_etp) ; % Rainfall data time in minutes
    ETP_Parameters.climatologic_spatial_duration = ETP_Parameters.climatologic_spatial_duration + ETP_Parameters.delta_ETP_date;
    [ETP_Parameters.climatologic_spatial_duration(1,1)] = 0;
    
    % Preallocating Arrays
    % Maps.Hydro.ETP_save = zeros(size(dem,1),size(dem,2),ETP_Parameters.n_obs_ETP); % Maps of ET
    Maps.Hydro.ETP_save = zeros(size(dem,1),size(dem,2),saver_memory_maps);
    Maps.Hydro.ETR_save = Maps.Hydro.ETP_save;
    % ---- Loop ETP --- %
    for i = 1:ETP_Parameters.n_max_etp_stations
        try
            ETP_Parameters.n_stations = 50;
            % Maximum Temperature
            ETP_Parameters.maxtemp_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 3));
            % Minimum Temperature
            ETP_Parameters.mintemp_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 4));
            % Average Temperature
            ETP_Parameters.avgtemp_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 5));
            % U2
            ETP_Parameters.u2_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 6));
            % UR
            ETP_Parameters.ur_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 7));
            % G
            ETP_Parameters.G_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 8));

            % Coordinates
            ETP_Parameters.coordinates_stations(i,1) = table2array(input_table(1,6*(i-1) + 6));
            ETP_Parameters.coordinates_stations(i,2) = table2array(input_table(1,6*(i-1) + 8));
        catch
            ETP_Parameters.n_stations = i-1;
            break
        end
        if i == ETP_Parameters.n_max_etp_stations
            ETP_Parameters.n_stations = ETP_Parameters.n_max_etp_stations;
        end
    end

    % DEM Raster Information for ETP Calcualtion

    ETP_Parameters.DEM_etp(idx) = nan;    
    ETP_Parameters.height = ETP_Parameters.info.Height; % Integer indicating the height of the image in pixels
    ETP_Parameters.width = ETP_Parameters.info.Width; % Integer indicating the width of the image in pixels
    [ETP_Parameters.cols_etp,ETP_Parameters.rows_etp] = meshgrid(1:ETP_Parameters.width,1:ETP_Parameters.height);
    [ETP_Parameters.x_etp,ETP_Parameters.y_etp] = pix2map(ETP_Parameters.info.RefMatrix, ETP_Parameters.rows_etp, ETP_Parameters.cols_etp); % Map Coordinates
    [ETP_Parameters.lat,ETP_Parameters.lon] = projinv(ETP_Parameters.info.SpatialRef.ProjectedCRS, ETP_Parameters.x_etp,ETP_Parameters.y_etp); % Latitude and Longitude
    ETP_Parameters.neg_DEM = ETP_Parameters.DEM_etp <= 0;
    ETP_Parameters.DEM_etp(ETP_Parameters.neg_DEM) = nan;
    ETP_Parameters.idx_cells = ETP_Parameters.DEM_etp >= min_dem_value;
    % ETP_Parameters.idx_cells = ETP_Parameters.double(idx_cells(:));
    ETP_Parameters.lat(ETP_Parameters.neg_DEM) = nan; % Latitude
    ETP_Parameters.lon(ETP_Parameters.neg_DEM) = nan; % Longitude
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
%% Fillsinks
% max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
% DEM_filled = fillsinks(DEM,max_depth);
if flags.flag_fill_DEM == 1
    DEM_filled = fillsinks(DEM_raster); % Filled DEM
    DIFFDEM = DEM_filled - DEM;
    DIFFDEM.Z(DIFFDEM.Z==0) = nan;
    DEM_raster = DEM_filled;
    % imageschs(DEM_raster,DIFFDEM.Z);
end

%% Smooth DEM
if flags.flag_smooth_cells == 1
    Vq = imgaussfilt(dem,'FilterSize',3);
    dem_diff_smooth = Vq;
    dem = Vq; % New dem
    pause(2)
    DEM_raster.Z = dem;
    close all
end

%% Decrese Elevations in Creeks
if flags.flag_reduce_DEM
    FD = FLOWobj(DEM_raster); % Flow direction
    As  = flowacc(FD); % Flow Accumulation
    Wshed_Properties.fac_area = As.Z*(Wshed_Properties.cell_area/1000/1000); % km2
    idx_facc = Wshed_Properties.fac_area <= GIS_data.min_area;

    % B
    B = GIS_data.beta_1*Wshed_Properties.fac_area.^(GIS_data.beta_2);
    H = GIS_data.alfa_1*Wshed_Properties.fac_area.^(GIS_data.alfa_2);
    Flow_Area = B.*H; % m2
    H_abg = Flow_Area/Wshed_Properties.Resolution;
    if max(max(H_abg)) > 10^2
        error('Be careful. The decreasing in the DEM to consider the water depths are larger than 100 m. You may change the parameters.')
    end
    H_abg(idx_facc) = 0;

    % Reducing Elevation in creeks
    DEM_raster.Z = DEM_raster.Z - H_abg; % [m]

    % New Data
    dem = DEM_raster.Z;
end

%% DEM Smoothening
if flags.flag_smoothening == 1
    [DEM_raster,DEM,S] = DEM_smoothening(DEM_raster,GIS_data.min_area,flags.flag_trunk,GIS_data.tau,GIS_data.K_value);
end

%% Imposing Minimum Slope - If Required
if flags.flag_diffusive ~= 1
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
SLP = arcslope(DEM_raster);
close all
imagesc(SLP); colorbar
close all

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

%% DEM with Streams and Obseved Points

FD = FLOWobj(DEM_raster);
area_km2 = GIS_data.min_area; % km2
area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
MS = STREAMobj2mapstruct(S);

FileName_String = 'streamnetwork.shp';
FileName = fullfile(folderName,strcat('\',FileName_String));
try
    shapewrite(MS,FileName);
catch ME
   error('The threshold for flow accumulation might be too high to create a stream network.')
end


if flags.flag_obs_gauges == 1
    % Example x and y coordinates
    gauges.x_coord_gauges = gauges.easting_obs_gauges_absolute;
    gauges.y_coord_gauges = gauges.northing_obs_gauges_absolute;
    labels_gauges = gauges.labels_observed_string;  % Labels for each point


    % Create a shapefile writer
    % shapefile = shapewrite('observed_gauges.shp');

    % Create a geoshape object with the points and labels
    points = mappoint(gauges.x_coord_gauges', gauges.y_coord_gauges', 'Label', labels_gauges');

    % Write the geoshape object to a shapefile
    FileName_String = 'observed_gauges.shp';
    FileName = fullfile(folderName,strcat('\',FileName_String));
    shapewrite(points, FileName); % This is in the same coordinate system of your DEM
end

close all

%% Drainage Basins
% D = drainagebasins(FD);
% imageschs(DEM_raster,shufflelabel(D))
% area_km2 = GIS_data.min_area;
% area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
% imagesc(DEM_raster);
% hold on
% S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
% if ~isempty(S)
%     plot(S,'linewidth',2,'color','red');
%     plot(trunk(S),'linewidth',2,'color','red');
%     D = drainagebasins(FD,S);
%     imageschs(DEM_raster,shufflelabel(D))
% end

% hold on
%
% plot(S,'linewidth',2,'color','red');% area_km2 = 500; % km2
% area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
% imagesc(DEM_raster);
% hold on
% S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
% plot(S,'linewidth',2,'color','red');
%
% D = drainagebasins(FD,S);
% imageschs(DEM_raster,shufflelabel(D))
% hold on
%
% plot(S,'linewidth',2,'color','red');


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
            if dem(i,j) > 0 % Boundary has values
                perimeter(i,j) = 1;
            else
                perimeter(i,j) = 0;
            end
        end
        if i ~= 1 && i ~= size(dem,1) && j ~= 1 && j ~= size(dem,2)
            if isnan(dem(i,j-1)) && isnan(dem(i,j+1)) && isnan(dem(i-1,j)) && isnan(dem(i+1,j))
                perimeter(i,j) = 0;
            elseif  dem(i,j) > 0 && (isnan(dem(i,j-1)) || isnan(dem(i,j+1)) || isnan(dem(i-1,j)) || isnan(dem(i+1,j)))
                perimeter(i,j) = 1;
            else
                perimeter(i,j) = 0;
            end
        end
    end
end
[row_boundary, col_boundary] = find(perimeter > 0); % Boundary Cells
max_length = (max(col_boundary) - min(col_boundary))*Wshed_Properties.Resolution - (max(row_boundary) - min(row_boundary))*Wshed_Properties.Resolution; % Length
max_width = (max(row_boundary) - min(row_boundary))*Wshed_Properties.Resolution - (max(row_boundary) - min(row_boundary))*Wshed_Properties.Resolution; % width

% Number of Outlets
idx_outlet = elevation_nan == min(min(elevation_nan)); % Finding cells with lowest elevations
[Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(idx_outlet == 1); % finding rows and cols of these cells
outlet_index(idx_outlet) = 1;

% Perimeter
Wshed_Properties.perimeter = perimeter;


%% Outlet Calculations - 1
% Here we divide the number of outlets symetrically from the prior outlet
% New way: Assuming that the outlets should be at the border of the domain, we can
% do:
for i = 1:length(Wshed_Properties.row_outlet)
    % Checking left, right up, and down
    row = Wshed_Properties.row_outlet(i); % Row
    col = Wshed_Properties.col_outlet(i); % Col
    % Left
    if  row - 1 == 0 || isnan(elevation_nan(row-1,col))
        outlet_index(Wshed_Properties.row_outlet(i),Wshed_Properties.col_outlet(i)) = 1; % This is an outlet
    end
    % Right
    if  row + 1 > size(elevation_nan,1) || isnan(elevation_nan(row+1,col))
        outlet_index(Wshed_Properties.row_outlet(i),Wshed_Properties.col_outlet(i)) = 1; % This is an outlet
    end
    % Up
    if col - 1 == 0 || isnan(elevation_nan(row,col-1))
        outlet_index(Wshed_Properties.row_outlet(i),Wshed_Properties.col_outlet(i)) = 1; % This is an outlet
    end
    % Down
    if  col + 1 > size(elevation_nan,2) || isnan(elevation_nan(row,col+1))
        outlet_index(Wshed_Properties.row_outlet(i),Wshed_Properties.col_outlet(i)) = 1; % This is an outlet
    end
end
[Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1);  % Refreshing outlet

% Minimum Row
min_Wshed_Properties.row_outlet = min(min(Wshed_Properties.row_outlet));
% Moving Up
pos = find(Wshed_Properties.row_outlet == min_Wshed_Properties.row_outlet);
col = Wshed_Properties.col_outlet(pos);
row = Wshed_Properties.row_outlet(pos);

% 0 , 45, 90, 135, 180, 225, 270, 315
check = zeros(8,1);
outlet_new = zeros(size(dem));
right_loop = round(n_outlets_data/2);
left_loop =  n_outlets_data - right_loop;
for i = 1:right_loop
    % 0, 45, 90, 135, 180
    if col + 1 <= size(dem,2) && perimeter(row,col + 1) > 0  % 0
        col = col+1;
    elseif row - 1 >= 1 && col + 1 <= size(dem,2) && perimeter(row - 1,col + 1) > 0  % 45
        row = row-1;
        col = col+1;
    elseif row - 1 >= 1 && perimeter(row-1,col) > 0  % 90
        row = row-1;
    elseif row - 1 >= 1 && col - 1 >= 1 && perimeter(row-1,col-1) > 0  % 135
        row = row-1;
        col = col-1;
    end
    outlet_new(row,col) = 1;
end

outlet_index = outlet_index + outlet_new; % Refreshing Outlet
% Making Sure to not have a raster larger than 1
outlet_index(outlet_index>0) = 1;

% Maximum Row
max_row_outlet = max(max(Wshed_Properties.row_outlet));
pos = find(Wshed_Properties.row_outlet == max_row_outlet);
col = Wshed_Properties.col_outlet(pos);
row = Wshed_Properties.row_outlet(pos);
for i = 1:left_loop
    % 180, 225, 270, 315, 360
    if  col - 1 >= 1 && perimeter(row,col-1) > 0  % 180
        col = col - 1;
    elseif row + 1 <= size(dem,1) && col - 1 >= 1 && perimeter(row+1,col - 1) > 0  % 225
        row = row + 1;
        col = col - 1;
    elseif row + 1 <= size(dem,1) && perimeter(row + 1,col) > 0  % 270
        row = row+1;
    elseif row + 1 <= size(dem,1) && col + 1 <= size(dem,2) && perimeter(row+1,col+1) > 0  % 315
        row = row+1;
        col = col+1;
    end
    outlet_new(row,col) = 1;
end

outlet_index = outlet_index + outlet_new; % Refreshing Outlet


% Making Sure to not have a raster larger than 1
outlet_index(outlet_index>0) = 1;

% Filling Inside Gaps
[row_outlets,col_outlets] = find(outlet_index == 1);
if n_outlets_data > 0
    for i = 1:sum(sum(outlet_index))
        row_out = row_outlets(i);
        col_out = col_outlets(i);
        % 0, 45, 90, 135, 180, 225, 270, 315, 360
        r = row_out; c = col_out + 1;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 0
            outlet_index(r,c) = 1;
        end
        r = row_out-1; c = col_out+1;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 45
            outlet_index(r,c) = 1;
        end
        r = row_out-1; c = col_out;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 90
            outlet_index(r,c) = 1;
        end
        r = row_out-1; c = col_out-1;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 135
            outlet_index(r,c) = 1;
        end
        r = row_out; c = col_out-1;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 180
            outlet_index(r,c) = 1;
        end
        r = row_out+1; c = col_out-1;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 225
            outlet_index(r,c) = 1;
        end
        r = row_out+1; c = col_out;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 270
            outlet_index(r,c) = 1;
        end
        r = row_out+1; c = col_out+1;
        if outlet_index(r,c) == 0 && perimeter(r,c) == 1 % 315
            outlet_index(r,c) = 1;
        end
    end
end

% Delete
% outlet_index = perimeter;
idx_outlet = logical(outlet_index); % Added refreshed idx_outlet


[Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1); % finding rows and cols of these cells


[Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1);  % Refreshing outlet
Wshed_Properties.el_outlet = dem;
Wshed_Properties.el_outlet(idx_outlet < 1) = nan;
Wshed_Properties.stage_min = min(min(Wshed_Properties.el_outlet));
[Wshed_Properties.row_min, Wshed_Properties.col_min] = find(dem == Wshed_Properties.stage_min);


% % Save Map
FileName = 'Outlet_Mask.tif';
FileName = fullfile(folderName,FileName);
% Exporting Outlet_Index as a Mask
raster_to_export = DEM_raster; % Just to get the properties
raster_to_export.Z = outlet_index; % Adding Outlets
raster_to_export.Z(isnan(DEM_raster.Z)) = nan;
try
    GRIDobj2geotiff(raster_to_export,FileName) % Exporting the Map
catch me
    warning('Outlet mask not exported.')
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
% In case flags.flag_abstraction == 1, we cut the domain assuming the cells
% entered in the input_data file
if flags.flag_abstraction ~= 1
    xmin = 1; % initial position x in the grid (collums)
    ymin = 1; % lines (up to down)
    xmax = zzz(2);
    ymax = zzz(1);
end
% ------------ Cutting Cells ------------
if flags.flag_inflow == 1
    Wshed_Properties.inflow_cells = Wshed_Properties.inflow_cells(ymin:ymax,xmin:xmax,:);
end
outlet_index = outlet_index(ymin:ymax,xmin:xmax);
elevation = dem(ymin:ymax,xmin:xmax); % using only specified grid
spatial_domain = zeros(size(dem));
dem = dem(ymin:ymax,xmin:xmax);
lulc_matrix = imp(ymin:ymax,xmin:xmax); % using only specified grid
soil_matrix  = soil(ymin:ymax,xmin:xmax); % using only specified grid;
Wshed_Properties.rainfall_matrix = double(elevation >= min_dem_value) ; % Elevation could be negative
Wshed_Properties.rainfall_matrix(Wshed_Properties.rainfall_matrix == 0) = nan;

% ----------- Struct Cells ------------- %
outflow_rate = struct('qout_left_t',spatial_domain,'qout_right_t',spatial_domain,'qout_up_t',spatial_domain,'qout_down_t',spatial_domain); % Struct array
% depths = struct('d_0',spatial_domain,'d_avg',spatial_domain,'d_left_cell',spatial_domain,'d_right_cell',spatial_domain,'d_up_cell',spatial_domain,'d_down_cell',spatial_domain,'d_tot',spatial_domain,'depth_wse',depth_wse); % struct
depths = struct('d_0',spatial_domain,'d_avg',spatial_domain,'d_tot',spatial_domain,'depth_wse',depth_wse); % struct
velocities = struct('vel_left',spatial_domain,'vel_right',spatial_domain,'vel_up',spatial_domain,'vel_down',spatial_domain,'velocity_raster',spatial_domain,'vmax_final',spatial_domain); % struct

Soil_Properties = struct('ksat',spatial_domain,'psi',spatial_domain,'teta_sat',spatial_domain,'teta_i',spatial_domain,'I_0',spatial_domain); % struct
LULC_Properties = struct('roughness',spatial_domain,'h_0',spatial_domain,'C_1',spatial_domain,'C_2',spatial_domain,'C_3',spatial_domain,'C_4',spatial_domain,'B_0',spatial_domain,'ADD',0); % struct

% Elevation_Properties = struct('elevation_cell',spatial_domain,'elevation_left_t',spatial_domain,'elevation_right_t',spatial_domain,'elevation_up_t',spatial_domain,'elevation_down_t',spatial_domain); % struct

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
for i = 1:LULC_Properties.n_lulc % 8 types of LULC
    % Only Roughness and h_0
    LULC_Properties.roughness(LULC_Properties.idx_lulc(:,:,i)) = lulc_parameters(i,1); % assigning values for roughness at impervious areas
    LULC_Properties.h_0(LULC_Properties.idx_lulc(:,:,i)) = lulc_parameters(i,2); % Initial Abstraction

    % ------------ Warm-up Data ------------:
    if flags.flag_warmup == 1
        if i == 1
            % Identifying positive values of the Warmup
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
        end
    else
        depths.d_0(LULC_Properties.idx_lulc(:,:,i)) = lulc_parameters(i,3);
    end
    if flags.flag_waterquality == 1 % Only if water quality is being modeled
        LULC_Properties.C_1(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,4);
        LULC_Properties.C_2(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,5);
        LULC_Properties.C_3(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,6);
        LULC_Properties.C_4(LULC_Properties.idx_lulc(:,:,i)) =  lulc_parameters(i,7);
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
            Initial_Buildup_Raster = GRIDobj(Warmup_Depth_path);
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

% Soil Parameter Assigning
Soil_Properties.Soil = SOIL;
for i = 1:Soil_Properties.n_soil
    Soil_Properties.ksat(idx_soil(:,:,i)) = soil_parameters(i,1); % similar things happening below    
    Soil_Properties.psi(idx_soil(:,:,i)) = soil_parameters(i,2);
    Soil_Properties.I_0(idx_soil(:,:,i)) = soil_parameters(i,3);
    Soil_Properties.teta_sat(idx_soil(:,:,i)) = soil_parameters(i,4);
    Soil_Properties.teta_i(idx_soil(:,:,i)) = soil_parameters(i,5);
end

% Mask in Impervious Areas. Cells that are built areas have no infiltration
Soil_Properties.ksat(LULC_Properties.idx_imp) = 0; % Impervious areas
Soil_Properties.soil_matrix = soil_matrix;
Soil_Properties.idx_soil = idx_soil;
if flags.flag_warmup == 1
    % Identifying positive values of the Warmup
    Warmup_Raster_Moisture = GRIDobj(Initial_Soil_Moisture_path);
    if sum(size(Warmup_Raster_Moisture.Z)) ~= sum(size(DEM))
        Warmup_Raster_Moisture = resample(Warmup_Raster_Moisture,DEM_raster);
    end
    Initial_Soil_Moisture = double(Warmup_Raster_Moisture.Z);
    %     Initial_Soil_Moisture = Initial_Soil_Moisture./Initial_Soil_Moisture*10; % DELETE
    idx_warmup_moisture = Warmup_Depths >= 0;
    Soil_Properties.I_0 =  idx_warmup_moisture.*Initial_Soil_Moisture; % Depths in m
    % Treat inf values
    Soil_Properties.I_0(idx_nan) = nan; % carefull here
end

% Initial Moisture
Soil_Properties.I_0(LULC_Properties.idx_imp) = 0; % Impervious areas

% Replenishing Coefficient
Soil_Properties.kr = (1/75*(sqrt(Soil_Properties.ksat/25.4))); % Replenishing rate (1/hr) (Check Rossman, pg. 110)
Soil_Properties.Tr = 4.5./sqrt(Soil_Properties.ksat/25.4); %  Recovery time (hr) Check rossman pg. 110
Soil_Properties.Lu = 4.*sqrt(Soil_Properties.ksat/25.4); % Inches - Uppermost layer of the soil
Soil_Properties.Lu = Soil_Properties.Lu*2.54/100; % Meters
Soil_Properties.k_out = (Soil_Properties.teta_sat - Soil_Properties.teta_i).*Soil_Properties.kr.*Soil_Properties.Lu*1000; % Rate of replenishing exfiltration from the saturated zone during recoverying times (mm/hr)
% clear idx % clear this matrices to avoid huge storage data

%%%% New xmax and ymax
xmin = 1; %initial position x in the grid (collums)
ymin = 1; % lines (up to down)
zzz = size(elevation);
xmax = zzz(2);
ymax = zzz(1);

%% Locally Reducing Manning of River Networks
if flags.flag_reduce_DEM == 1
    LULC_Properties.roughness(idx_facc) = 1*LULC_Properties.roughness(idx_facc);
    LULC_Properties.h_0(idx_facc) = 0;
end

%% Conversion of inflow into the time-step of calculations
if flags.flag_inflow == 1
    inflow_length = size(Inflow_Parameters.inflow_hydrograph_rate,1); % Number of interval
    inflow_discretized = zeros(size(Inflow_Parameters.inflow_hydrograph_rate,2),ceil(inflow_length*Inflow_Parameters.time_step_inflow/time_step_model)); % Preallocating
    for z = 1:Inflow_Parameters.n_stream_gauges
        for i =1:((inflow_length)*Inflow_Parameters.time_step_inflow/time_step_model)
            inflow_discretized(z,i) = Inflow_Parameters.inflow_hydrograph_rate(ceil((round(i*time_step_model/Inflow_Parameters.time_step_inflow,12))),z); % Discretized into moldel's time-step
        end
    end
end
% Each z row represent an stream gauge and each i col is a inflow value
%% Conversion of rainfall into the time-step of calculations for concentrated rainfall
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~=1  % Only for concentrated rainfall
    intensity_rainfall_length = length(Rainfall_Parameters.intensity_rainfall) ; % Number of intervals
    intensity_discretized = zeros(1,ceil(intensity_rainfall_length*Rainfall_Parameters.time_step_rainfall/time_step_model)); % Preallocating
    for i =1:(intensity_rainfall_length*Rainfall_Parameters.time_step_rainfall/time_step_model)
        intensity_discretized(i) = Rainfall_Parameters.intensity_rainfall(ceil((round(i*time_step_model/Rainfall_Parameters.time_step_rainfall,12))));  % Discretized into moldel's time-step
    end
end
%% Determination of grid parameters and outlet coordinates (Whole Domain)

% Calculates the number of non-inf cells to determine the watershed area
matrix_nan = isnan(elevation);
number_nan = sum(sum(matrix_nan));
Wshed_Properties.drainage_area = (nx*ny - number_nan)*Wshed_Properties.Resolution^2; % m2
Edges = bwperim(dem>=0,8);
Wshed_Properties.watershed_perimeter = sum(sum(Edges))*Wshed_Properties.Resolution/1000;  % km
Wshed_Properties.impervious_area = sum(sum(impervious_cells))*Wshed_Properties.cell_area;
Wshed_Properties.pervious_area = sum(sum(pervious_cells))*Wshed_Properties.cell_area;
Wshed_Properties.impervious_rate = Wshed_Properties.impervious_area/Wshed_Properties.drainage_area;
Wshed_Properties.pervious_rate = Wshed_Properties.impervious_area/Wshed_Properties.drainage_area;

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
%% Calculation of accumulated and incremental inflow on the inflow cells
if flags.flag_inflow == 1
    [~,BC_States.delta_inflow,inflow_intensity] = accumulated_incremental(running_control.steps,inflow_discretized,time_step_model);
    BC_States.time_deltainflow = cumsum(ones(size(BC_States.delta_inflow,2),1)*time_step_model); % Vector with indices
end
%% Calculation of accumulated and incremental precipitation in each cell
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~=1  % Only for concentrated rainfall
    [~,BC_States.delta_p,rainfall_intensity] = accumulated_incremental(running_control.steps,intensity_discretized,time_step_model);
    running_control.time_deltap = cumsum(ones(1,length(BC_States.delta_p))*time_step_model);
    outlet_runoff_volume = zeros(size(BC_States.delta_p));
end
%% Pre allocating more arrays
Total_Inflow = 0;
depths.d_t = spatial_domain;
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
outet_hydrograph = zeros(time_size,1);
time_hydrograph = zeros(time_size,1);
if flags.flag_waterquality == 1
    Maps.WQ_States.Pol_Conc_Map = zeros(ny,nx,time_size_2);
    Maps.WQ_States.Pol_mass_map = zeros(ny,nx,time_size_2);
    Maps.WQ_States.outet_pollutograph = zeros(time_size,1);
end
%% Clearing a few variables
clear accum_precipitation  precipitation imp time_size idx_soil  Land_Cover_Data Elevation_DATA intensity_discretized idx_ idx_1 matrix_nan idx accum_inflow col col_check d_0_imp d_0_per h_0_imp h_0_per I_0_per I_0_imp inflow_discretized

%% ---------------------- Main Routing (GA + 4D/8D CA + BW models) ------------------%%
tic % Start counting time
% ----------------- Initialize Variables -----------------
Soil_Properties.I_t = Soil_Properties.I_0;
Soil_Properties.I_p = Soil_Properties.I_0;
time_calculation_routing = zeros(running_control.steps,1);
k = 1; % Start Counter of the While Loop
recording_parameters.actual_record_state = 1;
recording_parameters.last_record_maps = 1;
recording_parameters.last_record_hydrograph = 1;
if flags.flag_satellite_rainfall == 1
    recording_parameters.last_record_maps_rainfall = 1;
end
if flags.flag_inflow == 1
    for i = 1:Inflow_Parameters.n_stream_gauges
        BC_States.delta_inflow_agg(i,1) = BC_States.delta_inflow(i,1);
    end
end
t = 0;
if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ~=1 && flags.flag_input_rainfall_map ~= 1 % Concentrated Rainfall
    BC_States.delta_p_agg = BC_States.delta_p(1,1)*flags.flag_rainfall;
end

if flags.flag_rainfall ~= 1  % No Rainfall
    BC_States.delta_p_agg = 0;
end

% Grid
Spatial_Rainfall_Parameters.x_grid = GIS_data.xulcorner + Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,2)]'; % Pixel easting coordinates
Spatial_Rainfall_Parameters.y_grid = GIS_data.yulcorner - Wshed_Properties.Resolution*[1:1:size(DEM_raster.Z,1)]'; % Pixel northing coordinates

if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall ==1 && flags.flag_input_rainfall_map ~= 1 && flags.flag_real_time_satellite_rainfall ~= 1  && flags.flag_input_rainfall_map ~= 1 && flags.flag_satellite_rainfall ~= 1  % Spatial Rainfall with Raingauges
    nsteps_spatial_rainfall = floor(running_control.routing_time/running_control.record_time_spatial_rainfall);
    % nsteps_spatial_rainfall = saver_memory_maps;
    if nsteps_spatial_rainfall < 1
        error('Please, enter a feasible time-step for the aggregation of spatial rainfall')
    end
    % Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1),size(dem,2),nsteps_spatial_rainfall); % Maps of Rainfall
    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1),size(dem,2),saver_memory_maps); % Maps of Rainfall
    % Spatial Rainfall
    if t == 0
        z = 1;
    else
        z = find(t <= Spatial_Rainfall_Parameters.rainfall_spatial_duration,1,'last'); % Duration
    end
    Spatial_Rainfall_Parameters.x_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,1); % Coordinates (easting) of each rain gauge
    Spatial_Rainfall_Parameters.y_coordinate = Spatial_Rainfall_Parameters.coordinates(1:Spatial_Rainfall_Parameters.n_raingauges,2); % Coordinates (northing) of each rain gauge
    rainfall = Spatial_Rainfall_Parameters.rainfall_raingauges(1,1:Spatial_Rainfall_Parameters.n_raingauges)'; % Values of rainfall at t for each rain gauge
    idx_rainfall = isnan(rainfall);
    % idx_rainfall = logical(isnan(rainfall) | rainfall == 0);    
    rainfall(idx_rainfall) = []; % Taking out nans
    Spatial_Rainfall_Parameters.x_coordinate(idx_rainfall) = []; % Taking out nans
    Spatial_Rainfall_Parameters.y_coordinate(idx_rainfall) = []; % Taking out nans
    [spatial_rainfall] = Rainfall_Interpolator(Spatial_Rainfall_Parameters.x_coordinate,Spatial_Rainfall_Parameters.y_coordinate,rainfall,Spatial_Rainfall_Parameters.x_grid,Spatial_Rainfall_Parameters.y_grid); % Interpolated Values
    spatial_rainfall(idx_nan) = nan;
    BC_States.delta_p_agg = spatial_rainfall/3600*time_step_model*60; % Matrix of delta P for each pixel
    Maps.Hydro.spatial_rainfall_maps(:,:,1) = spatial_rainfall;
    BC_States.average_spatial_rainfall(1,1) = mean(spatial_rainfall(spatial_rainfall>=0));
end

time_step = running_control.min_time_step*60;
t = time_step_model; % minutes

if flags.flag_rainfall == 1 && flags.flag_spatial_rainfall == 1 && flags.flag_input_rainfall_map == 1 && flags.flag_satellite_rainfall == 0 && flags.flag_real_time_satellite_rainfall == 0 % Input Rainfall Maps
    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1),size(dem,2),saver_memory_maps); % Maps of Rainfall
    % input_rainfall = GRIDobj(string(Input_Rainfall.labels_Directory{1,:}));
    try
        [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{1,:}),'CoordinateSystemType','geographic');
        % Reproject the coordinates from EPSG:4326 to EPSG:3857
        rR.GeographicCRS=geocrs(4326);
        if flags.flag_resample
            if rR.CellExtentInLatitude ~= GIS_data.resolution_resample
                input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,1);
            end
        else
            if rR.CellExtentInLatitude ~= DEM_raster.cellsize
                input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,1);
            end
        end
    catch
        [input_rainfall,rR] = readgeoraster(string(Input_Rainfall.labels_Directory{1,:}));
        if flags.flag_resample
            if rR.CellExtentInWorldX ~= GIS_data.resolution_resample
                input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,0);
            end
        else
            if rR.CellExtentInWorldX ~= DEM_raster.cellsize
                input_rainfall = raster_cutter(DEM_raster,rR,input_rainfall,0);
            end
        end
    end

    input_rainfall = input_rainfall.Z; % Only the values
    Maps.Hydro.spatial_rainfall_maps(:,:,1) = input_rainfall;
    BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
    % Spatial_Rainfall_Parameters.rainfall_spatial_duration = Input_Rainfall.time';
    BC_States.average_spatial_rainfall(1,1) = mean(input_rainfall(input_rainfall>=0));
end



% Satellite Rainfall
register = [];
if flags.flag_satellite_rainfall == 1
    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1),size(dem,2),saver_memory_maps); % Maps of Rainfall
    register = 0;
    product = 'PDIRNow1hourly';
    [rainfall_raster, register,register_data,~] = Satellite_rainfall_processing([],[],register,product,date_begin,date_end,flags.flag_satellite_rainfall,flags.flag_real_time_satellite_rainfall,DEM_raster);
  
    input_rainfall = rainfall_raster.Z; % Only the values
    Maps.Hydro.spatial_rainfall_maps(:,:,1) = input_rainfall;
    BC_States.delta_p_agg = input_rainfall*time_step/60; % mm
    Spatial_Rainfall_Parameters.rainfall_spatial_duration = Input_Rainfall.time;
    BC_States.average_spatial_rainfall(1,1) = mean(input_rainfall(input_rainfall>=0));
end

% Real-Time Rainfall
if flags.flag_real_time_satellite_rainfall == 1
    Maps.Hydro.spatial_rainfall_maps = zeros(size(dem,1),size(dem,2),saver_memory_maps); % Maps of Rainfall
    register = 0;
    product = 'PDIRNow1hourly';
    [rainfall_raster, register] = Satellite_rainfall_processing(register,product,date_begin,date_end,flags.flag_satellite_rainfall,flags.flag_real_time_satellite_rainfall,DEM_raster);
    input_rainfall = rainfall_raster.Z; % Only the values
    Maps.Hydro.spatial_rainfall_maps(:,:,1) = input_rainfall;
    BC_States.delta_p_agg = input_rainfall*time_step/60;
    Spatial_Rainfall_Parameters.rainfall_spatial_duration = Input_Rainfall.time;
    BC_States.average_spatial_rainfall(1,1) = mean(input_rainfall(input_rainfall>=0));
    % Block code for real time rainfall from Persiann
end

running_control.time_step_model = time_step_model;
running_control.time_save_previous = 0; % minutes
Previous_Volume = 0;
t_previous = 0;
step_error = zeros(1,running_control.steps);
% CA_States.I_cell = zeros(ny,nx,5);
WQ_States.P_conc = 0;
tmin_wq = running_control.max_time_step;
dmax_final = zeros(size(elevation));
WQ_States.mass_lost = 0; % initializing
WQ_States.mass_outlet = 0; % initializing
Out_Conc = 0;
WQ_States.vol_outlet = 0;
outlet_runoff_volume = 0; % initializing
WQ_States.Tot_Washed = spatial_domain;
BC_States.outflow_volume = 0;
BC_States.inflow_volume = 0;
CA_States.I_tot_end_cell = zeros(size(spatial_domain)); % Attention here
Hydro_States.ETP = zeros(size(DEM));

% Inflows Boundary Condition
if flags.flag_inflow == 1
    BC_States.inflow = spatial_domain; % This is to solve spatially, don't delete
    for i = 1:Inflow_Parameters.n_stream_gauges
        BC_States.inflow = BC_States.inflow + BC_States.delta_inflow_agg(i,1)*Wshed_Properties.inflow_cells(:,:,i); % mm
    end
else
   BC_States.inflow = 0*spatial_domain;
end

%%%% ELEVATIONS %%%
%%%% ASSIGNING VALUES %%%

Elevation_Properties.elevation_cell = elevation;
if flags.flag_D8 == 1
    zero_matrix = NaN(size(Elevation_Properties.elevation_cell));

    Elevation_Properties.elevation_left_t = [zeros(ny,1),elevation(:,1:(nx-1))];
    Elevation_Properties.elevation_right_t = [elevation(:,(2:(nx))) zeros(ny,1)];
    Elevation_Properties.elevation_up_t = [zeros(1,nx) ; elevation(1:(ny-1),:)];
    Elevation_Properties.elevation_down_t = [elevation(2:(ny),:) ; zeros(1,nx)];

    Elevation_Properties.elevation_NE_t = zero_matrix;
    Elevation_Properties.elevation_SE_t = zero_matrix;
    Elevation_Properties.elevation_SW_t = zero_matrix;
    Elevation_Properties.elevation_NW_t = zero_matrix;

    depths.d_NE_t = zero_matrix; depths.d_SE_t = zero_matrix; depths.d_SW_t = zero_matrix; depths.d_NW_t = zero_matrix;
    %
    Elevation_Properties.elevation_NE_t(2:(ny),1:(nx-1)) = Elevation_Properties.elevation_cell(1:(ny-1),2:nx); % OK
    Elevation_Properties.elevation_SE_t(1:(ny-1),1:(nx-1)) = Elevation_Properties.elevation_cell(2:ny,2:nx); % OK
    Elevation_Properties.elevation_SW_t(1:(ny-1),2:(nx)) = Elevation_Properties.elevation_cell(2:(ny),1:(nx-1)); % OK
    Elevation_Properties.elevation_NW_t(2:ny,2:nx) = Elevation_Properties.elevation_cell(1:(ny-1),1:(nx-1)); % OK
end


%% Checking if all cells of the domain have data
zzz = Wshed_Properties.rainfall_matrix; zzz = zzz > 0;
idx_cells = logical(zzz);
if min(min(LULC_Properties.roughness(idx_cells))) == 0
    warning('Cells with not associated LULC parameters')
    warning('Assuming n = 0.03 for these areas. Also assuming h0 = 0 for them.')
    idx_not_assigned = LULC_Properties.roughness == 0 & idx_cells == 1;
    LULC_Properties.roughness(idx_not_assigned) = 0.03;
    LULC_Properties.h_0(idx_not_assigned) = 0;
    pause(1)
end

if min(min(Soil_Properties.ksat(idx_cells))) == 0
    warning('Cells with not associated SOIL parameters or K = 0')
    warning('Assuming K = 0 for these areas.')
    idx_not_assigned = Soil_Properties.ksat == 0 & idx_cells == 1;
    Soil_Properties.ksat(idx_not_assigned) = 0; 
    pause(1)
end

%% Plotting Input Rasters
Plot_Initial_Maps; % Script to plot initial maps


%% CA-8D Matrices
% Slope
    if flags.flag_GPU
        dim1 = size(Elevation_Properties.elevation_cell,1);
        dim2 = size(Elevation_Properties.elevation_cell,2);
        if flags.flag_D8 == 1
            dim3 = 9;
        else
            dim3 = 5;
        end
        wse_slope_zeros = gpuArray(zeros(dim1,dim2,dim3));
        Distance_Matrix = gpuArray(zeros(dim1,dim2));
        outflow_bates = gpuArray(zeros(dim1,dim2,dim3));
    else
        dim1 = size(Elevation_Properties.elevation_cell,1);
        dim2 = size(Elevation_Properties.elevation_cell,2);
        if flags.flag_D8 == 1
            dim3 = 9;
        else
            dim3 = 5;
        end 
        wse_slope_zeros = (zeros(dim1,dim2,dim3));
        Distance_Matrix = (zeros(dim1,dim2));
        outflow_bates = (zeros(dim1,dim2,dim3));
    end


%% Reservoir Data
if flags.flag_reservoir == 1

end


%% Min Soil Moisture
min_soil_moisture = 5*(double(Soil_Properties.I_0>0));
min_soil_moisture(Soil_Properties.I_0==0) = 0.01;
min_soil_moisture(idx_nan) = nan;
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

%% Converting Arrays to Single Precision
if flags.flag_single == 1
    % Convert all Arrays to Single Arrays
    % Structure Arrays
    BC_States = structfun(@single, BC_States, 'UniformOutput', false);
    CA_States = structfun(@single, CA_States, 'UniformOutput', false);
    Courant_Parameters = structfun(@single, Courant_Parameters, 'UniformOutput', false);
    depths = structfun(@single, depths, 'UniformOutput', false);
    Elevation_Properties = structfun(@single, Elevation_Properties, 'UniformOutput', false);
    flags = structfun(@single, flags, 'UniformOutput', false);
    % Gauges Label
    % if flags.flag_obs_gauges == 1
    %     extra_parameters.gauges.label_observed_string = gauges.labels_observed_string;
    %     gauges.labels_observed_string = [];
    %     gauges = structfun(@gpuArray, gauges, 'UniformOutput', false);
    % end
    GIS_data = structfun(@single, GIS_data, 'UniformOutput', false);
    if flags.flag_human_instability > 0
        Human_Instability.slope = arcslope(DEM_raster,'degree');
        Human_Instability.slope = Human_Instability.slope.Z;
        Human_Instability = structfun(@single, Human_Instability, 'UniformOutput', false);
    end
    Hydro_States = structfun(@single, Hydro_States, 'UniformOutput', false);
    Inflow_Parameters = structfun(@single, Inflow_Parameters, 'UniformOutput', false);
    LULC_Properties = structfun(@single, LULC_Properties, 'UniformOutput', false);
    LULC_Properties.idx_lulc = logical(LULC_Properties.idx_lulc);
    LULC_Properties.idx_imp = logical(LULC_Properties.idx_imp);
    Rainfall_Parameters = structfun(@single, Rainfall_Parameters, 'UniformOutput', false);
    recording_parameters = structfun(@single, recording_parameters, 'UniformOutput', false);
    running_control = structfun(@single, running_control, 'UniformOutput', false);
    Soil_Properties = structfun(@single, Soil_Properties, 'UniformOutput', false);
    Soil_Properties.idx_soil = logical(Soil_Properties.idx_soil);
    if flags.flag_waterquality == 1
        WQ_States = structfun(@single, WQ_States, 'UniformOutput', false);
    end
    Wshed_Properties = structfun(@single, Wshed_Properties, 'UniformOutput', false);
    if flags.flag_reservoir == 1
        Reservoir_Data = structfun(@single, Reservoir_Data, 'UniformOutput', false);
    end
    if flags.flag_spatial_rainfall == 1
        Spatial_Rainfall_Parameters = structfun(@single, Spatial_Rainfall_Parameters, 'UniformOutput', false);
    end
    % ETP Data
    if flags.flag_ETP == 1
        extra_parameters.ETP.time_ETP_begin = ETP_Parameters.time_ETP_begin;
        ETP_Parameters.time_ETP_begin = [];
        extra_parameters.ETP.time_ETP = ETP_Parameters.time_ETP;
        extra_parameters.ETP.info = ETP_Parameters.info;
        ETP_Parameters.info = [];
        extra_parameters.ETP.time_ETP = ETP_Parameters.time_ETP;
        ETP_Parameters.time_ETP = [];
        ETP_Parameters = structfun(@single, ETP_Parameters, 'UniformOutput', false);
    end
    % Double Arrays
    % date_begin = gpuArray(date_begin);
    elevation = single(elevation);
    idx_nan = logical(idx_nan);
    if flags.flag_waterquality == 1
        idx_nan_5 = logical(idx_nan_5);
    else
        idx_nan_5 = [];
    end
    idx_outlet = logical(idx_outlet);
    k = single(k);
    nx = single(nx);
    ny = single(ny);
    if flags.flag_waterquality == 1
        Out_Conc = single(Out_Conc);
    end    
    outlet_index = single(outlet_index);
    outlet_runoff_volume = single(outlet_runoff_volume);
    outlet_type = single(outlet_type);
    slope_outlet = single(slope_outlet);
    spatial_domain = single(spatial_domain);
    t = single(t);
    t_previous = single(t_previous);
    time_calculation_routing = single(time_calculation_routing);
    time_step = single(time_step);
    time_step_model = single(time_step_model);
    tmin_wq = single(tmin_wq);
    C = single(C);
    min_soil_moisture = single(min_soil_moisture);
    k = single(1);

end

%% Calculating River Matrix to Estimate Groundwater
if flags.flag_groundwater_modeling == 1
    FD = FLOWobj(DEM_raster); % Flow direction
    As  = flowacc(FD); % Flow Accumulation
    Wshed_Properties.fac_area = As.Z*(Wshed_Properties.cell_area/1000/1000); % km2
    idx_rivers = Wshed_Properties.fac_area >= GIS_data.min_area;  % Logical Matrix with 1 being pixels with rivers   
else
    Lateral_Groundwater_Flux = 0; % m3/s/km of river
    idx_rivers = zeros(size(DEM_raster.Z));
end

%% Dam Break Parameters
% If dam breaking is simulated
if flags.flag_dam_break == 1
    input_table_breaker = readtable('Damns_breaks_points.xlsx'); 
    breakers.easting_absolute = table2array(input_table_breaker(:,2)); % Easting Coordinates
    breakers.northing_absolute = table2array(input_table_breaker(:,3)); % Northing Coordinates
    breakers.ID = table2array(input_table_breaker(:,1)); % ID
    % --- Converting coordinates to local coordinates in pixels
    breakers.easting= round((-GIS_data.xulcorner + breakers.easting_absolute)/Wshed_Properties.Resolution);
    breakers.northing = round((GIS_data.yulcorner - breakers.northing_absolute)/Wshed_Properties.Resolution);
    flag_break_1 = 1; flag_break_2 = 1;    
end

% We need to finish this part. Also, we need to fix it in the main while 

%% Clearing Variables

% clearvars  -except register saver_memory_maps idx_rivers rainfall_spatial_aggregation model_folder Input_Rainfall Reservoir_Data wse_slope_zeros Distance_Matrix depths Maps Spatial_Rainfall_Parameters GIS_data Inflow_Parameters ETP_Parameters Rainfall_Parameters CA_States BC_States Wshed_Properties Wshed_Properties Human_Instability gauges Hydro_States recording_parameters Courant_Parameters running_control Elevation_Properties inflow_volume idx_outlet outflow_volume outlet_runoff_volume I_t num_obs_gauges drainage_area northing_obs_gauges easting_obs_gauges depths time_record_hydrograph last_record_hydrograph initial_mass delta_p WQ_States routing_time flags LULC_Properties Soil_Properties topo_path idx_lulc idx_imp idx_soil d steps 	alfa_albedo_input 	alfa_max 	alfa_min 	alfa_save 	avgtemp_stations 	B_t   	C  	Cd 	cell_area 	climatologic_spatial_duration 	col_outlet 	coordinate_x 	coordinate_y 	coordinates_stations d_t  d_p 	date_begin  date_end	delta_p_agg  	DEM_etp 	DEM_raster 	depth_tolerance 	elevation    	ETP 	ETP_save 	factor_cells		flow_tolerance	flows_cells	G_stations	gravity	I_tot_end_cell	idx_nan	idx_nan_5	inflow	inflow_cells	k	k_out	Krs	ksat_fulldomain	last_record_maps	lat	mass_lost	mass_outlet	running_control.max_time_step	maxtemp_stations	min_time_step	mintemp_stations	mu	Inflow_Parameters.n_stream_gauges	nx	ny	Out_Conc	outlet_index	outlet_index_fulldomain	outlet_type	P_conc	psi_fulldomain	rainfall_matrix	rainfall_matrix_full_domain	Resolution	ro_water	roughness	roughness_fulldomain	row_outlet	slope_alfa	slope_outlet	spatial_domain	t	t_previous	teta_i_fulldomain	teta_sat	teta_sat_fulldomain	time_calculation_routing	time_change_matrices	time_change_records	time_deltap	time_ETP	time_records	time_save_previous	time_step	time_step_change	time_step_increments	time_step_model	time_step_save	tmin_wq	Tot_Washed	Tr	u2_stations	ur_stations	v_threshold	vel_down	vel_left	vel_right	vel_up	vol_outlet	weight_person	width1_person	width2_person
clearvars  -except register saver_memory_maps extra_parameters Lateral_Groundwater_Flux idx_rivers min_soil_moisture rainfall_spatial_aggregation model_folder Input_Rainfall Reservoir_Data wse_slope_zeros Distance_Matrix depths Maps Spatial_Rainfall_Parameters GIS_data Inflow_Parameters ETP_Parameters Rainfall_Parameters CA_States BC_States Wshed_Properties Wshed_Properties Human_Instability gauges Hydro_States recording_parameters Courant_Parameters running_control Elevation_Properties inflow_volume idx_outlet outflow_volume outlet_runoff_volume I_t num_obs_gauges drainage_area northing_obs_gauges easting_obs_gauges depths time_record_hydrograph last_record_hydrograph initial_mass delta_p WQ_States routing_time flags LULC_Properties Soil_Properties topo_path idx_lulc idx_imp idx_soil d steps 	alfa_albedo_input 	alfa_max 	alfa_min 	alfa_save 	avgtemp_stations 	B_t   	C  	Cd 	cell_area 	climatologic_spatial_duration 	col_outlet 	coordinate_x 	coordinate_y 	coordinates_stations d_t  d_p 	date_begin  date_end	delta_p_agg  	DEM_etp 	DEM_raster 	depth_tolerance 	elevation    	ETP 	ETP_save 	factor_cells		flow_tolerance	flows_cells	G_stations	gravity	I_tot_end_cell	idx_nan	idx_nan_5	inflow	inflow_cells	k	k_out	Krs	ksat_fulldomain	last_record_maps	lat	mass_lost	mass_outlet	running_control.max_time_step	maxtemp_stations	min_time_step	mintemp_stations	mu	Inflow_Parameters.n_stream_gauges	nx	ny	Out_Conc	outlet_index	outlet_index_fulldomain	outlet_type	P_conc	psi_fulldomain	rainfall_matrix	rainfall_matrix_full_domain	Resolution	ro_water	roughness	roughness_fulldomain	row_outlet	slope_alfa	slope_outlet	spatial_domain	t	t_previous	teta_i_fulldomain	teta_sat	teta_sat_fulldomain	time_calculation_routing	time_change_matrices	time_change_records	time_deltap	time_ETP	time_records	time_save_previous	time_step	time_step_change	time_step_increments	time_step_model	time_step_save	tmin_wq	Tot_Washed	Tr	u2_stations	ur_stations	v_threshold	vel_down	vel_left	vel_right	vel_up	vol_outlet	weight_person	width1_person	width2_person outflow_bates


%% Converting Arrays to GPU Arrays, if required
% Converting to GPU Arrays
if flags.flag_GPU == 1
    % Convert all Arrays to GPU Arrays
    % Structure Arrays
    BC_States = structfun(@gpuArray, BC_States, 'UniformOutput', false);
    CA_States = structfun(@gpuArray, CA_States, 'UniformOutput', false);
    Courant_Parameters = structfun(@gpuArray, Courant_Parameters, 'UniformOutput', false);
    depths = structfun(@gpuArray, depths, 'UniformOutput', false);
    Elevation_Properties = structfun(@gpuArray, Elevation_Properties, 'UniformOutput', false);
    flags = structfun(@gpuArray, flags, 'UniformOutput', false);
    % Gauges Label
    if flags.flag_obs_gauges == 1
        extra_parameters.gauges.label_observed_string = gauges.labels_observed_string;
        gauges.labels_observed_string = [];
        gauges = structfun(@gpuArray, gauges, 'UniformOutput', false);
    end
    GIS_data = structfun(@gpuArray, GIS_data, 'UniformOutput', false);
    if flags.flag_human_instability > 0
        Human_Instability.slope = arcslope(DEM_raster,'degree');
        Human_Instability.slope = Human_Instability.slope.Z;
        Human_Instability_text.list = Human_Instability.list;
        Human_Instability.list = [];
        Human_Instability_text.names = Human_Instability.names;
        Human_Instability.names = [];
        Human_Instability = structfun(@gpuArray, Human_Instability, 'UniformOutput', false);
    end
    Hydro_States = structfun(@gpuArray, Hydro_States, 'UniformOutput', false);
    Inflow_Parameters = structfun(@gpuArray, Inflow_Parameters, 'UniformOutput', false);
    LULC_Properties = structfun(@gpuArray, LULC_Properties, 'UniformOutput', false);
    Rainfall_Parameters = structfun(@gpuArray, Rainfall_Parameters, 'UniformOutput', false);
    recording_parameters = structfun(@gpuArray, recording_parameters, 'UniformOutput', false);
    running_control = structfun(@gpuArray, running_control, 'UniformOutput', false);
    Soil_Properties = structfun(@gpuArray, Soil_Properties, 'UniformOutput', false);
    if flags.flag_waterquality == 1
        WQ_States = structfun(@gpuArray, WQ_States, 'UniformOutput', false);
    else
        WQ_States = [];
    end
    Wshed_Properties = structfun(@gpuArray, Wshed_Properties, 'UniformOutput', false);
    if flags.flag_reservoir == 1
        Reservoir_Data = structfun(@gpuArray, Reservoir_Data, 'UniformOutput', false);
    end
    if flags.flag_spatial_rainfall == 1
        Spatial_Rainfall_Parameters = structfun(@gpuArray, Spatial_Rainfall_Parameters, 'UniformOutput', false);
    end
    % ETP Data
    if flags.flag_ETP == 1
        if flags.flag_single ~= 1
            extra_parameters.ETP.time_ETP_begin = ETP_Parameters.time_ETP_begin;
            ETP_Parameters.time_ETP_begin = [];
            extra_parameters.ETP.time_ETP = ETP_Parameters.time_ETP;
            ETP_Parameters.time_ETP = [];
            extra_parameters.ETP.info = ETP_Parameters.info;
            ETP_Parameters.info = [];
        end
        ETP_Parameters = structfun(@gpuArray, ETP_Parameters, 'UniformOutput', false);
    end
    % Double Arrays
    % date_begin = gpuArray(date_begin);
    elevation = gpuArray(elevation);
    idx_nan = gpuArray(idx_nan);
    if flags.flag_waterquality == 1
        idx_nan_5 = gpuArray(idx_nan_5);
    end
    idx_outlet = gpuArray(idx_outlet);
    k = gpuArray(k);
    nx = gpuArray(nx);
    ny = gpuArray(ny);
    Out_Conc = gpuArray(Out_Conc);
    outlet_index = gpuArray(outlet_index);
    outlet_runoff_volume = gpuArray(outlet_runoff_volume);
    outlet_type = gpuArray(outlet_type);
    slope_outlet = gpuArray(slope_outlet);
    spatial_domain = gpuArray(spatial_domain);
    t = gpuArray(t);
    t_previous = gpuArray(t_previous);
    time_calculation_routing = gpuArray(time_calculation_routing);
    time_step = gpuArray(time_step);
    time_step_model = gpuArray(time_step_model);
    tmin_wq = gpuArray(tmin_wq);
    C = gpuArray(C);
    min_soil_moisture = gpuArray(min_soil_moisture);
    k = gpuArray(1);
end


clear spatial_domain





