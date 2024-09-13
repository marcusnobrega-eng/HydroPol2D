%% 2-D Overland Flow Model
% % % % % % % % % % % % % Model Status % % % % % % % % % % % % %
% ---------- Version 3.0 -------------
% Last Update - 9/21/2022
% Next step: Read TIFF Files
%            Read Spatio-Time-Varying Rainfall

% ----------------- Initial Data - General -----------------
clc
clear all
tic

%% Calibration Results
observed_concentration = [28607.59494	27468.35443	25822.78481	24303.79747	21518.98734	17848.10127	15443.03797	13037.97468	11772.1519	11265.82278	10759.49367	10253.16456	9240.506329	7848.101266	6708.860759	5443.037975	4303.797468	3670.886076	2658.227848	1898.734177	1645.56962	1139.240506	1139.240506]'; % mg/L
time_observed_concentration = [0.956989247	1.204301075	1.47311828	1.720430108	1.967741935	2.204301075	2.462365591	2.720430108	2.967741935	3.47311828	3.978494624	4.47311828	4.978494624	5.483870968	5.989247312	6.494623656	7	7.483870968	8.010752688	8.516129032	9	9.505376344	9.989247312]'; % min
C3_test = [960000:100000:1560000];
C4_test = [1.2:0.01:1.5]';
n_tests = length(C3_test)*length(C4_test);
%% Convert Rasters and Pre-Processing

input_table = readtable('general_data.xlsx');
% Load TopoToolBox Tools
topo_path = table2cell(input_table(1,31));
addpath(genpath(char(topo_path)));

% Read Plane Watershed Data
% [~,~,~] = plane_watershed(0.02,0.01,1.48,2.96,0);
[~,~,~] = plane_watershed(0.15,0.01,1.5,3,0);


% Rasters
% fname_LULC = 'LULC.tif'; fname_DEM = 'DEM.tif';
% fname_SOIL = 'SOIL.tif';

fname_LULC = 'LULC_Plane.asc'; fname_DEM = 'DEM_Plane.asc';
fname_SOIL = 'SOIL_Plane.asc';


LULC_raster = GRIDobj(fname_LULC); % Land Use and Land Cover Classification
DEM_raster = GRIDobj(fname_DEM); % Digital Elevation Model (m)
SOIL_raster = GRIDobj(fname_SOIL); % Soil Map

% Extent Problem
if sum(size(DEM_raster.Z)) > sum(size(LULC_raster.Z)) && sum(size(DEM_raster.Z)) > sum(size(SOIL_raster.Z)) % DEM is larger
    raster_resample = DEM_raster;
    % Resample other two rasters
    % ---- Constraint at LULC Raster
    LULC_raster = resample(LULC_raster,raster_resample);
    LULC_raster.Z = round(LULC_raster.Z); % Only Integers
    % ---- Constraint at SOIL Raster
    SOIL_raster = resample(SOIL_raster,raster_resample);
    SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
end

if sum(size(SOIL_raster.Z)) > sum(size(DEM_raster.Z)) && sum(size(SOIL_raster.Z)) > sum(size(LULC_raster.Z))  % SOIL is larger
    raster_resample = SOIL_raster;
    % Resample other two rasters
    LULC_raster = resample(LULC_raster,raster_resample);
    % ---- Constraint at LULC Raster
    LULC_raster = resample(LULC_raster,raster_resample);
    LULC_raster.Z = round(LULC_raster.Z); % Only Integers
    DEM_raster = resample(DEM_raster,raster_resample);
end

if sum(size(LULC_raster.Z)) > sum(size(DEM_raster.Z)) && sum(size(LULC_raster.Z)) > sum(size(SOIL_raster.Z))  % SOIL is larger
    raster_resample = DEM_raster;
    % Resample other two rasters
    % ---- Constraint at SOIL Raster
    SOIL_raster = resample(SOIL_raster,raster_resample);
    SOIL_raster.Z =  round(SOIL_raster.Z); % Only Integers
    DEM_raster = resample(DEM_raster,raster_resample);
end

% Raster Extent
xulcorner = DEM_raster.refmat(3,1); % Up Left Corner
yulcorner = DEM_raster.refmat(3,2);
% - Extent is already solved, we can login the input data
Resolution = DEM_raster.cellsize; % m

%% Drainage Basins


% %% Fillsinks
% % max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
% % DEM_filled = fillsinks(DEM,max_depth);
%     DEM_filled = fillsinks(DEM_raster); % Filled DEM
%     DIFFDEM = DEM_filled - DEM_raster;
%     DIFFDEM.Z(DIFFDEM.Z==0) = nan;
%     DEM_raster = DEM_filled;
%     imageschs(DEM_raster,DIFFDEM.Z);

FD = FLOWobj(DEM_raster);

% area_km2 = 500; % km2
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

%% Load General Data Input
input_data_script;  % Load general data, soil, and LULC parameters

%% Resample DEM, if Required
% In case we want to resample the DEM
resolution = resolution_resample; % m
if flag_resample == 1
    % DEM
    DEM_raster = resample(DEM_raster,resolution);
    % LULC
    LULC_raster = resample(LULC_raster,resolution);
    % SOIL
    SOIL_raster = resample(SOIL_raster,resolution);
end

% ----- Transforming Raster into Matrix with Values ----- %

LULC = double(LULC_raster.Z);
DEM = double(DEM_raster.Z);
SOIL = double(SOIL_raster.Z);

neg_DEM = DEM <= 0;
neg_LULC = LULC < 0;
neg_SOIL = SOIL < 0;
inf_nan_MAPS = isinf(DEM) + isnan(DEM) + neg_DEM + isnan(LULC) + isnan(SOIL) + neg_LULC + neg_SOIL + isinf(LULC) + isinf(SOIL); % Logical array
idx = inf_nan_MAPS > 0;

% Rebuilding Rasters to the Lowest Extent
LULC_raster.Z = LULC;% Land Use and Land Cover Classification
DEM_raster.Z = DEM; % Digital Elevation Model
SOIL_raster.Z = SOIL; % Soil Map

% Replacing Values with Issues
DEM_raster.Z(idx) = nan;

dem = DEM; % Further used in the elevation data
imp = LULC;
soil = SOIL;

cell_area = Resolution^2; % cell area in square meters
zzz = size(dem); % Dimensions of DEM matrix

% ---- DEM Dimensions --- %
[ny_max,nx_max] = size(dem);

% ------------ Inflow Cells  ------------ %
% Here we read where the stream gauges are located
inflow_cells = zeros(ny_max,nx_max,n_stream_gauges);
if flag_inflow == 1
    for i = 1:n_stream_gauges
        for z = 1:n_inlets(i)
            x = easting_inlet_cells(z,i);
            y = northing_inlet_cells(z,i);
            inflow_cells(y,x,i) = 1; % Coordinates of each stream gauge
        end
    end
end

% ------------ Rainfall Matrix ------------ %
if flag_rainfall == 0 % No rainfall
    rainfall_matrix = flag_rainfall*zeros(size(dem));
elseif flag_rainfall == 1 && flag_spatial_rainfall == 1
    % Spatial Rainfall Case
    input_table = readtable('Rainfall_Spatial_Input.xlsx');
    % Observations
    n_obs = sum((table2array(input_table(:,2))>=0)); % Number of observations
    n_max_raingauges = 50;
    time_step_spatial = table2array(input_table(7,2)) - table2array(input_table(6,2)); % min
    end_rain = (n_obs-1)*time_step_spatial;
    rainfall_spatial_duration = 0:time_step_spatial:(end_rain); % Rainfall data time in minutes

    % Rainfall Data
    for i = 1:n_max_raingauges
        rainfall_raingauges(:,i) = table2array(input_table(6:end,3 + (i-1)));
        coordinates(i,1) = table2array(input_table(3,3 + (i-1)));
        coordinates(i,2) = table2array(input_table(4,3 + (i-1)));
    end
    n_raingauges = sum(rainfall_raingauges(1,:) >= 0); % Number of raingauges
elseif flag_rainfall == 1 && flag_spatial_rainfall ~= 1
    % Lumped Rainfall Case
    rainfall_matrix = flag_rainfall*ones(size(dem));
end


%% ------------ ETP Matrices ------------
if flag_ETP == 1
    % We are running ETP
    Krs = 0.16; % Parameter
    alfa_albedo_input = 0.23; % Parameter
    input_table = readtable('ETP_input_data.xlsx');
    % Observations
    n_obs_ETP = sum((table2array(input_table(:,3))>=0)); % Number of observations
    n_max_etp_stations = 50;
    time_step_etp = minutes(table2array(input_table(4,2)) - table2array(input_table(3,2))); % min
    end_etp = (n_obs_ETP-1)*time_step_etp;
    time_ETP = table2array(input_table(3:end,2));
    time_ETP_begin = time_ETP(1);

    % Check Initial ETP Dates
    delta_ETP_date = minutes(time_ETP_begin - date_begin); % minutes positive or negative
    climatologic_spatial_duration = 0:time_step_etp:(end_etp) ; % Rainfall data time in minutes
    climatologic_spatial_duration = climatologic_spatial_duration + + double(delta_ETP_date);

    % Preallocating Arrays
    ETP_save = zeros(size(dem,1),size(dem,2),n_obs_ETP); % Maps of ET

    % ---- Loop ETP --- %
    for i = 1:n_max_etp_stations
        try
            n_stations = 50;
            % Maximum Temperature
            maxtemp_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 3));
            % Minimum Temperature
            mintemp_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 4));
            % Average Temperature
            avgtemp_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 5));
            % U2
            u2_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 6));
            % UR
            ur_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 7));
            % G
            G_stations(:,i) = table2array(input_table(3:end,6*(i-1) + 8));

            % Coordinates
            coordinates_stations(i,1) = table2array(input_table(1,6*(i-1) + 6));
            coordinates_stations(i,2) = table2array(input_table(1,6*(i-1) + 8));
        catch
            n_stations = i-1;
            break
        end
        if i == n_max_etp_stations
            n_stations = n_max_etp_stations;
        end
    end

    % DEM Raster Information for ETP Calcualtion

    [DEM_etp,R_etp] = readgeoraster(fname_DEM); % Getting Raster Information
    info = geotiffinfo(fname_DEM);
    height = info.Height; % Integer indicating the height of the image in pixels
    width = info.Width; % Integer indicating the width of the image in pixels
    [cols_etp,rows_etp] = meshgrid(1:width,1:height);
    [x_etp,y_etp] = pix2map(info.RefMatrix, rows_etp, cols_etp); % Map Coordinates
    [lat,lon] = projinv(info, x_etp,y_etp); % Latitude and Longitude
    neg_DEM = DEM_etp <= 0;
    DEM_etp(neg_DEM) = nan;
    idx_cells = DEM_etp >= 0;
    idx_cells = double(idx_cells(:));
    lat(neg_DEM) = nan; % Latitude
    lon(neg_DEM) = nan; % Longitude

end
%% ------------ Recording time of outputs (i.e., flows, concentrations ...) ------------
% Calculations
steps = routing_time/time_step_model; % number of calculation steps
number_of_records = floor(steps*time_step_model/record_time_maps); % number of stored data (size of the vector)
% Checking if the recording time is working
if number_of_records == 0
    error('The recording time is larger than the routing time, please change it')
end
time_records = (0:record_time_maps:routing_time); % time in minutes
time_record_hydrograph = (0:record_time_hydrographs:routing_time); % time in minutes
time_change_records = (0:time_step_change/60:routing_time); % time in minutes
time_change_matrices = (0:time_step_matrices/60:routing_time); % time in minutes
% vector to store data
time_store = time_records./time_step_model; % number of steps necessary to reach the record vector
time_store(1) = 1; % the zero is the firt time step
time_step_save = zeros(steps,1);
flow_acceleration = zeros(steps,1);
alfa_save = zeros(steps,1);
%% Inflow and Precipitation data
inflow_hydrograph_rate = inflow_hydrograph_rate*flag_inflow;
if flag_rainfall == 1 && flag_spatial_rainfall ~=1 % Only for concentrated rainfall
    intensity_rainfall = intensity_rainfall*flag_rainfall; % If flag_rainfall is zero, no rainfall is considered
end

%% Fillsinks
% max_depth = 0.1; % meters, positive value. If you don't want to use, delete it from the function
% DEM_filled = fillsinks(DEM,max_depth);
if flag_fill_DEM == 1
    DEM_filled = fillsinks(DEM_raster); % Filled DEM
    DIFFDEM = DEM_filled - DEM;
    DIFFDEM.Z(DIFFDEM.Z==0) = nan;
    DEM_raster = DEM_filled;
    % imageschs(DEM_raster,DIFFDEM.Z);
end

%% DEM Smoothening
if flag_smoothening == 1
    [DEM_raster,DEM,S] = DEM_smoothening(DEM_raster,min_area,flag_trunk,tau,K_value);
end

%% Imposing Minimum Slope - If Required
if flag_diffusive ~= 1
    % Impose Mininum Slope
    if flag_smoothening ~=1 % We don't have a S, so we need to calculate it
        FD = FLOWobj(DEM_raster);
        area_km2 = min_area; % km2
        area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
        S = STREAMobj(FD,'minarea',area_cells); % Flow accumulation
    end
    DEM_new = imposemin(S,DEM_raster,sl);
    DEM_DIFF = DEM_new - DEM_raster;
    DEM_raster = DEM_new;
    %     imagesc(DEM_DIFF); colorbar; % if you want to plot
end

%% Slope Map
% SLP = arcslope(DEM_raster);
% imagesc(SLP); colorbar

%% Smooth DEM
if flag_smooth_cells == 1
    Vq = imgaussfilt(dem,'FilterSize',3);
    dem_diff_smooth = Vq;
    %     subplot(2,1,1)
    %     surf(dem); view(0,90); shading interp
    %     subplot(2,1,2)
    %     surf(Vq); view(0,90); shading interp
    dem = Vq; % New dem
    pause(2)
    close all
end

%% Decrese Elevations in Creeks
if flag_reduce_DEM
    FD = FLOWobj(DEM_raster); % Flow direction
    As  = flowacc(FD); % Flow Accumulation
    fac_area = As.Z*(cell_area/1000/1000); % km2
    idx_facc = fac_area <= min_area;

    % B
    B = beta_1*fac_area.^(beta_2);
    H = alfa_1*fac_area.^(alfa_2);
    Flow_Area = B.*H; % m2
    H_abg = Flow_Area/Resolution;
    H_abg(idx_facc) = 0;

    % Reducing Elevation in creeks
    DEM_raster.Z = DEM_raster.Z - H_abg; % [m]

    % New Data
    dem = DEM_raster.Z;
end

%% Drainage Basins
D = drainagebasins(FD);
imageschs(DEM_raster,shufflelabel(D))
area_km2 = min_area;
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
outlet_index_fulldomain = zeros(zzz);
% The following for-loop is used when we enter the boundaries of the domain
% as the outlet
% if sum(sum(isnan(dem))) + sum(sum(isinf(dem))) + sum(sum(dem(dem<0))) > 0
%     for i = x_outlet_begin:x_outlet_end
%         for j = y_outlet_begin:y_outlet_end
%             outlet_index_fulldomain(j,i) = 1;
%         end
%     end
% end

elevation_nan = dem;
elevation_nan(elevation_nan < 0) = nan; % Replacing no-info to nan

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
max_length = (max(col_boundary) - min(col_boundary))*Resolution - (max(row_boundary) - min(row_boundary))*Resolution; % Length
max_width = (max(row_boundary) - min(row_boundary))*Resolution - (max(row_boundary) - min(row_boundary))*Resolution; % width
avg_river_length = sqrt(max_length^2 + max_width^2); % Average river length

% Number of Outlets
idx_outlet = elevation_nan == min(min(elevation_nan)); % Finding cells with lowest elevations
[row_outlet,col_outlet] = find(idx_outlet == 1); % finding rows and cols of these cells
outlet_index = outlet_index_fulldomain;

% New way: Assuming that the outlets should be at the border of the domain, we can
% do:
for i = 1:length(row_outlet)
    % Checking left, right up, and down
    row = row_outlet(i); % Row
    col = col_outlet(i); % Col
    % Left
    if  row - 1 == 0 || isnan(elevation_nan(row-1,col))
        outlet_index(row_outlet(i),col_outlet(i)) = 1; % This is an outlet
    end
    % Right
    if  row + 1 > size(elevation_nan,1) || isnan(elevation_nan(row+1,col))
        outlet_index(row_outlet(i),col_outlet(i)) = 1; % This is an outlet
    end
    % Up
    if col - 1 == 0 || isnan(elevation_nan(row,col-1))
        outlet_index(row_outlet(i),col_outlet(i)) = 1; % This is an outlet
    end
    % Down
    if  col + 1 > size(elevation_nan,2) || isnan(elevation_nan(row,col+1))
        outlet_index(row_outlet(i),col_outlet(i)) = 1; % This is an outlet
    end
end
[row_outlet,col_outlet] = find(outlet_index == 1);  % Refreshing outlet

% Minimum Row
min_row_outlet = min(min(row_outlet));
% Moving Up
pos = find(row_outlet == min_row_outlet);
col = col_outlet(pos);
row = row_outlet(pos);

% 0 , 45, 90, 135, 180, 225, 270, 315
check = zeros(8,1);
i = 1;
k = 1;
outlet_new = zeros(size(dem));
while i <= round(n_outlets_data/2)
    % 0 , 45, 90, 135
    if col + 1 <= size(dem,2) && perimeter(row,col + 1) > 0  % 0
        col = col+1;
        k = k + 1;
    elseif row - 1 >= 1 && col + 1 <= size(dem,2) && perimeter(row - 1,col + 1) > 0  % 45
        row = row-1;
        col = col+1;
        k = k + 1;
    elseif row - 1 >= 1 && perimeter(row-1,col) > 0  % 90
        row = row-1;
        k = k + 1;
    elseif row - 1 >= 1 && col - 1 >= 1 && perimeter(row-1,col-1) > 0  % 135
        row = row-1;
        col = col-1;
        k = k + 1;
    end
    outlet_new(row,col) = 1;
    i = i + 1;
end
% Maximum Row
max_row_outlet = max(max(row_outlet));
pos = find(row_outlet == max_row_outlet);
col = col_outlet(pos);
row = row_outlet(pos);
i = 1;
while i <= n_outlets_data - k
    % 180, 225, 270, 315
    if  col - 1 >= 1 && perimeter(row,col-1) > 0  % 180
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
    i = i + 1;
end

outlet_index = outlet_index + outlet_new; % Refreshing Outlet
[row_outlet,col_outlet] = find(outlet_index == 1); % finding rows and cols of these cells
el_outlet = dem;
el_outlet(idx_outlet < 1) = nan;
stage_min = min(min(el_outlet));

[row_min, col_min] = find(dem == stage_min);

%% Grid Domain
%%%%%% ORIGINAL GRID %%%%%%
% In case flag_abstraction == 1, we cut the domain assuming the cells
% entered in the input_data file
if flag_abstraction ~= 1
    xmin = 1; % initial position x in the grid (collums)
    ymin = 1; % lines (up to down)
    xmax = zzz(2);
    ymax = zzz(1);
end
% ------------ Cutting Cells ------------
inflow_cells = inflow_cells(ymin:ymax,xmin:xmax,:);
outlet_index = outlet_index(ymin:ymax,xmin:xmax);
elevation = dem(ymin:ymax,xmin:xmax); % using only specified grid
spatial_domain = zeros(size(dem));
dem = dem(ymin:ymax,xmin:xmax);
lulc_matrix = imp(ymin:ymax,xmin:xmax); % using only specified grid
soil_matrix  = soil(ymin:ymax,xmin:xmax); % using only specified grid;
rainfall_matrix = double(elevation >= 0) ; % Elevation could be 0
rainfall_matrix(rainfall_matrix == 0) = nan;

% ------------ GRID ------------- %
x_grid = xulcorner + Resolution*[1:1:size(dem,2)];
y_grid = yulcorner - Resolution*[1:1:size(dem,1)];

% Lower Left Corner
xllcorner = xulcorner; % Same
yllcorner = yulcorner - Resolution*(size(dem,1)); % Subtract Y distance

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
idx_nan_5 = zeros(zzz(1),zzz(2),5);
for i = 1:5 % Left, right, up, down, outlet
    idx_nan_5(:,:,i) = idx_nan; % This is used in CA model
end
idx_nan_5 = logical(idx_nan_5); % Converting to logical

elevation(idx) = nan; % change those values for nan
lulc_matrix(idx) = nan; % change those values for nan
soil_matrix(idx) = nan; % change those values for nan

spatial_domain(idx) = nan; % change those values for inf
% ------------- Preallocating arrays to avoid large computational efforts
y_1 = length(ymin:ymax);
x_1 = length(xmin:xmax);
roughness = spatial_domain;
d_0 = spatial_domain; h_0 = spatial_domain ;psi = spatial_domain;
ksat = spatial_domain; teta_i = spatial_domain;
teta_sat = spatial_domain; I_0 = spatial_domain;
if flag_waterquality == 1 % Build-up and Wash-off matrices
    C_1 = spatial_domain;
    C_2  = spatial_domain;
    C_3 = spatial_domain;
    C_4 = spatial_domain;
    B_0 = spatial_domain;
end
%% ----------------- Fill Variables ----------------- %%
% Correcting LULC_index
lulc_matrix(idx_nan) = -9999;
idx_lulc = zeros(size(DEM,1),size(dem,2),n_lulc);
% Assuming that we know where is the impervious areas
for i = 1:n_lulc
    index = lulc_matrix == LULC_index(i,1);
    idx_lulc(:,:,i) = index;
end
idx_imp = idx_lulc(:,:,imp_index);
% idx = cat(3,idx_1,idx_2,idx_3,idx_4,idx_5,idx_6); % Concatenating all of them
impervious_cells = sum(sum(idx_imp));
pervious_cells = sum(sum(sum(idx_lulc))) - impervious_cells;
% Converting to Logical Values
idx_lulc = logical(idx_lulc);
idx_imp = logical(idx_imp);
for i = 1:n_lulc % 8 types of LULC
    % Only Roughness and h_0
    roughness(idx_lulc(:,:,i)) = lulc_parameters(i,1); % assigning values for roughness at impervious areas
    h_0(idx_lulc(:,:,i)) = lulc_parameters(i,2); % Initial Abstraction

    % ------------ Warm-up Data ------------:
    if flag_warmup == 1
        % Identifying positive values of the Warmup
        Warmup_Raster = GRIDobj('Warmup_Depths.asc');
        Warmup_Depths = Warmup_Raster.Z;
        idx_warmup = Warmup_Depths >= 0;
        d_0 =  idx_warmup.*Warmup_Depths; % Depths in m
        % Treat inf values
        d_0(idx_nan) = nan; % carefull here
        % Since d_0 is in meters, we have to convert it to mm
        d_0 = d_0*1000; % That means we have to input d_0 in meters
    else
        d_0(idx_lulc(:,:,i)) = lulc_parameters(i,3);
    end
    if flag_waterquality == 1 % Only if water quality is being modeled
        C_1(idx_lulc(:,:,i)) =  lulc_parameters(i,4);
        C_2(idx_lulc(:,:,i)) =  lulc_parameters(i,5);
        C_3(idx_lulc(:,:,i)) =  lulc_parameters(i,6);
        C_4(idx_lulc(:,:,i)) =  lulc_parameters(i,7);
    end
    if i == n_lulc % Last land use
        % Do we add another mass in the initial buildup or not?
        flag_mass_sum = 1;
        if flag_initial_buildup == 1 && flag_mass_sum == 1 && flag_waterquality == 1
            load Initial_Buildup.txt
            B_0 = Initial_Buildup; % kg/cell
            % Treat inf values
            B_0(idx_nan) = inf; % carefull here
            % Another filter
            B_0(B_0 < 0) = 0; % Attention here
            B_0_add = C_1.*(1 - exp(1).^(-C_2*ADD)); % kg/ha
            B_t = B_0 + B_0_add*cell_area/10^4; % kg
            initial_buildup_map = B_t; % kg

        elseif flag_initial_buildup == 1 && flag_mass_sum ~= 1 && flag_waterquality == 1
            load Initial_Buildup.txt
            B_0 = Initial_Buildup; % kg/cell
            % Treat inf values
            B_0(idx_nan) = inf; % carefull here
            % Another filter
            B_0(B_0 < 0) = 0; % Attention here
            B_t = B_0;
            initial_buildup_map = B_t; % kg

        elseif flag_waterquality == 1 % Let's calculate iT
            % Calculate Build-up using C1 and C2
            B_0 = C_1.*(1 - exp(1).^(-C_2*ADD)); % kg/ha
            B_t = B_0*cell_area/10^4; % kg
            % Entering the Value of B1 manually
            B_0 = 125; % g/m2
            B_t = (B_0/1000)*cell_area.*idx_not_nan; % kg in each pixel
            initial_buildup_map = B_t; % kg
        end
    end
end

idx_soil = zeros(size(dem,1),size(dem,2),n_soil);
for i = 1:n_soil
    index = soil_matrix == SOIL_index(i,1);
    idx_soil(:,:,i) = index;
end

idx_soil = logical(idx_soil);

% Soil Parameter Assigning
for i = 1:n_soil
    ksat(idx_soil(:,:,i)) = soil_parameters(i,1); % similar things happening below
    psi(idx_soil(:,:,i)) = soil_parameters(i,2);
    I_0(idx_soil(:,:,i)) = soil_parameters(i,3);
    teta_sat(idx_soil(:,:,i)) = soil_parameters(i,4);
    teta_i(idx_soil(:,:,i)) = soil_parameters(i,5);
end

% Mask in Impervious Areas
ksat(idx_imp) = 0;

% Replenishing Coefficient
kr = (1/75*(sqrt(ksat/25.4))); % Replenishing rate (1/hr) (Check Rossman, pg. 110)
Tr = 4.5./sqrt(ksat/25.4); %  Recovery time (hr) Check rossman pg. 110
Lu = 4.*sqrt(ksat/25.4); % Inches - Uppermost layer of the soil
Lu = Lu*2.54/100; % Meters
k_out = (teta_sat - teta_i).*kr.*Lu*1000; % Rate of replenishing exfiltration from the saturated zone during recoverying times (mm/hr)
% clear idx % clear this matrices to avoid huge storage data
% ----------------- Full Domain Cells -----------------
rainfall_matrix_full_domain = rainfall_matrix;
teta_sat_fulldomain = teta_sat;
teta_i_fulldomain = teta_i;
psi_fulldomain = psi;
roughness_fulldomain = roughness;
ksat_fulldomain = ksat;
h_0_fulldomain = h_0;
I_p_fulldomain = I_0;
I_t_fulldomain = I_0;
outlet_index_fulldomain = outlet_index_fulldomain(ymin:ymax,xmin:xmax);
%%%% New xmax and ymax
xmin = 1; %initial position x in the grid (collums)
ymin = 1; % lines (up to down)
zzz = size(elevation);
xmax = zzz(2);
ymax = zzz(1);
%% Conversion of inflow into the time-step of calculations
if flag_inflow == 1
    inflow_length = size(inflow_hydrograph_rate,1)- 1; % Number of interval
    inflow_discretized = zeros(size(inflow_hydrograph_rate,2),ceil(inflow_length*time_step_inflow/time_step_model)); % Preallocating
    for z = 1:n_stream_gauges
        for i =1:((inflow_length-1)*time_step_inflow/time_step_model)
            inflow_discretized(z,i) = inflow_hydrograph_rate(ceil((round(i*time_step_model/time_step_inflow,12))),z); % Discretized into moldel's time-step
        end
    end
end
% Each z row represent an stream gauge and each i col is a inflow value
%% Conversion of rainfall into the time-step of calculations for concentrated rainfall
if flag_rainfall == 1 && flag_spatial_rainfall ~=1  % Only for concentrated rainfall
    intensity_rainfall_length = length(intensity_rainfall) - 1; % Number of intervals
    intensity_discretized = zeros(1,ceil(intensity_rainfall_length*time_step_rainfall/time_step_model)); % Preallocating
    for i =1:(intensity_rainfall_length*time_step_rainfall/time_step_model)
        intensity_discretized(i) = intensity_rainfall(ceil((round(i*time_step_model/time_step_rainfall,12))));  % Discretized into moldel's time-step
    end
end
%% Determination of grid parameters and outlet coordinates (Whole Domain)
nx_max = length(elevation(1,:)); % number of x cells
ny_max = length(elevation(:,1)); % number of y cells

% Calculates the number of non-inf cells to determine the watershed area
matrix_nan = isnan(elevation);
number_nan = sum(sum(matrix_nan));
drainage_area = (nx_max*ny_max - number_nan)*Resolution^2; % m2
watershed_perimeter = sum(sum(perimeter))*12/1000; % kilometers
impervious_area = sum(sum(impervious_cells))*cell_area;
pervious_area = sum(sum(pervious_cells))*cell_area;
impervious_rate = impervious_area/drainage_area;
pervious_rate = pervious_area/drainage_area;

% Geometrical Properties
compactness_coefficient = 0.28*watershed_perimeter*1000/sqrt(drainage_area);
form_factor = drainage_area/(avg_river_length^2); % A / L^2
circularity_index = 12.57*drainage_area/((watershed_perimeter*1000)^2); % more close to 1, closer to a circle
% Pollutant Mass
if flag_waterquality == 1
    initial_mass = sum(sum(B_t(~isinf(B_t)))); % kg of pollutant
end
%% Calculation of accumulated and incremental inflow on the inflow cells
if flag_inflow == 1
    [~,delta_inflow,inflow_intensity] = accumulated_incremental(steps,inflow_discretized,time_step_model);
    time_deltainflow = cumsum(ones(size(delta_inflow,2),1)*time_step_inflow); % Vector with indices
end
%% Calculation of accumulated and incremental precipitation in each cell
if flag_rainfall == 1 && flag_spatial_rainfall ~=1  % Only for concentrated rainfall
    [~,delta_p,rainfall_intensity] = accumulated_incremental(steps,intensity_discretized,time_step_model);
    time_deltap = cumsum(ones(1,length(delta_p))*time_step_model);
    outlet_runoff_volume = zeros(size(delta_p));
end
%% Pre allocating more arrays
Total_Inflow = 0;
d_t = spatial_domain;
outflow_cms_t = spatial_domain;
inf_m3s_t = spatial_domain;
q_exut_t = spatial_domain;
% ----------------- Space dependent arrays -----------------
qout_left_t = spatial_domain;
qout_right_t = spatial_domain;
qout_up_t = spatial_domain;
qout_down_t = spatial_domain;
d_p = spatial_domain;
d_avg = spatial_domain;
qin_t = spatial_domain;
vel_left = spatial_domain;
vel_right = spatial_domain;
vel_up = spatial_domain;
vel_down = spatial_domain;
I_tot_end = zeros(size(spatial_domain)); %%% ATTENTION HERE
% inundated_cells = spatial_domain;
I_tot_end_cell_fulldomain = spatial_domain;
% ----------------- Time dependent arrays -----------------
time_size = length(time_store);
EMC_outlet = zeros(time_size,1);
d = zeros(ny_max,nx_max,time_size);
risk = zeros(ny_max,nx_max,time_size);
outet_hydrograph = zeros(time_size,1);
time_hydrograph = zeros(time_size,1);
if flag_waterquality == 1
    Pol_Conc_Map = zeros(ny_max,nx_max,time_size);
    Pol_mass_map = zeros(ny_max,nx_max,time_size);
    outet_pollutograph = zeros(time_size,1);
end
%% Clearing a few variables
clear accum_precipitation precipitation imp Land_Cover_Data Elevation_DATA intensity_discretized idx_ idx_1 matrix_nan idx accum_inflow col col_check d_0_imp d_0_per h_0_imp h_0_per I_0_per I_0_imp inflow_discretized

%% ---------------------- Main Routing (GA + 4D/8D CA + BW models) ------------------%%
tic % Start counting time
% ----------------- Find Initial Domain -----------------
[y_depth_stored, x_depth_stored] = find(d_0 > 0);
[row_check, col_check] = find(inflow_cells > 0);
[row_check_rainfall, col_check_rainfall] = find(rainfall_matrix > 0);
cells_check = [y_depth_stored x_depth_stored; row_check col_check; row_check_rainfall, col_check_rainfall];
cells_check = unique(cells_check,'rows');
flows_cells = zeros(ny_max,nx_max); %%% ATTENTION HERE
% ----------------- Initialize Variables -----------------
I_t = I_0;
I_p = I_0;
time_calculation_routing = zeros(steps,1);
k = 1; % Start Counter of the While Loop
actual_record_state = 1;
last_record_maps = 1;
last_record_hydrograph = 1;
if flag_inflow == 1
    for i = 1:n_stream_gauges
        delta_inflow_agg(i) = delta_inflow(i,1);
    end
end
t = 0;
if flag_rainfall == 1 && flag_spatial_rainfall ~=1
    delta_p_agg = delta_p(1,1)*flag_rainfall;
end

if flag_rainfall ~= 1
    delta_p_agg = 0;
end

if flag_rainfall == 1 && flag_spatial_rainfall ==1
    spatial_rainfall_maps = zeros(size(dem,1),size(dem,2),size(d,3)); % Maps of Rainfall
    % Spatial Rainfall
    if t == 0
        z = 1;
    else
        z = find(t <= rainfall_spatial_duration,1,'last'); % Duration
    end
    rainfall = rainfall_raingauges(z,1:n_raingauges)'; % Values of rainfall at t for each rain gauge
    x_grid = xulcorner + Resolution*[1:1:size(DEM_raster.Z,2)]'; % Pixel easting coordinates
    y_grid = yulcorner - Resolution*[1:1:size(DEM_raster.Z,1)]'; % Pixel northing coordinates
    x_coordinate = coordinates(1:n_raingauges,1); % Coordinates (easting) of each rain gauge
    y_coordinate = coordinates(1:n_raingauges,2); % Coordinates (northing) of each rain gauge
    [spatial_rainfall] = Rainfall_Interpolator(x_coordinate,y_coordinate,rainfall,x_grid,y_grid); % Interpolated Values
    spatial_rainfall_maps(:,:,1) = spatial_rainfall;
    delta_p_agg = spatial_rainfall/3600*time_step_model*60; % Matrix of delta P for each pixel
end
time_step = time_step_model;
I = zeros(4,1); % Volumes in each direction
t = time_step_model; % minutes
time_save_previous = 0; % minutes
Previous_Volume = 0;
t_previous = 0;
coordinate_x = 1;
coordinate_y = 1;
time_matrices_previous = 0;
step_error = zeros(1,steps);
time_step_cell = max_time_step*ones(1,steps);
relative_vol_error = zeros(1,steps);
I_cell = zeros(ny_max,nx_max,5);
C = 0;
P_conc = 0;
tmin_wq = max_time_step;
dmax_final = zeros(size(elevation));
mass_lost = 0; % initializing
mass_outlet = 0; % initializing
Out_Conc = 0;
vol_outlet = 0;
outlet_runoff_volume = 0; % initializing
Tot_Washed = spatial_domain;
outflow_volume = 0;
inflow_volume = 0;
I_tot_end_cell = zeros(size(spatial_domain)); % Attention here
inflow = zeros(size(dem,1),size(dem,2));
ETP = 0;

%%%% ELEVATIONS %%%
elevation_left_t = [zeros(ny_max,1),elevation(:,1:(nx_max-1))];
elevation_right_t = [elevation(:,(2:(nx_max))) zeros(ny_max,1)];
elevation_up_t = [zeros(1,nx_max) ; elevation(1:(ny_max-1),:)];
elevation_down_t = [elevation(2:(ny_max),:) ; zeros(1,nx_max)];
i_save = 1;
data_water_quality = zeros(n_tests,length(observed_concentration) + 3);
%%
for nc3 = 1:length(C3_test)
    for nc4 = 1:length(C4_test)
        % Values Changed for Calibration
        C_3 = double(idx_not_nan)*C3_test(nc3);
        C_4 = double(idx_not_nan)*C4_test(nc4);        
        while t <= (routing_time + min_time_step/60)
            % Infiltration and Available Depth
            % Show stats
            %perc____duremain______tsec_______dtmm______infmmhr____CmgL_____dtmWQ = [(t)/routing_time*100, (toc/((t)/routing_time) - toc)/3600,time_step*60,max(max(d_t(~isinf(d_t)))),max(max(C)), max(max((P_conc))), tmin_wq]
            if tmin_wq < 0 || isnan(tmin_wq) || isinf(tmin_wq)
                error('Instability')
            end
            % Inflows
            if flag_inflow == 1
                inflow = zeros(size(dem,1),size(dem,2)); % This is to solve spatially, don't delete
                for i = 1:n_stream_gauges
                    inflow = inflow + delta_inflow_agg(i)*inflow_cells(:,:,i);
                end
            end

            if k == 1
                if flag_infiltration == 1
                    i_a = (delta_p_agg.*rainfall_matrix + inflow + d_0 - ETP/24)./(time_step/60);
                    C = ksat.*(1 + ((d_0 + psi).*(teta_sat - teta_i))./I_0); % matrix form
                    f = min(C,i_a);
                    I_t = I_0 + (f - ETP)*time_step/60; % Extracting ETP from soil
                    pef = delta_p_agg.*rainfall_matrix + inflow - f*time_step/60;
                    d_t = d_0 + pef;
                    inf_m3s_t = f/1000.*cell_area/3600;
                    T = Tr; % Beginning to track the replenishing time
                else
                    %i_a = (delta_p_agg*rainfall_matrix + delta_inflow_agg*inflow_cells + d_0)./(time_step/60);
                    pef = delta_p_agg.*rainfall_matrix + inflow;
                    d_t = d_0 + pef;
                end
            else
                if flag_infiltration == 1
                    % Effective precipitation - Green-Ampt(1911)
                    i_a = (delta_p_agg.*rainfall_matrix + inflow + d_p - ETP/24)./(time_step/60);
                    C = ksat.*(1 + ((d_p + psi).*(teta_sat - teta_i))./I_p); % matrix form
                    idx_C = (i_a <= C); % Values where i_a is below C
                    idx_i_a = (i_a > C); % Values where i_a is larger than C
                    T = T - time_step/60/60; % Recoverying time (hours)
                    T(idx_i_a) = Tr(idx_i_a); % Saturated Areas we begin the process again
                    idx_T = T < 0; % Cells where the replenishing time is reached
                    f = min(C,i_a); % Infiltration rate (mm/hr)
                    I_t = max(I_p + f*time_step/60 - ETP/24 - k_out.*double(idx_C)*time_step/60,0);
                    % Refreshing I_t to I_0 for cases where idx_T > 0
                    I_t(idx_T) = I_0(idx_T);
                    pef = delta_p_agg.*rainfall_matrix + inflow - f*time_step/60;
                    d_t = d_p + pef; %% ATTENTION HERE
                    inf_m3s_t = f/1000.*cell_area/3600;
                else
                    i_a = (delta_p_agg.*rainfall_matrix + inflow + d_p - ETP/24)./(time_step/60);
                    pef = delta_p_agg.*rainfall_matrix + inflow - ETP/24;
                    d_t = d_p + pef;
                end
            end
            d_tot = d_t;

            %%%% WATER DEPTHS %%%
            d_left_cell = [zeros(ny_max,1),d_tot(:,1:(nx_max-1))];
            d_right_cell = [d_tot(:,(2:(nx_max))) zeros(ny_max,1)];
            d_up_cell = [zeros(1,nx_max) ; d_tot(1:(ny_max-1),:)];
            d_down_cell = [d_tot(2:(ny_max),:) ; zeros(1,nx_max)];
            if flag_D8 == 1
                % Simulate with D-8 flow direction
                % --- Adding NE, SE, SW, NW --- %
                zero_matrix = NaN(size(spatial_domain));
                elevation_NE_t = zero_matrix; elevation_SE_t = zero_matrix;
                elevation_SW_t = zero_matrix; elevation_NW_t = zero_matrix;
                d_NE_t = zero_matrix; d_SE_t = zero_matrix; d_SW_t = zero_matrix; d_NW_t = zero_matrix;

                elevation_NE_t(2:(ny_max),1:(nx_max-1)) = elevation(1:(ny_max-1),2:nx_max); % OK
                elevation_SE_t(1:(ny_max-1),1:(nx_max-1)) = elevation(2:ny_max,2:nx_max); % OK
                elevation_SW_t(1:(ny_max-1),2:(nx_max)) = elevation(2:(ny_max),1:(nx_max-1)); % OK
                elevation_NW_t(2:ny_max,2:nx_max) = elevation(1:(ny_max-1),1:(nx_max-1)); % OK

                d_NE_t(2:(ny_max),1:(nx_max-1)) = d_tot(1:(ny_max-1),2:nx_max); % OK
                d_SE_t(1:(ny_max-1),1:(nx_max-1)) = d_tot(2:ny_max,2:nx_max); % OK
                d_SW_t(1:(ny_max-1),2:(nx_max)) = d_tot(2:(ny_max),1:(nx_max-1)); % OK
                d_NW_t(2:ny_max,2:nx_max) = d_tot(1:(ny_max-1),1:(nx_max-1)); % OK
            else
                % Simulate with D-4 flow direction
                % Everything already calculated
            end
            %%%% ASSIGNING VALUES %%%
            elevation_cell = elevation;

            %  Flood routing for each cell
            % Check if Nan or Inf occured in d_t
            if max(max(isinf(d_t))) == 1
                ttt = 1;
            end
            if flag_D8 == 1 % D-8 C-A Model
                [qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,qout_ne_t,qout_se_t,qout_sw_t,qout_nw_t,d_t,I_tot_end_cell,I_cell] = CA_Routing_8D(elevation_cell,elevation_left_t,elevation_right_t,elevation_up_t,elevation_down_t,d_tot,d_left_cell,d_right_cell,d_up_cell,d_down_cell,roughness,cell_area,time_step,h_0,Resolution,I_tot_end_cell,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,idx_nan,flag_critical,elevation_NE_t, elevation_SE_t, elevation_SW_t, elevation_NW_t, d_NE_t, d_SE_t, d_SW_t, d_NW_t);
            else % 4-D CA Model
                [qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,d_t,I_tot_end_cell,I_cell] = CA_Routing(elevation_cell,elevation_left_t,elevation_right_t,elevation_up_t,elevation_down_t,d_tot,d_left_cell,d_right_cell,d_up_cell,d_down_cell,roughness,cell_area,time_step,h_0,Resolution,I_tot_end_cell,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,idx_nan,flag_critical);
            end
            %     qout_t = qout_left + qout_right + qout_up + qout_down + qout_ne + qout_se + qout_sw + qout_ne;
            % Inflows - Von-Neuman
            qin_left_t = [zeros(ny_max,1),qout_right_t(:,1:(nx_max-1))];
            qin_right_t = [qout_left_t(:,(2:(nx_max))) zeros(ny_max,1)];
            qin_up_t = [zeros(1,nx_max) ; qout_down_t(1:(ny_max-1),:)];
            qin_down_t = [qout_up_t(2:(ny_max),:) ; zeros(1,nx_max)];
            % Inflows - Inclined Directions
            if flag_D8 == 1
                zero_matrix = zeros(size(spatial_domain));
                qin_ne_t = zero_matrix; qin_se_t = zero_matrix; qin_sw_t = zero_matrix; qin_nw_t = zero_matrix;
                qin_ne_t(2:(ny_max),1:(nx_max-1)) = qout_sw_t(1:(ny_max-1),2:(nx_max)); % OK
                qin_se_t(1:(ny_max-1),1:(nx_max-1)) = qout_nw_t(2:ny_max,2:nx_max); % OK
                qin_sw_t(1:(ny_max-1),2:(nx_max)) = qout_ne_t(2:(ny_max),1:(nx_max-1)); % OK
                qin_nw_t(2:ny_max,2:nx_max) = qout_se_t(1:(ny_max-1),1:(nx_max-1)); % OK
            end
            if flag_D8 == 1
                qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t + qin_ne_t + qin_se_t + qin_sw_t + qin_nw_t;
                %         qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t;
            else
                qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t;
            end
            idx = qin_t > flow_tolerance;
            idx2 = qin_t <= flow_tolerance;
            idx3 = logical(isnan(qin_t) + isinf(qin_t));
            qin_t(idx3) = 0;
            flows_cells(idx) = 1;
            flows_cells(idx2) = 0;

            if max(max(qout_up_t)) > 0
                ttt = 1;
            end
            if t >= routing_time/2
                ttt = 1;
            end

            % Water Quality Parameters for f(B(t))
            if flag_waterquality == 1
                if flag_D8 ~= 1
                    [B_t,P_conc,Out_Conc,tmin_wq,tot_W_out,mass_lost,Tot_Washed] = build_up_wash_off(C_3,C_4,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,B_t,time_step,nx_max,ny_max,cell_area,outlet_index,idx_nan_5,flag_wq_model,mass_lost,Tot_Washed,Bmin,Bmax,min_Bt);
                else
                    [B_t,P_conc,Out_Conc,tmin_wq,tot_W_out,mass_lost,Tot_Washed] = build_up_wash_off_8D(C_3,C_4,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,qout_ne_t,qout_se_t,qout_sw_t,qout_nw_t,B_t,time_step,nx_max,ny_max,cell_area,outlet_index,idx_nan_5,flag_wq_model,mass_lost,Tot_Washed,Bmin,Bmax,min_Bt);
                end
            end
            %%%% Checking Mass Balance
            if flag_waterquality == 1
                if sum(sum(B_t(~isinf(B_t)))) > 1.2*initial_mass  % More than 5%
                    error('Brutal instability in B(t). More than 20% difference')
                end
            end
            % New Time-step Calculation
            pos_save = find(time_change_records < t,1,'last');
            time_save = time_change_records(pos_save); % min
            delta_time_save = time_save - time_save_previous;
            time_save_previous = time_save;
            actual_record_timestep = find(time_change_records < t,1,'last');
            if delta_time_save > 0 || k == 1 % First time-step
                if flag_timestep == 0
                    %%% SOLUTION FOR COURANT METHOD %%%
                    d_left_cell = [zeros(ny_max,1),d_t(:,1:(nx_max-1))];
                    d_right_cell = [d_t(:,(2:(nx_max))) zeros(ny_max,1)];
                    d_up_cell = [zeros(1,nx_max) ; d_t(1:(ny_max-1),:)];%
                    d_down_cell = [d_t(2:(ny_max),:) ; zeros(1,nx_max)];
                    if flag_D8 == 1
                        d_NE_cell = zero_matrix; d_SE_cell = zero_matrix; d_SW_cell = zero_matrix; d_NW_cell = zero_matrix;
                        d_NE_cell(2:(ny_max),1:(nx_max-1)) = d_t(1:(ny_max-1),2:nx_max); % OK
                        d_SE_cell(1:(ny_max-1),1:(nx_max-1)) = d_t(2:ny_max,2:nx_max); % OK
                        d_SW_cell(1:(ny_max-1),2:(nx_max)) = d_t(2:(ny_max),1:(nx_max-1)); % OK
                        d_NW_cell(2:ny_max,2:nx_max) = d_t(1:(ny_max-1),1:(nx_max-1)); % OK
                    end
                    % Considering Velocity as Celerity
                    %             vel_left = (I_cell(:,:,1)./(0.5/1000.*(d_t + d_left_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_left_cell)/2/1000); % m/s
                    %             vel_right = (I_cell(:,:,2)./(0.5/1000.*(d_t + d_right_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_right_cell)/2/1000);
                    %             vel_up = (I_cell(:,:,3)./(0.5/1000.*(d_t + d_up_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_up_cell)/2/1000); % m/s;;
                    %             vel_down = (I_cell(:,:,4)./(0.5/1000.*(d_t + d_down_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_down_cell)/2/1000); % m/s;

                    % Using CA Velocities
                    %             vel_left = (qout_left_t/1000/3600); % m/s
                    %             vel_right = (qout_right_t/1000/3600); % m/s
                    %             vel_up = (qout_up_t/1000/3600); % m/s
                    %             vel_down = (qout_down_t/1000/3600); % m/s

                    z = d_t;
                    z(d_t < 10) = 1e12;
                    vel_left = (qout_left_t/1000/3600)*Resolution^2./(Resolution*z/1000); % m/s
                    vel_right = (qout_right_t/1000/3600)*Resolution./(z/1000); % m/s
                    vel_up = (qout_up_t/1000/3600)*Resolution./(z/1000); % m/s
                    vel_down = (qout_down_t/1000/3600)*Resolution./(z/1000); % m/s

                    if flag_D8 == 1
                        %                 vel_ne = (I_cell(:,:,6)./(0.5/1000.*(d_t + d_NE_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_NE_cell)/2/1000); % m/s
                        %                 vel_se = (I_cell(:,:,7)./(0.5/1000.*(d_t + d_SE_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_SE_cell)/2/1000);
                        %                 vel_sw = (I_cell(:,:,8)./(0.5/1000.*(d_t + d_SW_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_SW_cell)/2/1000); % m/s;;
                        %                 vel_nw = (I_cell(:,:,9)./(0.5/1000.*(d_t + d_NW_cell).*Resolution))/(time_step*60) + sqrt(9.81*(d_t + d_NW_cell)/2/1000); % m/s;
                        vel_ne = (qout_ne_t/1000/3600); % m/s
                        vel_se = (qout_se_t/1000/3600); % m/s
                        vel_sw = (qout_sw_t/1000/3600); % m/s
                        vel_nw = (qout_nw_t/1000/3600); % m/s
                    end
                    %%%%%%%%%%%%%% Find the Maximum Velocity
                    max_velocity_left = max(max(vel_left));
                    max_velocity_right = max(max(vel_right));
                    max_velocity_up = max(max(vel_up));
                    max_velocity_down = max(max(vel_down));
                    if flag_D8 == 1
                        max_velocity_ne = max(max(vel_ne));
                        max_velocity_se = max(max(vel_se));
                        max_velocity_sw = max(max(vel_sw));
                        max_velocity_nw = max(max(vel_nw));
                    end
                    % - Velocit Raster - %
                    %%% Maximum of All of Them %%%
                    velocity_raster = max(vel_left,vel_right);
                    velocity_raster = max(velocity_raster,vel_up);
                    velocity_raster = max(velocity_raster,vel_down);
                    if flag_D8 == 1
                        velocity_raster = max(velocity_raster,vel_ne);
                        velocity_raster = max(velocity_raster,vel_se);
                        velocity_raster = max(velocity_raster,vel_sw);
                        velocity_raster = max(velocity_raster,vel_nw);
                        velocity_vector = [max_velocity_left, max_velocity_right, max_velocity_up, max_velocity_down, max_velocity_ne, max_velocity_se, max_velocity_sw, max_velocity_nw];
                    else
                        velocity_vector = [max_velocity_left, max_velocity_right, max_velocity_up, max_velocity_down];
                    end

                    max_velocity = max(velocity_vector);
                    if k == 1
                        dvdt = 0;
                    else
                        dvdt = max((max_velocity - max_velocity_previous)/(time_step_change),0); % m/s2
                    end
                    flow_velocity(pos_save,1) = max_velocity;
                    flow_acceleration(pos_save,1) =  dvdt; % m/s2
                    if flag_D8 == 1
                        factor_grid = sqrt(1/2);
                    else
                        factor_grid = 1;
                    end
                    if max_velocity == 0
                        new_timestep = max_time_step;
                    else
                        new_timestep = (factor_grid*Resolution/max_velocity); % seconds
                    end
                    time_step_factor = max(alfa_max - slope_alfa*(max(max_velocity - v_threshold,0)),alfa_min);
                    alfa_save(pos_save,1) = time_step_factor;
                    new_timestep = new_timestep*alfa_save(pos_save,1);

                    % ---- Calculation of Stability --- %
                    F_person = weight_person*gravity; % N
                    % Buyoance
                    F_buoy = (width1_person*width2_person)*d_t/1000*ro_water*gravity; % N
                    % Available Friction
                    available_friction = mu*(max(F_person - F_buoy,0));
                    % Hydrodynamic Force
                    hydro_force = max(1/2*(Cd*width1_person*d_t/1000.*velocity_raster.^2),0);
                    % Risk Factor
                    risk_t = min(hydro_force./available_friction,1);
                    risk_t(idx_nan) = nan;
                else             % % % % % % % %  Solution for the Stable Method - Time-step refreshment
                    % Water Slopes Calculation (THERE IS A MISTAKE HERE)
                    error('This method is currenly not working, please choose the courant method')
                    %             wse_cell = elevation_cell + d_t_cell/1000;
                    %             slope(:,:,1) = 1/Resolution*(wse_cell - elevation_left_t - d_left_cell/1000);
                    %             slope(:,:,2) = 1/Resolution*(wse_cell - elevation_right_t - d_right_cell/1000);
                    %             slope(:,:,3) = 1/Resolution*(wse_cell - elevation_up_t - d_up_cell/1000);
                    %             slope(:,:,4) = 1/Resolution*(wse_cell - elevation_down_t - d_down_cell/1000);
                    %             slope_cell = max(min(slope,[],3),slope_min); % Minimum assumed slope
                    %             time_step_cell = max((Resolution^2/4)*(2*roughness./((d_t_cell/1000).^(5/3)).*sqrt(slope_cell)),min_time_step);
                    %             time_step_cell(idx_nan_5(:,:,1)) = inf;
                    %             new_timestep = min(min((time_step_cell)));
                end
                if flag_waterquality == 1
                    % Adding Water Quality Minimum Time-step
                    alfa_wq = 1; % Reduction factor of water quality
                    new_timestep = min(alfa_wq*tmin_wq,new_timestep);
                else
                    new_timestep = min(time_step_factor*new_timestep);
                end
                % Rounding time-step to min_timestep or max_timestep with the required
                % precision. This is not very interesting, maybe we should delete
                % it
                time_calculation_routing(k,1) = new_timestep;
                time_calculation_routing(k,1) = max(time_step_increments*floor(time_calculation_routing(k,1)/time_step_increments),min_time_step);
                time_calculation_routing(k,1) = min(time_calculation_routing(k,1),max_time_step);
                time_step = time_calculation_routing(k,1)/60; % new time-step for the next run
                if time_calculation_routing(k,1) == min_time_step % Check if we reached the minimum time-step
                    unstable = 1; % If this is 1, it means we possibly had an unstable period at least
                end
            elseif  k == 1
                time_calculation_routing(k,1) = time_step*60; % sec
            else
                time_calculation_routing(k,1) = time_calculation_routing(k-1,1);
            end
            max_velocity_previous = max_velocity; % Assigning right velocities
            % Agregating Inflows to the New Time-step
            if flag_inflow > 0
                for z = 1:n_stream_gauges
                    z1 = find(time_deltainflow(z,:) > t_previous,1,'first'); % begin of the time-step
                    z2 = find(time_deltainflow(z,:) <= t,1,'last'); % end of the time-step
                    if isempty(z1)
                        z1 = 1;
                    end
                    if isempty(z2) || z2 < z1
                        z2 = z1;
                    end
                    if time_step >= time_step_model
                        delta_inflow_agg(z,1) = mean(delta_inflow(z,z1:z2))/(time_step_model*60)*time_step*60;
                    else
                        delta_inflow_agg(z,1) = delta_inflow(z,z1)/(time_step_model*60)*time_step*60;
                    end
                end
            end
            % Agregating Precipitation to the New Time-step
            if flag_rainfall > 0
                if flag_spatial_rainfall ~= 1
                    z1 = find(time_deltap > t_previous,1,'first'); % begin of the time-step
                    z2 = find(time_deltap <= t,1,'last'); % end of the time-step
                    if z2 < z1
                        z2 = z1;
                    end
                    if time_step >= time_step_model
                        delta_p_agg = mean(delta_p(1,z1:z2))/(time_step_model*60)*time_step*60;
                    else
                        delta_p_agg = delta_p(1,z1)/(time_step_model*60)*time_step*60;
                    end
                else
                    % Spatial Rainfall
                    % Code for Spatial-Varying Rainfall
                    z1 = find(rainfall_spatial_duration <= t_previous,1,'last');
                    z2 = find(rainfall_spatial_duration <= t,1,'last');
                    if z1 ~= z2 || z2 == length(rainfall_spatial_duration)
                        x_coordinate = coordinates(1:n_raingauges,1); % Coordinates (easting) of each rain gauge
                        y_coordinate = coordinates(1:n_raingauges,2); % Coordinates (northing) of each rain gauge
                        x_grid = xulcorner + Resolution*[1:1:size(DEM_raster.Z,2)]'; % Pixel easting coordinates
                        y_grid = yulcorner - Resolution*[1:1:size(DEM_raster.Z,1)]'; % Pixel northing coordinates

                        % Spatial Rainfall
                        if z2 == length(rainfall_spatial_duration)
                            spatial_rainfall = zeros(size(dem,1),size(dem,2));

                        else
                            rainfall = rainfall_raingauges(z2,1:n_raingauges)'; % Values of rainfall at t for each rain gauge
                            [spatial_rainfall] = Rainfall_Interpolator(x_coordinate,y_coordinate,rainfall,x_grid,y_grid); % Interpolated Values
                        end
                        delta_p_agg = spatial_rainfall/3600*time_step_model*60; % Matrix of delta P for each pixel
                        spatial_rainfall_maps(:,:,z2) = spatial_rainfall;
                        average_spatial_rainfall(z2,1) = mean(spatial_rainfall(spatial_rainfall>=0));
                    end
                end
            end

            % Aggregating ETP for next time-step

            if flag_ETP == 1
                z1 = find(climatologic_spatial_duration <= t_previous,1,'last');
                z2 = find(climatologic_spatial_duration <= t,1,'last');
                if isempty(z1) && isempty(z2)
                    % No data, No ETP.
                    ETP = zeros(size(dem));
                else
                    if ~isempty(z1) && z1 == z2 && z2 == length(climatologic_spatial_duration)
                        % Outside of ETP data maximum duration
                        ETP = zeros(size(dem));
                        ETP_save(:,:,z2) = ETP; % Saving ETP Maps
                    elseif  (isempty(z1) && z2 > 0) || z2 > z1 && z2 < length(climatologic_spatial_duration)
                        day_of_year = day(time_ETP(z2,1),'dayofyear');
                        [ETP] = ETP_model(z2,day_of_year,coordinates_stations(:,1),coordinates_stations(:,2),x_grid',y_grid',maxtemp_stations,mintemp_stations,avgtemp_stations,u2_stations,ur_stations,G_stations,DEM_etp,lat,Krs,alfa_albedo_input,idx_nan);
                        ETP_save(:,:,z2) = ETP; % Saving ETP Maps
                    elseif z1 == 1 && z2 == 1
                        % First data. We assume a constant ETP using
                        % 1st data
                        day_of_year = day(time_ETP(z2,1),'dayofyear');
                        [ETP] = ETP_model(z2,day_of_year,coordinates_stations(:,1),coordinates_stations(:,2),x_grid',y_grid',maxtemp_stations,mintemp_stations,avgtemp_stations,u2_stations,ur_stations,G_stations,DEM_etp,lat,Krs,alfa_albedo_input,idx_nan);
                        ETP_save(:,:,z2) = ETP; % Saving ETP Maps
                    end
                end
            end

            % Previous Time-step
            if k == 1
                t_previous = time_calculation_routing(k,1)/60;
                t_previous_date = date_begin + time_calculation_routing(k,1)/60/60/24; % Datetime
            else
                t_previous = t;
                t_previous_date = t_previous_date + t/60/24; % Datetime
            end
            % New Domain - Matrices
            % We are temporaly deactivating this code due to lack of checking
            flag_out = 0;
            if (flag_rainfall == 0 && flag_inflow == 1 && flag_out == 1) % We are simulating only inflow hydrograph
                pos_matrices = find(time_change_matrices < t,1,'last');
                time_matrices = time_change_matrices(pos_matrices); % min
                if k == 1
                    time_matrices = 1;
                end
                delta_time_matrices = time_matrices - time_matrices_previous;
                time_matrices_previous = time_matrices;

                if delta_time_matrices > 0
                    % Finding cells receiving water %%
                    %%% Local Coordinate System from the begining of the Domain
                    [y_receiving, x_receiving] = find(flows_cells>0);
                    %%%% Finding cells with stored water %%%%
                    [y_depth_stored, x_depth_stored] = find(d_t > depth_tolerance);
                    [row_check, col_check] = find(inflow_cells > 0);
                    cells_check = [y_receiving x_receiving ; y_depth_stored x_depth_stored; row_check col_check];
                    cells_check = unique(cells_check,'rows');    % Domain Fronteirs - Global Coordinate System

                    xmin_cells= (min(cells_check(:,2)) -1 ) + coordinate_x;
                    xmin_cells_t = xmin_cells; % Refresh
                    xmax_cells = (max(cells_check(:,2)) -1 ) + coordinate_x;

                    ymin_cells = (min(cells_check(:,1)) - 1 ) + coordinate_y;
                    ymin_cells_t = ymin_cells;
                    ymax_cells = (max(cells_check(:,1)) - 1 ) + coordinate_y;

                    %%% Min Values (Only works in the absolute coordinate reference systen)
                    xmin_domain_new = max(xmin_cells - factor_cells,1);
                    ymin_domain_new = max(ymin_cells - factor_cells,1);

                    %%% Max Values
                    xmax_domain_new = min(xmax_cells + factor_cells,xmax); %xmax is wrong here
                    ymax_domain_new = min(ymax_cells + factor_cells,ymax);
                    domain = zeros(ymax_domain_new - ymin_domain_new + 1,xmax_domain_new - xmin_domain_new  + 1);
                    %%%%% Domain is in the local coordinate system
                    % New Coordinates - Local of the 1st cell in the Global Coordinate System
                    coordinate_x_previous = coordinate_x;
                    coordinate_y_previous = coordinate_y;
                    if coordinate_x == 1 % First cut in the domain (outside border)
                        coordinate_x = xmin_domain_new; % Using the outside reference (new reference)
                        coordinate_y = ymin_domain_new;
                        begin = 1; % Initial Abstraction of the Matrices
                    else
                        begin = 0;
                        coordinate_x = xmin_domain_new;
                        coordinate_y = ymin_domain_new;
                    end

                    % Matrices Sizes - Deslocation (Local Coordinate System)
                    dymax_matrix = ymax_cells - coordinate_y + 1; % Top of the matrix
                    dymin_matrix = ymin_cells - coordinate_y + 1; % Bottom
                    dxmax_matrix = xmax_cells - coordinate_x + 1; % Right side
                    dxmin_matrix = xmin_cells - coordinate_x + 1; % Left side
                    [ny_max,nx_max] = size(domain);

                    %%%% Abstracted Matrices
                    if begin == 1
                        ymin_domain_abstracted = ymin_cells - 0; % Already in the right coordinate system
                        ymax_domain_abstracted = ymax_cells - 0;
                        xmin_domain_abstracted = xmin_cells - 0;
                        xmax_domain_abstracted = xmax_cells - 0;
                    else
                        ymin_domain_abstracted = ymin_cells - coordinate_y_previous + 1;
                        ymax_domain_abstracted = ymax_cells - coordinate_y_previous + 1;
                        xmin_domain_abstracted = xmin_cells - coordinate_x_previous + 1;
                        xmax_domain_abstracted = xmax_cells - coordinate_x_previous + 1;
                    end

                    % New Matrices
                    %%% Depths
                    d_t = resize_matrix(d_t,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %     surf(d_t); view(0,90); xlabel('x'); ylabel('y')
                    %%% Inflow Cells
                    inflow_cells = resize_matrix(inflow_cells,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %%% Rainfall Matrix
                    rainfall_matrix = resize_matrix(rainfall_matrix,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %%% Flows Cells
                    flows_cells = resize_matrix(flows_cells,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %%% Inundation Cells
                    %     inundated_cells = resize_matrix(inundated_cells,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %%% Qins
                    qin_t = resize_matrix(qin_t,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %     qout_left_t = resize_matrix(qout_left_t,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %     qout_right_t = resize_matrix(qout_right_t,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %     qout_up_t = resize_matrix(qout_up_t,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    %     qout_down_t = resize_matrix(qout_down_t,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    % Velocities
                    vel_down = resize_matrix(vel_down,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    vel_up = resize_matrix(vel_up,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    vel_left = resize_matrix(vel_left,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    vel_right = resize_matrix(vel_right,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    % I_tot_end_cell
                    I_tot_end = resize_matrix(I_tot_end,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    I_tot_end_cell = I_tot_end;
                    % I_cell
                    I_cell = resize_matrix(I_cell,domain,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,5);

                    %%%% All Domain Cells %%%%
                    %%% Elevation
                    elevation = dem(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    % Infiltration Parameters
                    teta_sat = teta_sat_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    teta_i = teta_i_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    psi = psi_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    ksat = ksat_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    %  Roughness
                    roughness = roughness_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    % Initial Abstraction
                    h_0 = h_0_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    % Accumulated Infiltration
                    domain_infiltration = I_0(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    I_t = resize_matrix(I_t,domain_infiltration,ymin_domain_abstracted,ymax_domain_abstracted,xmin_domain_abstracted,xmax_domain_abstracted,dymin_matrix,dymax_matrix,dxmin_matrix,dxmax_matrix,1);
                    I_p = I_t;
                    % Outlet Cells
                    outlet_index = outlet_index_fulldomain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);
                    % Rainfall Cells
                    rainfall_matrix = rainfall_matrix_full_domain(ymin_domain_new:ymax_domain_new,xmin_domain_new:xmax_domain_new);

                    % Cells Check
                    [y_receiving, x_receiving] = find(flows_cells>0);
                    [y_depth_stored, x_depth_stored] = find(d_t > depth_tolerance);
                    [row_check, col_check] = find(inflow_cells > 0);
                    cells_check = [y_receiving x_receiving ; y_depth_stored x_depth_stored; row_check col_check];
                    cells_check = unique(cells_check,'rows');
                end
            end
            % Inflows and Depth Refreshments
            d_t = d_t + qin_t*time_step/60;
            % Clearing stored values
            qout_left_t = zeros(ny_max,nx_max);
            qout_right_t = zeros(ny_max,nx_max);
            qout_up_t = zeros(ny_max,nx_max);
            qout_down_t = zeros(ny_max,nx_max);
            qout_ne_t = zeros(ny_max,nx_max);
            qout_se_t = zeros(ny_max,nx_max);
            qout_sw_t = zeros(ny_max,nx_max);
            qout_nw_t = zeros(ny_max,nx_max);
            qin_t = zeros(ny_max,nx_max);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Checking Water Balance and Saving State Variables %%%%
            % FIX IT!!!
            if flag_waterbalance  == 1
                n_rain_cells = sum(sum(idx_not_nan));
                if k == 1
                    storage_begin = sum(sum(cell_area/1000.*d_0(idx_not_nan))); % m3
                else
                    storage_begin = storage_depth;
                end
                inflows = delta_inflow_agg/1000*n_inlets*cell_area + delta_p_agg/1000*(n_rain_cells)*cell_area;
                outflows = sum(sum(outlet_flow))/1000/3600*cell_area*time_step/60 + sum(sum((I_t(idx_not_nan) - I_p(idx_not_nan))/1000*cell_area));
                delta_flows = inflows - outflows ; % m3
                storage_depth = sum(sum(cell_area/1000.*d_t(idx_not_nan)));
                delta_storage = storage_depth - storage_begin;
                relative_vol_error(k) = 1 - delta_flows/delta_storage;
                extra_depth = (delta_flows - delta_storage)/((sum(sum(idx_not_nan)))*cell_area)*1000;
                d_t = d_t + idx_not_nan*extra_depth;

                % Water Balance Error
                relative_vol_error(k) = 1 - delta_flows/delta_storage;
            end
            % Saving Plotting Values - Recording Time
            % Maps of Flood Depths, WSE and Pollutant Concentrations
            % --- Calculating EMC --- %
            mass_outlet = max(mass_outlet + Out_Conc*((nansum(nansum(outlet_flow)/1000/3600*1000)))*(time_step*60),0); % mg
            vol_outlet = max((nansum(nansum(outlet_flow))/1000/3600*1000)*(time_step*60) + vol_outlet,0);
            % --- Calculating M(V) Curve --- %
            t_save = t + time_calculation_routing(k,1)/60;
            actual_record_state = find(time_records < t_save,1,'last');
            delta_record = actual_record_state - last_record_maps;
            last_record_maps = actual_record_state;
            if k == 1
                d((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),1) = d_t;
                risk((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),1) = risk_t;
                if flag_waterquality(1,1) == 1
                    Pol_Conc_Map((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),1) = P_conc;
                    Pol_mass_map((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),1) = B_t;
                    EMC_outlet(1,1) = mass_outlet/vol_outlet; % mg/L
                    mass_outlet_save(1,1) = mass_outlet;
                    vol_outlet_save(1,1) = vol_outlet;
                end
            elseif delta_record > 0
                t_store = actual_record_state;
                d((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),t_store) = d_t;
                risk((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),t_store) = risk_t;
                if flag_waterquality == 1
                    Pol_Conc_Map((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),t_store) = P_conc;
                    Pol_mass_map((coordinate_y:coordinate_y + ny_max - 1),(coordinate_x:coordinate_x + nx_max - 1),t_store) = B_t;
                    EMC_outlet(t_store,1) = mass_outlet/vol_outlet; % mg/L
                    mass_outlet_save(t_store,1) = mass_outlet;
                    vol_outlet_save(t_store,1) = vol_outlet;
                end
            end % Calls the sub
            % Hydrographs
            actual_record_hydrograph = find(time_record_hydrograph < t_save,1,'last');
            delta_record_hydrograph = actual_record_hydrograph - last_record_hydrograph;
            last_record_hydrograph = actual_record_hydrograph;

            if k == 1
                outet_hydrograph(1,1) = nansum(nansum(outlet_flow))/1000*cell_area/3600; % m3/s
                time_hydrograph(1,1) = time_calculation_routing(k,1)/60;
                depth_outlet(1,1) = mean(d_t(idx_outlet));
                % Saving Data of Input Gauges
                if flag_obs_gauges == 1
                    for i = 1:num_obs_gauges
                        x_cell = easting_obs_gauges(i); y_cell = northing_obs_gauges(i);
                        wse_cell(1,1) = d_t(y_cell,x_cell)/1000 + elevation(y_cell,x_cell); % m
                        depth_cell(1,1) = d_t(y_cell,x_cell)/1000; % m
                        hydrograph_cell(1,i) = I_tot_end_cell(y_cell,x_cell)/(time_calculation_routing(k)); % m3/s
                    end
                end
                if flag_waterquality == 1
                    outet_pollutograph(1,1) = Out_Conc; % Already averaged for all outlet cells
                end
            elseif delta_record_hydrograph > 0
                t_store = actual_record_hydrograph;
                outet_hydrograph(t_store,1) = nansum(nansum(outlet_flow))/1000*cell_area/3600; % m3/s
                time_hydrograph(t_store,1) = t;
                depth_outlet(t_store,1) = mean(d_t(idx_outlet));
                if flag_obs_gauges == 1
                    % Saving Data of Input Gauges
                    for i = 1:num_obs_gauges
                        x_cell = easting_obs_gauges(i); y_cell = northing_obs_gauges(i);
                        wse_cell(t_store,i)= d_t(y_cell,x_cell)/1000 + elevation(y_cell,x_cell); % m
                        depth_cell(t_store,i) = d_t(y_cell,x_cell)/1000; % m
                        hydrograph_cell(t_store,i) = I_tot_end_cell(y_cell,x_cell)/(time_calculation_routing(k)); % m3/s
                    end
                end

                if flag_waterquality == 1
                    outet_pollutograph(t_store,1) = Out_Conc;
                    % Saving Data of Input Gauges
                    for i = 1:num_obs_gauges
                        x_cell = easting_obs_gauges(i); y_cell = northing_obs_gauges(i);
                        pollutograph_cell(t_store,i)= P_conc(y_cell,x_cell); % mg/L
                    end
                end
            end % Calls the sub
            % Saving Maximum Depths and Outlet Flows
            if k == 1
                dmax_final = d_t;
                vmax_final = velocity_raster;
            else
                dmax_final = max(d_t,dmax_final);
                vmax_final = max(velocity_raster,vmax_final);
            end
            outlet_runoff_volume = nansum(nansum(outlet_flow))*time_step/60*cell_area/drainage_area + outlet_runoff_volume; % mm

            % Previous Depths
            d_p = d_t;
            I_p = I_t;

            % Check if Nan or Inf occured in d_t
            if max(max(isinf(d_t))) == 1
                unstable_dt = 1;
            end

            % Runoff Coefficient Calculation
            outflow_volume = nansum(nansum(outlet_flow))/1000*cell_area/3600*time_step*60 + outflow_volume;
            if flag_spatial_rainfall == 1
                inflow_volume = sum(sum(delta_inflow_agg'/1000.*n_inlets*cell_area)) + nansum(nansum(delta_p_agg))/1000*cell_area + inflow_volume; % check future
            elseif flag_spatial_rainfall ~= 1 && flag_inflow == 1
                inflow_volume = sum(sum(delta_inflow_agg'/1000.*n_inlets*cell_area)) + delta_p_agg/1000*drainage_area + inflow_volume; % check future
            else
                inflow_volume = delta_p_agg/1000*drainage_area + inflow_volume; % check future
            end
            % Saving Modeled Values
            obs_index = find(t >= time_observed_concentration,1,'last');
            modeled_concentration(1,obs_index) = Out_Conc;
            
            % Increase the counter
            t = time_calculation_routing(k,1)/60 + t;
            time_step_save(k,2) = time_calculation_routing(k,1);
            time_step_save(k,1) = t;
            k = k + 1;        
        end
        % Comparing Pollutographs
        perc_tot = i_save/(n_tests)*100
        error = sum((observed_concentration' - modeled_concentration).^2);
        % Save Results
        data_water_quality(i_save,:) = [C3_test(nc3),C4_test(nc4),error,modeled_concentration];
        i_save = i_save + 1;    
        t = time_step;
        B_t = initial_mass/(nx_max*ny_max);
    end
end

%% Post-Processing Results
    close all
    post_processing