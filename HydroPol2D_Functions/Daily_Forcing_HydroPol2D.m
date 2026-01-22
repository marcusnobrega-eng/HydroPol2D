%% ========================================================================
% SCRIPT: Generate Rainfall and ETP Spreadsheets from GEE Output
% AUTHOR: (Marcus Nobrega, Ph.D)
% DATE: (Jun 2025)
%
% PURPOSE:
%   This script takes a CSV file exported from Google Earth Engine (GEE)
%   containing gridded daily climate data (rainfall, temperature, wind).
%   It organizes and outputs two Excel spreadsheets:
%     1. 'Rainfall_Structured.xlsx' with rainfall and spatial metadata
%     2. 'ETP_Structured.xlsx' with ETP-related parameters and spatial info
%
% GEE OUTPUT NEEDED:
%   Use a script (see below) in Google Earth Engine to export a CSV with:
%     - date
%     - x_mercator (easting in EPSG:3857)
%     - y_mercator (northing in EPSG:3857)
%     - rainfall
%     - temperature_mean, temperature_max, temperature_min
%     - wind_speed
%
% AFTER GENERATING THE SPREADSHEETS
%  Fill your Spatial_Rainfall_Input and ETP_Input with the data from the
%  spreadsheet generated
%
% GEE TOOL REFERENCE:
%   • CHIRPS:  'UCSB-CHG/CHIRPS/DAILY'
%   • ERA5-Land: 'ECMWF/ERA5_LAND/DAILY_AGGR'
%   • See the GEE script block provided in your project notes
%% ========================================================================

%% === LOAD DATA ==========================================================

% Load climate data (from GEE export)
climate_data = readtable('Gauge_Daily_Rainfall_Temperature_Wind.csv');
climate_data.date = datetime(climate_data.date);  % Convert date column

% Extract unique spatial coordinates for gauges
[unique_coords, ~, ~] = unique([climate_data.x_mercator, climate_data.y_mercator], 'rows');
gauge_easting = unique_coords(:,1);
gauge_northing = unique_coords(:,2);
num_gauges = size(unique_coords, 1);

% Extract unique dates
unique_dates = unique(climate_data.date);
num_days = length(unique_dates);
time_minutes = (0:num_days-1)' * 60;  % Assume each record is 1-day step

%% === RAINFALL MATRIX ====================================================

rainfall_matrix = nan(num_days, num_gauges);
for g = 1:num_gauges
    gx = gauge_easting(g);
    gy = gauge_northing(g);
    gauge_data = climate_data(climate_data.x_mercator == gx & climate_data.y_mercator == gy, :);
    [~, loc] = ismember(gauge_data.date, unique_dates);
    rainfall_matrix(loc, g) = gauge_data.rainfall;
end

% Build rainfall headers
gauge_names = strcat("Gauge ", string(1:num_gauges));
header1 = ['Index', num2cell(1:num_gauges)];
header2 = ['Rain Gauge', cellstr(gauge_names)];
header3 = ['Easting [m]', num2cell(gauge_easting')];
header4 = ['Northing [m]', num2cell(gauge_northing')];
header5 = ['Date (min)', repmat({'Rainfall Intensity (mm/h)'}, 1, num_gauges)];

% Combine header and data
rain_table = [num2cell(time_minutes), num2cell(rainfall_matrix)];
rain_output = [header1; header2; header3; header4; header5; rain_table];

% Write rainfall sheet
writecell(rain_output, 'Rainfall_Structured.xlsx', 'Sheet', 'Rainfall');

%% === ETP MATRIX =========================================================

% Columns: Tmax, Tmin, Tmed, Wind, UR, G
etp_matrix = nan(num_days, num_gauges * 6);
param_labels = {'Tmax (°C)', 'Tmin (°C)', 'Tmed (°C)', 'U2 [m/s]', 'UR (%)', 'G (MJ/(m2.dia))'};

for g = 1:num_gauges
    gx = gauge_easting(g);
    gy = gauge_northing(g);
    gauge_data = climate_data(climate_data.x_mercator == gx & climate_data.y_mercator == gy, :);
    [~, loc] = ismember(gauge_data.date, unique_dates);

    Tmax = nan(num_days,1);
    Tmin = nan(num_days,1);
    Tmed = nan(num_days,1);
    U2   = nan(num_days,1);
    UR   = 50 * ones(num_days,1);  % Constant humidity
    G    = zeros(num_days,1);      % Constant soil heat flux

    Tmax(loc) = gauge_data.temperature_max;
    Tmin(loc) = gauge_data.temperature_min;
    Tmed(loc) = gauge_data.temperature_mean;
    U2(loc)   = gauge_data.wind_speed;

    block = [Tmax, Tmin, Tmed, U2, UR, G];
    etp_matrix(:, (g-1)*6+1 : g*6) = block;
end

% Build ETP header rows
etp_header1 = {'Date'};
etp_header2 = {};

for g = 1:num_gauges
    etp_header1 = [etp_header1, {'Index'}, {['A', num2str(100 + g)]}, {'x easting (m)'}, gauge_easting(g), {'y northing (m)'}, gauge_northing(g)];
    etp_header2 = [etp_header2, param_labels];
end

etp_table = [cellstr(datestr(unique_dates)), num2cell(etp_matrix)];
etp_output = [etp_header1; ['Date', etp_header2]; etp_table];

% Write ETP sheet
writecell(etp_output, 'ETP_Structured.xlsx', 'Sheet', 'ETP');

%% === DONE ===============================================================
disp('Rainfall and ETP sheets successfully written.');


%% GEE CODE TO GENERATE DATA (JAVASCRIPT)
% // Scale of image
% var scale_of_image = 1000;
% 
% // Set start year and end year
% var startyear = 2015;
% var endyear = 2020;
% 
% // Define Projections
% var crs = 'EPSG:3857';  // 3857 is pseudomercartor
% 
% // Define Geometry
% var basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_5");
% var Folder_Name = 'Amazon_5km';
% var hybas = basins;
% var catchmentId = 7050585120;
% var crs = 'EPSG:3857';
% var roi = hybas.filter(ee.Filter.eq('HYBAS_ID', catchmentId));
% 
% // Input ROI
% var roi = amazon
% 
% // Visualize the catchment on the map
% Map.addLayer(roi, {color: 'red'}, 'Catchment Area');
% Map.centerObject(roi, 9);
% 
% // Load CHIRPS daily precipitation dataset
% var rainfall_dataset = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY');
% 
% // Load ERA5-Land dataset for temperature and fluxes
% var era5_land = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
%   .filterDate(ee.Date.fromYMD(startyear, 1, 1), ee.Date.fromYMD(endyear, 12, 31))
%   .select([
%     'temperature_2m', 'temperature_2m_min', 'temperature_2m_max', 
%     'surface_net_solar_radiation_sum', 'surface_sensible_heat_flux_sum', 'surface_latent_heat_flux_sum',
%     'u_component_of_wind_10m', 'v_component_of_wind_10m' // Added wind components
%   ]);
% 
% // Generate equally spaced gauge points within the catchment
% var num_x = 20;
% var num_y = 20;
% // Define the region of interest (roi) and bounds
% var bounds = roi.geometry().bounds().coordinates().flatten();
% var xmin = ee.Number(ee.List(bounds).get(0));
% var ymin = ee.Number(ee.List(bounds).get(1));
% var xmax = ee.Number(ee.List(bounds).get(4));
% var ymax = ee.Number(ee.List(bounds).get(5));
% 
% // Calculate xStep and yStep in degrees
% var xStep = xmax.subtract(xmin).divide(ee.Number(num_x - 1));
% var yStep = ymax.subtract(ymin).divide(ee.Number(num_y - 1));
% 
% // Convert xStep and yStep from degrees to meters
% var lat = ee.Number(roi.geometry().centroid().coordinates().get(1)); // Latitude (y)
% var latDeg = lat.multiply(Math.PI).divide(180);
% var degreesToMeters = ee.Number(111320).multiply(ee.Number(latDeg).cos()); // Cosine correction for longitude
% var xStepMeters = xStep.multiply(degreesToMeters); // Convert longitude to meters
% var yStepMeters = yStep.multiply(111320); // Convert latitude to meters
% 
% // Assuming scale_of_image is in meters
% var scale_of_image = ee.Number(scale_of_image);
% 
% // Find the minimum between scale_of_image, xStepMeters, and yStepMeters
% var minStep = xStepMeters.min(yStepMeters);
% 
% // Generate points for gauging locations
% var points = ee.List.sequence(0, num_x - 1).map(function(i) {
%   return ee.List.sequence(0, num_y - 1).map(function(j) {
%     var x = xmin.add(xStep.multiply(i));
%     var y = ymin.add(yStep.multiply(j));
%     var point = ee.Geometry.Point([x, y]);
%     return ee.Feature(point);
%   });
% }).flatten();
% 
% var gauge_points = ee.FeatureCollection(points).filterBounds(roi);
% Map.addLayer(gauge_points, {color: 'blue'}, 'Gauge Locations');
% 
% // Extract all daily images within the given years
% var years = ee.List.sequence(startyear, endyear);
% var days = ee.List.sequence(0, 364);
% 
% var dailyRain = ee.ImageCollection.fromImages(
%   years.map(function(y) {
%     return days.map(function(d) {
%       var date = ee.Date.fromYMD(y, 1, 1).advance(ee.Number(d), 'day');
%       var rain = rainfall_dataset.filterDate(date, date.advance(1, 'day'))
%                     .sum()
%                     .clip(roi);
%       return rain.set('year', y)
%                  .set('day', ee.Number(d).add(1))
%                  .set('system:time_start', date.millis());
%     });
%   }).flatten()
% );
% 
% // Optimized function to extract daily data at all gauge points
% var extractDailyData = function(image) {
%   var date = image.date();
%   var era5_image = era5_land.filterDate(date, date.advance(1, 'day')).first();
% 
%   // Combine all reduceRegion calls into one
%   var extractedData = image.addBands(era5_image).reduceRegions({
%     collection: gauge_points,
%     reducer: ee.Reducer.mean(),
%     scale: scale_of_image
%   }).map(function(feature) {
%     var projectedGeom = feature.geometry().transform(crs); // Reproject to EPSG:3857
%     var coords = projectedGeom.coordinates(); // Get projected coordinates
% 
%     return feature.set({
%       'date': image.date().format('yyyy-MM-dd HH:mm'),
%       'rainfall': feature.get('precipitation'), // Extract rainfall data
%       'x_mercator': coords.get(0),  // X-coordinate (Eastings)
%       'y_mercator': coords.get(1)   // Y-coordinate (Northings)
%     });
%   });
% 
% 
%   return extractedData.map(function(feature) {
%     return feature.set({
%       date: date.format('yyyy-MM-dd'),
%       rainfall: ee.Algorithms.If(feature.get('precipitation'), feature.get('precipitation'), 0),
%       temperature_mean: ee.Algorithms.If(feature.get('temperature_2m'),
%                                          ee.Number(feature.get('temperature_2m')).subtract(273.15),
%                                          null),
%       temperature_max: ee.Algorithms.If(feature.get('temperature_2m_max'),
%                                         ee.Number(feature.get('temperature_2m_max')).subtract(273.15),
%                                         null),
%       temperature_min: ee.Algorithms.If(feature.get('temperature_2m_min'),
%                                         ee.Number(feature.get('temperature_2m_min')).subtract(273.15),
%                                         null),
%       wind_speed: ee.Algorithms.If(
%         feature.get('u_component_of_wind_10m'),
%         ee.Number(feature.get('u_component_of_wind_10m')).pow(2)
%           .add(ee.Number(feature.get('v_component_of_wind_10m')).pow(2)).sqrt(),
%         null
%       )
%     });
%   });
% };
% 
% 
% // Create daily data (rainfall, temperature, wind speed)
% var dailyData = dailyRain.map(extractDailyData).flatten();
% 
% // Export daily data (rainfall, temperature, wind speed) to CSV
% Export.table.toDrive({
%   collection: dailyData,
%   description: 'Gauge_Daily_Rainfall_Temperature_Wind',
%   folder: Folder_Name,
%   selectors: ['date', 'x_mercator', 'y_mercator', 'rainfall', 'temperature_mean', 'temperature_max', 'temperature_min', 'wind_speed'],
%   fileFormat: 'CSV'
% });
% 
% 
% // Function to create time series chart for Rainfall and Temperature (Daily)
% var plotDailyData = function(dailyData, parameter, title) {
%   var chart = ui.Chart.feature.byFeature(dailyData, 'date')
%     .setChartType('LineChart')
%     .setOptions({
%       title: title,
%       vAxis: {title: parameter},
%       hAxis: {title: 'Date'},
%       lineWidth: 2,
%       pointSize: 4,
%       series: {
%         0: {color: 'green'}
%       }
%     });
%   print(chart);
% };
% 
% // Plot Daily Rainfall and Mean Temperature
% plotDailyData(dailyData.filter(ee.Filter.notNull(['rainfall'])), 'Rainfall (mm)', 'Daily Rainfall Timeseries');
% plotDailyData(dailyData.filter(ee.Filter.notNull(['temperature_mean'])), 'Mean Temperature (°C)', 'Daily Mean Temperature Timeseries');
% plotDailyData(dailyData.filter(ee.Filter.notNull(['wind_speed'])), 'Wind Speed (m/s)', 'Daily Wind Speed Timeseries');
