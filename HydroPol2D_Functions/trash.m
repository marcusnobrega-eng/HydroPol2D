% Define the bounding box for Arizona (xmin, ymin, xmax, ymax)
bbox = [113.5, 31.3, 115.0, 37.0];  % Approximate coordinates for Arizona

% USGS parameter code for river stage (gage height)
param_code = '00065';  % Gage height

% Construct the URL for the USGS API to get active stations within the bounding box
url = sprintf(['https://waterservices.usgs.gov/nwis/site/?format=rdb', ...
               '&bBox=%.4f,%.4f,%.4f,%.4f', ...
               '&parameterCd=%s&siteStatus=active'], ...
               bbox(1), bbox(2), bbox(3), bbox(4), param_code);

% Fetch data from USGS API
options = weboptions('Timeout', 30);
dataText = webread(url, options);

% Split the data into lines and remove comments
lines = splitlines(dataText);
lines = lines(~startsWith(lines, '#') & ~startsWith(lines, '5s'));

% Check if we have valid data
if numel(lines) < 3
    error('No valid data found from USGS API response.');
end

% Extract the header line and split by tab to get the column names
headerLine = lines{1};
headers = split(headerLine, '\t');

% Display headers for debugging
disp('Headers in the response:');
disp(headers);

% Extract the data rows (station information)
dataLines = lines(3:end);
dataFields = cellfun(@(l) split(l, '\t'), dataLines, 'UniformOutput', false);

% Initialize arrays to store station data
nStations = numel(dataFields);
site_no = strings(nStations, 1);
station_name = strings(nStations, 1);
lat = zeros(nStations, 1);
lon = zeros(nStations, 1);

% Find the indices of the relevant columns
latIndex = find(strcmp(headers, 'dec_lat_va'), 1);
lonIndex = find(strcmp(headers, 'dec_long_va'), 1);
siteIndex = find(strcmp(headers, 'site_no'), 1);
nameIndex = find(strcmp(headers, 'station_nm'), 1);

% Check if we have all necessary columns
if isempty(latIndex) || isempty(lonIndex) || isempty(siteIndex) || isempty(nameIndex)
    error('Could not find necessary columns in the header.');
end

% Extract the station data
for i = 1:nStations
    row = dataFields{i};
    site_no(i) = row{siteIndex};
    station_name(i) = row{nameIndex};
    lat(i) = str2double(row{latIndex});
    lon(i) = str2double(row{lonIndex});
end

% Store the data in a table
stationsTable = table(site_no, station_name, lat, lon, 'VariableNames', {'Site_No', 'Station_Name', 'Latitude', 'Longitude'});

% Save the result to a CSV file
writetable(stationsTable, 'Arizona_river_stage_stations.csv');

disp(['Found ', num2str(height(stationsTable)), ' USGS river stage stations in Arizona.']);
