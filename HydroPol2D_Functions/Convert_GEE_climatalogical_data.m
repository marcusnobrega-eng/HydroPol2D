%% Input data

% Rainfall directory
rainfall_path = 'G:\My Drive\Beaver_Creek_30m_correct\Gauge_Rainfall_IMERG_Optimized.csv';

% Climatological directory
climate_path = 'G:\My Drive\Beaver_Creek_30m_correct\Gauge_Daily_Rainfall_Temperature_Wind.csv';

%% Rainfall Routine

% Assuming your data is in a CSV file named 'rainfall_data.csv'
% with columns: date, x_mercator, y_mercator, mean

% Step 1: Load the data
data = readtable(rainfall_path);

% Step 2: Extract unique stations based on coordinates
stations = unique(data(:, {'x_mercator', 'y_mercator'}), 'rows');
station_names = strcat('Gauge ', string(1:height(stations)));

% Step 3: Create the first table for station coordinates
station_table = table(station_names', stations.x_mercator, stations.y_mercator, ...
    'VariableNames', {'Rain Gauge', 'Easting [m]', 'Northing [m]'});

% Reshape the station table so that each station gets its own column (e.g., Gauge 1, Gauge 2...)
station_table_reshaped = table('Size', [2, height(stations)], ...
    'VariableTypes', repmat({'double'}, 1, height(stations)), ...
    'VariableNames', station_names');

station_table_reshaped(1, :) = num2cell(stations.x_mercator');
station_table_reshaped(2, :) = num2cell(stations.y_mercator');

% Step 4: Organize the data by each unique station (Optimized)
timestamps = unique(data.date);
num_timestamps = length(timestamps);
num_stations = height(stations);

% Preallocate the rainfall data matrix (faster than dynamic appending)
rainfall_data = NaN(num_timestamps, num_stations);  % Use NaN to represent missing values

% Set up the status display
fprintf('Processing %d timestamps across %d stations...\n', num_timestamps, num_stations);

% Get the x and y coordinates of all stations for efficient matching
station_coords = [stations.x_mercator, stations.y_mercator];

% Step 5: Process the data
for t = 1:num_timestamps
    % Status update: Print progress every 10%
    if mod(t, ceil(num_timestamps / 10)) == 0
        fprintf('Progress: %.0f%%\n', (t / num_timestamps) * 100);
    end
    
    % Get the current time slice of data
    current_time_data = data(data.date == timestamps(t), :);
    
    % Use ismember to find matching stations for current timestamp
    [~, idx_stations] = ismember([current_time_data.x_mercator, current_time_data.y_mercator], station_coords, 'rows');
    
    % Fill in the rainfall data for the stations that match this timestamp
    for idx = 1:length(idx_stations)
        if idx_stations(idx) > 0
            rainfall_data(t, idx_stations(idx)) = current_time_data.mean(idx);
        end
    end
end

% Step 6: Create the second table with rainfall intensity over time for each station
rainfall_intensity_table = array2table(rainfall_data, ...
    'VariableNames', station_names);

% % Display the station table and rainfall intensity table
% disp('Station Coordinates Table (Reshaped):');
% disp(station_table_reshaped);
% 
% disp('Rainfall Intensity Table:');
% disp(rainfall_intensity_table);

% Step 7: Statistical metrics

% Convert rainfall intensity from mm/h to total volume per time step (volume = intensity * time)
% Since we assume the time step is 1 hour, volume per hour = rainfall intensity in mm
% We will also compute yearly metrics assuming that the data spans over multiple years.

% For simplicity, assume that the data is equally distributed across the years
% (this can be adjusted based on your dataset if necessary).

% Assuming a time step of 360 minutes (1 hour)
time_step_hours = 1;

% Preallocate arrays to store the calculated statistics for each station
mean_rainfall = NaN(1, num_stations);
std_rainfall = NaN(1, num_stations);
total_volume = NaN(1, num_stations);
avg_yearly_volume = NaN(1, num_stations);
max_intensity = NaN(1, num_stations);

% Calculate statistics for each station
for s = 1:num_stations
    % Get the rainfall intensity data for the current station
    station_rainfall = rainfall_data(:, s);
    
    % Remove NaN values (missing data)
    station_rainfall = station_rainfall(~isnan(station_rainfall));
    
    % Mean rainfall intensity (mm/h)
    mean_rainfall(s) = mean(station_rainfall);
    
    % Standard deviation of rainfall intensity (mm/h)
    std_rainfall(s) = std(station_rainfall);
    
    % Total rainfall volume (mm)
    total_volume(s) = sum(station_rainfall) * time_step_hours;  % Sum of intensities 
    
    % Assuming 365 days per year and that the dataset spans a full year
    total_volume_per_year = total_volume(s);  % Total volume in one year (mm)
    avg_yearly_volume(s) = total_volume_per_year / length(unique(year(timestamps)));
    
    % Maximum rainfall intensity (mm/h)
    max_intensity(s) = max(station_rainfall);
end

% Step 8: Display the statistical metrics for each station
metrics_table = table(station_names', mean_rainfall', std_rainfall', total_volume', avg_yearly_volume', max_intensity', ...
    'VariableNames', {'Rain Gauge', 'Mean Rainfall (mm/h)', 'Std Rainfall (mm/h)', 'Total Volume (mm)', 'Avg Yearly Volume (mm)', 'Max Intensity (mm/h)'});

disp('Statistical Metrics for Each Station:');
disp(metrics_table);

% Optional: Save the metrics to CSV file
writetable(metrics_table, 'rainfall_metrics.csv');

% Optional: Save the tables to CSV files
writetable(station_table_reshaped, 'station_coordinates_rainfall.csv');
writetable(rainfall_intensity_table, 'rainfall_data_organized.csv');

% Create plots
figure;
tiledlayout(2, 2);  % Creating a 2x2 grid of subplots

% Plot 1: Rainfall intensity over time for the first station
nexttile;
plot(timestamps, rainfall_data(:,1), 'LineWidth', 2, 'Color', [0.2 0.6 1]);  % First station's rainfall intensity
xlabel('Timestamp');
ylabel('Rainfall Intensity (mm/h)');
title(['Rainfall Intensity: ', station_names{1}]);
grid on;
set(gca, 'FontName', 'Montserrat', 'TickDir', 'out');  % Set font and ticks outside

% Plot 2: Rainfall intensity over time for the second station
nexttile;
plot(timestamps, rainfall_data(:,2), 'LineWidth', 2, 'Color', [1 0.6 0]);  % Second station's rainfall intensity
xlabel('Timestamp');
ylabel('Rainfall Intensity (mm/h)');
title(['Rainfall Intensity: ', station_names{2}]);
grid on;
set(gca, 'FontName', 'Montserrat', 'TickDir', 'out');  % Set font and ticks outside

% Plot 3: Total volume per station
nexttile;
bar(avg_yearly_volume, 'FaceColor', [0.1 0.8 0.3]);
xlabel('Station');
ylabel('Yearly Average Volume (mm)');
title('Mean Rainfall for Each Station');
xticklabels(station_names);
grid on;
set(gca, 'FontName', 'Montserrat', 'TickDir', 'out');  % Set font and ticks outside

% Plot 4: Max intensity per station
nexttile;
bar(max_intensity, 'FaceColor', [0.8 0.1 0.1]);
xlabel('Station');
ylabel('Max Intensity (mm/h)');
title('Max Rainfall Intensity for Each Station');
xticklabels(station_names);
grid on;
set(gca, 'FontName', 'Montserrat', 'TickDir', 'out');  % Set font and ticks outside

% Optional: Save the figure as a PNG file
saveas(gcf, 'rainfall_plots.png');


%% Climatological Data

% Step 1: Load the data (assuming the CSV file is named 'climate_data.csv')
data = readtable(climate_path);  % Update the file path if necessary

% Step 2: Extract unique stations based on latitude and longitude
stations = unique(data(:, {'x_mercator', 'y_mercator'}), 'rows');

% Step 3: Create the table for station coordinates in the desired format
station_indices = (1:height(stations))';  % Create an Index column starting from 1
station_coordinates = table(station_indices, stations.x_mercator, stations.y_mercator);

% Convert latitude/longitude to 'easting' and 'northing'
% Assuming you have a function to convert lat/lon to Mercator coordinates (you may need to replace this if you don't)
% If you already have the Mercator coordinates in the data, you can use those directly.
% For now, I'll assume that 'stations.x_mercator' and 'stations.y_mercator' are the 'easting' and 'northing' values.
station_coordinates.Properties.VariableNames = {'Index', 'x easting (m)', 'y northing (m)'};

% Step 4: Prepare the new table structure for climatological variables
climate_data = [];

% Get all unique dates in the data
dates = unique(data.date);

% Define the column names for the climatological variables
var_names = {'Tmax (°C)', 'Tmin (°C)', 'Tmed (°C)', 'U2 [m/s]', 'UR [%]', 'G (MJ/(m2.dia))'};

% Step 5: Initialize a cell array for unique column names
column_names = {'Date'};  % First column will be Date

% Generate the column names for each station and climatological variable
for s = 1:height(stations)
    for v = 1:length(var_names)
        column_names{end+1} = strcat(var_names{v}, ' ', 'Index ', num2str(s));  % Concatenate variable with station index
    end
end

% Step 6: Loop through each date and gather the data for each station
disp('Processing climatological data...');
for d = 1:length(dates)
    % Show progress (optional)
    if mod(d, 100) == 0  % Display progress every 100 dates
        fprintf('Processing date %d of %d...\n', d, length(dates));
    end
    
    % Get the current time slice of data (for each date)
    current_date_data = data(data.date == dates(d), :);
    
    % Initialize a row for the current date
    row_data = {dates(d)};  % The first element in the row will be the date
    
    % Loop through each station and get the corresponding climatological data
    for s = 1:height(stations)
        % Get the station's lat and lon (Mercator coordinates)
        x_mercator = stations.x_mercator(s);
        y_mercator = stations.y_mercator(s);
        
        % Find the data for the current station and date
        station_data = current_date_data(current_date_data.x_mercator == x_mercator & ...
                                         current_date_data.y_mercator == y_mercator, :);
        
        % If data exists for this station, extract the climatological variables
        if ~isempty(station_data)
            % Extract Tmax, Tmin, Tmed, U2, UR, G (assuming they are in the columns as per the input table)
            Tmax = station_data.temperature_max(1);
            Tmin = station_data.temperature_min(1);
            Tmed = station_data.temperature_mean(1);
            U2 = station_data.wind_speed(1);
            UR = 0;  % Assuming rainfall is used as UR (%)
            G = 0;  % No data provided for G, assuming it as 0 for now

            % Append the climatological data for this station to the row
            row_data = [row_data, Tmax, Tmin, Tmed, U2, UR, G];
        else
            % If no data, use NaN or 0 (depending on your preference)
            row_data = [row_data, NaN, NaN, NaN, NaN, NaN, NaN];  % Use NaN for missing data
        end
    end
    
    % Append the row data to the climatological data array
    climate_data = [climate_data; row_data];
end
column_names = string(column_names);

% Step 7: Create the table for climatological data with stations as columns
climate_table = array2table(climate_data, 'VariableNames', column_names);

% Step 8: Compute statistics for each climatological variable
climate_stats = table;

for v = 2:(size(climate_table,2))  % Skip the first column (Date)
    % Get the column data for the current variable (across all stations and dates)
    data_column = cell2mat(climate_table{:, v});
    
    % Calculate statistics
    climate_stats.(column_names{v}) = struct( ...
        'Mean', nanmean(data_column), ...
        'Std', nanstd(data_column), ...
        'Max', nanmax(data_column), ...
        'Yearly_Avg', nansum(data_column) / length(dates));
end

% Display the station coordinates and climatological data table
% disp('Station Coordinates Table:');
% disp(station_coordinates);
% disp('Climatological Data Table:');
% disp(climate_table);
% disp('Climatological Statistics:');
% disp(climate_stats);

% Optional: Save the tables to CSV files
writetable(station_coordinates, 'station_coordinates_ETP.csv');
writetable(climate_table, 'climatological_data_organized.csv');
save('climate_statistics.mat', 'climate_stats');  % Save statistics as a .mat file

% --- Visualization (Plots) ---
figure;

% Set figure size to A4 dimensions (in inches)
% A4 size in inches is 8.27 x 11.69
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 8.27, 11.69]);

tiledlayout(2, 2);  % Creating a 2x2 grid of subplots

% Check column names to debug
disp(climate_table.Properties.VariableNames);

% Plot 1: Max Temperature (Tmax) for each station over time
nexttile;
hold on;
for s = 1:height(stations)
    % Find the correct column for each station's Tmax
    Tmax_col_name = strcat('Tmax (°C)Index ', num2str(s));
    
    % Check if the column exists in the climate_table (debugging step)
    if ismember(Tmax_col_name, climate_table.Properties.VariableNames)
        Tmax_col = find(strcmp(climate_table.Properties.VariableNames, Tmax_col_name));
        plot(dates, cell2mat(climate_table{:, Tmax_col}), 'LineWidth', 1.5);
    else
        fprintf('Column %s not found.\n', Tmax_col_name);  % Print a message if the column doesn't exist
    end
end
xlabel('Date', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Max Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
title('Max Temperature (Tmax) Over Time for Each Station', 'FontSize', 14);
legend(string(station_indices(:)), 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Montserrat', 'FontSize', 14, 'TickDir', 'out');  % Set font and ticks outside

% Plot 2: Wind Speed (U2) for each station over time
nexttile;
hold on;
for s = 1:height(stations)
    U2_col = find(contains(climate_table.Properties.VariableNames, strcat('U2 [m/s]Index ', num2str(s))));
    plot(dates, cell2mat(climate_table{:, U2_col}), 'LineWidth', 1.5);
end
xlabel('Date', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Wind Speed (m/s)', 'FontSize', 14, 'FontWeight', 'bold');
title('Wind Speed (U2) Over Time for Each Station', 'FontSize', 14);
legend(string(station_indices(:)), 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Montserrat', 'FontSize', 14, 'TickDir', 'out');  % Set font and ticks outside

% Plot 3: Mean Temperature (Tmed) Distribution across stations
nexttile;
hold on;
for s = 1:height(stations)
    Tmed_col = find(contains(climate_table.Properties.VariableNames, strcat('Tmed (°C)Index ', num2str(s))));
    plot(dates, cell2mat(climate_table{:, Tmed_col}), 'LineWidth', 1.5);
end
xlabel('Date', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Mean Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
title('Mean Temperature (Tmed) Over Time for Each Station', 'FontSize', 14);
legend(string(station_indices(:)), 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Montserrat', 'FontSize', 14, 'TickDir', 'out');  % Set font and ticks outside

% Plot 4: Temperature Min (Tmin) for each station over time
nexttile;
hold on;
for s = 1:height(stations)
    Tmin_col = find(contains(climate_table.Properties.VariableNames, strcat('Tmin (°C)Index ', num2str(s))));
    plot(dates, cell2mat(climate_table{:, Tmin_col}), 'LineWidth', 1.5);
end
xlabel('Date', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Min Temperature (°C)', 'FontSize', 14, 'FontWeight', 'bold');
title('Min Temperature (Tmin) Over Time for Each Station', 'FontSize', 14);
legend(string(station_indices(:)), 'Location', 'best', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Montserrat', 'FontSize', 14, 'TickDir', 'out');  % Set font and ticks outside

% Optional: Save the figure as a PNG file
saveas(gcf, 'climatological_plots.png');

