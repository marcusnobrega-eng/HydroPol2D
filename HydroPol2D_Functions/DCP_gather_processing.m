function [DCP_data] = DCP_gather_processing(labels,date_start)
%Assessing which case is to download the data
    if  isempty(labels) == 1 %% HADS data collector
            % Define the URL
            url = 'https://hads.ncep.noaa.gov/nexhads2/servlet/DecodedData?sinceday=-6&hsa=nil&state=HN&nesdis_ids=nil&of=1';
            
            % Fetch the data from the URL
            trier = 1;
            while trier==1
                try
                    data = webread(url);
                    trier=0;
                catch
                    fprintf('Lost connection with HADS server, retrying...');
                    pause(5);
                end
            end
            % Split the data into rows based on newline characters
            rows = strsplit(data, '\n');
            % Initialize an empty cell array to store the data
            dataCell = cell(length(rows), 6); % Assuming you have 6 columns
        
            % Split each row into columns based on the "|" delimiter
            for i = 1:length(rows)
                columns = strsplit(rows{i}, '|');
                % Check if the row has the expected number of columns
                if numel(columns) >= 6
                    dataCell(i, :) = columns(1:6); % Store the first 6 columns
                end
            end
            %Droping the last row because is empty
            dataCell(end,:) = [];
        
            % Create a MATLAB table from the cell array
            temp = cell2table(dataCell, 'VariableNames', {'NESDIS ID', 'NWSLI', 'SHEF Code', 'Time', 'Value', '-'}); % Assuming you have 6 columns
            % Convert the 'Value' column from string to numeric
            temp.Value = str2double(temp.Value);
            % Convert the 'Time' column from string to datetime
            temp.Time = datetime(temp.Time, 'InputFormat', 'yyyy-MM-dd HH:mm');
            %Filtering the data
            logicalMask = strcmp(temp{:, 'SHEF Code'}, 'HG');
            temp = temp(logicalMask, {'Time', 'Value', 'NWSLI'});
            % Identify unique NWSLI values
            uniqueNWSLI = unique(temp.NWSLI);
            % Initialize a new table for the pivoted data
            % Initialize a new table for the pivoted data
            numCols = numel(uniqueNWSLI) + 1; % Include one column for 'Time' and 'Value'
            varNames = ['Time'; uniqueNWSLI];
            varTypes = [{'datetime'}, repmat({'double'}, 1, numel(uniqueNWSLI))]; % Use 'cell' data type for all columns
            % Get unique times directly from the 'Time' column
            uniqueTimes = unique(temp.Time);
            % Create an array of NaN values with the same size as the table
            % Create the pivotedTable specifying variable names
            pivotedTable = table('Size', [height(uniqueTimes), numCols], 'VariableNames', varNames, 'VariableType', varTypes);
            % Loop through the numeric columns (excluding the 'Time' column) and set them to NaN
            for i = 2:numCols
                pivotedTable{:, i} = NaN;
            end
        
            for i = 1:numel(uniqueNWSLI)
                currentNWSLI = uniqueNWSLI{i};
                idx = strcmp(temp.NWSLI, currentNWSLI); % Use strcmp to compare cell values
                nwsliTable = temp(idx, {'Time', 'Value'});
                for j = 1:length(nwsliTable.Time)
                    % Find the row with the matching time
                    timeMatch_index = find(ismember(uniqueTimes,nwsliTable.Time(j))); % Use strcmpi for case-insensitive comparison
                    pivotedTable{timeMatch_index,"Time"} = nwsliTable{j,"Time"};
                    pivotedTable{timeMatch_index,currentNWSLI} = nwsliTable{j,"Value"};
                end
            end
            DCP_data = pivotedTable;
    else  %% ANA data collector
        T = table();
        temp_names={};
        for i =1:length(labels)
            % Define the URL
            url = strcat('http://telemetriaws1.ana.gov.br/ServiceANA.asmx/DadosHidrometeorologicos?CodEstacao=',labels{i}{1}(1:8),'&DataInicio=',string(date_start),'&DataFim=',string(datetime()));
    
            % Fetch the data from the URL
            % request data (html)
            trier = 1;
            while trier==1
                try
                    html = webread(url);
                    trier=0;
                catch
                    fprintf('Lost connection with the ANA server, retrying...');
                    pause(5);
                end
            end
            % Gather each row of codes
            code = regexp(html,'<CodEstacao>(.*?)</CodEstacao>','tokens');
            code = cellfun(@(x) x{1}, code, 'UniformOutput', false);
            code = cellfun(@str2num,code');
            % Gather each row of dates
            date = regexp(html,'<DataHora>(.*?)</DataHora>','tokens');
            date = cellfun(@(x) x{1}, date, 'UniformOutput', false);
            date = cellfun(@datetime, date');
            % Filling flow missing data, and gather each row
            html = regexprep(html,'<Vazao />','<Vazao> </Vazao>');
            flow = regexp(html,'<Vazao>(.*?)</Vazao>','tokens');
            flow = cellfun(@(x) x{1}, flow, 'UniformOutput', false);
            flow = cellfun(@str2double,flow');
            % Filling stage missing data, and gather each row
            html = regexprep(html,'<Nivel />','<Nivel> </Nivel>');
            stage = regexp(html,'<Nivel>(.*?)</Nivel>','tokens');
            stage = cellfun(@(x) x{1}, stage, 'UniformOutput', false);
            stage = cellfun(@str2double,stage')./100; %data comes in cm
            % Filling chuva missing data, and gather each row
            html = regexprep(html,'<Chuva />','<Chuva> </Chuva>');
            rainfall = regexp(html,'<Chuva>(.*?)</Chuva>','tokens');
            rainfall = cellfun(@(x) x{1}, rainfall, 'UniformOutput', false);
            rainfall = cellfun(@str2double,rainfall');
            % Create table
            T = [T;table(date, stage,repmat(labels{i}, numel(date),1))];
        end
        T.Properties.VariableNames{'date'} = 'Time';
        % Initialize a new table for the pivoted data
        % numCols = numel(labels) + 1 - missing; % Include one column for 'Time' and 'Value'
        for k = 1:length(labels)
            string_to_add = labels{k}{1};
            temp_names{end+1} = string_to_add;
        end
        temp_names=temp_names';
        varNames = ['Time';temp_names];
        varTypes = [{'datetime'}, repmat({'double'}, 1, numel(labels))]; % Use 'cell' data type for all columns
        % Get unique times directly from the 'Time' column
        uniqueTimes = unique(T.Time);
        % Create an array of NaN values with the same size as the table
        % Create the pivotedTable specifying variable names
        pivotedTable = table('Size', [height(uniqueTimes), numel(labels) + 1], 'VariableNames', varNames, 'VariableType', varTypes);
        % Loop through the numeric columns (excluding the 'Time' column) and set them to NaN
        for i = 2:numel(labels) + 1
            pivotedTable{:, i} = NaN;
        end

        for i = 1:numel(labels)
            currentNWSLI = labels{i};
            idx = strcmp(T.Var3, currentNWSLI); % Use strcmp to compare cell values
            nwsliTable = T(idx, {'Time', 'stage'});
            for j = 1:length(nwsliTable.Time)
                % Find the row with the matching time
                timeMatch_index = find(ismember(uniqueTimes,nwsliTable.Time(j))); % Use strcmpi for case-insensitive comparison
                pivotedTable{timeMatch_index,"Time"} = nwsliTable{j,"Time"};
                pivotedTable{timeMatch_index,currentNWSLI} = nwsliTable{j,"stage"};
            end
        end
        DCP_data = pivotedTable; 
    end
end























