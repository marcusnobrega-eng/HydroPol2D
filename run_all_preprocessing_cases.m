% =========================================================================
% 🧪 Run HydroPol2D Preprocessing for All Valid Flag Combinations
% =========================================================================
clear; clc;

% === Set output folder to save .mat workspaces ===
outputFolder = 'C:\Users\marcu\Documents\GitHub\HydroPol2D\Synthetic_Workspaces';  % ← Customize this!
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% === Define all flag names ===
flagNames = {'Concentrated', 'Spatial_Gauges', 'Spatial_Map', 'Inflow', ...
             'Alternated_Blocks', 'Huff', 'Stage', 'WaterQuality'};

nFlags = numel(flagNames);

% === Manually define compatibility matrix based on final table ===
compatMatrix = [
    1 0 0 1 0 0 1 1;  % Concentrated
    0 1 0 1 0 0 1 1;  % Spatial_Gauges
    0 0 1 1 0 0 1 1;  % Spatial_Map
    1 1 1 1 1 1 1 1;  % Inflow
    0 0 0 1 1 0 1 1;  % Alternated_Blocks
    0 0 0 1 0 1 1 1;  % Huff
    1 1 1 0 1 1 1 1;  % Stage
    1 1 1 1 1 1 1 0]; % WaterQuality

% === Generate valid combinations (1, 2, and 3 flags) ===
flagCombos = {};

% 1-flag
for i = 1:nFlags
    flagCombos{end+1} = {flagNames{i}}; %#ok<SAGROW>
end

% 2-flags
for i = 1:nFlags
    for j = i+1:nFlags
        if compatMatrix(i,j) && compatMatrix(j,i)
            flagCombos{end+1} = sort({flagNames{i}, flagNames{j}}); %#ok<SAGROW>
        end
    end
end

% 3-flags
for i = 1:nFlags
    for j = i+1:nFlags
        for k = j+1:nFlags
            if all([compatMatrix(i,j), compatMatrix(j,i), ...
                    compatMatrix(i,k), compatMatrix(k,i), ...
                    compatMatrix(j,k), compatMatrix(k,j)])
                flagCombos{end+1} = sort({flagNames{i}, flagNames{j}, flagNames{k}}); %#ok<SAGROW>
            end
        end
    end
end

% === Remove duplicates ===
flagCombos = unique(cellfun(@(c) sort(c), flagCombos, 'UniformOutput', false), 'stable');

% === Loop through combinations and run preprocessing ===
for i = 1:length(flagCombos)
    flags = flagCombos{i};

    % Define current flag struct
    f = struct();
    for j = 1:nFlags
        f.(flagNames{j}) = ismember(flagNames{j}, flags);
    end

    % Define case name
    caseName = strjoin(flags, '_');

    fprintf('🧪 [%02d/%02d] Running case: %s\n', i, numel(flagCombos), caseName);

    % Call main preprocessing function
    Validation_Domains(f, caseName);

    % Save .mat workspace
    save(fullfile(outputFolder, ['WS_' caseName '.mat']));
end
