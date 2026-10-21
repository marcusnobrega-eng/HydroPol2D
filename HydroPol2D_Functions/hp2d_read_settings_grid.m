function GD = hp2d_read_settings_grid(workbookPath)
%HP2D_READ_SETTINGS_GRID The settings tabs of a case workbook stacked into one grid.
%   HydroPol2D looks inputs up by label (xlget), not by cell address, so a workbook may
%   spread them over tabs. The grid stacks, in this order and when present: the settings
%   tabs of the single case workbook (tools/workbook_builder), the older themed sheets,
%   and General_Data. The class-table and time-series tabs are not part of it
%   (hp2d_read_case_table reads them). A workbook with none of these sheets is read
%   as General_Data, as before.
settingsSheets = { ...
    'Run control', 'Input maps', 'Soil and groundwater', 'Surface flow', ...
    'Rain and ET', 'Design storm', 'Water quality', 'Risk and gauges', ...
    'Run_Control', 'Static_Maps', 'Hydrology_Groundwater', 'Hydrodynamics_Boundaries', ...
    'Water_Quality', 'Human_Risk', 'Observation_Points', 'Spatial_Forcing', ...
    'Design_Storm', 'General_Data'};
try
    available = cellstr(string(sheetnames(workbookPath)));
catch
    available = {'General_Data'};
end
parts = {};
maxCols = 0;
for i = 1:numel(settingsSheets)
    hit = find(strcmpi(available, settingsSheets{i}), 1);
    if isempty(hit)
        continue
    end
    C = readcell(workbookPath, 'Sheet', available{hit});
    if isempty(C)
        continue
    end
    parts{end + 1} = C; %#ok<AGROW>
    maxCols = max(maxCols, size(C, 2));
end
if isempty(parts)
    GD = readcell(workbookPath, 'Sheet', 'General_Data');
    return
end
if isscalar(parts)
    GD = parts{1};   % one sheet (an older General_Data workbook): exactly what readcell gives
    return
end
GD = {};
for i = 1:numel(parts)
    C = parts{i};
    if size(C, 2) < maxCols
        C(:, end + 1:maxCols) = {[]};
    end
    GD = [GD; C; cell(2, maxCols)]; %#ok<AGROW>
end
end
