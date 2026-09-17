function [settings, source_file] = HydroPol2D_Read_Voronoi_Settings(input_excel_file, legacy_general)
%HYDROPOL2D_READ_VORONOI_SETTINGS Load dedicated or legacy Excel settings.

arguments
    input_excel_file (1,:) char
    legacy_general cell = {}
end

source_file = fullfile(fileparts(input_excel_file), 'Voronoi_Settings.xlsx');
if isfile(source_file)
    settings = readcell(source_file, 'Sheet', 'Voronoi_Settings');
    return
end

% Existing projects remain valid while users migrate the Voronoi block out
% of General_Data.xlsx.
source_file = input_excel_file;
if isempty(legacy_general)
    legacy_general = readcell(input_excel_file, 'Sheet', 'General_Data');
end
settings = legacy_general;
end
