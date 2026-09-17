example_root = fileparts(mfilename('fullpath'));
model_root = fileparts(fileparts(example_root));
workbook = fullfile(example_root,'Input_Data_Sheets','General_Data.xlsx');

flag_data = readcell(workbook,'Sheet','Flags');
[flag_row,flag_col] = find(strcmpi(string(flag_data),'flag_voronoi'),1);
if isempty(flag_row) || flag_col >= size(flag_data,2)
    error('HydroPol2D:MissingVoronoiFlag', ...
        'The Flags sheet must contain flag_voronoi and its value to the right.');
end
flag_value = flag_data{flag_row,flag_col + 1};
if isnumeric(flag_value)
    flag_voronoi = double(flag_value);
else
    flag_voronoi = str2double(string(flag_value));
end
if ~ismember(flag_voronoi,[0 1])
    error('HydroPol2D:InvalidVoronoiFlag', ...
        'flag_voronoi must be 0 (regular grid) or 1 (Voronoi).');
end

if flag_voronoi == 1
    output_name = 'voronoi';
else
    output_name = 'regular';
end

setenv('HYDROPOL_RUN_MODE','excel');
setenv('HYDROPOL_INPUT_EXCEL_FILE',workbook);
setenv('HYDROPOL_EXPORT_ROOT_DIR',fullfile(example_root,'Outputs',output_name));
setenv('HYDROPOL_SKIP_POSTPROCESS','');

run(fullfile(model_root,'HydroPol2D_V115.m'));
