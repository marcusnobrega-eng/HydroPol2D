function selection = HydroPol2D_Voronoi_Selection(run_mode,input_excel_file,bypass_script,InputPaths)
%HYDROPOL2D_VORONOI_SELECTION Read the compact raster/Voronoi switch.

selection = struct('enabled',false,'case_file','','options',struct());
if strcmpi(run_mode,'bypass')
    run(bypass_script);
    assert(exist('InputData_Bypass','var')==1 && isstruct(InputData_Bypass), ...
        'HydroPol2D:InvalidBypass','Bypass script must create InputData_Bypass.');
    flags = InputData_Bypass.flags;
    if isfield(flags,'flag_voronoi'), selection.enabled=logical(flags.flag_voronoi); end
    if isfield(InputData_Bypass,'Voronoi')
        selection.options=InputData_Bypass.Voronoi;
        if isfield(InputData_Bypass.Voronoi,'case_file')
            selection.case_file=char(InputData_Bypass.Voronoi.case_file);
        end
    end
else
    flags = readcell(input_excel_file,'Sheet','Flags');
    selection.enabled = logical(cell_right(flags,'flag_voronoi',0));
    general = readcell(input_excel_file,'Sheet','General_Data');
    selection.case_file=char(string(cell_right(general,'Voronoi case file','')));
    selection.options=struct( ...
        'background_target_width_m',cell_right(general,'Background width',2000), ...
        'minimum_cell_width_m',cell_right(general,'Minimum cell width',100), ...
        'maximum_adjacent_size_ratio',cell_right(general,'Adjacent size ratio',2), ...
        'urban_target_width_m',cell_right(general,'Urban target width',200), ...
        'urban_transition_buffer_m',cell_right(general,'Urban buffer',1000), ...
        'river_preferred_cells_across',cell_right(general,'River cells across',3), ...
        'unresolved_river_policy',char(string(cell_right(general,'Unresolved rivers','neal_subgrid'))));
end
selection.options=HydroPol2D_Voronoi_Options(double(selection.enabled),selection.options);
if selection.enabled && isempty(strtrim(selection.case_file))
    error('HydroPol2D:MissingVoronoiCase', ...
        'flag_voronoi=1 requires a prepared Voronoi case MAT file.');
end
end

function value = cell_right(cells,label,default_value)
value=default_value;
for row=1:size(cells,1)
    for column=1:size(cells,2)-1
        entry=cells{row,column};
        if (ischar(entry) || isstring(entry)) && strcmpi(strtrim(string(entry)),label)
            candidate=cells{row,column+1};
            if ~(isempty(candidate) || (isstring(candidate) && ismissing(candidate)))
                value=candidate;
            end
            return
        end
    end
end
end
