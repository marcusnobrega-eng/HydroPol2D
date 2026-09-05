function selection = HydroPol2D_Voronoi_Selection(run_mode,input_excel_file,bypass_script,InputPaths)
%HYDROPOL2D_VORONOI_SELECTION Read the compact raster/Voronoi switch.

selection = struct('enabled',false,'case_file','','options',struct(), ...
    'output_controls',struct(),'subgrid_choice_explicit',false);
if strcmpi(run_mode,'bypass')
    run(bypass_script);
    assert(exist('InputData_Bypass','var')==1 && isstruct(InputData_Bypass), ...
        'HydroPol2D:InvalidBypass','Bypass script must create InputData_Bypass.');
    flags = InputData_Bypass.flags;
    selection.subgrid_choice_explicit=isfield(flags,'flag_voronoi_subgrid');
    flags=HydroPol2D_Normalize_Subgrid_Flags(flags);
    if isfield(flags,'flag_voronoi'), selection.enabled=logical(flags.flag_voronoi); end
    if isfield(InputData_Bypass,'Voronoi')
        selection.options=InputData_Bypass.Voronoi;
        if isfield(InputData_Bypass.Voronoi,'case_file')
            selection.case_file=char(InputData_Bypass.Voronoi.case_file);
        end
    end
    if isfield(InputData_Bypass,'VoronoiOutput')
        selection.output_controls=InputData_Bypass.VoronoiOutput;
    end
else
    flag_cells = readcell(input_excel_file,'Sheet','Flags');
    [found_voronoi,value_voronoi]=cell_right_found(flag_cells,'flag_voronoi',0);
    [found_voronoi_subgrid,value_voronoi_subgrid]=cell_right_found( ...
        flag_cells,'flag_voronoi_subgrid',0);
    raw_flags=struct();
    if found_voronoi, raw_flags.flag_voronoi=value_voronoi; end
    if found_voronoi_subgrid
        raw_flags.flag_voronoi_subgrid=value_voronoi_subgrid;
        selection.subgrid_choice_explicit=true;
    end
    flags=HydroPol2D_Normalize_Subgrid_Flags(raw_flags);
    selection.enabled = logical(flags.flag_voronoi);
    general = readcell(input_excel_file,'Sheet','General_Data');
    selection.case_file=char(string(cell_right(general,'Voronoi case file','')));
    selection.options=struct( ...
        'background_target_width_m',cell_right(general,'Background width',2000), ...
        'minimum_cell_width_m',cell_right(general,'Minimum cell width',100), ...
        'maximum_adjacent_size_ratio',cell_right(general,'Adjacent size ratio',2), ...
        'urban_target_width_m',cell_right(general,'Urban target width',200), ...
        'urban_transition_buffer_m',cell_right(general,'Urban buffer',1000), ...
        'unresolved_river_policy',char(string(cell_right(general,'Unresolved rivers','neal_subgrid'))), ...
        'subgrid_table_path',char(string(cell_right(general,'Voronoi subgrid table',''))));
    selection.output_controls=struct( ...
        'write_native_archive',logical(cell_right(general,'Write native archive',1)), ...
        'native_output_interval_s',cell_right(general,'Native output interval s',300), ...
        'write_final_geotiffs',logical(cell_right(general,'Write final GeoTIFFs',1)), ...
        'write_temporal_geotiffs',logical(cell_right(general,'Write temporal GeoTIFFs',1)), ...
        'raster_stack_interval_s',cell_right(general,'Raster stack interval s',3600), ...
        'write_figures',logical(cell_right(general,'Write figures',1)), ...
        'write_videos',logical(cell_right(general,'Write videos',1)), ...
        'write_gauge_hydrographs',logical(cell_right(general,'Write gauge hydrographs',1)));
end
if selection.subgrid_choice_explicit
    selection.options.voronoi_subgrid_enabled=logical(flags.flag_voronoi_subgrid);
end
selection.options=HydroPol2D_Voronoi_Options(double(selection.enabled),selection.options);
if selection.enabled && isempty(strtrim(selection.case_file))
    error('HydroPol2D:MissingVoronoiCase', ...
        'flag_voronoi=1 requires a prepared Voronoi case MAT file.');
end
end

function [found,value] = cell_right_found(cells,label,default_value)
found=false;
value=default_value;
for row=1:size(cells,1)
    for column=1:size(cells,2)-1
        entry=cells{row,column};
        if (ischar(entry) || isstring(entry)) && strcmpi(strtrim(string(entry)),label)
            found=true;
            candidate=cells{row,column+1};
            if ~(isempty(candidate) || (isstring(candidate) && ismissing(candidate)))
                value=candidate;
            end
            return
        end
    end
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
