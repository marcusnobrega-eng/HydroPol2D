function results = HydroPol2D_Run_Voronoi_Case(case_definition, output_directory, write_rasters)
%HYDROPOL2D_RUN_VORONOI_CASE Execute one prepared UGRID case end to end.

arguments
    case_definition struct
    output_directory (1,:) char
    write_rasters (1,1) logical = true
end
required = {'mesh_file','config','forcing'};
for k = 1:numel(required)
    assert(isfield(case_definition,required{k}), 'HydroPol2D:InvalidVoronoiCase', ...
        'Prepared Voronoi case is missing %s.', required{k});
end
if exist(output_directory,'dir') ~= 7, mkdir(output_directory); end
prefix = optional(case_definition,'output_prefix','voronoi');
raw_options=optional(case_definition,'options',struct());
options = HydroPol2D_Voronoi_Options(1,raw_options);
validate_mesh_metadata(case_definition.mesh_file,options);
output_controls = output_defaults(optional(case_definition,'output_controls',struct()));
if ~write_rasters
    output_controls.write_final_geotiffs=false;
    output_controls.write_temporal_geotiffs=false;
    output_controls.write_figures=false;
    output_controls.write_videos=false;
    output_controls.write_gauge_hydrographs=false;
end
overlap_file = optional(case_definition,'overlap_file','');
postprocess_requested=output_controls.write_final_geotiffs || ...
    output_controls.write_temporal_geotiffs || output_controls.write_figures || ...
    output_controls.write_videos || output_controls.write_gauge_hydrographs;
if postprocess_requested
    preflight_postprocessing(case_definition.mesh_file,overlap_file,output_controls);
end

config = case_definition.config;
if isfield(raw_options,'voronoi_subgrid_enabled') || ...
        strlength(strtrim(string(options.subgrid_table_path))) > 0
    config.voronoi_subgrid_enabled=options.voronoi_subgrid_enabled;
end
if strlength(strtrim(string(options.subgrid_table_path))) > 0
    config.subgrid_table_path=options.subgrid_table_path;
end
if ~isfield(config,'routing_solver') || isempty(config.routing_solver)
    config.routing_solver = 'local_inertial';
end
assert(any(strcmpi(config.routing_solver, {'local_inertial','kinematic','diffusive','full_momentum'})), ...
    'HydroPol2D:VoronoiSolverUnavailable', ...
    'Voronoi routing_solver must be local_inertial, kinematic, diffusive, or full_momentum.');
config.compute_backend = 'cpu';
config.minimum_cell_width_m = options.minimum_cell_width_m;
config.maximum_adjacent_size_ratio = options.maximum_adjacent_size_ratio;
native_directory=fullfile(output_directory,'Modeling_Results','Native');
if exist(native_directory,'dir')~=7, mkdir(native_directory); end
config.output_netcdf = fullfile(native_directory,'unstructured_results.nc');
config.output_interval_s=output_controls.native_output_interval_s;
config.raster_output_interval_s=output_controls.raster_stack_interval_s;
config.overwrite_output = true;

forcing = case_definition.forcing;
if ~isempty(overlap_file), forcing.overlap_file = overlap_file; end
results = HydroPol2D_Voronoi_Run(case_definition.mesh_file,config,forcing);

assert(all(isfinite(results.surface_depth_m),'all') && all(results.surface_depth_m >= 0,'all'), ...
    'HydroPol2D:VoronoiAcceptance','Surface depth contains invalid values.');
assert(all(isfinite(results.surface_velocity_m_s),'all') && all(results.surface_velocity_m_s >= 0,'all'), ...
    'HydroPol2D:VoronoiAcceptance','Surface velocity contains invalid values.');
d = results.diagnostics;
recorded_input = sum([d.surface_source_volume_m3]) + ...
    sum(max([d.boundary_net_inflow_volume_m3],0));
if isfield(case_definition,'expected_input_volume_m3')
    expected = double(case_definition.expected_input_volume_m3);
    input_error = abs(recorded_input-expected)/max(abs(expected),eps);
    assert(input_error <= optional(case_definition,'maximum_input_volume_error_fraction',1e-10), ...
        'HydroPol2D:VoronoiAcceptance','Recorded forcing volume differs from the prepared case.');
else
    expected = NaN; input_error = NaN;
end
mass_scale = max(recorded_input,1);
if isfinite(expected), mass_scale=max(mass_scale,abs(expected)); end
mass_residual_fraction = max(abs([d.step_mass_residual_m3]))/mass_scale;
assert(mass_residual_fraction <= optional(case_definition,'maximum_mass_residual_fraction',1e-3), ...
    'HydroPol2D:VoronoiAcceptance','Voronoi water-balance residual exceeds the case tolerance.');

summary = struct( ...
    'mesh_file',string(case_definition.mesh_file), ...
    'output_netcdf',string(config.output_netcdf), ...
    'cell_count',results.preflight.n_cells, ...
    'edge_count',results.preflight.n_edges, ...
    'recorded_input_volume_m3',recorded_input, ...
    'expected_input_volume_m3',expected, ...
    'input_volume_error_fraction',input_error, ...
    'maximum_mass_residual_fraction',mass_residual_fraction, ...
    'minimum_dt_s',min([d.dt_s]), ...
    'maximum_dt_s',max([d.dt_s]), ...
    'maximum_surface_depth_m',max([d.max_surface_depth_m]), ...
    'maximum_surface_velocity_m_s',max([d.max_surface_velocity_m_s]), ...
    'maximum_channel_depth_m',max([d.max_channel_depth_m]), ...
    'maximum_channel_velocity_m_s',max([d.max_channel_velocity_m_s]), ...
    'passed',true);
writetable(struct2table(summary),fullfile(output_directory,[prefix '-summary.csv']));
writetable(struct2table(d),fullfile(output_directory,[prefix '-diagnostics.csv']));
save(fullfile(output_directory,[prefix '-prepared-case.mat']),'case_definition');

results.summary = summary;
results.output_netcdf = config.output_netcdf;
results.raster_files = strings(0,1);
if postprocess_requested
    results.postprocessing = HydroPol2D_Postprocess_Voronoi(config.output_netcdf, ...
        case_definition.mesh_file,overlap_file,fullfile(output_directory,'Modeling_Results'), ...
        output_controls);
    results.raster_files=results.postprocessing.raster_files;
end
if ~output_controls.write_native_archive && exist(config.output_netcdf,'file')==2
    delete(config.output_netcdf);
end

function preflight_postprocessing(mesh_file,overlap_file,controls)
assert(~isempty(overlap_file) && exist(overlap_file,'file')==2, ...
    'HydroPol2D:InvalidVoronoiCase','Voronoi post-processing requires overlap_file.');
assert(exist(mesh_file,'file')==2,'HydroPol2D:InvalidVoronoiCase', ...
    'Voronoi post-processing mesh file is unavailable.');
if controls.write_final_geotiffs || controls.write_temporal_geotiffs
    assert(exist('geotiffwrite','file')==2 && exist('maprefcells','file')==2 && ...
        exist('Tiff','class')==8,'HydroPol2D:MissingMappingToolbox', ...
        'Voronoi raster output requires Mapping Toolbox and TIFF support.');
end
if controls.write_figures
    assert(exist('exportgraphics','file')~=0,'HydroPol2D:MissingGraphicsSupport', ...
        'Voronoi figure output requires exportgraphics.');
end
if controls.write_videos
    assert(exist('VideoWriter','class')==8,'HydroPol2D:MissingVideoSupport', ...
        'Voronoi video output requires VideoWriter.');
end
mesh_info=ncinfo(mesh_file,'cell_area_m2');
mesh_cells=mesh_info.Size;
overlap_cells=numel(ncread(overlap_file,'mesh_area_m2'));
assert(mesh_cells==overlap_cells,'HydroPol2D:InvalidVoronoiOutput', ...
    'Overlap weights do not match the selected UGRID mesh.');
end
end

function value = optional(source,name,default_value)
if isfield(source,name) && ~isempty(source.(name)), value=source.(name); else, value=default_value; end
end

function controls=output_defaults(controls)
values=struct('write_native_archive',true,'native_output_interval_s',300, ...
    'write_final_geotiffs',true,'write_temporal_geotiffs',true, ...
    'raster_stack_interval_s',3600,'write_figures',true,'write_videos',true, ...
    'write_gauge_hydrographs',true,'video_fps',8);
names=fieldnames(values);
for k=1:numel(names), if ~isfield(controls,names{k}) || isempty(controls.(names{k})), controls.(names{k})=values.(names{k}); end, end
ratio=controls.raster_stack_interval_s/controls.native_output_interval_s;
assert(ratio>=1 && abs(ratio-round(ratio))<=1e-9*max(ratio,1), ...
    'HydroPol2D:InvalidVoronoiConfiguration', ...
    'raster_stack_interval_s must be an integer multiple of native_output_interval_s.');
end

function validate_mesh_metadata(mesh_file,options)
info = ncinfo(mesh_file);
attribute_names = string({info.Attributes.Name});
numeric = {'background_target_width_m','minimum_cell_width_m', ...
    'maximum_adjacent_size_ratio','urban_target_width_m','urban_transition_buffer_m'};
for k=1:numel(numeric)
    if ~any(attribute_names == numeric{k}), continue; end
    saved=double(ncreadatt(mesh_file,'/',numeric{k}));
    assert(abs(saved-options.(numeric{k})) <= 1e-9*max(abs(saved),1), ...
        'HydroPol2D:VoronoiMeshConfigurationMismatch', ...
        'Prepared mesh metadata does not match %s.',numeric{k});
end
if any(attribute_names == "unresolved_river_policy")
    saved_policy=char(ncreadatt(mesh_file,'/','unresolved_river_policy'));
    assert(strcmp(saved_policy,options.unresolved_river_policy), ...
        'HydroPol2D:VoronoiMeshConfigurationMismatch', ...
        'Prepared mesh unresolved-river policy does not match the selected policy.');
end
end
