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
options = HydroPol2D_Voronoi_Options(1,optional(case_definition,'options',struct()));
validate_mesh_metadata(case_definition.mesh_file,options);

config = case_definition.config;
if ~isfield(config,'routing_solver') || isempty(config.routing_solver)
    config.routing_solver = 'local_inertial';
end
assert(any(strcmpi(config.routing_solver, {'local_inertial','kinematic','diffusive','full_momentum'})), ...
    'HydroPol2D:VoronoiSolverUnavailable', ...
    'Voronoi routing_solver must be local_inertial, kinematic, diffusive, or full_momentum.');
config.compute_backend = 'cpu';
config.minimum_cell_width_m = options.minimum_cell_width_m;
config.maximum_adjacent_size_ratio = options.maximum_adjacent_size_ratio;
if ~isfield(config,'output_netcdf') || isempty(config.output_netcdf)
    config.output_netcdf = fullfile(output_directory,[prefix '-results.nc']);
end
config.overwrite_output = true;

forcing = case_definition.forcing;
overlap_file = optional(case_definition,'overlap_file','');
if ~isempty(overlap_file), forcing.overlap_file = overlap_file; end
results = HydroPol2D_Voronoi_Run(case_definition.mesh_file,config,forcing);

assert(all(isfinite(results.surface_depth_m),'all') && all(results.surface_depth_m >= 0,'all'), ...
    'HydroPol2D:VoronoiAcceptance','Surface depth contains invalid values.');
assert(all(isfinite(results.surface_velocity_m_s),'all') && all(results.surface_velocity_m_s >= 0,'all'), ...
    'HydroPol2D:VoronoiAcceptance','Surface velocity contains invalid values.');
d = results.diagnostics;
recorded_input = sum([d.precipitation_volume_m3]) + ...
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
if write_rasters
    assert(~isempty(overlap_file), 'HydroPol2D:InvalidVoronoiCase', ...
        'Raster export requires overlap_file.');
    results.raster_files = HydroPol2D_Export_Voronoi_Rasters( ...
        case_definition.mesh_file,overlap_file,results,output_directory,prefix);
end
end

function value = optional(source,name,default_value)
if isfield(source,name) && ~isempty(source.(name)), value=source.(name); else, value=default_value; end
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
if any(attribute_names == "preferred_cells_across")
    saved_cells=double(ncreadatt(mesh_file,'/','preferred_cells_across'));
    assert(saved_cells==options.river_preferred_cells_across, ...
        'HydroPol2D:VoronoiMeshConfigurationMismatch', ...
        'Prepared mesh metadata does not match river_preferred_cells_across.');
end
if any(attribute_names == "unresolved_river_policy")
    saved_policy=char(ncreadatt(mesh_file,'/','unresolved_river_policy'));
    assert(strcmp(saved_policy,options.unresolved_river_policy), ...
        'HydroPol2D:VoronoiMeshConfigurationMismatch', ...
        'Prepared mesh unresolved-river policy does not match the selected policy.');
end
end
