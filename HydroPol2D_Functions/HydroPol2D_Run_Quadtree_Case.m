function results = HydroPol2D_Run_Quadtree_Case(case_definition, output_directory, write_outputs)
%HYDROPOL2D_RUN_QUADTREE_CASE Run a square 2:1 UGRID quadtree.
%
% The hydraulic equations remain those of the validated unstructured
% finite-volume solver. This wrapper validates the restricted quadtree
% contract and supplies the matching mesh options.

arguments
    case_definition struct
    output_directory (1,:) char
    write_outputs (1,1) logical = true
end
required = {'mesh_file','fine_resolution_m','coarse_resolution_m', ...
    'urban_refinement_buffer_m','config','forcing'};
for k = 1:numel(required)
    assert(isfield(case_definition,required{k}), ...
        'HydroPol2D:InvalidQuadtreeCase', ...
        'Prepared quadtree case is missing %s.',required{k});
end
fine = double(case_definition.fine_resolution_m);
coarse = double(case_definition.coarse_resolution_m);
ratio = coarse/fine;
assert(fine>0 && coarse>fine && abs(log2(ratio)-round(log2(ratio)))<1e-10, ...
    'HydroPol2D:InvalidQuadtreeCase', ...
    'coarse_resolution_m must be a power-of-two multiple of fine_resolution_m.');
width = double(ncread(case_definition.mesh_file,'cell_target_width_m'));
allowed = fine.*2.^(0:round(log2(ratio)));
valid = arrayfun(@(value) any(abs(value-allowed)<=1e-9*max(value,1)),width);
assert(all(valid),'HydroPol2D:InvalidQuadtreeCase', ...
    'The UGRID mesh contains a cell width outside the configured quadtree levels.');

options = struct( ...
    'background_target_width_m',coarse, ...
    'minimum_cell_width_m',fine, ...
    'maximum_adjacent_size_ratio',2, ...
    'urban_target_width_m',fine, ...
    'urban_transition_buffer_m',double(case_definition.urban_refinement_buffer_m), ...
    'unresolved_river_policy','none');
if isfield(case_definition,'options')
    names=fieldnames(case_definition.options);
    for k=1:numel(names), options.(names{k})=case_definition.options.(names{k}); end
end
case_definition.options=options;
results=HydroPol2D_Run_Voronoi_Case(case_definition,output_directory,write_outputs);
end
