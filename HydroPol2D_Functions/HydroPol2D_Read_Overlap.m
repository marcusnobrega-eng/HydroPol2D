function mapping = HydroPol2D_Read_Overlap(path, mesh_file, options)
%HYDROPOL2D_READ_OVERLAP Read conservative raster/polygon overlap weights.
%
%   mapping.mesh_to_raster is CONSERVATIVE: its weights are
%   overlap_area/raster_area, so for a raster cell only partly covered by the
%   mesh they sum to the coverage fraction rather than to one.  That is correct
%   for an EXTENSIVE field and wrong for an INTENSIVE one -- remapping a constant
%   returns constant*coverage.  Use hp2d_remap_intensive for depth, stage,
%   velocity and any per-area quantity.
%
%   mapping.raster_coverage gives that fraction per raster cell, flattened
%   south-up to match mesh_to_raster.

arguments
    path (1,:) char
    mesh_file (1,:) char
    options.allow_legacy (1,1) logical = false
end

mapping.contract_report = HydroPol2D_Validate_Mesh_Bundle( ...
    mesh_file, path, '', allow_legacy=options.allow_legacy);

mesh_index = double(ncread(path, 'overlap_mesh_index')) + 1;
raster_index = double(ncread(path, 'overlap_raster_index')) + 1;
area = double(ncread(path, 'overlap_area_m2'));
mesh_area = double(ncread(path, 'mesh_area_m2'));
raster_area = double(ncread(path, 'raster_area_m2'));
mapping.raster_to_mesh = sparse(mesh_index, raster_index, area ./ mesh_area(mesh_index), numel(mesh_area), numel(raster_area));
mapping.mesh_to_raster = sparse(raster_index, mesh_index, area ./ raster_area(raster_index), numel(raster_area), numel(mesh_area));
mapping.mesh_area = mesh_area;
mapping.raster_area = raster_area;
mapping.raster_coverage = accumarray(raster_index, area, [numel(raster_area) 1], @sum, 0) ...
    ./ max(raster_area, realmin);
mapping.raster_rows = double(ncreadatt(path, '/', 'raster_rows'));
mapping.raster_cols = double(ncreadatt(path, '/', 'raster_cols'));
mapping.x_edges = double(ncread(path, 'x_edges'));
mapping.y_edges = double(ncread(path, 'y_edges'));
end
