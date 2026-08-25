function mapping = HydroPol2D_Read_Overlap(path)
%HYDROPOL2D_READ_OVERLAP Read conservative raster/polygon overlap weights.

mesh_index = double(ncread(path, 'overlap_mesh_index')) + 1;
raster_index = double(ncread(path, 'overlap_raster_index')) + 1;
area = double(ncread(path, 'overlap_area_m2'));
mesh_area = double(ncread(path, 'mesh_area_m2'));
raster_area = double(ncread(path, 'raster_area_m2'));
mapping.raster_to_mesh = sparse(mesh_index, raster_index, area ./ mesh_area(mesh_index), numel(mesh_area), numel(raster_area));
mapping.mesh_to_raster = sparse(raster_index, mesh_index, area ./ raster_area(raster_index), numel(raster_area), numel(mesh_area));
mapping.mesh_area = mesh_area;
mapping.raster_area = raster_area;
mapping.raster_rows = double(ncreadatt(path, '/', 'raster_rows'));
mapping.raster_cols = double(ncreadatt(path, '/', 'raster_cols'));
mapping.x_edges = double(ncread(path, 'x_edges'));
mapping.y_edges = double(ncread(path, 'y_edges'));
end
