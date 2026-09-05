function tables = HydroPol2D_Read_Subgrid_Tables(path, n_cells, n_edges)
%HYDROPOL2D_READ_SUBGRID_TABLES Read HEC-RAS style sub-grid property tables.
%
%   These are the tables HydroBathyDEM writes next to an unstructured mesh: per
%   cell an elevation-volume curve, and per face elevation versus flow area,
%   wetted perimeter and conveyance.  Method follows Casulli (2009) as implemented
%   in the HEC-RAS 2D engine.
%
%   This is NOT the same method as Subgrid_Properties_Lookup.m, which builds
%   SFINCS-style grid-averaged havg/nrep tables on a structured coarse raster.
%   That path cannot be applied to an unstructured mesh.  The two coexist.
%
%   tables.cell_is_subgrid selects which cells use the level-pool closure.  A FACE
%   may use the tables only where BOTH of its cells do -- switching the cell
%   closure alone leaves a flat-prism cell referenced to its area-mean bed while
%   the face table measures depth from the face minimum, about 1.9 m lower on a
%   90 m mesh, and a nearly dry cell then hands its faces metres of phantom head.

assert(isfile(path), 'HydroPol2D:MissingSubgridTables', ...
    'Sub-grid table file %s does not exist.', path);

required = {'cell_datum_m','cell_zeta_m','cell_volume_m3','cell_wet_area_m2', ...
    'cell_point_count','cell_plan_area_m2','face_datum_m','face_zeta_m', ...
    'face_flow_area_m2','face_perimeter_m','face_conveyance','face_point_count', ...
    'face_length_m'};
info = ncinfo(path);
present = {info.Variables.Name};
missing = setdiff(required, present);
assert(isempty(missing), 'HydroPol2D:IncompleteSubgridTables', ...
    'Sub-grid tables in %s are missing: %s', path, strjoin(missing, ', '));

for k = 1:numel(required)
    tables.(required{k}) = double(ncread(path, required{k}));
end
% netCDF is written cell-major; MATLAB reads the fastest dimension first, so the
% (cell, point) tables arrive transposed.
for name = ["cell_zeta_m","cell_volume_m3","cell_wet_area_m2", ...
            "face_zeta_m","face_flow_area_m2","face_perimeter_m","face_conveyance"]
    if size(tables.(name), 1) ~= n_cells && size(tables.(name), 2) == n_cells
        tables.(name) = tables.(name).';
    elseif size(tables.(name), 1) ~= n_edges && size(tables.(name), 2) == n_edges
        tables.(name) = tables.(name).';
    end
end
tables.cell_point_count = round(tables.cell_point_count(:));
tables.face_point_count = round(tables.face_point_count(:));
tables.cell_datum_m = tables.cell_datum_m(:);
tables.cell_plan_area_m2 = tables.cell_plan_area_m2(:);
tables.face_datum_m = tables.face_datum_m(:);
tables.face_length_m = tables.face_length_m(:);

if any(strcmp(present, 'cell_is_subgrid'))
    tables.cell_is_subgrid = logical(double(ncread(path, 'cell_is_subgrid')));
    tables.cell_is_subgrid = tables.cell_is_subgrid(:);
else
    tables.cell_is_subgrid = true(n_cells, 1);   % absent means every cell
end

assert(numel(tables.cell_datum_m) == n_cells, 'HydroPol2D:SubgridCellCount', ...
    'Sub-grid tables cover %d cells but the mesh has %d. Rebuild the tables whenever the mesh changes.', ...
    numel(tables.cell_datum_m), n_cells);
assert(numel(tables.face_datum_m) == n_edges, 'HydroPol2D:SubgridFaceCount', ...
    'Sub-grid tables cover %d faces but the mesh has %d. Rebuild the tables whenever the mesh changes.', ...
    numel(tables.face_datum_m), n_edges);
assert(numel(tables.cell_is_subgrid) == n_cells, 'HydroPol2D:SubgridFlagCount', ...
    'cell_is_subgrid covers %d cells but the mesh has %d.', ...
    numel(tables.cell_is_subgrid), n_cells);

summary = hp2d_validate_voronoi_subgrid_tables(tables);
tables.validation = summary;
end
