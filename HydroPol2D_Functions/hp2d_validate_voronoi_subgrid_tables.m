function summary = hp2d_validate_voronoi_subgrid_tables(tables)
%HP2D_VALIDATE_VORONOI_SUBGRID_TABLES Structural checks on HEC-RAS style tables.
%
%   Separate from hp2d_validate_subgrid_tables, which checks the SFINCS-style
%   structured-grid tables and expects entirely different field names.

n_cells = numel(tables.cell_datum_m);
n_edges = numel(tables.face_datum_m);
summary = struct();

rows = (1:n_cells).';
last = max(tables.cell_point_count, 1);
lin = sub2ind(size(tables.cell_volume_m3), rows, last);
summary.cell_volume_nonmonotonic = 0;
for i = 1:n_cells
    k = max(tables.cell_point_count(i), 1);
    v = tables.cell_volume_m3(i, 1:k);
    summary.cell_volume_nonmonotonic = summary.cell_volume_nonmonotonic + ...
        nnz(diff(v) < -1e-9);
end
summary.cell_negative_area = nnz(tables.cell_wet_area_m2 < -eps);
summary.face_negative_area = nnz(tables.face_flow_area_m2 < -eps);
summary.face_negative_conveyance = nnz(tables.face_conveyance < -eps);
summary.cell_first_volume_nonzero = nnz(abs(tables.cell_volume_m3(:,1)) > 1e-9);
summary.subgrid_cells = nnz(tables.cell_is_subgrid);
summary.n_cells = n_cells;
summary.n_edges = n_edges;

if summary.cell_volume_nonmonotonic > 0
    error('HydroPol2D:NonMonotonicSubgridVolume', ...
        'Cell elevation-volume curves have %d monotonicity violations.', ...
        summary.cell_volume_nonmonotonic);
end
if summary.cell_negative_area > 0 || summary.face_negative_area > 0 || ...
        summary.face_negative_conveyance > 0
    error('HydroPol2D:NegativeSubgridGeometry', ...
        'Sub-grid tables contain negative area or conveyance.');
end
if summary.cell_first_volume_nonzero > 0
    error('HydroPol2D:SubgridVolumeOrigin', ...
        'Cell volume curves must start at zero volume; %d do not.', ...
        summary.cell_first_volume_nonzero);
end
end
