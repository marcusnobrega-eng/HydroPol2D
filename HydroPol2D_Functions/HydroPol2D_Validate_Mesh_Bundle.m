function report = HydroPol2D_Validate_Mesh_Bundle(mesh_file, overlap_file, subgrid_file, options)
%HYDROPOL2D_VALIDATE_MESH_BUNDLE Check a HydroBathyDEM bundle against the contract.
%
%   HydroBathyDEM is the only mesh generator and this model only consumes, so the
%   integration IS the file contract. Until now neither side validated it, and
%   every bug in the 2026-09 sub-grid work lived in that gap: a south-up overlap
%   index read as north-up, boundary edge ids off by one base, face tables whose
%   datum sat below their own cells, tables whose counts silently disagreed with
%   the mesh.
%
%   Each check exists because the corresponding breach happened, and the message
%   says what it cost. Companion validator: hydropol2d.mesh_contract in Python.
%   Specification: UNSTRUCTURED_MESH_CONTRACT.md.
%
%   IMPORTANT on indexing. In the FILE cells are numbered from ZERO and a boundary
%   face carries edge_neighbor = -1. HydroPol2D_Read_UGRID adds one, so cells are
%   1-based here and a boundary face becomes 0. This function reads the RAW file,
%   so it works in file convention throughout.

arguments
    mesh_file (1,:) char
    overlap_file (1,:) char = ''
    subgrid_file (1,:) char = ''
    options.area_tolerance_m2 (1,1) double = 1e-3
    options.allow_legacy (1,1) logical = false
end

FILE_BOUNDARY_NEIGHBOR = -1;
MESH_CONTRACT_VERSION = '1.0';
MESH_PRODUCT_SCHEMA = 'hydrobathydem-mesh-1.0';
COORDINATE_UNITS = 'm';
EDGE_NORMAL_CONVENTION = 'unit normal points outward from edge_owner';
OVERLAP_RASTER_ORDER = 'south_up_row_major';
SUBGRID_CONVEYANCE_CONVENTION = 'K=sum((A/n)*(A/P)^(2/3))';
report = struct();

if exist(mesh_file, 'file') ~= 2
    error('HydroPol2D:ContractMissingMesh', 'Mesh file not found: %s', mesh_file);
end

mesh_required = {'cell_area_m2','cell_bed_elevation_m','cell_hydraulic_roughness', ...
    'edge_owner','edge_neighbor','edge_length_m','edge_center_distance_m', ...
    'edge_midpoint_x','edge_midpoint_y','edge_normal_x','edge_normal_y'};
require_variables(mesh_file, mesh_required, 'Mesh');

strict_contract = has_attribute(mesh_file, 'mesh_contract_version');
if ~strict_contract && ~options.allow_legacy
    error('HydroPol2D:LegacyMeshContractRejected', ...
        ['Mesh %s does not declare mesh_contract_version. Rebuild it with ' ...
         'HydroBathyDEM, or pass allow_legacy=true only for an explicit migration.'], ...
        mesh_file);
end
if strict_contract
    require_attribute(mesh_file, 'mesh_contract_version', MESH_CONTRACT_VERSION, 'Mesh');
    require_attribute(mesh_file, 'schema_version', MESH_PRODUCT_SCHEMA, 'Mesh');
    require_attribute(mesh_file, 'product_type', 'hydraulic_mesh', 'Mesh');
    require_attribute(mesh_file, 'file_index_base', 0, 'Mesh');
    require_attribute(mesh_file, 'boundary_neighbor_sentinel', FILE_BOUNDARY_NEIGHBOR, 'Mesh');
    require_attribute(mesh_file, 'coordinate_units', COORDINATE_UNITS, 'Mesh');
    require_attribute(mesh_file, 'edge_normal_convention', EDGE_NORMAL_CONVENTION, 'Mesh');
    require_nonempty_attribute(mesh_file, 'crs_wkt', 'Mesh');
    require_variables(mesh_file, ...
        {'cell_cfl_width_m','mesh2d_face_x','mesh2d_face_y','mesh2d_face_nodes'}, 'Mesh');
    require_units(mesh_file, 'cell_area_m2', 'm2', 'Mesh');
    require_units(mesh_file, 'mesh2d_face_x', 'm', 'Mesh');
    require_units(mesh_file, 'mesh2d_face_y', 'm', 'Mesh');
    require_units(mesh_file, 'cell_bed_elevation_m', 'm', 'Mesh');
    require_units(mesh_file, 'cell_cfl_width_m', 'm', 'Mesh');
    require_units(mesh_file, 'edge_length_m', 'm', 'Mesh');
    require_units(mesh_file, 'edge_center_distance_m', 'm', 'Mesh');
    report.contract_version = MESH_CONTRACT_VERSION;
else
    report.contract_version = 'legacy';
    report.contract_note = ['Legacy mesh accepted for compatibility. Rebuild it with ' ...
        'HydroBathyDEM before release or parity validation.'];
end

area      = double(ncread(mesh_file,'cell_area_m2'));      area = area(:);
bed       = double(ncread(mesh_file,'cell_bed_elevation_m')); bed = bed(:);
cell_x = []; cell_y = [];
if strict_contract
    cell_x = double(ncread(mesh_file,'mesh2d_face_x')); cell_x = cell_x(:);
    cell_y = double(ncread(mesh_file,'mesh2d_face_y')); cell_y = cell_y(:);
end
owner     = double(ncread(mesh_file,'edge_owner'));        owner = owner(:);
neighbor  = double(ncread(mesh_file,'edge_neighbor'));     neighbor = neighbor(:);
len       = double(ncread(mesh_file,'edge_length_m'));     len = len(:);
distance  = double(ncread(mesh_file,'edge_center_distance_m')); distance = distance(:);
midpoint_x = double(ncread(mesh_file,'edge_midpoint_x')); midpoint_x = midpoint_x(:);
midpoint_y = double(ncread(mesh_file,'edge_midpoint_y')); midpoint_y = midpoint_y(:);
normal_x = double(ncread(mesh_file,'edge_normal_x')); normal_x = normal_x(:);
normal_y = double(ncread(mesh_file,'edge_normal_y')); normal_y = normal_y(:);

n_cells = numel(area); n_faces = numel(owner);
report.n_cells = n_cells; report.n_faces = n_faces;

if any([numel(neighbor), numel(len), numel(distance), numel(midpoint_x), ...
        numel(midpoint_y), numel(normal_x), numel(normal_y)] ~= n_faces)
    error('HydroPol2D:ContractFaceArrays', ...
        ['Mesh %s: edge_owner, edge_neighbor and edge_length_m must have equal ' ...
         'length; got %d, %d, %d.'], mesh_file, n_faces, numel(neighbor), numel(len));
end
if strict_contract && (numel(cell_x) ~= n_cells || numel(cell_y) ~= n_cells || ...
        any(~isfinite(cell_x)) || any(~isfinite(cell_y)) ...
    )
    error('HydroPol2D:ContractCellCenters', ...
        'Mesh %s: mesh2d_face_x/y must contain one finite center per cell.', mesh_file);
end
if strict_contract
    face_info = ncinfo(mesh_file, 'mesh2d_face_nodes');
    attr_names = string({face_info.Attributes.Name});
    start_position = find(attr_names == "start_index", 1);
    if isempty(start_position) || double(face_info.Attributes(start_position).Value) ~= 0
        error('HydroPol2D:ContractConnectivityBase', ...
            'Mesh %s: mesh2d_face_nodes start_index must be 0.', mesh_file);
    end
end
if min(owner) < 0 || max(owner) >= n_cells
    error('HydroPol2D:ContractOwnerRange', ...
        ['Mesh %s: edge_owner out of range [0, %d] -- the FILE numbers cells from ' ...
         'ZERO. Got %g..%g.'], mesh_file, n_cells-1, min(owner), max(owner));
end

internal = neighbor >= 0;
report.n_internal_faces = sum(internal);
report.n_boundary_faces = sum(~internal);
bad = neighbor(~internal);
if ~isempty(bad) && any(bad ~= FILE_BOUNDARY_NEIGHBOR)
    error('HydroPol2D:ContractBoundarySentinel', ...
        ['Mesh %s: a boundary face must carry edge_neighbor = %d in the file; ' ...
         'found %g. This reader adds one, so any other sentinel becomes a valid ' ...
         'cell index here -- an off-by-one inside the boundary block still lands ' ...
         'on a boundary face and passes every other check, which moved the Pune ' ...
         'outlet 528 m silently.'], mesh_file, FILE_BOUNDARY_NEIGHBOR, bad(find(bad ~= FILE_BOUNDARY_NEIGHBOR,1)));
end
if any(internal) && max(neighbor(internal)) >= n_cells
    error('HydroPol2D:ContractNeighborRange', ...
        'Mesh %s: edge_neighbor out of range [0, %d].', mesh_file, n_cells-1);
end

normal_norm = hypot(normal_x, normal_y);
if any(~isfinite(normal_norm)) || any(abs(normal_norm - 1) > 1e-6)
    error('HydroPol2D:ContractNormalLength', ...
        'Mesh %s: edge normals must be finite unit vectors.', mesh_file);
end
if strict_contract && any(internal)
    direction_x = cell_x(neighbor(internal)+1) - cell_x(owner(internal)+1);
    direction_y = cell_y(neighbor(internal)+1) - cell_y(owner(internal)+1);
    outward_dot = normal_x(internal) .* direction_x + normal_y(internal) .* direction_y;
    bad_direction = find(~isfinite(outward_dot) | outward_dot <= 0, 1);
    if ~isempty(bad_direction)
        internal_indices = find(internal);
        error('HydroPol2D:ContractNormalOrientation', ...
            ['Mesh %s: edge normal %d does not point outward from edge_owner ' ...
             'toward its neighbor.'], mesh_file, internal_indices(bad_direction)-1);
    end
end

zero_faces = find(len <= 0);
if ~isempty(zero_faces)
    error('HydroPol2D:ContractZeroFace', ...
        ['Mesh %s: %d face(s) have edge_length_m <= 0 (first at index %d). A ' ...
         'zero-length face divides by zero in the momentum step and produced NaN ' ...
         'at t = 3 s; two of them were why no Pune run would start.'], ...
        mesh_file, numel(zero_faces), zero_faces(1)-1);
end
if any(distance(internal) <= 0)
    error('HydroPol2D:ContractZeroDistance', ...
        ['Mesh %s: edge_center_distance_m must be positive on internal faces; it ' ...
         'is the denominator of the water-surface gradient.'], mesh_file);
end
if ~all(isfinite(bed))
    error('HydroPol2D:ContractNonFiniteBed', ...
        'Mesh %s: cell_bed_elevation_m has %d non-finite entries.', ...
        mesh_file, sum(~isfinite(bed)));
end
if any(area <= 0)
    error('HydroPol2D:ContractCellArea', 'Mesh %s: cell_area_m2 must be positive.', mesh_file);
end
if strict_contract
    cfl_width = double(ncread(mesh_file,'cell_cfl_width_m')); cfl_width = cfl_width(:);
    if numel(cfl_width) ~= n_cells || any(~isfinite(cfl_width) | cfl_width <= 0)
        error('HydroPol2D:ContractCFLWidth', ...
            'Mesh %s: cell_cfl_width_m must contain one positive finite value per cell.', mesh_file);
    end
end

perimeter = accumarray(owner+1, len, [n_cells 1]);
perimeter = perimeter + accumarray(neighbor(internal)+1, len(internal), [n_cells 1]);
ap = inf(n_cells,1);
positive = perimeter > 0;
ap(positive) = area(positive) ./ perimeter(positive);
report.minimum_area_over_perimeter_m = min(ap);

% ---- overlap -------------------------------------------------------------
if ~isempty(overlap_file)
    if exist(overlap_file,'file') ~= 2
        error('HydroPol2D:ContractMissingOverlap', 'Overlap file not found: %s', overlap_file);
    end
    require_variables(overlap_file, {'overlap_mesh_index','overlap_raster_index', ...
        'overlap_area_m2','mesh_area_m2','raster_area_m2','x_edges','y_edges'}, 'Overlap');
    if strict_contract
        require_attribute(overlap_file, 'mesh_contract_version', MESH_CONTRACT_VERSION, 'Overlap');
        require_attribute(overlap_file, 'schema_version', MESH_PRODUCT_SCHEMA, 'Overlap');
        require_attribute(overlap_file, 'product_type', 'conservative_raster_overlap', 'Overlap');
        require_attribute(overlap_file, 'file_index_base', 0, 'Overlap');
        require_attribute(overlap_file, 'raster_index_order', OVERLAP_RASTER_ORDER, 'Overlap');
        require_attribute(overlap_file, 'coordinate_units', COORDINATE_UNITS, 'Overlap');
        require_matching_crs(mesh_file, overlap_file, 'Overlap');
        require_units(overlap_file, 'overlap_area_m2', 'm2', 'Overlap');
        require_units(overlap_file, 'mesh_area_m2', 'm2', 'Overlap');
        require_units(overlap_file, 'raster_area_m2', 'm2', 'Overlap');
        require_units(overlap_file, 'x_edges', 'm', 'Overlap');
        require_units(overlap_file, 'y_edges', 'm', 'Overlap');
    end
    mi = double(ncread(overlap_file,'overlap_mesh_index')); mi = mi(:);
    ri = double(ncread(overlap_file,'overlap_raster_index')); ri = ri(:);
    oa = double(ncread(overlap_file,'overlap_area_m2')); oa = oa(:);
    ra = double(ncread(overlap_file,'raster_area_m2')); ra = ra(:);
    if numel(mi) ~= numel(ri) || numel(mi) ~= numel(oa)
        error('HydroPol2D:ContractOverlapArrays', ...
            'Overlap %s: the three overlap_* arrays must have equal length.', overlap_file);
    end
    if ~isempty(mi) && (min(mi) < 0 || max(mi) >= n_cells)
        error('HydroPol2D:ContractOverlapMeshIndex', ...
            'Overlap %s: overlap_mesh_index out of range [0, %d].', overlap_file, n_cells-1);
    end
    if ~isempty(ri) && (min(ri) < 0 || max(ri) >= numel(ra))
        error('HydroPol2D:ContractOverlapRasterIndex', ...
            'Overlap %s: overlap_raster_index out of range [0, %d].', overlap_file, numel(ra)-1);
    end
    covered = accumarray(mi+1, oa, [n_cells 1]);
    err = abs(covered - area);
    [worst_err, worst] = max(err);
    if worst_err > options.area_tolerance_m2
        error('HydroPol2D:ContractOverlapCoverage', ...
            ['Overlap %s: cell %d is covered by %.6f m2 of raster against its own ' ...
             'area %.6f m2 (error %.3e > %g). Every cell must be fully covered or ' ...
             'the raster-to-mesh aggregation loses water.'], overlap_file, worst-1, ...
            covered(worst), area(worst), worst_err, options.area_tolerance_m2);
    end
    report.overlap_pairs = numel(oa);
    report.maximum_cell_area_error_m2 = worst_err;
    % overlap_raster_index is SOUTH-UP: index 0 is the BOTTOM row of a north-up
    % raster. Reading it north-up silently samples the wrong terrain.
    report.overlap_note = ['overlap_raster_index is SOUTH-UP -- flip a north-up ' ...
        'raster before indexing it.'];
end

% ---- sub-grid tables -----------------------------------------------------
if ~isempty(subgrid_file)
    if exist(subgrid_file,'file') ~= 2
        error('HydroPol2D:ContractMissingSubgrid', ...
            'Sub-grid table file not found: %s', subgrid_file);
    end
    require_variables(subgrid_file, {'cell_datum_m','cell_zeta_m','cell_volume_m3', ...
        'cell_wet_area_m2','cell_point_count','cell_plan_area_m2','face_datum_m', ...
        'face_zeta_m','face_flow_area_m2','face_perimeter_m','face_conveyance', ...
        'face_point_count','face_length_m'}, 'Sub-grid table');
    if strict_contract
        require_attribute(subgrid_file, 'mesh_contract_version', MESH_CONTRACT_VERSION, 'Sub-grid table');
        require_attribute(subgrid_file, 'schema_version', MESH_PRODUCT_SCHEMA, 'Sub-grid table');
        require_attribute(subgrid_file, 'product_type', 'hydraulic_subgrid_tables', 'Sub-grid table');
        require_attribute(subgrid_file, 'file_index_base', 0, 'Sub-grid table');
        require_attribute(subgrid_file, 'coordinate_units', COORDINATE_UNITS, 'Sub-grid table');
        require_attribute(subgrid_file, 'conveyance_convention', SUBGRID_CONVEYANCE_CONVENTION, 'Sub-grid table');
        require_matching_crs(mesh_file, subgrid_file, 'Sub-grid table');
        require_units(subgrid_file, 'cell_datum_m', 'm', 'Sub-grid table');
        require_units(subgrid_file, 'cell_zeta_m', 'm', 'Sub-grid table');
        require_units(subgrid_file, 'cell_volume_m3', 'm3', 'Sub-grid table');
        require_units(subgrid_file, 'cell_wet_area_m2', 'm2', 'Sub-grid table');
        require_units(subgrid_file, 'face_datum_m', 'm', 'Sub-grid table');
        require_units(subgrid_file, 'face_zeta_m', 'm', 'Sub-grid table');
        require_units(subgrid_file, 'face_flow_area_m2', 'm2', 'Sub-grid table');
        require_units(subgrid_file, 'face_perimeter_m', 'm', 'Sub-grid table');
        require_units(subgrid_file, 'face_length_m', 'm', 'Sub-grid table');
    end
    cell_datum = double(ncread(subgrid_file,'cell_datum_m')); cell_datum = cell_datum(:);
    cell_zeta  = double(ncread(subgrid_file,'cell_zeta_m'));
    cell_vol   = double(ncread(subgrid_file,'cell_volume_m3'));
    cell_cnt   = double(ncread(subgrid_file,'cell_point_count')); cell_cnt = cell_cnt(:);
    plan       = double(ncread(subgrid_file,'cell_plan_area_m2')); plan = plan(:);
    face_datum = double(ncread(subgrid_file,'face_datum_m')); face_datum = face_datum(:);
    face_cnt   = double(ncread(subgrid_file,'face_point_count')); face_cnt = face_cnt(:);

    if numel(cell_datum) ~= n_cells
        error('HydroPol2D:ContractSubgridCells', ...
            ['Sub-grid tables %s cover %d cells but the mesh %s has %d. The tables ' ...
             'must be rebuilt whenever the mesh changes.'], ...
            subgrid_file, numel(cell_datum), mesh_file, n_cells);
    end
    if numel(face_datum) ~= n_faces
        error('HydroPol2D:ContractSubgridFaces', ...
            'Sub-grid tables %s cover %d faces but the mesh has %d.', ...
            subgrid_file, numel(face_datum), n_faces);
    end
    % ncread returns (cell, point); MATLAB reads NetCDF dimensions reversed only
    % when the file declares them in C order, so orient by the cell count.
    if size(cell_zeta,1) ~= n_cells, cell_zeta = cell_zeta.'; cell_vol = cell_vol.'; end
    width = size(cell_zeta,2);
    if any(cell_cnt < 1) || any(cell_cnt > width)
        error('HydroPol2D:ContractCellPointCount', ...
            'Sub-grid tables %s: cell_point_count must be in [1, %d].', subgrid_file, width);
    end
    if any(face_cnt < 1)
        error('HydroPol2D:ContractFacePointCount', ...
            'Sub-grid tables %s: face_point_count must be >= 1.', subgrid_file);
    end
    valid = (1:width) <= cell_cnt;
    for pair = {{'cell_zeta_m', cell_zeta}, {'cell_volume_m3', cell_vol}}
        name = pair{1}{1}; table = pair{1}{2};
        step = diff(table, 1, 2);
        bad_rows = any(valid(:,2:end) & (step < -1e-9), 2);
        if any(bad_rows)
            error('HydroPol2D:ContractMonotonic', ...
                ['Sub-grid tables %s: %s decreases inside the valid prefix of row ' ...
                 '%d. Rows must be non-decreasing across the FULL width -- the ' ...
                 'padding repeats the last value so one vectorised lookup can ' ...
                 'serve every row.'], subgrid_file, name, find(bad_rows,1));
        end
    end
    if any(plan <= 0)
        error('HydroPol2D:ContractPlanArea', ...
            ['Sub-grid tables %s: cell_plan_area_m2 must be positive; it is the ' ...
             'curve slope above the top breakpoint.'], subgrid_file);
    end

    % A face lies ON the cell boundary, so it belongs to both cells' footprints and
    % cannot dip below either. Violating this put 3,920 of 7,969 faces below their
    % own cells by up to 1.0 m, and the sub-grid arm then failed to converge to the
    % flat-prism answer even at 20 m cells where the two are identical.
    safe = max(neighbor, 0);
    lower = cell_datum(owner+1);
    lower(internal) = min(cell_datum(owner(internal)+1), cell_datum(safe(internal)+1));
    gap = lower - face_datum;
    offenders = find(gap > 1e-6);
    if ~isempty(offenders)
        [worst_gap, k] = max(gap(offenders));
        error('HydroPol2D:ContractFaceDatum', ...
            ['Sub-grid tables %s: %d face(s) have a datum BELOW both of their ' ...
             'cells, worst face %d by %.4f m. cell_datum must be <= the datum of ' ...
             'each of its faces, or a nearly dry cell hands its faces metres of ' ...
             'phantom head.'], subgrid_file, numel(offenders), offenders(k)-1, worst_gap);
    end

    info = ncinfo(subgrid_file);
    names = string({info.Variables.Name});
    if any(names == "cell_is_subgrid")
        flag = double(ncread(subgrid_file,'cell_is_subgrid'));
        if any(flag ~= 0 & flag ~= 1)
            error('HydroPol2D:ContractSubgridFlag', ...
                ['Sub-grid tables %s: cell_is_subgrid must be 0 or 1. It gates the ' ...
                 'FACE closure too, so any other value selects a closure for faces ' ...
                 'it was never meant to.'], subgrid_file);
        end
        report.subgrid_cells = sum(flag == 1);
    else
        report.subgrid_cells = n_cells;
        report.subgrid_note = 'no cell_is_subgrid variable; every cell uses the tables.';
    end
    report.subgrid_cell_points = width;
end
end

function require_variables(file, names, kind)
info = ncinfo(file);
present = string({info.Variables.Name});
missing = setdiff(string(names), present);
if ~isempty(missing)
    error('HydroPol2D:ContractMissingVariables', ...
        ['%s %s does not satisfy the mesh contract: missing %s. See ' ...
         'UNSTRUCTURED_MESH_CONTRACT.md.'], kind, file, strjoin(cellstr(missing), ', '));
end
end

function present = has_attribute(file, name)
info = ncinfo(file);
present = ~isempty(info.Attributes) && ...
    any(string({info.Attributes.Name}) == string(name));
end

function require_attribute(file, name, expected, kind)
if ~has_attribute(file, name)
    error('HydroPol2D:ContractMissingAttribute', ...
        '%s %s does not satisfy the mesh contract: missing attribute %s.', kind, file, name);
end
actual = ncreadatt(file, '/', name);
if isnumeric(expected)
    matches = isnumeric(actual) && isscalar(actual) && double(actual) == double(expected);
else
    matches = strcmp(char(string(actual)), char(string(expected)));
end
if ~matches
    error('HydroPol2D:ContractAttribute', ...
        '%s %s: attribute %s is %s; expected %s.', kind, file, name, ...
        char(string(actual)), char(string(expected)));
end
end

function require_nonempty_attribute(file, name, kind)
if ~has_attribute(file, name) || strlength(strtrim(string(ncreadatt(file, '/', name)))) == 0
    error('HydroPol2D:ContractMissingAttribute', ...
        '%s %s does not satisfy the mesh contract: missing nonempty attribute %s.', ...
        kind, file, name);
end
end

function require_units(file, variable, expected, kind)
info = ncinfo(file, variable);
names = string({info.Attributes.Name});
position = find(names == "units", 1);
if isempty(position)
    error('HydroPol2D:ContractUnits', ...
        '%s %s: variable %s has no units attribute.', kind, file, variable);
end
actual = string(info.Attributes(position).Value);
if actual ~= string(expected)
    error('HydroPol2D:ContractUnits', ...
        '%s %s: variable %s uses units %s; expected %s.', ...
        kind, file, variable, actual, string(expected));
end
end

function require_matching_crs(mesh_file, child_file, kind)
require_nonempty_attribute(child_file, 'crs_wkt', kind);
mesh_crs = string(ncreadatt(mesh_file, '/', 'crs_wkt'));
child_crs = string(ncreadatt(child_file, '/', 'crs_wkt'));
if mesh_crs ~= child_crs
    error('HydroPol2D:ContractCRS', ...
        '%s %s uses a CRS different from mesh %s.', kind, child_file, mesh_file);
end
end
