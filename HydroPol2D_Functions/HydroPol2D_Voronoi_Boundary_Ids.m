function edge_id = HydroPol2D_Voronoi_Boundary_Ids(boundary, n_edges)
%HYDROPOL2D_VORONOI_BOUNDARY_IDS Normalise boundary face ids to 1-based MATLAB indices.
%
%   The UGRID file numbers faces from zero, and so does HydroPol2D-Python.  MATLAB
%   numbers them from one.  The same wall therefore carries two labels: the Pune
%   outlet is face 831522 in the file and 831523 here.  Copying a list between the
%   two without saying which base it is written in moves the outlet to a different
%   wall hundreds of metres away, and it does so silently, because an off-by-one
%   inside the boundary block still lands on a valid boundary face and passes every
%   other check.
%
%   boundary.edge_id_base is 0 or 1.  It defaults to 1 so existing MATLAB cases are
%   unaffected; a case exported from Python must set it to 0.

arguments
    boundary struct
    n_edges (1,1) double
end

if ~isfield(boundary, 'edge_id') || isempty(boundary.edge_id)
    edge_id = zeros(0, 1);
    return
end
edge_id = double(boundary.edge_id(:));

base = 1;
if isfield(boundary, 'edge_id_base') && ~isempty(boundary.edge_id_base)
    base = double(boundary.edge_id_base);
end
assert(isscalar(base) && (base == 0 || base == 1), ...
    'HydroPol2D:InvalidBoundaryEdgeBase', 'boundary.edge_id_base must be 0 or 1.');

if base == 0
    assert(all(edge_id >= 0), 'HydroPol2D:InvalidBoundaryEdgeBase', ...
        'boundary.edge_id_base is 0 but an id is negative.');
    edge_id = edge_id + 1;
else
    assert(all(edge_id >= 1), 'HydroPol2D:InvalidBoundaryEdgeBase', ...
        ['boundary.edge_id_base is 1 but an id is below 1; if these ids came ' ...
         'from a Python case set boundary.edge_id_base = 0.']);
end
assert(all(edge_id <= n_edges), 'HydroPol2D:InvalidBoundaryEdge', ...
    'A boundary edge id is beyond the number of faces in this mesh.');
end
