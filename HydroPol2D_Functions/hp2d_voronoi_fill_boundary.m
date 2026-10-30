function boundary = hp2d_voronoi_fill_boundary(mesh, boundary)
%HP2D_VORONOI_FILL_BOUNDARY Give every unlisted boundary edge one condition.
%   With boundary.fill_unassigned_type (for example "normal_flow") and
%   boundary.fill_unassigned_value (its slope), every exterior edge that the
%   boundary does not list, and whose cell is not an inflow cell, gets that
%   condition. This is normal_flow_all_unassigned_boundary_edges in
%   HydroPol2D-Python (voronoi.py). Without the fields the boundary is left
%   as it is, and unlisted edges stay walls.

if ~isstruct(boundary) || ~isfield(boundary, 'fill_unassigned_type') || ...
        isempty(boundary.fill_unassigned_type)
    return
end
listed = zeros(0, 1); types = strings(0, 1); values = zeros(0, 1);
if isfield(boundary, 'edge_id') && ~isempty(boundary.edge_id)
    listed = double(boundary.edge_id(:));
    types = string(boundary.type(:));
    if isscalar(types), types = repmat(types, numel(listed), 1); end
    values = double(boundary.value(:));
    if isscalar(values), values = repmat(values, numel(listed), 1); end
end
inflow_cells = mesh.edge_owner(listed(types == "inflow"));
exterior = find(mesh.edge_neighbor(:) == 0);
fill = exterior(~ismember(exterior, listed) & ~ismember(mesh.edge_owner(exterior), inflow_cells));
fill_value = 0;
if isfield(boundary, 'fill_unassigned_value'), fill_value = double(boundary.fill_unassigned_value); end
boundary.edge_id = [listed; fill];
boundary.type = [types; repmat(string(boundary.fill_unassigned_type), numel(fill), 1)];
boundary.value = [values; repmat(fill_value, numel(fill), 1)];
boundary = rmfield(boundary, intersect({'fill_unassigned_type','fill_unassigned_value'}, fieldnames(boundary)));
end
