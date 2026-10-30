function [boundary, inflow_volume_m3] = hp2d_voronoi_split_inflow(mesh, boundary, dt_s)
%HP2D_VORONOI_SPLIT_INFLOW Take inflow edges off the routing boundary as cell volume.
%   HydroPol2D-Python adds an inflow hydrograph to its host cells as volume
%   before the vertical step (voronoi.py), so the inflow can infiltrate in
%   the step it arrives. Each "inflow" edge's discharge [m3/s] times dt goes
%   to the edge's cell (inflow_volume_m3, n_cells x 1), and the edge becomes
%   a wall for routing, so the water enters once.

inflow_volume_m3 = zeros(mesh.n_cells, 1);
if ~isstruct(boundary) || ~isfield(boundary, 'edge_id') || isempty(boundary.edge_id)
    return
end
edge_id = double(boundary.edge_id(:));
types = string(boundary.type(:));
if isscalar(types), types = repmat(types, numel(edge_id), 1); end
values = zeros(numel(edge_id), 1);
if isfield(boundary, 'value')
    values = double(boundary.value(:));
    if isscalar(values), values = repmat(values, numel(edge_id), 1); end
end
inflow = types == "inflow";
if ~any(inflow), return; end
inflow_volume_m3 = accumarray(mesh.edge_owner(edge_id(inflow)), values(inflow) .* dt_s, ...
    [mesh.n_cells 1], @sum, 0);
types(inflow) = "wall";
values(inflow) = 0;
boundary.type = types;
boundary.value = values;
end
