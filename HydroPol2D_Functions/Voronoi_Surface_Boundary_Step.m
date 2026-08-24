function [surface_volume, edge_q, diagnostics] = Voronoi_Surface_Boundary_Step(mesh, surface_volume, edge_q, roughness, boundary, dt, options)
%VORONOI_SURFACE_BOUNDARY_STEP Apply signed fluxes on selected boundary edges.
%
% Positive discharge is out of the modeled domain. Supported types are
% wall, inflow (value is m3/s into domain), stage (value is WSE m),
% normal_flow (value is slope), and critical_flow.

arguments
    mesh struct
    surface_volume (:,1) double
    edge_q (:,1) double
    roughness (:,1) double
    boundary struct
    dt (1,1) double {mustBePositive}
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.evaluation_volume (:,1) double = surface_volume
end
if isempty(boundary) || ~isfield(boundary, 'edge_id') || isempty(boundary.edge_id)
    diagnostics = struct('net_inflow_volume_m3', 0, 'max_discharge_m3_s', 0); return
end
edge_id = double(boundary.edge_id(:));
assert(all(edge_id >= 1 & edge_id <= mesh.n_edges & mesh.edge_neighbor(edge_id) == 0));
types = string(boundary.type(:));
if isscalar(types), types = repmat(types, numel(edge_id), 1); end
values = zeros(numel(edge_id),1);
if isfield(boundary, 'value'), values = expand(boundary.value, numel(edge_id)); end
owner = mesh.edge_owner(edge_id);
assert(numel(options.evaluation_volume) == mesh.n_cells);
depth = max(options.evaluation_volume(owner) ./ mesh.surface_area(owner), 0);
wse = mesh.surface_bed(owner) + depth;
width = mesh.edge_length(edge_id);
qout = zeros(numel(edge_id),1); % m3/s
for k = 1:numel(edge_id)
    switch types(k)
        case "wall"
            qout(k) = 0;
        case "inflow"
            qout(k) = -values(k);
        case "stage"
            hflow = max(max(wse(k), values(k)) - mesh.surface_bed(owner(k)), 0);
            slope = (values(k) - wse(k)) / mesh.edge_distance(edge_id(k));
            qold = edge_q(edge_id(k));
            denominator = 1 + options.gravity * dt * roughness(owner(k))^2 * abs(qold) / max(hflow^(7/3), eps);
            qout(k) = width(k) * (qold - options.gravity * hflow * dt * slope) / denominator;
        case "normal_flow"
            qout(k) = width(k) / roughness(owner(k)) * depth(k)^(5/3) * sqrt(max(values(k), 0));
        case "critical_flow"
            qout(k) = width(k) * depth(k) * sqrt(options.gravity * depth(k));
        otherwise
            error('HydroPol2D:UnknownBoundaryType', 'Unknown Voronoi boundary type %s.', types(k));
    end
end
outgoing = zeros(mesh.n_cells, 1);
if any(qout > 0)
    outgoing = accumarray(owner(qout > 0), qout(qout > 0) .* dt, [mesh.n_cells 1], @sum, 0);
end
scale = min(1, surface_volume ./ max(outgoing, eps));
qout(qout > 0) = qout(qout > 0) .* scale(owner(qout > 0));
delta = accumarray(owner, -qout .* dt, [mesh.n_cells 1], @sum, 0);
surface_volume = max(surface_volume + delta, 0);
edge_q(edge_id) = qout ./ width;
diagnostics.net_inflow_volume_m3 = -sum(qout) * dt;
diagnostics.max_discharge_m3_s = max(abs(qout), [], 'omitnan');
end

function value = expand(value, count)
if isscalar(value), value = repmat(double(value), count, 1); else, value = double(value(:)); end
assert(numel(value) == count && all(isfinite(value)));
end
