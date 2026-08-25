function [edge_q, face_depth, edge_conductance] = Voronoi_Diffusive_Face_Flux(mesh, surface_volume, roughness, options)
%VORONOI_DIFFUSIVE_FACE_FLUX Explicit Manning diffusive-wave face fluxes.
%
% edge_conductance is a regularized dQ/d(eta_owner-eta_neighbor) used
% only for the explicit parabolic timestep restriction.

arguments
    mesh struct
    surface_volume (:,1) double
    roughness double
    options.dry_tolerance_m (1,1) double = 1e-6
    options.slope_regularization (1,1) double {mustBePositive} = 1e-4
end

n = mesh.n_cells;
if isscalar(roughness), roughness = repmat(roughness, n, 1); end
assert(numel(surface_volume) == n && numel(roughness) == n);
owner = mesh.edge_owner(:); neighbor = mesh.edge_neighbor(:);
internal = neighbor > 0;
o = owner(internal); d = neighbor(internal);
depth = max(surface_volume(:) ./ mesh.surface_area(:), 0);
eta = mesh.surface_bed(:) + depth;
bed = max(mesh.surface_bed(o), mesh.surface_bed(d));
h = max(max(eta(o), eta(d)) - bed, 0);
slope = (eta(o) - eta(d)) ./ mesh.edge_distance(internal);
nface = 0.5 .* (roughness(o) + roughness(d));
nface(~isfinite(nface) | nface <= 0) = 1e-6;
coefficient = h.^(5/3) ./ nface;
q = sign(slope) .* coefficient .* sqrt(abs(slope));
q(h <= options.dry_tolerance_m | abs(slope) <= 0) = 0;

% The first term is the regularized slope derivative. The second safely
% accounts for the changing hydrostatic face depth when either stage rises.
slope_abs = max(abs(slope), options.slope_regularization);
dq_deta = 0.5 .* coefficient ./ sqrt(slope_abs) ./ mesh.edge_distance(internal) + ...
    (5/3) .* h.^(2/3) ./ nface .* sqrt(abs(slope));
conductance = mesh.edge_length(internal) .* dq_deta;
conductance(~isfinite(conductance) | h <= options.dry_tolerance_m) = 0;

edge_q = zeros(mesh.n_edges, 1);
face_depth = zeros(mesh.n_edges, 1);
edge_conductance = zeros(mesh.n_edges, 1);
edge_q(internal) = q;
face_depth(internal) = h;
edge_conductance(internal) = conductance;
end
