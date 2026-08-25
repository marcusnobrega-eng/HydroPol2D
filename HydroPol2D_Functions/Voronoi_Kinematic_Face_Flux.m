function [edge_q, face_depth] = Voronoi_Kinematic_Face_Flux(mesh, surface_volume, roughness, options)
%VORONOI_KINEMATIC_FACE_FLUX Upwind Manning fluxes on an arbitrary mesh.
%
% The kinematic approximation follows bed slope only; it cannot convey
% backwater or flow reversal caused by the water-surface gradient.

arguments
    mesh struct
    surface_volume (:,1) double
    roughness double
    options.gravity (1,1) double = 9.81
    options.dry_tolerance_m (1,1) double = 1e-6
    options.critical_flow (1,1) logical = false
end

n = mesh.n_cells;
if isscalar(roughness), roughness = repmat(roughness, n, 1); end
assert(numel(surface_volume) == n && numel(roughness) == n);
owner = mesh.edge_owner(:); neighbor = mesh.edge_neighbor(:);
internal = neighbor > 0;
o = owner(internal); d = neighbor(internal);
bed_slope = (mesh.surface_bed(o) - mesh.surface_bed(d)) ./ mesh.edge_distance(internal);
donor = o;
donor(bed_slope < 0) = d(bed_slope < 0);
depth = max(surface_volume(:) ./ mesh.surface_area(:), 0);
h = depth(donor);
nface = roughness(donor);
nface(~isfinite(nface) | nface <= 0) = 1e-6;
q = sign(bed_slope) .* h.^(5/3) ./ nface .* sqrt(abs(bed_slope));
q(h <= options.dry_tolerance_m | abs(bed_slope) <= 0) = 0;
if options.critical_flow
    qcrit = h .* sqrt(options.gravity .* h);
    q = sign(q) .* min(abs(q), qcrit);
end

edge_q = zeros(mesh.n_edges, 1);
face_depth = zeros(mesh.n_edges, 1);
edge_q(internal) = q;
face_depth(internal) = h;
end
