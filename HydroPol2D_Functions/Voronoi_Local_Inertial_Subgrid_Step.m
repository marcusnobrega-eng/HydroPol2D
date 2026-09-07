function [surface_volume, edge_q, diagnostics] = Voronoi_Local_Inertial_Subgrid_Step( ...
    mesh, surface_volume, edge_q, roughness, tables, dt, options)
%VORONOI_LOCAL_INERTIAL_SUBGRID_STEP Local-inertial update on HEC-RAS sub-grid tables.
%
%   Two substitutions relative to the flat-prism step:
%
%     * stage inverts the cell elevation-volume curve instead of dividing volume
%       by plan area, so water held in a channel narrower than the cell reports
%       the channel water surface;
%     * the friction denominator uses |Q|*A/K^2 in place of n^2|q|/hflow^(7/3).
%       These are the same number when the face is a flat rectangle of uniform
%       roughness -- with K = A^(5/3)/(n P^(2/3)), n^2 P^(4/3)/A^(7/3) reduces to
%       A/K^2 -- and the K form is correct when the face is neither flat nor
%       uniform.
%
%   Working in TOTAL discharge is what makes this possible: A and K are face
%   integrals with no per-unit-width equivalent once the face stops being a
%   rectangle.
%
%   SELECTIVE SUB-GRID.  A face uses the tables only where BOTH its cells do; the
%   rest keep the flat-prism closure, evaluated here so the two share one state.
%   Switching the cell closure alone leaves a flat-prism cell referenced to its
%   area-mean bed while the face table measures from the face minimum -- about
%   1.9 m lower on a 90 m mesh -- and a nearly dry cell then presents its faces
%   metres of phantom head.
%
%   Verified against the flat-prism step: on perfectly flat terrain, where level
%   pool is exact and every face profile is flat, the two closures agree to
%   round-off (discharge ratio 1.0000 at three timesteps).

arguments
    mesh struct
    surface_volume (:,1) double
    edge_q (:,1) double
    roughness double
    tables struct
    dt (1,1) double {mustBePositive}
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
g = options.gravity;

% ---- stage from the cell tables ------------------------------------------
[wse, ~] = hp2d_voronoi_subgrid_cell_stage(tables, surface_volume);

% ---- selective closure ---------------------------------------------------
face_subgrid = tables.cell_is_subgrid(o) & tables.cell_is_subgrid(d);
length_i = mesh.edge_length(internal);
internal_idx = find(internal);
sub_idx = internal_idx(face_subgrid);

% Look the face tables up ONLY on faces that actually use them. Every other face
% takes the flat-prism formula below regardless, so evaluating all of them was
% wasted: on the Pune mesh 674 of 28,055 faces are sub-grid faces, so this is 42x
% less work in the hottest routine of the step. A sub-table is built so the row ids
% still line up -- the lookup is indexed by face row, so a gathered subset must
% carry its own datum, zeta and count.
area = zeros(numel(o),1); perim = zeros(numel(o),1); conveyance = zeros(numel(o),1);
if ~isempty(sub_idx)
    sub = struct( ...
        'face_datum_m', tables.face_datum_m(sub_idx), ...
        'face_zeta_m', tables.face_zeta_m(sub_idx,:), ...
        'face_flow_area_m2', tables.face_flow_area_m2(sub_idx,:), ...
        'face_perimeter_m', tables.face_perimeter_m(sub_idx,:), ...
        'face_conveyance', tables.face_conveyance(sub_idx,:), ...
        'face_point_count', tables.face_point_count(sub_idx));
    sub_stage = max(wse(o(face_subgrid)), wse(d(face_subgrid)));
    [a_s, p_s, c_s] = hp2d_voronoi_subgrid_face_state(sub, sub_stage);
    area(face_subgrid) = a_s;
    perim(face_subgrid) = p_s;
    conveyance(face_subgrid) = c_s;
end
hflat = max(max(wse(o), wse(d)) - max(mesh.surface_bed(o), mesh.surface_bed(d)), 0);
nflat = 0.5 .* (roughness(o) + roughness(d));
nflat(~isfinite(nflat) | nflat <= 0) = 1e-6;
aflat = hflat .* length_i;
kflat = aflat ./ nflat .* max(hflat, 0).^(2/3);
area(~face_subgrid) = aflat(~face_subgrid);
perim(~face_subgrid) = length_i(~face_subgrid);
conveyance(~face_subgrid) = kflat(~face_subgrid);

% hydraulic depth A/P: the depth that sets the wave speed through this face
hflow = zeros(size(area));
positive = perim > 0;
hflow(positive) = area(positive) ./ perim(positive);
slope = (wse(d) - wse(o)) ./ mesh.edge_distance(internal);
Qold = edge_q(internal) .* length_i;

% A dry cell's water surface sits at its own invert, which on sloping ground is
% ABOVE the low point of the face it shares downslope, so the face tables report
% area on dry ground.  Require water to exist here rather than clipping the face
% profile up to a controlling invert: clipping flattens the face into a weir and
% inflated conveyance 57x in the Python implementation.
has_water = surface_volume(:) > 0;
wet = (has_water(o) | has_water(d)) & (hflow > options.dry_tolerance_m) & ...
    (area > 0) & (conveyance > 0) & isfinite(slope);

Q = zeros(size(Qold));
den = 1 + g .* dt .* abs(Qold(wet)) .* area(wet) ./ max(conveyance(wet).^2, eps);
Q(wet) = (Qold(wet) - g .* area(wet) .* dt .* slope(wet)) ./ den;
if options.critical_flow
    qcrit = area .* sqrt(g .* max(hflow, 0));
    Q = sign(Q) .* min(abs(Q), qcrit);
end

out_cell = o; out_cell(Q < 0) = d(Q < 0);
requested = accumarray(out_cell, abs(Q) .* dt, [n 1], @sum, 0);
scale = min(1, surface_volume(:) ./ max(requested, eps));
Q = Q .* scale(out_cell);

delta = accumarray(o, -Q .* dt, [n 1], @sum, 0) + accumarray(d, Q .* dt, [n 1], @sum, 0);
surface_volume = surface_volume(:) + delta;
surface_volume(abs(surface_volume) < 10 * eps(max(1, max(surface_volume)))) = 0;
if any(surface_volume < -1e-10)
    error('HydroPol2D:NegativeSurfaceVolume', 'Sub-grid surface draining limiter failed.');
end
surface_volume = max(surface_volume, 0);
edge_q(internal) = Q ./ length_i;
edge_q(~internal) = 0;   % boundary fluxes are supplied separately by the runner

velocity = zeros(size(Q));
velocity(area > 0) = abs(Q(area > 0)) ./ area(area > 0);
diagnostics.face_flow_depth_m = zeros(mesh.n_edges,1);
diagnostics.face_flow_depth_m(internal) = hflow;
% Reported depth stays volume/plan-area so it is the same quantity the flat-prism
% run reports; the sub-grid water surface is recovered from the stored volumes and
% the cell table when mapping depths.
diagnostics.max_depth_m = max(max(surface_volume ./ mesh.surface_area(:), 0));
diagnostics.max_velocity_m_s = max(velocity, [], 'omitnan');
diagnostics.internal_flux_volume_m3 = sum(abs(Q)) * dt;
diagnostics.mass_change_m3 = sum(delta);
diagnostics.subgrid_faces = nnz(face_subgrid);
end
