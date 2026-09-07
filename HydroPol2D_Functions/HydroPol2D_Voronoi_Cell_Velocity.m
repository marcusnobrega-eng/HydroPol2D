function velocity = HydroPol2D_Voronoi_Cell_Velocity(mesh, depth_m, edge_q_m2_s, dry_tolerance_m, face_depth_m)
%HYDROPOL2D_VORONOI_CELL_VELOCITY Maximum adjacent-face speed per polygon.
%
%   FACE_DEPTH_M (optional, n_edges x n_times) is the face flow depth that the
%   momentum step PAIRED with EDGE_Q_M2_S.  Pass it whenever it is available.
%
%   Recomputing the depth from the post-step state instead is wrong, and
%   measurably so: critical_flow caps q at h*sqrt(g*h) during the step, so
%   v <= sqrt(g*h) for the h used then, but the volume update that follows changes
%   h.  Measured on a 90 m Pune mesh, the recomputed pair put 31.56% of wet faces
%   above critical, up to 182x, and reported a maximum velocity of 34.996 m/s
%   where the deepest cell's critical velocity is 1.69 m/s.  The limiter had held;
%   only the diagnostic was inconsistent.  With the paired depth the same run
%   reports 1.020 m/s and no cell above 5 m/s.

if nargin < 4 || isempty(dry_tolerance_m), dry_tolerance_m = 1e-6; end
if nargin < 5, face_depth_m = []; end
assert(size(depth_m,1) == mesh.n_cells && size(edge_q_m2_s,1) == mesh.n_edges);
assert(size(depth_m,2) == size(edge_q_m2_s,2));
paired = ~isempty(face_depth_m);
if paired
    assert(size(face_depth_m,1) == mesh.n_edges && ...
        size(face_depth_m,2) == size(edge_q_m2_s,2), ...
        'HydroPol2D:PairedFaceDepthShape', ...
        'face_depth_m must be n_edges by n_times to match edge_q_m2_s.');
end

velocity = zeros(size(depth_m));
internal = mesh.edge_neighbor > 0;
for k = 1:size(depth_m,2)
    if paired
        face_depth = face_depth_m(:,k);
    else
        stage = mesh.surface_bed(:) + depth_m(:,k);
        face_depth = max(stage(mesh.edge_owner) - mesh.surface_bed(mesh.edge_owner), 0);
        face_depth(internal) = max( ...
            max(stage(mesh.edge_owner(internal)), stage(mesh.edge_neighbor(internal))) - ...
            max(mesh.surface_bed(mesh.edge_owner(internal)), mesh.surface_bed(mesh.edge_neighbor(internal))), 0);
    end
    face_speed = zeros(mesh.n_edges,1);
    wet = face_depth > dry_tolerance_m;
    face_speed(wet) = abs(edge_q_m2_s(wet,k)) ./ face_depth(wet);
    owner_speed = accumarray(mesh.edge_owner, face_speed, [mesh.n_cells 1], @max, 0);
    neighbor_speed = accumarray(mesh.edge_neighbor(internal), face_speed(internal), ...
        [mesh.n_cells 1], @max, 0);
    velocity(:,k) = max(owner_speed, neighbor_speed);
end
end
