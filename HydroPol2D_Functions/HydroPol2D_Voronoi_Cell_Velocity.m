function velocity = HydroPol2D_Voronoi_Cell_Velocity(mesh, depth_m, edge_q_m2_s, dry_tolerance_m)
%HYDROPOL2D_VORONOI_CELL_VELOCITY Maximum adjacent-face speed per polygon.

if nargin < 4, dry_tolerance_m = 1e-6; end
assert(size(depth_m,1) == mesh.n_cells && size(edge_q_m2_s,1) == mesh.n_edges);
assert(size(depth_m,2) == size(edge_q_m2_s,2));

velocity = zeros(size(depth_m));
internal = mesh.edge_neighbor > 0;
for k = 1:size(depth_m,2)
    stage = mesh.surface_bed(:) + depth_m(:,k);
    face_depth = max(stage(mesh.edge_owner) - mesh.surface_bed(mesh.edge_owner), 0);
    face_depth(internal) = max( ...
        max(stage(mesh.edge_owner(internal)), stage(mesh.edge_neighbor(internal))) - ...
        max(mesh.surface_bed(mesh.edge_owner(internal)), mesh.surface_bed(mesh.edge_neighbor(internal))), 0);
    face_speed = zeros(mesh.n_edges,1);
    wet = face_depth > dry_tolerance_m;
    face_speed(wet) = abs(edge_q_m2_s(wet,k)) ./ face_depth(wet);
    owner_speed = accumarray(mesh.edge_owner, face_speed, [mesh.n_cells 1], @max, 0);
    neighbor_speed = accumarray(mesh.edge_neighbor(internal), face_speed(internal), ...
        [mesh.n_cells 1], @max, 0);
    velocity(:,k) = max(owner_speed, neighbor_speed);
end
end
