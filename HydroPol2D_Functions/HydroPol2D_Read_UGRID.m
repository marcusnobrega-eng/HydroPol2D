function mesh = HydroPol2D_Read_UGRID(mesh_file)
%HYDROPOL2D_READ_UGRID Read the HydroPolMesh UGRID exchange contract.

arguments
    mesh_file (1,:) char
end
if exist(mesh_file, 'file') ~= 2
    error('HydroPol2D:MissingMesh', 'Mesh file not found: %s', mesh_file);
end
info = ncinfo(mesh_file);
variables = string({info.Variables.Name});
if isempty(info.Attributes)
    attributes = strings(0,1);
else
    attributes = string({info.Attributes.Name});
end
mesh.crs_wkt = '';
if any(attributes == "crs_wkt"), mesh.crs_wkt=char(ncreadatt(mesh_file,'/','crs_wkt')); end

mesh.cell_area = double(ncread(mesh_file, 'cell_area_m2'));
mesh.cell_bed = double(ncread(mesh_file, 'cell_bed_elevation_m'));
mesh.cell_x = double(ncread(mesh_file, 'mesh2d_face_x'));
mesh.cell_y = double(ncread(mesh_file, 'mesh2d_face_y'));
mesh.cell_target_width = double(ncread(mesh_file, 'cell_target_width_m'));
mesh.cell_refinement_source = double(ncread(mesh_file, 'cell_refinement_source'));
if any(variables == "cell_hydraulic_roughness")
    mesh.cell_roughness = double(ncread(mesh_file, 'cell_hydraulic_roughness'));
else
    mesh.cell_roughness = nan(size(mesh.cell_area));
end
mesh.edge_owner = double(ncread(mesh_file, 'edge_owner')) + 1;
raw_neighbor = double(ncread(mesh_file, 'edge_neighbor'));
mesh.edge_neighbor = raw_neighbor + 1; % NetCDF -1 becomes MATLAB boundary index 0.
mesh.edge_length = double(ncread(mesh_file, 'edge_length_m'));
mesh.edge_distance = double(ncread(mesh_file, 'edge_center_distance_m'));
mesh.edge_midpoint_x = double(ncread(mesh_file, 'edge_midpoint_x'));
mesh.edge_midpoint_y = double(ncread(mesh_file, 'edge_midpoint_y'));
mesh.edge_normal_x = double(ncread(mesh_file, 'edge_normal_x'));
mesh.edge_normal_y = double(ncread(mesh_file, 'edge_normal_y'));
mesh.edge_boundary_type = double(ncread(mesh_file, 'edge_boundary_type'));
mesh.n_cells = numel(mesh.cell_area);
mesh.n_edges = numel(mesh.edge_owner);

% The hydrographically relevant length scale is the minimum face-normal
% centre separation.  It is exported by HydroPolMesh as cell_cfl_width_m.
% Do not use 2*A/P here: that compactness measure is deliberately smaller
% for a valid elongated, flow-aligned floodplain/river polygon.
if any(variables == "cell_cfl_width_m")
    mesh.cell_cfl_width = double(ncread(mesh_file, 'cell_cfl_width_m'));
else
    edge_cfl_width = mesh.edge_distance(:);
    boundary = mesh.edge_neighbor == 0;
    edge_cfl_width(boundary) = 2 .* edge_cfl_width(boundary);
    owner_width = accumarray(mesh.edge_owner(:), edge_cfl_width, [mesh.n_cells 1], @min, inf);
    internal = ~boundary;
    neighbor_width = accumarray(mesh.edge_neighbor(internal), edge_cfl_width(internal), ...
        [mesh.n_cells 1], @min, inf);
    mesh.cell_cfl_width = min(owner_width, neighbor_width);
end

assert(all(mesh.cell_area > 0 & isfinite(mesh.cell_area)), 'HydroPol2D:InvalidMesh', 'Cell areas must be positive and finite.');
assert(all(mesh.edge_owner >= 1 & mesh.edge_owner <= mesh.n_cells), 'HydroPol2D:InvalidMesh', 'Invalid edge owner.');
assert(all(mesh.edge_neighbor >= 0 & mesh.edge_neighbor <= mesh.n_cells), 'HydroPol2D:InvalidMesh', 'Invalid edge neighbor.');
assert(all(mesh.edge_length > 0 & mesh.edge_distance > 0), 'HydroPol2D:InvalidMesh', 'Edge geometry must be positive.');
assert(numel(mesh.cell_cfl_width) == mesh.n_cells && ...
    all(mesh.cell_cfl_width > 0 & isfinite(mesh.cell_cfl_width)), ...
    'HydroPol2D:InvalidMesh', 'Cell CFL widths must be positive and finite.');

mesh.channel = struct('n_nodes', 0, 'n_links', 0, 'n_transitions', 0);
if any(variables == "channel_node_host_cell")
    mesh.channel.x = double(ncread(mesh_file, 'channel_node_x'));
    mesh.channel.y = double(ncread(mesh_file, 'channel_node_y'));
    mesh.channel.host_cell = double(ncread(mesh_file, 'channel_node_host_cell')) + 1;
    mesh.channel.bed = double(ncread(mesh_file, 'channel_node_bed_elevation_m'));
    mesh.channel.bank = double(ncread(mesh_file, 'channel_node_bank_elevation_m'));
    mesh.channel.plan_area = double(ncread(mesh_file, 'channel_node_plan_area_m2'));
    mesh.channel.n_nodes = numel(mesh.channel.host_cell);
end
if any(variables == "channel_link_upstream_node")
    mesh.channel.link_up = double(ncread(mesh_file, 'channel_link_upstream_node')) + 1;
    mesh.channel.link_down = double(ncread(mesh_file, 'channel_link_downstream_node')) + 1;
    mesh.channel.link_length = double(ncread(mesh_file, 'channel_link_length_m'));
    mesh.channel.link_width = double(ncread(mesh_file, 'channel_link_width_m'));
    mesh.channel.link_depth = double(ncread(mesh_file, 'channel_link_bankfull_depth_m'));
    mesh.channel.link_roughness = double(ncread(mesh_file, 'channel_link_roughness'));
    mesh.channel.link_reach_id = double(ncread(mesh_file, 'channel_link_reach_id'));
    mesh.channel.n_links = numel(mesh.channel.link_up);
end
if any(variables == "channel_transition_node")
    mesh.channel.transition_node = double(ncread(mesh_file, 'channel_transition_node')) + 1;
    mesh.channel.transition_cell = double(ncread(mesh_file, 'channel_transition_resolved_cell')) + 1;
    mesh.channel.transition_positive_node_to_cell = logical(ncread(mesh_file, 'channel_transition_positive_node_to_cell'));
    mesh.channel.transition_length = double(ncread(mesh_file, 'channel_transition_length_m'));
    mesh.channel.transition_width = double(ncread(mesh_file, 'channel_transition_width_m'));
    mesh.channel.transition_depth = double(ncread(mesh_file, 'channel_transition_bankfull_depth_m'));
    mesh.channel.transition_roughness = double(ncread(mesh_file, 'channel_transition_roughness'));
    mesh.channel.n_transitions = numel(mesh.channel.transition_node);
end

if mesh.channel.n_nodes > 0
    assert(numel(unique(mesh.channel.host_cell)) == mesh.channel.n_nodes, ...
        'HydroPol2D:InvalidChannelGraph', 'Each host polygon may own at most one aggregated channel node.');
    assert(all(mesh.channel.plan_area > 0), 'HydroPol2D:InvalidChannelGraph', 'Channel storage area must be positive.');
end
if mesh.channel.n_links > 0
    assert(all(mesh.channel.link_up >= 1 & mesh.channel.link_up <= mesh.channel.n_nodes));
    assert(all(mesh.channel.link_down >= 1 & mesh.channel.link_down <= mesh.channel.n_nodes));
    assert(all(mesh.channel.link_length > 0 & mesh.channel.link_width > 0));
end
if mesh.channel.n_transitions > 0
    assert(all(mesh.channel.transition_node >= 1 & mesh.channel.transition_node <= mesh.channel.n_nodes));
    assert(all(mesh.channel.transition_cell >= 1 & mesh.channel.transition_cell <= mesh.n_cells));
    assert(all(mesh.channel.transition_length > 0 & mesh.channel.transition_width > 0));
end

% Hybrid cells retain separate channel and floodplain plan areas. The
% surface hydraulic bed at a subgrid host is the mapped channel bank, not
% the channel bed; this keeps one stage without lowering the whole polygon.
mesh.surface_area = mesh.cell_area(:);
mesh.surface_bed = mesh.cell_bed(:);
if mesh.channel.n_nodes > 0
    host = mesh.channel.host_cell(:);
    assert(all(mesh.channel.plan_area(:) < mesh.cell_area(host)), ...
        'HydroPol2D:InvalidChannelGraph', ...
        'Subgrid channel plan area must be smaller than its host polygon.');
    mesh.surface_area(host) = mesh.cell_area(host) - mesh.channel.plan_area(:);
    mesh.surface_bed(host) = mesh.channel.bank(:);
end
end
