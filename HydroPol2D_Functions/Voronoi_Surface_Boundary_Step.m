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
    options.dry_tolerance_m (1,1) double = 1e-3
    options.evaluation_volume (:,1) double = surface_volume
    options.routing_solver (1,:) char = 'local_inertial'
    % Edge discharge per width before this step's surface step. The surface step sets
    % boundary edges to zero, so the stage boundary must read q_old from here (Python
    % 12f494b, 2026-10-10). Empty keeps edge_q.
    options.previous_edge_q (:,1) double = []
    % Hydraulic tables (HydroPol2D_Read_Subgrid_Tables). On a boundary edge whose owner
    % uses them, the water surface comes from the cell volume curve and the face table
    % gives A, P and K at the boundary water surface (Python 12f494b): stage
    % Q = (Q_old - g A dt S)/(1 + g dt |Q_old| A/K^2), normal Q = K sqrt(S), critical
    % Q = A sqrt(g A/P). On flat tables these equal the flat-prism formulas.
    options.subgrid_tables struct = struct()
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
            if strcmpi(options.routing_solver, 'kinematic')
                error('HydroPol2D:KinematicStageBoundaryUnavailable', ...
                    'A kinematic-wave boundary cannot prescribe stage; use inflow, normal_flow, critical_flow, or wall.');
            end
            hflow = max(max(wse(k), values(k)) - mesh.surface_bed(owner(k)), 0);
            slope = (values(k) - wse(k)) / mesh.edge_distance(edge_id(k));
            if strcmpi(options.routing_solver, 'diffusive')
                qout(k) = width(k) / roughness(owner(k)) * hflow^(5/3) * sign(-slope) * sqrt(abs(slope));
            else
                if isempty(options.previous_edge_q)
                    qold = edge_q(edge_id(k));
                else
                    qold = options.previous_edge_q(edge_id(k));
                end
                denominator = 1 + options.gravity * dt * roughness(owner(k))^2 * abs(qold) / max(hflow^(7/3), eps);
                qout(k) = width(k) * (qold - options.gravity * hflow * dt * slope) / denominator;
            end
        case "normal_flow"
            qout(k) = width(k) / roughness(owner(k)) * depth(k)^(5/3) * sqrt(max(values(k), 0));
        case "critical_flow"
            qout(k) = width(k) * depth(k) * sqrt(options.gravity * depth(k));
        otherwise
            error('HydroPol2D:UnknownBoundaryType', 'Unknown Voronoi boundary type %s.', types(k));
    end
end
tables = options.subgrid_tables;
if isfield(tables, 'face_datum_m')
    if isfield(tables, 'cell_is_subgrid') && ~isempty(tables.cell_is_subgrid)
        owner_uses_tables = logical(tables.cell_is_subgrid(owner));
    else
        owner_uses_tables = true(size(owner));
    end
    rows = find(owner_uses_tables & ismember(types, ["stage", "normal_flow", "critical_flow"]));
    if ~isempty(rows) && ~strcmpi(options.routing_solver, 'kinematic')
        faces = edge_id(rows);
        sub_tables = struct('face_datum_m', tables.face_datum_m(faces), ...
            'face_zeta_m', tables.face_zeta_m(faces, :), ...
            'face_flow_area_m2', tables.face_flow_area_m2(faces, :), ...
            'face_perimeter_m', tables.face_perimeter_m(faces, :), ...
            'face_conveyance', tables.face_conveyance(faces, :), ...
            'face_point_count', tables.face_point_count(faces));
        [cell_wse, ~] = hp2d_voronoi_subgrid_cell_stage(tables, options.evaluation_volume);
        table_wse = cell_wse(owner(rows));
        is_stage = types(rows) == "stage";
        face_stage = table_wse;
        face_stage(is_stage) = max(table_wse(is_stage), values(rows(is_stage)));
        [area, perimeter, conveyance] = hp2d_voronoi_subgrid_face_state(sub_tables, face_stage);
        hydraulic_depth = zeros(size(area));
        wet = perimeter > 0;
        hydraulic_depth(wet) = area(wet)./perimeter(wet);
        slope = (values(rows) - table_wse)./mesh.edge_distance(faces);
        if isempty(options.previous_edge_q)
            q_unit_old = edge_q(faces);
        else
            q_unit_old = options.previous_edge_q(faces);
        end
        discharge_old = q_unit_old.*width(rows);
        stage_q = (discharge_old - options.gravity.*area.*dt.*slope)./ ...
            (1 + options.gravity.*dt.*abs(discharge_old).*area./max(conveyance.^2, eps));
        stage_q(~(area > 0 & conveyance > 0)) = 0;
        normal_q = conveyance.*sqrt(max(values(rows), 0));
        critical_q = area.*sqrt(options.gravity.*hydraulic_depth);
        table_q = critical_q;
        table_q(types(rows) == "normal_flow") = normal_q(types(rows) == "normal_flow");
        table_q(is_stage) = stage_q(is_stage);
        qout(rows) = table_q;
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
