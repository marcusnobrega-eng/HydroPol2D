function results = HydroPol2D_Voronoi_Run(mesh_file, config, forcing)
%HYDROPOL2D_VORONOI_RUN Separate conservative Voronoi/local-inertial runner.
%
% forcing.surface_source_m_s may be a scalar, n-cell vector, or function
% handle source = f(time_s, state, mesh). Positive values add water.

arguments
    mesh_file (1,:) char
    config struct
    forcing struct = struct()
end
config = defaults(config);
if ~strcmpi(config.routing_solver, 'local_inertial')
    error('HydroPol2D:VoronoiSolverUnavailable', ...
        'Voronoi hybrid rivers currently require routing_solver=local_inertial.');
end
[mesh, preflight] = HydroPol2D_Voronoi_Preflight(mesh_file, config);
n = mesh.n_cells;
if isfield(forcing, 'overlap_file')
    forcing.mapping = HydroPol2D_Read_Overlap(forcing.overlap_file);
    assert(size(forcing.mapping.raster_to_mesh, 1) == n, ...
        'HydroPol2D:InvalidVoronoiForcing', 'Overlap file does not match mesh.');
end
surface_volume = initial_vector(config.initial_surface_depth_m, n) .* mesh.surface_area(:);
edge_q = signed_vector(config.initial_surface_discharge_per_width_m2_s, mesh.n_edges);
channel_volume = zeros(mesh.channel.n_nodes, 1);
link_q = signed_vector(config.initial_channel_discharge_m3_s, mesh.channel.n_links);
transition_q = zeros(mesh.channel.n_transitions, 1);
channel_boundary_count = 0;
if isfield(forcing, 'channel_boundary') && isfield(forcing.channel_boundary, 'node_id')
    channel_boundary_count = numel(forcing.channel_boundary.node_id);
end
channel_boundary_q = zeros(channel_boundary_count, 1);
if mesh.channel.n_nodes > 0
    channel_volume = initial_vector(config.initial_channel_depth_m, mesh.channel.n_nodes) .* mesh.channel.plan_area(:);
    [surface_volume, channel_volume] = Voronoi_Equilibrate_Channel_Storage(mesh, surface_volume, channel_volume);
end

roughness = initial_vector(config.surface_roughness, n);
roughness(isfinite(mesh.cell_roughness)) = mesh.cell_roughness(isfinite(mesh.cell_roughness));
groundwater = initialize_groundwater(config, mesh);
t = 0; next_output = 0; step = 0; initial_mass = sum(surface_volume) + sum(channel_volume);
if groundwater.enabled, initial_mass = initial_mass + groundwater_mass(groundwater, mesh); end
times = []; depths = {}; channel_depths = {}; groundwater_heads = {}; diagnostics = struct([]);
while t < config.duration_s - eps(config.duration_s)
    dt = stable_timestep(mesh, surface_volume, channel_volume, groundwater, config, t);
    dt = min([dt, config.max_dt_s, config.duration_s - t]);
    if dt < config.min_dt_s && config.duration_s - t > config.min_dt_s
        error('HydroPol2D:VoronoiTimestepFloor', 'Required timestep %.6g s is below minimum %.6g s.', dt, config.min_dt_s);
    end

    state = struct('surface_volume_m3', surface_volume, 'channel_volume_m3', channel_volume, 'groundwater_head_m', groundwater.head);
    source = source_at(forcing, t, state, mesh, n);
    source_volume = source .* mesh.cell_area(:) .* dt;
    if any(source_volume < 0)
        source_volume = max(source_volume, -surface_volume);
    end
    surface_volume = surface_volume + source_volume;
    groundwater_source_volume = zeros(n,1);
    if groundwater.enabled
        recharge = groundwater_source_at(forcing, t, state, mesh, n);
        groundwater_source_volume = recharge .* mesh.cell_area(:) .* dt;
        groundwater_source_volume = max(groundwater_source_volume, -groundwater.storage);
        groundwater.storage = max(groundwater.storage + groundwater_source_volume, 0);
        groundwater.head = groundwater.bottom + groundwater.storage ./ (groundwater.specific_yield .* mesh.cell_area(:));
        [groundwater.head, groundwater_diag] = Voronoi_Boussinesq_Step(mesh, groundwater.head, groundwater.bottom, ...
            groundwater.hydraulic_conductivity, groundwater.specific_yield, dt);
        groundwater.storage = groundwater.specific_yield .* max(groundwater.head - groundwater.bottom, 0) .* mesh.cell_area(:);
        excess = groundwater.specific_yield .* max(groundwater.head - mesh.cell_bed(:), 0) .* mesh.cell_area(:);
        groundwater.storage = groundwater.storage - excess;
        groundwater.head = min(groundwater.head, mesh.cell_bed(:));
        surface_volume = surface_volume + excess;
    else
        groundwater_diag = struct('max_flux_m3_s', 0, 'mass_change_m3', 0);
    end
    surface_channel_exchange = 0;
    if mesh.channel.n_nodes > 0
        [surface_volume, channel_volume, ~, exchange] = Voronoi_Equilibrate_Channel_Storage(mesh, surface_volume, channel_volume);
        surface_channel_exchange = surface_channel_exchange + exchange;
    end
    surface_flux_evaluation_volume = surface_volume;
    [surface_volume, edge_q, surface_diag] = Voronoi_Local_Inertial_Step( ...
        mesh, surface_volume, edge_q, roughness, dt, gravity=config.gravity, ...
        dry_tolerance_m=config.dry_tolerance_m, critical_flow=config.critical_flow);
    surface_boundary = boundary_at(forcing, 'surface_boundary', t, state, mesh);
    [surface_volume, edge_q, surface_boundary_diag] = Voronoi_Surface_Boundary_Step( ...
        mesh, surface_volume, edge_q, roughness, surface_boundary, dt, gravity=config.gravity, ...
        dry_tolerance_m=config.dry_tolerance_m, evaluation_volume=surface_flux_evaluation_volume);
    if mesh.channel.n_nodes > 0
        [surface_volume, channel_volume, ~, exchange] = Voronoi_Equilibrate_Channel_Storage(mesh, surface_volume, channel_volume);
        surface_channel_exchange = surface_channel_exchange + exchange;
        channel_flux_evaluation_volume = channel_volume;
        [channel_volume, link_q, channel_diag] = Voronoi_Neal_Channel_Step( ...
            mesh.channel, channel_volume, link_q, dt, gravity=config.gravity, ...
            dry_tolerance_m=config.dry_tolerance_m, critical_flow=config.critical_flow);
        channel_boundary = boundary_at(forcing, 'channel_boundary', t, state, mesh);
        [channel_volume, channel_boundary_q, channel_boundary_diag] = Voronoi_Channel_Boundary_Step( ...
            mesh.channel, channel_volume, channel_boundary_q, channel_boundary, dt, ...
            gravity=config.gravity, evaluation_volume=channel_flux_evaluation_volume);
        [surface_volume, channel_volume, transition_q, transition_diag] = Voronoi_Channel_Transition_Step( ...
            mesh, surface_volume, channel_volume, transition_q, dt, gravity=config.gravity, ...
            dry_tolerance_m=config.dry_tolerance_m, critical_flow=config.critical_flow);
        [surface_volume, channel_volume, ~, exchange] = Voronoi_Equilibrate_Channel_Storage(mesh, surface_volume, channel_volume);
        surface_channel_exchange = surface_channel_exchange + exchange;
    else
        channel_diag = struct('max_depth_m', 0, 'max_velocity_m_s', 0, 'mass_change_m3', 0);
        transition_diag = struct('max_velocity_m_s', 0, 'mass_change_m3', 0);
        channel_boundary_diag = struct('net_inflow_volume_m3', 0, 'max_discharge_m3_s', 0);
    end
    t = t + dt; step = step + 1;
    mass = sum(surface_volume) + sum(channel_volume);
    if groundwater.enabled, mass = mass + sum(groundwater.storage); end
    boundary_volume = surface_boundary_diag.net_inflow_volume_m3 + channel_boundary_diag.net_inflow_volume_m3;
    expected = initial_mass + sum(source_volume) + sum(groundwater_source_volume) + boundary_volume;
    diagnostics(step).time_s = t;
    diagnostics(step).dt_s = dt;
    diagnostics(step).mass_m3 = mass;
    diagnostics(step).step_mass_residual_m3 = mass - expected;
    diagnostics(step).max_surface_depth_m = surface_diag.max_depth_m;
    diagnostics(step).max_surface_velocity_m_s = surface_diag.max_velocity_m_s;
    diagnostics(step).max_channel_depth_m = channel_diag.max_depth_m;
    diagnostics(step).max_channel_velocity_m_s = channel_diag.max_velocity_m_s;
    diagnostics(step).max_transition_velocity_m_s = transition_diag.max_velocity_m_s;
    diagnostics(step).max_groundwater_flux_m3_s = groundwater_diag.max_flux_m3_s;
    diagnostics(step).surface_channel_exchange_m3 = surface_channel_exchange;
    diagnostics(step).boundary_net_inflow_volume_m3 = boundary_volume;
    initial_mass = expected;
    if t + eps(t) >= next_output
        times(end+1,1) = t; %#ok<AGROW>
        depths{end+1,1} = surface_volume ./ mesh.surface_area(:); %#ok<AGROW>
        if mesh.channel.n_nodes > 0
            channel_depths{end+1,1} = channel_volume ./ mesh.channel.plan_area(:); %#ok<AGROW>
        end
        if groundwater.enabled, groundwater_heads{end+1,1} = groundwater.head; end %#ok<AGROW>
        next_output = next_output + config.output_interval_s;
    end
end
results.mesh_file = mesh_file;
results.preflight = preflight;
results.time_s = times;
results.surface_depth_m = cat(2, depths{:});
results.channel_depth_m = cat_or_empty(channel_depths, mesh.channel.n_nodes);
results.groundwater_head_m = cat_or_empty(groundwater_heads, n);
results.final_surface_volume_m3 = surface_volume;
results.final_channel_volume_m3 = channel_volume;
results.final_groundwater_head_m = groundwater.head;
results.edge_discharge_per_width_m2_s = edge_q;
results.channel_discharge_m3_s = link_q;
results.channel_transition_discharge_m3_s = transition_q;
results.channel_boundary_discharge_m3_s = channel_boundary_q;
results.diagnostics = diagnostics;
if ~isempty(config.output_netcdf)
    HydroPol2D_Write_Voronoi_Output(config.output_netcdf, mesh, results, config.overwrite_output);
end
end

function config = defaults(config)
values = struct('routing_solver','local_inertial','duration_s',3600,'initial_surface_depth_m',0, ...
    'initial_channel_depth_m',0,'surface_roughness',0.05,'min_dt_s',0.01,'max_dt_s',30, ...
    'output_interval_s',300,'courant',0.6,'gravity',9.81, ...
    'dry_tolerance_m',1e-6,'critical_flow',false,'groundwater_enabled',false, ...
    'initial_surface_discharge_per_width_m2_s',0,'initial_channel_discharge_m3_s',0, ...
    'initial_groundwater_head_m',0,'aquifer_bottom_m',-10,'hydraulic_conductivity_m_s',1e-5, ...
    'specific_yield',0.2,'compute_backend','cpu','output_netcdf','', ...
    'overwrite_output',false,'forcing_interval_s',inf);
names = fieldnames(values);
for k = 1:numel(names)
    if ~isfield(config, names{k}), config.(names{k}) = values.(names{k}); end
end
end

function value = initial_vector(value, count)
if isscalar(value), value = repmat(double(value), count, 1); else, value = double(value(:)); end
assert(numel(value) == count && all(isfinite(value) & value >= 0));
end

function source = source_at(forcing, time_s, state, mesh, n)
if isfield(forcing, 'surface_source_m_s')
    source = forcing.surface_source_m_s;
elseif isfield(forcing, 'raster_source_m_s') && isfield(forcing, 'mapping')
    source = forcing.raster_source_m_s;
    if isa(source, 'function_handle'), source = source(time_s, state, mesh); end
    source = forcing.mapping.raster_to_mesh * double(source(:));
    return
else
    source = zeros(n,1); return
end
if isa(source, 'function_handle'), source = source(time_s, state, mesh); end
if isscalar(source), source = repmat(double(source), n, 1); else, source = double(source(:)); end
assert(numel(source) == n && all(isfinite(source)), 'HydroPol2D:InvalidVoronoiForcing', 'Invalid surface source.');
end

function source = groundwater_source_at(forcing, time_s, state, mesh, n)
if ~isfield(forcing, 'groundwater_recharge_m_s'), source = zeros(n,1); return; end
source = forcing.groundwater_recharge_m_s;
if isa(source, 'function_handle'), source = source(time_s, state, mesh); end
if isscalar(source), source = repmat(double(source), n, 1); else, source = double(source(:)); end
assert(numel(source) == n && all(isfinite(source)));
end

function boundary = boundary_at(forcing, name, time_s, state, mesh)
if ~isfield(forcing, name), boundary = struct(); return; end
boundary = forcing.(name);
if isa(boundary, 'function_handle'), boundary = boundary(time_s, state, mesh); end
end

function groundwater = initialize_groundwater(config, mesh)
groundwater.enabled = logical(config.groundwater_enabled);
if ~groundwater.enabled
    groundwater.head = zeros(0,1); groundwater.storage = zeros(0,1); return
end
n = mesh.n_cells;
groundwater.head = signed_vector(config.initial_groundwater_head_m, n);
groundwater.bottom = signed_vector(config.aquifer_bottom_m, n);
groundwater.hydraulic_conductivity = initial_vector(config.hydraulic_conductivity_m_s, n);
groundwater.specific_yield = initial_vector(config.specific_yield, n);
assert(all(groundwater.head >= groundwater.bottom));
groundwater.storage = groundwater.specific_yield .* (groundwater.head - groundwater.bottom) .* mesh.cell_area(:);
end

function value = signed_vector(value, count)
if isscalar(value), value = repmat(double(value), count, 1); else, value = double(value(:)); end
assert(numel(value) == count && all(isfinite(value)));
end

function mass = groundwater_mass(groundwater, mesh)
mass = sum(groundwater.specific_yield .* max(groundwater.head - groundwater.bottom, 0) .* mesh.cell_area(:));
end

function dt = stable_timestep(mesh, surface_volume, channel_volume, groundwater, config, time_s)
depth = max(surface_volume ./ mesh.surface_area(:), 0);
perimeter = accumarray(mesh.edge_owner(:), mesh.edge_length(:), [mesh.n_cells 1], @sum, 0);
internal = mesh.edge_neighbor > 0;
perimeter = perimeter + accumarray(mesh.edge_neighbor(internal), mesh.edge_length(internal), [mesh.n_cells 1], @sum, 0);
length_scale = 2 .* mesh.cell_area(:) ./ max(perimeter, eps);
celerity = sqrt(config.gravity .* max(depth, config.dry_tolerance_m));
dt = config.courant * min(length_scale ./ celerity);
if mesh.channel.n_links > 0
    channel_depth = max(channel_volume ./ mesh.channel.plan_area(:), config.dry_tolerance_m);
    link_depth = max(channel_depth(mesh.channel.link_up), channel_depth(mesh.channel.link_down));
    dt = min(dt, config.courant * min(mesh.channel.link_length ./ sqrt(config.gravity .* link_depth)));
end
if mesh.channel.n_transitions > 0
    transition_depth = max(channel_volume(mesh.channel.transition_node) ./ mesh.channel.plan_area(mesh.channel.transition_node), config.dry_tolerance_m);
    dt = min(dt, config.courant * min(mesh.channel.transition_length ./ sqrt(config.gravity .* transition_depth)));
end
if groundwater.enabled
    thickness = max(groundwater.head - groundwater.bottom, 0);
    transmissivity = groundwater.hydraulic_conductivity .* thickness;
    internal = mesh.edge_neighbor > 0;
    o = mesh.edge_owner(internal); d = mesh.edge_neighbor(internal);
    conductance = 2 .* transmissivity(o) .* transmissivity(d) ./ max(transmissivity(o) + transmissivity(d), eps) ...
        .* mesh.edge_length(internal) ./ mesh.edge_distance(internal);
    sum_conductance = accumarray(o,conductance,[mesh.n_cells 1],@sum,0) + ...
        accumarray(d,conductance,[mesh.n_cells 1],@sum,0);
    capacity = groundwater.specific_yield .* mesh.cell_area(:);
    active = sum_conductance > 0;
    if any(active)
        dt = min(dt, config.courant * min(capacity(active) ./ sum_conductance(active)));
    end
end
if isfinite(config.forcing_interval_s)
    remaining = config.forcing_interval_s - mod(time_s, config.forcing_interval_s);
    if remaining <= 10 * eps(max(time_s,1)), remaining = config.forcing_interval_s; end
    dt = min(dt, remaining);
end
end

function out = cat_or_empty(values, rows)
if isempty(values), out = zeros(rows,0); else, out = cat(2, values{:}); end
end
