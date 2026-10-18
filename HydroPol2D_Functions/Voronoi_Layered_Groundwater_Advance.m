function [surface_volume, channel_volume, groundwater, state, diagnostics] = Voronoi_Layered_Groundwater_Advance( ...
    mesh, surface_volume, channel_volume, groundwater, state, duration_s, config)
%VORONOI_LAYERED_GROUNDWATER_ADVANCE Pending exchange, lateral flow, river exchange, fixed heads.
%   Layered counterpart of Voronoi_Groundwater_Advance (HydroPol2D-Python
%   voronoi_layered.layered_advance_groundwater). The pending vertical exchange and
%   every lateral/river volume go through hp2d_solve_shared_water_table, so the HRUs of a
%   cell share one table and the MMF rules hold. Lateral flow: face T = mean of the cell
%   area-weighted layered T, conductance = T*edge_length/edge_distance, explicit
%   sub-steps limited by courant*Sy*A/sum(conductance), outflow limited by drainable water.
diagnostics = struct('max_flux_m3_s', 0, 'substep_count', 0, 'river_exchange_m3', 0, ...
    'seepage_volume_m3', 0, 'fixed_head_exchange_m3', 0, 'mass_change_m3', 0, ...
    'rejected_steps', 0, 'suggested_dt_s', inf);
if ~groundwater.enabled || duration_s <= 0
    return
end
if ~isfield(config, 'groundwater_courant'), config.groundwater_courant = 0.25; end
L = state.layered;
n = mesh.n_cells;
host = L.host;
area = L.area_m2;
cell_area = mesh.cell_area(:);
seepage = zeros(n, 1);
[groundwater, L, exfiltration, unmet] = apply_cell_volume(groundwater, L, groundwater.pending_exchange_m3);
tolerance = 1e-9*max(max(abs(groundwater.pending_exchange_m3), [], 'all'), 1);
if max(unmet) > tolerance
    error('HydroPol2D:GroundwaterMass', ...
        'Vadose exchange requested more groundwater than is stored (deficit %.6g m3).', max(unmet));
end
seepage = seepage + exfiltration;
groundwater.pending_exchange_m3(:) = 0;

internal = mesh.edge_neighbor(:) > 0;
owner = mesh.edge_owner(internal); owner = owner(:);
neighbor = mesh.edge_neighbor(internal); neighbor = neighbor(:);
lateral = ~(isfield(config, 'groundwater_lateral_flow') && ~config.groundwater_lateral_flow);
river = isfield(config, 'river_groundwater_exchange') && logical(config.river_groundwater_exchange);
fixed_mask = false(n, 1);
if isfield(config, 'groundwater_fixed_head_mask') && ~isempty(config.groundwater_fixed_head_mask)
    fixed_mask = logical(config.groundwater_fixed_head_mask(:));
end
remaining = duration_s;
maximum_flux = 0; substeps = 0; fixed_exchange = 0; river_total = 0;
exchange_by_cell = zeros(n, 1);
while (lateral || river) && remaining > 10*eps(max(duration_s, 1))
    head = groundwater.head(:);
    drainable = accumarray(host, hp2d_drainable_storage_mm(L.column, head(host))/1000.*area, [n 1]);
    delta = zeros(n, 1);
    out_volume = zeros(n, 1);
    stable = inf;
    if lateral
        T = accumarray(host, hp2d_layered_transmissivity(L.column, head(host)).*area, [n 1])./cell_area;
        face_t = 0.5*(T(owner) + T(neighbor));
        conductance = face_t.*mesh.edge_length(internal)./mesh.edge_distance(internal);
        conductance = conductance(:);
        total = accumarray(owner, conductance, [n 1]) + accumarray(neighbor, conductance, [n 1]);
        sy = max(accumarray(host, hp2d_effective_specific_yield(L.column, head(host)).*area, [n 1])./cell_area, 1e-3);
        active = total > 0;
        if any(active)
            stable = config.groundwater_courant*min(sy(active).*cell_area(active)./total(active));
        end
    end
    dt = min(remaining, stable);
    if lateral
        q = -conductance.*(head(neighbor) - head(owner));
        donor = owner;
        donor(q < 0) = neighbor(q < 0);
        requested = accumarray(donor, abs(q)*dt, [n 1]);
        scale = min(1, drainable./max(requested, eps));
        q = q.*scale(donor);
        out_volume = accumarray(donor, abs(q)*dt, [n 1]);
        delta = accumarray(owner, -q*dt, [n 1]) + accumarray(neighbor, q*dt, [n 1]);
        maximum_flux = max(maximum_flux, max(abs(q), [], 'all'));
    end
    if river
        [surface_volume, channel_volume, by_cell] = layered_river_exchange(mesh, ...
            surface_volume, channel_volume, head, dt, drainable - out_volume, config);
        delta = delta - by_cell;
        river_total = river_total + sum(by_cell);
        exchange_by_cell = exchange_by_cell + by_cell;
        maximum_flux = max(maximum_flux, max(abs(by_cell/dt), [], 'all'));
    end
    [groundwater, L, exfiltration] = apply_cell_volume(groundwater, L, delta);
    seepage = seepage + exfiltration;
    if any(fixed_mask)
        [groundwater, L, added] = hold_fixed_heads(groundwater, L, fixed_mask, config.groundwater_fixed_head_m(:));
        fixed_exchange = fixed_exchange + added;
    end
    remaining = remaining - dt;
    substeps = substeps + 1;
    if substeps > 1000000
        error('HydroPol2D:GroundwaterSubcycling', 'Groundwater subcycling did not converge.');
    end
end
if ~(lateral || river) && any(fixed_mask)
    [groundwater, L, added] = hold_fixed_heads(groundwater, L, fixed_mask, config.groundwater_fixed_head_m(:));
    fixed_exchange = fixed_exchange + added;
end
surface_volume = surface_volume + seepage;
sat = hp2d_saturated_thickness(L.column, groundwater.head(host));
saturated_mm = sum((L.column.theta_s - L.column.theta_r).*sat*1000, 2);
groundwater.storage = accumarray(host, saturated_mm/1000.*area, [n 1]);
L.head_m = groundwater.head(:);
state.layered = L;
shape = size(state.fraction);
state.near_storage_m = reshape(L.water_mm(:, 1)/1000, shape);
state.root_storage_m = reshape(L.water_mm(:, 2)/1000, shape);
state.transmission_storage_m = reshape(L.water_mm(:, 3)/1000, shape);
if isfield(groundwater, 'cumulative_river_exchange_m')
    groundwater.cumulative_river_exchange_m = groundwater.cumulative_river_exchange_m + exchange_by_cell./cell_area;
end
if isfield(groundwater, 'cumulative_seepage_m')
    groundwater.cumulative_seepage_m = groundwater.cumulative_seepage_m + seepage./cell_area;
end
groundwater.last_river_exchange_rate_m_s = exchange_by_cell./cell_area/max(duration_s, eps);
groundwater.last_seepage_rate_m_s = seepage./cell_area/max(duration_s, eps);
diagnostics.max_flux_m3_s = maximum_flux;
diagnostics.substep_count = substeps;
diagnostics.river_exchange_m3 = river_total;
diagnostics.seepage_volume_m3 = sum(seepage);
diagnostics.fixed_head_exchange_m3 = fixed_exchange;
end

function [groundwater, L, exfiltration, unmet] = apply_cell_volume(groundwater, L, volume_m3)
R = hp2d_solve_shared_water_table(L.column, groundwater.head(:), L.water_mm, volume_m3(:), ...
    L.host, L.area_m2, L.ground_m, L.base_depth_m);
groundwater.head = R.head_m;
L.water_mm = R.water_mm;
exfiltration = R.exfiltration_m3;
unmet = R.unmet_removal_m3;
end

function [groundwater, L, added_m3] = hold_fixed_heads(groundwater, L, mask, fixed_head)
target = min(max(L.ground_m - groundwater.head(:), 0), L.base_depth_m);
target(mask) = min(max(L.ground_m(mask) - fixed_head(mask), 0), L.base_depth_m(mask));
[delta_mm, L.water_mm] = hp2d_storage_change_to_depth(L.column, groundwater.head(L.host), ...
    L.water_mm, target(L.host));
added = accumarray(L.host, delta_mm/1000.*L.area_m2, [numel(mask) 1]);
groundwater.head = L.ground_m - target;
added_m3 = sum(added(mask));
end

function [surface_volume, channel_volume, by_cell] = layered_river_exchange(mesh, surface_volume, ...
    channel_volume, head, dt, available_m3, config)
%Riverbed exchange per cell (m3, positive from the aquifer to surface water); same
%explicit conductance form as HydroPol2D-Python voronoi_layered._layered_river_exchange.
n = mesh.n_cells;
by_cell = zeros(n, 1);
thickness = max(double(config.riverbed_thickness_m), eps);
kbed = mesh.cell_riverbed_hydraulic_conductivity(:);
exchange_area = mesh.cell_riverbed_area(:);
available = max(available_m3(:), 0);
has_channel = isfield(mesh, 'channel') && mesh.channel.n_nodes > 0;
if has_channel, exchange_area(mesh.channel.host_cell(:)) = 0; end
resolved = find(exchange_area > 0);
if ~isempty(resolved)
    stage_all = mesh.surface_bed(:) + surface_volume./mesh.surface_area(:);
    if isfield(config, 'subgrid_tables') && ~isempty(fieldnames(config.subgrid_tables))
        [subgrid_stage, ~] = hp2d_voronoi_subgrid_cell_stage(config.subgrid_tables, surface_volume);
        selected = logical(config.subgrid_tables.cell_is_subgrid(:));
        stage_all(selected) = subgrid_stage(selected);
    end
    q = kbed(resolved).*exchange_area(resolved)/thickness.*(head(resolved) - stage_all(resolved));
    volume = q*dt;
    positive = volume > 0;
    volume(positive) = min(volume(positive), available(resolved(positive)));
    negative = volume < 0;
    volume(negative) = -min(-volume(negative), surface_volume(resolved(negative)));
    surface_volume(resolved) = surface_volume(resolved) + volume;
    by_cell(resolved) = by_cell(resolved) + volume;
    available(resolved) = available(resolved) - max(volume, 0);
end
if has_channel
    cells = mesh.channel.host_cell(:);
    stage = mesh.channel.bed(:) + channel_volume./mesh.channel.plan_area(:);
    q = kbed(cells).*mesh.cell_riverbed_area(cells)/thickness.*(head(cells) - stage);
    volume = q*dt;
    positive = volume > 0;
    volume(positive) = min(volume(positive), available(cells(positive)));
    negative = volume < 0;
    volume(negative) = -min(-volume(negative), channel_volume(negative));
    channel_volume = channel_volume + volume;
    by_cell = by_cell + accumarray(cells, volume, [n 1]);
end
end
