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
hydrology = struct('enabled',false);
if config.hydrology_enabled
    if groundwater.enabled
        hydrology_head = groundwater.head;
    else
        hydrology_head = mesh.surface_bed(:) - 2;
    end
    hydrology = Voronoi_Hydrology_Initialize(mesh,config.hydrology,hydrology_head);
end
t = 0; next_output = 0; step = 0; initial_mass = sum(surface_volume) + sum(channel_volume);
if groundwater.enabled, initial_mass = initial_mass + groundwater_mass(groundwater, mesh) + sum(groundwater.pending_exchange_m3); end
if hydrology.enabled, initial_mass = initial_mass + hydrology_mass(hydrology); end
times = []; depths = {}; channel_depths = {}; groundwater_heads = {}; groundwater_exchange = {}; groundwater_seepage = {}; hydrology_outputs = {}; diagnostics = struct([]);
cached_et=[]; cached_open_water_et=[]; next_meteorology_update=0;
time_tolerance=max(1e-9,100*eps(max(config.duration_s,1)));
while t < config.duration_s-time_tolerance
    dt = stable_timestep(mesh, surface_volume, channel_volume, groundwater, config, t);
    dt = min([dt, config.max_dt_s, config.duration_s - t]);
    if dt < config.min_dt_s && config.duration_s - t > config.min_dt_s
        error('HydroPol2D:VoronoiTimestepFloor', 'Required timestep %.6g s is below minimum %.6g s.', dt, config.min_dt_s);
    end

    state = struct('surface_volume_m3', surface_volume, 'channel_volume_m3', channel_volume, 'groundwater_head_m', groundwater.head);
    surface_channel_exchange=0;
    source = source_at(forcing, t, state, mesh, n);
    source_volume = source .* mesh.cell_area(:) .* dt;
    hydrology_diag = empty_hydrology_diagnostics();
    if hydrology.enabled
        assert(all(source >= 0), 'HydroPol2D:InvalidVoronoiForcing', ...
            'Hydrology-enabled surface forcing is precipitation and cannot be negative.');
        has_meteorology=isfield(forcing,'meteorology') || isfield(forcing,'raster_meteorology');
        if has_meteorology && ~isfield(forcing,'potential_et_m_s') && ...
                ~isfield(forcing,'raster_potential_et_m_s')
            if isempty(cached_et) || t+time_tolerance>=next_meteorology_update
                if isfield(forcing,'meteorology')
                    meteorology=forcing.meteorology; mapping=[];
                else
                    assert(isfield(forcing,'mapping'),'HydroPol2D:InvalidVoronoiForcing', ...
                        'raster_meteorology requires an overlap mapping.');
                    meteorology=forcing.raster_meteorology; mapping=forcing.mapping;
                end
                [cached_et,cached_open_water_et]=meteorological_et_at(meteorology,t,state,mesh,n,mapping);
                next_meteorology_update=(floor(t/config.meteorology_interval_s)+1)*config.meteorology_interval_s;
            end
            potential_et=cached_et; open_water_evaporation=cached_open_water_et;
        else
            potential_et = forcing_at(forcing,'potential_et_m_s','raster_potential_et_m_s',t,state,mesh,n);
            open_water_evaporation = forcing_at(forcing,'open_water_evaporation_m_s', ...
                'raster_open_water_evaporation_m_s',t,state,mesh,n,potential_et);
        end
        [surface_volume,hydrology,groundwater,hydrology_diag] = Voronoi_Hydrology_Step( ...
            mesh,surface_volume,hydrology,groundwater,source,potential_et,open_water_evaporation,dt);
    else
        if any(source_volume < 0), source_volume=max(source_volume,-surface_volume); end
        surface_volume=surface_volume+source_volume;
    end
    if mesh.channel.n_nodes>0
        [surface_volume,channel_volume,~,exchange]=Voronoi_Equilibrate_Channel_Storage(mesh,surface_volume,channel_volume);
        surface_channel_exchange=surface_channel_exchange+exchange;
        if hydrology.enabled
            host=mesh.channel.host_cell(:);
            channel_evaporation=min(channel_volume,open_water_evaporation(host).*dt.*mesh.channel.plan_area(:));
            channel_volume=channel_volume-channel_evaporation;
            evaporation_by_cell=accumarray(host,channel_evaporation,[n 1],@sum,0);
            hydrology.cumulative_surface_evaporation_m=hydrology.cumulative_surface_evaporation_m+evaporation_by_cell./mesh.cell_area(:);
            hydrology.last_surface_evaporation_rate_m_s=hydrology.last_surface_evaporation_rate_m_s+evaporation_by_cell./mesh.cell_area(:)./dt;
            hydrology_diag.actual_et_volume_m3=hydrology_diag.actual_et_volume_m3+sum(channel_evaporation);
        end
    end
    groundwater_source_volume = zeros(n,1);
    if groundwater.enabled
        recharge = groundwater_source_at(forcing, t, state, mesh, n);
        groundwater_source_volume = recharge .* mesh.cell_area(:) .* dt;
        available = groundwater.storage + groundwater.pending_exchange_m3;
        groundwater_source_volume = max(groundwater_source_volume,-available);
        groundwater.pending_exchange_m3 = groundwater.pending_exchange_m3 + groundwater_source_volume;
        groundwater.elapsed_s = groundwater.elapsed_s + dt;
        update_due = groundwater.elapsed_s >= config.groundwater_update_interval_s-10*eps(max(t+dt,1)) || ...
            config.duration_s-(t+dt) <= 10*eps(max(config.duration_s,1));
        if update_due
            [surface_volume,channel_volume,groundwater,groundwater_diag] = Voronoi_Groundwater_Advance( ...
                mesh,surface_volume,channel_volume,groundwater,groundwater.elapsed_s,config);
            groundwater.elapsed_s=0;
        else
            groundwater_diag=empty_groundwater_diagnostics();
        end
    else
        groundwater_diag=empty_groundwater_diagnostics();
    end
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
    if groundwater.enabled, mass = mass + sum(groundwater.storage) + sum(groundwater.pending_exchange_m3); end
    if hydrology.enabled, mass = mass + hydrology_mass(hydrology); end
    boundary_volume = surface_boundary_diag.net_inflow_volume_m3 + channel_boundary_diag.net_inflow_volume_m3;
    expected = initial_mass + sum(source_volume) + sum(groundwater_source_volume) + boundary_volume ...
        - hydrology_diag.actual_et_volume_m3 - hydrology_diag.deep_drainage_volume_m3;
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
    diagnostics(step).groundwater_substep_count = groundwater_diag.substep_count;
    diagnostics(step).groundwater_river_exchange_m3 = groundwater_diag.river_exchange_m3;
    diagnostics(step).groundwater_seepage_volume_m3 = groundwater_diag.seepage_volume_m3;
    diagnostics(step).precipitation_volume_m3 = hydrology_diag.precipitation_volume_m3;
    diagnostics(step).infiltration_volume_m3 = hydrology_diag.infiltration_volume_m3;
    diagnostics(step).actual_et_volume_m3 = hydrology_diag.actual_et_volume_m3;
    diagnostics(step).recharge_volume_m3 = hydrology_diag.recharge_volume_m3;
    diagnostics(step).capillary_volume_m3 = hydrology_diag.capillary_volume_m3;
    diagnostics(step).hydrology_mass_residual_m3 = hydrology_diag.mass_residual_m3;
    diagnostics(step).surface_channel_exchange_m3 = surface_channel_exchange;
    diagnostics(step).boundary_net_inflow_volume_m3 = boundary_volume;
    initial_mass = expected;
    if t + eps(t) >= next_output
        times(end+1,1) = t; %#ok<AGROW>
        depths{end+1,1} = surface_volume ./ mesh.surface_area(:); %#ok<AGROW>
        if mesh.channel.n_nodes > 0
            channel_depths{end+1,1} = channel_volume ./ mesh.channel.plan_area(:); %#ok<AGROW>
        end
        if groundwater.enabled
            groundwater_heads{end+1,1}=groundwater.head; %#ok<AGROW>
            groundwater_exchange{end+1,1}=groundwater.last_river_exchange_rate_m_s; %#ok<AGROW>
            groundwater_seepage{end+1,1}=groundwater.last_seepage_rate_m_s; %#ok<AGROW>
        end
        if hydrology.enabled, hydrology_outputs{end+1,1}=aggregate_hydrology(hydrology); end %#ok<AGROW>
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
results.groundwater_river_exchange_m_s=cat_or_empty(groundwater_exchange,n);
results.groundwater_seepage_m_s=cat_or_empty(groundwater_seepage,n);
results.hydrology = pack_hydrology_outputs(hydrology_outputs,hydrology,n);
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
    'overwrite_output',false,'forcing_interval_s',inf,'hydrology_enabled',false, ...
    'hydrology',struct(),'groundwater_update_interval_s',3600, ...
    'meteorology_interval_s',86400, ...
    'riverbed_hydraulic_conductivity_m_s',0,'riverbed_thickness_m',0.5, ...
    'resolved_river_mask',[]);
names = fieldnames(values);
for k = 1:numel(names)
    if ~isfield(config, names{k}), config.(names{k}) = values.(names{k}); end
end
assert(config.groundwater_update_interval_s>0 && config.meteorology_interval_s>0 && config.riverbed_thickness_m>0, ...
    'HydroPol2D:InvalidVoronoiConfiguration','Groundwater interval and riverbed thickness must be positive.');
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

function value = forcing_at(forcing,cell_name,raster_name,time_s,state,mesh,n,default_value)
if nargin<8, default_value=0; end
if isfield(forcing,cell_name)
    value=forcing.(cell_name);
elseif isfield(forcing,raster_name) && isfield(forcing,'mapping')
    value=forcing.(raster_name);
    if isa(value,'function_handle'), value=value(time_s,state,mesh); end
    value=forcing.mapping.raster_to_mesh*double(value(:));
    return
else
    value=default_value;
end
if isa(value,'function_handle'), value=value(time_s,state,mesh); end
if isscalar(value), value=repmat(double(value),n,1); else, value=double(value(:)); end
assert(numel(value)==n && all(isfinite(value) & value>=0), ...
    'HydroPol2D:InvalidVoronoiForcing','Invalid %s forcing.',cell_name);
end

function [et_m_s,open_water_m_s] = meteorological_et_at(meteorology,time_s,state,mesh,n,mapping)
if nargin<6, mapping=[]; end
if isa(meteorology,'function_handle'), meteorology=meteorology(time_s,state,mesh); end
if ~isempty(mapping)
    names=fieldnames(meteorology);
    for k=1:numel(names)
        value=meteorology.(names{k});
        if ~isscalar(value), meteorology.(names{k})=mapping.raster_to_mesh*double(value(:)); end
    end
end
temperature=meteorology_field(meteorology,'temperature_c',n);
maximum=meteorology_field(meteorology,'maximum_temperature_c',n);
minimum=meteorology_field(meteorology,'minimum_temperature_c',n);
day=meteorology_field(meteorology,'day_of_year',n);
latitude=meteorology_field(meteorology,'latitude_deg',n);
wind=meteorology_field(meteorology,'wind_speed_m_s',n);
humidity=meteorology_field(meteorology,'relative_humidity_pct',n,50);
krs=meteorology_field(meteorology,'krs',n,0.16);
albedo=meteorology_field(meteorology,'albedo',n,0.23);
ground_heat=meteorology_field(meteorology,'ground_heat_flux_mj_m2_day',n,0);
assert(all(day==day(1)),'HydroPol2D:InvalidVoronoiForcing','day_of_year must be spatially uniform.');
[et_mm_d,open_mm_d]=Evapotranspiration(mesh.surface_bed(:),temperature,maximum,minimum, ...
    day(1),latitude,wind,humidity,krs,albedo,ground_heat);
et_m_s=max(et_mm_d(:),0)./1000/86400;
open_water_m_s=max(open_mm_d(:),0)./1000/86400;
assert(all(isfinite(et_m_s) & isfinite(open_water_m_s)), ...
    'HydroPol2D:InvalidVoronoiForcing','Penman-Monteith returned non-finite ET.');
end

function value = meteorology_field(meteorology,name,n,default_value)
if nargin<4
    assert(isfield(meteorology,name),'HydroPol2D:InvalidVoronoiForcing', ...
        'Meteorological forcing is missing %s.',name);
    value=meteorology.(name);
elseif isfield(meteorology,name)
    value=meteorology.(name);
else
    value=default_value;
end
if isscalar(value), value=repmat(double(value),n,1); else, value=double(value(:)); end
assert(numel(value)==n && all(isfinite(value)), ...
    'HydroPol2D:InvalidVoronoiForcing','Invalid meteorological field %s.',name);
end

function boundary = boundary_at(forcing, name, time_s, state, mesh)
if ~isfield(forcing, name), boundary = struct(); return; end
boundary = forcing.(name);
if isa(boundary, 'function_handle'), boundary = boundary(time_s, state, mesh); end
end

function groundwater = initialize_groundwater(config, mesh)
groundwater.enabled = logical(config.groundwater_enabled);
if ~groundwater.enabled
    groundwater.head = zeros(0,1); groundwater.storage = zeros(0,1);
    groundwater.pending_exchange_m3=zeros(0,1); groundwater.elapsed_s=0; return
end
n = mesh.n_cells;
groundwater.head = signed_vector(config.initial_groundwater_head_m, n);
groundwater.bottom = signed_vector(config.aquifer_bottom_m, n);
groundwater.hydraulic_conductivity = initial_vector(config.hydraulic_conductivity_m_s, n);
groundwater.specific_yield = initial_vector(config.specific_yield, n);
assert(all(groundwater.head >= groundwater.bottom));
groundwater.storage = groundwater.specific_yield .* (groundwater.head - groundwater.bottom) .* mesh.cell_area(:);
groundwater.pending_exchange_m3=zeros(n,1);
groundwater.elapsed_s=0;
groundwater.last_river_exchange_rate_m_s=zeros(n,1);
groundwater.last_seepage_rate_m_s=zeros(n,1);
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
    remaining_groundwater=config.groundwater_update_interval_s-groundwater.elapsed_s;
    if remaining_groundwater>10*eps(max(time_s,1)), dt=min(dt,remaining_groundwater); end
end
if config.hydrology_enabled && isfinite(config.meteorology_interval_s)
    remaining_meteorology=config.meteorology_interval_s-mod(time_s,config.meteorology_interval_s);
    if remaining_meteorology>10*eps(max(time_s,1)), dt=min(dt,remaining_meteorology); end
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

function total = hydrology_mass(state)
total=sum((state.canopy_storage_m+state.near_storage_m+state.root_storage_m+ ...
    state.transmission_storage_m).*state.area_m2,'all');
end

function output = aggregate_hydrology(state)
weighted=@(value) sum(value.*state.fraction,2);
output.soil_water_m=weighted(state.near_storage_m+state.root_storage_m+state.transmission_storage_m);
output.canopy_storage_m=weighted(state.canopy_storage_m);
output.cumulative_infiltration_m=weighted(state.cumulative_infiltration_m);
output.cumulative_recharge_m=weighted(state.cumulative_recharge_m);
output.cumulative_actual_et_m=weighted(state.cumulative_actual_et_m)+state.cumulative_surface_evaporation_m;
output.infiltration_rate_m_s=weighted(state.last_infiltration_rate_m_s);
output.recharge_rate_m_s=weighted(state.last_recharge_rate_m_s);
output.capillary_rate_m_s=weighted(state.last_capillary_rate_m_s);
output.actual_et_rate_m_s=weighted(state.last_actual_et_rate_m_s+state.last_canopy_evaporation_rate_m_s) ...
    + state.last_surface_evaporation_rate_m_s;
end

function output = pack_hydrology_outputs(saved,state,n)
output=struct();
if ~state.enabled, return; end
names=fieldnames(saved{1});
for k=1:numel(names)
    output.(names{k})=cell2mat(cellfun(@(item) item.(names{k}),saved,'UniformOutput',false)');
end
output.hru_fraction=state.fraction;
output.final_canopy_storage_m=state.canopy_storage_m;
output.final_near_surface_storage_m=state.near_storage_m;
output.final_root_zone_storage_m=state.root_storage_m;
output.final_transmission_storage_m=state.transmission_storage_m;
assert(size(output.soil_water_m,1)==n);
end

function diagnostics = empty_hydrology_diagnostics()
diagnostics=struct('precipitation_volume_m3',0,'infiltration_volume_m3',0, ...
    'actual_et_volume_m3',0,'recharge_volume_m3',0,'capillary_volume_m3',0, ...
    'deep_drainage_volume_m3',0,'saturation_excess_volume_m3',0,'mass_residual_m3',0);
end

function diagnostics = empty_groundwater_diagnostics()
diagnostics=struct('max_flux_m3_s',0,'substep_count',0,'river_exchange_m3',0, ...
    'seepage_volume_m3',0,'mass_change_m3',0);
end
