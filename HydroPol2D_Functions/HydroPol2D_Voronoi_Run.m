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
solver = lower(string(config.routing_solver));
if ~any(solver == ["local_inertial" "kinematic" "diffusive" "full_momentum"])
    error('HydroPol2D:VoronoiSolverUnavailable', ...
        'Voronoi routing_solver must be local_inertial, kinematic, diffusive, or full_momentum.');
end
[mesh, preflight] = HydroPol2D_Voronoi_Preflight(mesh_file, config);
n = mesh.n_cells;
if isfield(forcing, 'overlap_file')
    forcing.mapping = HydroPol2D_Read_Overlap(forcing.overlap_file, mesh_file);
    assert(size(forcing.mapping.raster_to_mesh, 1) == n, ...
        'HydroPol2D:InvalidVoronoiForcing', 'Overlap file does not match mesh.');
end
surface_volume = initial_vector(config.initial_surface_depth_m, n) .* mesh.surface_area(:);
edge_q = signed_vector(config.initial_surface_discharge_per_width_m2_s, mesh.n_edges);
hu = signed_vector(config.initial_surface_momentum_x_m2_s, n);
hv = signed_vector(config.initial_surface_momentum_y_m2_s, n);
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
t = 0; next_output = config.output_interval_s; next_progress = config.progress_interval_s;
step = 0; running_min_dt = inf; initial_mass = sum(surface_volume) + sum(channel_volume);
if groundwater.enabled, initial_mass = initial_mass + groundwater_mass(groundwater, mesh) + sum(groundwater.pending_exchange_m3); end
if hydrology.enabled, initial_mass = initial_mass + hydrology_mass(hydrology); end
% ---- optional HEC-RAS style sub-grid property tables --------------------
use_subgrid = logical(config.voronoi_subgrid_enabled);
subgrid_tables = struct();
if use_subgrid
    assert(solver == "local_inertial", 'HydroPol2D:SubgridSolver', ...
        ['Sub-grid property tables are implemented for the local-inertial solver ' ...
         'only; got %s.'], solver);
    subgrid_tables = HydroPol2D_Read_Subgrid_Tables(config.subgrid_table_path, ...
        mesh.n_cells, mesh.n_edges);
    fprintf('Sub-grid tables: %d of %d cells use the level-pool closure.\n', ...
        subgrid_tables.validation.subgrid_cells, mesh.n_cells);
end

times = []; depths = {}; edge_discharges = {}; momentum_x = {}; momentum_y = {};
face_flow_depths = {};   % the face depth PAIRED with each stored edge_q
channel_depths = {}; channel_discharges = {}; transition_discharges = {};
groundwater_heads = {}; groundwater_exchange = {}; groundwater_seepage = {};
hydrology_outputs = {}; surface_source_rates = {}; cumulative_surface_source = zeros(n,1);
cumulative_surface_sources = {}; potential_et_rates = {}; outlet_discharges = [];
[gauge_names,gauge_cells,gauge_edges,gauge_signs] = gauge_configuration(config,mesh);
gauge_discharges = {}; gauge_depths = {}; gauge_wse = {}; diagnostics = struct([]);
cached_et=[]; cached_open_water_et=[]; next_meteorology_update=0;
[initial_gauge_q,initial_gauge_depth,initial_gauge_wse] = ...
    gauge_values(gauge_cells,gauge_edges,gauge_signs,edge_q,mesh,surface_volume);
times=0; depths={surface_volume./mesh.surface_area(:)}; edge_discharges={edge_q};
face_flow_depths={zeros(mesh.n_edges,1)};
surface_source_rates={zeros(n,1)}; cumulative_surface_sources={zeros(n,1)};
potential_et_rates={zeros(n,1)}; outlet_discharges= ...
    sum(max(edge_q(mesh.edge_neighbor==0).*mesh.edge_length(mesh.edge_neighbor==0),0));
gauge_discharges={initial_gauge_q}; gauge_depths={initial_gauge_depth}; gauge_wse={initial_gauge_wse};
if solver=="full_momentum", momentum_x={hu}; momentum_y={hv}; end
if mesh.channel.n_nodes>0
    channel_depths={channel_volume./mesh.channel.plan_area(:)};
    channel_discharges={link_q}; transition_discharges={transition_q};
end
if groundwater.enabled
    groundwater_heads={groundwater.head};
    groundwater_exchange={groundwater.last_river_exchange_rate_m_s};
    groundwater_seepage={groundwater.last_seepage_rate_m_s};
end
if hydrology.enabled, hydrology_outputs={aggregate_hydrology(hydrology)}; end
time_tolerance=max(1e-9,100*eps(max(config.duration_s,1)));
while t < config.duration_s-time_tolerance
    cfl_state=struct('surface_volume_m3',surface_volume,'channel_volume_m3',channel_volume, ...
        'groundwater_head_m',groundwater.head,'surface_momentum_x_m2_s',hu,'surface_momentum_y_m2_s',hv);
    cfl_boundary=boundary_at(forcing,'surface_boundary',t,cfl_state,mesh);
    [dt,cfl] = stable_timestep(mesh, surface_volume, edge_q, channel_volume, link_q, transition_q, ...
        groundwater, roughness, hu, hv, cfl_boundary, config, t);
    stability_dt = min(dt,config.max_dt_s);
    if stability_dt < config.min_dt_s && config.duration_s - t > config.min_dt_s
        error('HydroPol2D:VoronoiTimestepFloor', 'Required timestep %.6g s is below minimum %.6g s.', stability_dt, config.min_dt_s);
    end
    dt = min([stability_dt, config.duration_s - t, next_output-t]);
    if groundwater.enabled
        remaining_groundwater=config.groundwater_update_interval_s-groundwater.elapsed_s;
        if remaining_groundwater>time_tolerance, dt=min(dt,remaining_groundwater); end
    end
    if config.hydrology_enabled && isfinite(config.meteorology_interval_s)
        remaining_meteorology=config.meteorology_interval_s-mod(t,config.meteorology_interval_s);
        if remaining_meteorology>time_tolerance, dt=min(dt,remaining_meteorology); end
    end
    if isfinite(config.forcing_interval_s)
        remaining_forcing=config.forcing_interval_s-mod(t,config.forcing_interval_s);
        if remaining_forcing<=time_tolerance, remaining_forcing=config.forcing_interval_s; end
        dt=min(dt,remaining_forcing);
    end

    state = struct('surface_volume_m3', surface_volume, 'channel_volume_m3', channel_volume, 'groundwater_head_m', groundwater.head, 'surface_momentum_x_m2_s',hu,'surface_momentum_y_m2_s',hv);
    surface_before_vertical = surface_volume;
    surface_channel_exchange=0;
    source = source_at(forcing, t, state, mesh, n);
    source_volume = source .* mesh.cell_area(:) .* dt;
    cumulative_surface_source = cumulative_surface_source + max(source,0).*dt;
    potential_et = zeros(n,1); open_water_evaporation = zeros(n,1);
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
        event_tolerance=max(1e-9,100*eps(max(t+dt,1)));
        update_due = groundwater.elapsed_s >= config.groundwater_update_interval_s-event_tolerance || ...
            config.duration_s-(t+dt) <= max(1e-9,100*eps(max(config.duration_s,1)));
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
    % Rainfall adds zero horizontal momentum.  Vertical withdrawals remove
    % water with its local horizontal velocity; this scaling is needed only
    % by the resolved full-momentum state.
    if solver == "full_momentum"
        withdrawal_scale = min(1, surface_volume ./ max(surface_before_vertical, config.dry_tolerance_m .* mesh.surface_area(:)));
        hu = hu .* withdrawal_scale; hv = hv .* withdrawal_scale;
    end
    surface_flux_evaluation_volume = surface_volume;
    surface_boundary = boundary_at(forcing, 'surface_boundary', t, state, mesh);
    if solver == "full_momentum"
        [surface_volume,hu,hv,edge_q,surface_diag,surface_boundary_diag] = Voronoi_Full_Momentum_Step( ...
            mesh,surface_volume,hu,hv,roughness,surface_boundary,dt,gravity=config.gravity, ...
            dry_tolerance_m=config.dry_tolerance_m,maximum_velocity_m_s=config.full_momentum_maximum_velocity_m_s);
    elseif solver == "local_inertial"
        if use_subgrid
            [surface_volume, edge_q, surface_diag] = Voronoi_Local_Inertial_Subgrid_Step( ...
                mesh, surface_volume, edge_q, roughness, subgrid_tables, dt, ...
                gravity=config.gravity, dry_tolerance_m=config.dry_tolerance_m, ...
                critical_flow=config.critical_flow);
        else
            [surface_volume, edge_q, surface_diag] = Voronoi_Local_Inertial_Step( ...
                mesh, surface_volume, edge_q, roughness, dt, gravity=config.gravity, ...
                dry_tolerance_m=config.dry_tolerance_m, critical_flow=config.critical_flow);
        end
    elseif solver == "kinematic"
        [surface_volume, edge_q, surface_diag] = Voronoi_Kinematic_Step( ...
            mesh, surface_volume, roughness, dt, gravity=config.gravity, ...
            dry_tolerance_m=config.dry_tolerance_m, critical_flow=config.critical_flow);
    else
        [surface_volume, edge_q, surface_diag] = Voronoi_Diffusive_Step( ...
            mesh, surface_volume, roughness, dt, dry_tolerance_m=config.dry_tolerance_m, ...
            slope_regularization=config.diffusive_slope_regularization);
    end
    if solver ~= "full_momentum"
        [surface_volume, edge_q, surface_boundary_diag] = Voronoi_Surface_Boundary_Step( ...
            mesh, surface_volume, edge_q, roughness, surface_boundary, dt, gravity=config.gravity, ...
            dry_tolerance_m=config.dry_tolerance_m, evaluation_volume=surface_flux_evaluation_volume, ...
            routing_solver=char(solver));
    end
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
    t = t + dt;
    if abs(t-next_output)<=time_tolerance, t=next_output; end
    step = step + 1;
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
    diagnostics(step).surface_source_volume_m3 = sum(source_volume);
    diagnostics(step).infiltration_volume_m3 = hydrology_diag.infiltration_volume_m3;
    diagnostics(step).actual_et_volume_m3 = hydrology_diag.actual_et_volume_m3;
    diagnostics(step).recharge_volume_m3 = hydrology_diag.recharge_volume_m3;
    diagnostics(step).capillary_volume_m3 = hydrology_diag.capillary_volume_m3;
    diagnostics(step).hydrology_mass_residual_m3 = hydrology_diag.mass_residual_m3;
    diagnostics(step).surface_channel_exchange_m3 = surface_channel_exchange;
    diagnostics(step).boundary_net_inflow_volume_m3 = boundary_volume;
    running_min_dt = min(running_min_dt,stability_dt);
    if isfinite(config.progress_interval_s) && t + time_tolerance >= next_progress
        live_depth = surface_volume ./ mesh.surface_area(:);
        [live_max_depth,live_max_cell] = max(live_depth);
        fprintf(['[Voronoi live] t=%.3f h step=%d dt=%.6g s min_dt=%.6g s ' ...
            'stable_dt=%.6g s CFLcell=%d (%.1f, %.1f) CFLdepth=%.6g m max_depth=%.6g m ' ...
            'max_depth_cell=%d max_velocity=%.6g m/s mass_residual=%.6g m3\n'], ...
            t/3600,step,dt,running_min_dt,stability_dt,cfl.cell_id,cfl.x_m,cfl.y_m,cfl.depth_m, ...
            live_max_depth,live_max_cell,surface_diag.max_velocity_m_s,mass-expected);
        if ~isempty(config.progress_checkpoint_file)
            live = struct('time_s',t,'step',step,'dt_s',dt,'stability_dt_s',stability_dt, ...
                'minimum_stability_dt_s',running_min_dt, ...
                'surface_depth_m',live_depth,'maximum_depth_m',live_max_depth, ...
                'maximum_depth_cell',live_max_cell,'maximum_velocity_m_s',surface_diag.max_velocity_m_s, ...
                'mass_residual_m3',mass-expected,'cfl',cfl); %#ok<NASGU>
            save(config.progress_checkpoint_file,'live','-v7.3');
        end
        next_progress = next_progress + config.progress_interval_s;
    end
    initial_mass = expected;
    if t + eps(t) >= next_output
        times(end+1,1) = t; %#ok<AGROW>
        depths{end+1,1} = surface_volume ./ mesh.surface_area(:); %#ok<AGROW>
        edge_discharges{end+1,1} = edge_q; %#ok<AGROW>
        if isfield(surface_diag,'face_flow_depth_m')
            face_flow_depths{end+1,1} = surface_diag.face_flow_depth_m; %#ok<AGROW>
        else
            face_flow_depths{end+1,1} = nan(mesh.n_edges,1); %#ok<AGROW>
        end
        surface_source_rates{end+1,1} = source; %#ok<AGROW>
        cumulative_surface_sources{end+1,1} = cumulative_surface_source; %#ok<AGROW>
        potential_et_rates{end+1,1} = potential_et; %#ok<AGROW>
        boundary_edges = mesh.edge_neighbor == 0;
        outlet_discharges(end+1,1) = sum(max(edge_q(boundary_edges).*mesh.edge_length(boundary_edges),0)); %#ok<AGROW>
        if solver == "full_momentum"
            momentum_x{end+1,1}=hu; momentum_y{end+1,1}=hv; %#ok<AGROW>
        end
        if mesh.channel.n_nodes > 0
            channel_depths{end+1,1} = channel_volume ./ mesh.channel.plan_area(:); %#ok<AGROW>
            channel_discharges{end+1,1} = link_q; %#ok<AGROW>
            transition_discharges{end+1,1} = transition_q; %#ok<AGROW>
        end
        if groundwater.enabled
            groundwater_heads{end+1,1}=groundwater.head; %#ok<AGROW>
            groundwater_exchange{end+1,1}=groundwater.last_river_exchange_rate_m_s; %#ok<AGROW>
            groundwater_seepage{end+1,1}=groundwater.last_seepage_rate_m_s; %#ok<AGROW>
        end
        if hydrology.enabled, hydrology_outputs{end+1,1}=aggregate_hydrology(hydrology); end %#ok<AGROW>
        [gauge_discharges{end+1,1},gauge_depths{end+1,1},gauge_wse{end+1,1}] = ... %#ok<AGROW>
            gauge_values(gauge_cells,gauge_edges,gauge_signs,edge_q,mesh,surface_volume);
        next_output = next_output + config.output_interval_s;
    end
end
results.mesh_file = mesh_file;
results.preflight = preflight;
results.time_s = times;
results.surface_depth_m = cat(2, depths{:});
results.edge_discharge_per_width_history_m2_s = cat_or_empty(edge_discharges, mesh.n_edges);
results.face_flow_depth_history_m = cat_or_empty(face_flow_depths, mesh.n_edges);
paired_depth = results.face_flow_depth_history_m;
if isempty(paired_depth) || any(~isfinite(paired_depth(:)))
    paired_depth = [];   % a solver that does not report it falls back
end
results.surface_velocity_m_s = HydroPol2D_Voronoi_Cell_Velocity(mesh, ...
    results.surface_depth_m, results.edge_discharge_per_width_history_m2_s, ...
    config.dry_tolerance_m, paired_depth);
results.surface_momentum_x_m2_s = cat_or_empty(momentum_x,n);
results.surface_momentum_y_m2_s = cat_or_empty(momentum_y,n);
if solver == "full_momentum"
    results.surface_velocity_x_m_s = results.surface_momentum_x_m2_s ./ max(results.surface_depth_m, config.dry_tolerance_m);
    results.surface_velocity_y_m_s = results.surface_momentum_y_m2_s ./ max(results.surface_depth_m, config.dry_tolerance_m);
    results.surface_velocity_m_s = sqrt(results.surface_velocity_x_m_s.^2 + results.surface_velocity_y_m_s.^2);
else
    results.surface_velocity_x_m_s = zeros(n,0); results.surface_velocity_y_m_s = zeros(n,0);
end
results.channel_depth_m = cat_or_empty(channel_depths, mesh.channel.n_nodes);
results.channel_discharge_history_m3_s = cat_or_empty(channel_discharges, mesh.channel.n_links);
results.transition_discharge_history_m3_s = cat_or_empty(transition_discharges, mesh.channel.n_transitions);
results.surface_source_rate_m_s = cat_or_empty(surface_source_rates,n);
results.cumulative_surface_source_m = cat_or_empty(cumulative_surface_sources,n);
results.potential_et_rate_m_s = cat_or_empty(potential_et_rates,n);
results.outlet_discharge_m3_s = outlet_discharges;
results.gauge_names = gauge_names;
results.gauge_discharge_m3_s = cat_or_empty(gauge_discharges,numel(gauge_names));
results.gauge_surface_depth_m = cat_or_empty(gauge_depths,numel(gauge_names));
results.gauge_water_surface_elevation_m = cat_or_empty(gauge_wse,numel(gauge_names));
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
results.map_time_s = selected_output_times(results.time_s,config.raster_output_interval_s);
if isfield(forcing,'overlap_file'), results.overlap_file=forcing.overlap_file; else, results.overlap_file=''; end
if ~isempty(config.output_netcdf)
    HydroPol2D_Write_Voronoi_Output(config.output_netcdf, mesh, results, config.overwrite_output);
end
end

function config = defaults(config)
values = struct('routing_solver','local_inertial','duration_s',3600,'initial_surface_depth_m',0, ...
    'initial_channel_depth_m',0,'surface_roughness',0.05,'min_dt_s',0.01,'max_dt_s',30, ...
    'output_interval_s',300,'courant',0.2,'gravity',9.81, ...
    'dry_tolerance_m',1e-6,'critical_flow',false,'groundwater_enabled',false, ...
    'voronoi_subgrid_enabled',[],'subgrid_table_path','', ...
    'initial_surface_discharge_per_width_m2_s',0,'initial_surface_momentum_x_m2_s',0, ...
    'initial_surface_momentum_y_m2_s',0,'full_momentum_maximum_velocity_m_s',10, ...
    'initial_channel_discharge_m3_s',0, ...
    'initial_groundwater_head_m',0,'aquifer_bottom_m',-10,'hydraulic_conductivity_m_s',1e-5, ...
    'specific_yield',0.2,'compute_backend','cpu','output_netcdf','', ...
    'overwrite_output',false,'forcing_interval_s',inf,'hydrology_enabled',false, ...
    'allow_legacy_mesh',false, ...
    'hydrology',struct(),'groundwater_update_interval_s',3600, ...
    'meteorology_interval_s',86400, ...
    'diffusive_slope_regularization',1e-4, ...
    'riverbed_hydraulic_conductivity_m_s',0,'riverbed_thickness_m',0.5, ...
    'resolved_river_mask',[],'raster_output_interval_s',3600,'gauges',struct([]), ...
    'progress_interval_s',inf,'progress_checkpoint_file','');
names = fieldnames(values);
for k = 1:numel(names)
    if ~isfield(config, names{k}), config.(names{k}) = values.(names{k}); end
end
subgrid_path_present = strlength(strtrim(string(config.subgrid_table_path))) > 0;
if isempty(config.voronoi_subgrid_enabled)
    % Migration path for prepared cases created before the explicit flag.
    config.voronoi_subgrid_enabled = subgrid_path_present;
else
    value=config.voronoi_subgrid_enabled;
    if islogical(value), value=double(value); end
    assert(isscalar(value) && isnumeric(value) && isfinite(value) && ismember(double(value),[0 1]), ...
        'HydroPol2D:InvalidVoronoiConfiguration', ...
        'voronoi_subgrid_enabled must be 0 or 1.');
    config.voronoi_subgrid_enabled=logical(value);
end
assert(~config.voronoi_subgrid_enabled || subgrid_path_present, ...
    'HydroPol2D:MissingVoronoiSubgridTables', ...
    'Voronoi subgrid is enabled but subgrid_table_path is empty.');
assert(config.voronoi_subgrid_enabled || ~subgrid_path_present, ...
    'HydroPol2D:VoronoiSubgridConfigurationConflict', ...
    'subgrid_table_path is set while Voronoi subgrid is disabled.');
assert(config.groundwater_update_interval_s>0 && config.meteorology_interval_s>0 && config.riverbed_thickness_m>0, ...
    'HydroPol2D:InvalidVoronoiConfiguration','Groundwater interval and riverbed thickness must be positive.');
assert(config.diffusive_slope_regularization>0 && isfinite(config.diffusive_slope_regularization), ...
    'HydroPol2D:InvalidVoronoiConfiguration','diffusive_slope_regularization must be positive and finite.');
assert(config.full_momentum_maximum_velocity_m_s>0 && isfinite(config.full_momentum_maximum_velocity_m_s), ...
    'HydroPol2D:InvalidVoronoiConfiguration','full_momentum_maximum_velocity_m_s must be positive and finite.');
assert(isfinite(config.courant) && config.courant>0 && config.courant<=0.3, ...
    'HydroPol2D:InvalidVoronoiCourant', ...
    'Voronoi Courant must be in (0, 0.3]; use 0.2 unless a validated case requires less.');
ratio=config.raster_output_interval_s/config.output_interval_s;
assert(ratio>=1 && abs(ratio-round(ratio))<=1e-9*max(ratio,1), ...
    'HydroPol2D:InvalidVoronoiConfiguration', ...
    'raster_output_interval_s must be an integer multiple of output_interval_s.');
end

function selected = selected_output_times(times,interval)
if isempty(times), selected=times; return; end
keep=abs(times./interval-round(times./interval))<=1e-9;
keep(end)=true; selected=times(keep);
end

function [names,cells,edges,signs] = gauge_configuration(config,mesh)
gauges=config.gauges;
if isempty(gauges), names=strings(0,1); cells=zeros(0,1); edges=zeros(0,1); signs=zeros(0,1); return; end
n_gauges=numel(gauges); names=strings(n_gauges,1); cells=nan(n_gauges,1);
edges=nan(n_gauges,1); signs=ones(n_gauges,1);
for k=1:n_gauges
    if isfield(gauges(k),'name'), names(k)=string(gauges(k).name); else, names(k)="gauge_"+k; end
    if isfield(gauges(k),'cell_id'), cells(k)=double(gauges(k).cell_id); end
    if isfield(gauges(k),'edge_id'), edges(k)=double(gauges(k).edge_id); end
    if isfield(gauges(k),'sign'), signs(k)=double(gauges(k).sign); end
end
assert(all(isnan(cells) | (cells>=1 & cells<=mesh.n_cells & cells==fix(cells))) && ...
    all(isnan(edges) | (edges>=1 & edges<=mesh.n_edges & edges==fix(edges))) && ...
    all(abs(signs)==1) && all(~isnan(cells) | ~isnan(edges)), ...
    'HydroPol2D:InvalidVoronoiConfiguration','Each gauge needs a valid cell_id or edge_id.');
end

function [discharge,depth,wse] = gauge_values(cells,edges,signs,edge_q,mesh,surface_volume)
n=numel(cells); discharge=nan(n,1); depth=nan(n,1); wse=nan(n,1);
for k=1:n
    if ~isnan(edges(k)), discharge(k)=signs(k)*edge_q(edges(k))*mesh.edge_length(edges(k)); end
    if ~isnan(cells(k))
        depth(k)=surface_volume(cells(k))/mesh.surface_area(cells(k));
        wse(k)=mesh.surface_bed(cells(k))+depth(k);
    end
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

function [dt,diagnostic] = stable_timestep(mesh, surface_volume, edge_q, channel_volume, link_q, transition_q, groundwater, roughness, hu, hv, boundary, config, time_s)
depth = max(surface_volume ./ mesh.surface_area(:), 0);
internal = mesh.edge_neighbor > 0;
if strcmpi(config.routing_solver, 'kinematic')
    [kinematic_q, face_depth] = Voronoi_Kinematic_Face_Flux(mesh, surface_volume, roughness, ...
        gravity=config.gravity, dry_tolerance_m=config.dry_tolerance_m, critical_flow=config.critical_flow);
    out_cell = mesh.edge_owner(internal);
    neighbor = mesh.edge_neighbor(internal);
    negative = kinematic_q(internal) < 0;
    out_cell(negative) = neighbor(negative);
    face_signal = mesh.edge_length(internal) .* (5/3) .* abs(kinematic_q(internal)) ./ ...
        max(face_depth(internal), config.dry_tolerance_m);
    cell_signal = accumarray(out_cell, face_signal, [mesh.n_cells 1], @sum, 0);
    cell_signal = cell_signal + kinematic_boundary_signal(mesh,depth,roughness,boundary,config.gravity,config.dry_tolerance_m);
elseif strcmpi(config.routing_solver, 'diffusive')
    [~, ~, conductance] = Voronoi_Diffusive_Face_Flux(mesh, surface_volume, roughness, ...
        dry_tolerance_m=config.dry_tolerance_m, slope_regularization=config.diffusive_slope_regularization);
    cell_signal = accumarray(mesh.edge_owner(internal), conductance(internal), [mesh.n_cells 1], @sum, 0);
    cell_signal = cell_signal + accumarray(mesh.edge_neighbor(internal), conductance(internal), [mesh.n_cells 1], @sum, 0);
    cell_signal = cell_signal + diffusive_boundary_signal(mesh,depth,roughness,boundary,config.gravity,config.dry_tolerance_m,config.diffusive_slope_regularization);
elseif strcmpi(config.routing_solver, 'full_momentum')
    speed=sqrt(hu(:).^2+hv(:).^2)./max(depth,config.dry_tolerance_m) + sqrt(config.gravity.*max(depth,config.dry_tolerance_m));
    face_signal=mesh.edge_length(:).*speed(mesh.edge_owner);
    face_signal(internal)=mesh.edge_length(internal).*max(speed(mesh.edge_owner(internal)),speed(mesh.edge_neighbor(internal)));
    cell_signal=accumarray(mesh.edge_owner(:),face_signal,[mesh.n_cells 1],@sum,0);
    cell_signal=cell_signal+accumarray(mesh.edge_neighbor(internal),face_signal(internal),[mesh.n_cells 1],@sum,0);
else
    % Finite-volume CFL: each control volume is constrained by the summed
    % signal capacity of every incident face.
    stage = mesh.surface_bed(:) + depth;
    face_depth = depth(mesh.edge_owner);
    face_depth(internal) = max( ...
        max(stage(mesh.edge_owner(internal)), stage(mesh.edge_neighbor(internal))) - ...
        max(mesh.surface_bed(mesh.edge_owner(internal)), mesh.surface_bed(mesh.edge_neighbor(internal))), 0);
    % Use the discharge that the upcoming local-inertial update can
    % actually apply at the current face depth.  Reusing a wet-step flux
    % after its face has drained creates a fictitious q/h velocity and can
    % collapse the CFL timestep before the solver zeros/clips that flux.
    internal_depth = face_depth(internal);
    cfl_q = abs(edge_q(internal));
    cfl_q(internal_depth <= config.dry_tolerance_m) = 0;
    if config.critical_flow
        cfl_q = min(cfl_q, internal_depth .* sqrt(config.gravity .* internal_depth));
    end
    face_speed = cfl_q ./ max(internal_depth, config.dry_tolerance_m);
    face_signal = mesh.edge_length(internal) .* ...
        (face_speed + sqrt(config.gravity .* max(internal_depth, config.dry_tolerance_m)));
    cell_signal = accumarray(mesh.edge_owner(internal), face_signal, [mesh.n_cells 1], @sum, 0);
    cell_signal = cell_signal + accumarray(mesh.edge_neighbor(internal), face_signal, [mesh.n_cells 1], @sum, 0);

    % Closed exterior faces do not exchange water and therefore do not
    % constrain the explicit update. Add only configured open boundaries.
    [boundary_edge,~,~] = surface_boundary_parameters(mesh,boundary);
    if ~isempty(boundary_edge)
        boundary_owner = mesh.edge_owner(boundary_edge);
        boundary_depth = depth(boundary_owner);
        boundary_q = abs(edge_q(boundary_edge));
        boundary_q(boundary_depth <= config.dry_tolerance_m) = 0;
        boundary_depth = max(boundary_depth, config.dry_tolerance_m);
        boundary_signal = mesh.edge_length(boundary_edge) .* ...
            (boundary_q ./ boundary_depth + sqrt(config.gravity .* boundary_depth));
        cell_signal = cell_signal + accumarray(boundary_owner, boundary_signal, ...
            [mesh.n_cells 1], @sum, 0);
    end
end

% Neal links are a one-dimensional finite-volume graph. Their timestep is
% constrained by every incident link and transition, not only the shortest
% link's gravity-wave speed. This branch is local-inertial only by preflight.
channel_dt=inf;
if mesh.channel.n_nodes > 0
    node_signal=zeros(mesh.channel.n_nodes,1);
    if mesh.channel.n_links > 0
        channel_depth=max(channel_volume./mesh.channel.plan_area(:),config.dry_tolerance_m);
        up=mesh.channel.link_up(:); down=mesh.channel.link_down(:);
        link_depth=max(channel_depth(up),channel_depth(down));
        area=mesh.channel.link_width(:).*link_depth;
        signal_speed=abs(link_q(:))./max(area,eps)+sqrt(config.gravity.*link_depth);
        signal=mesh.channel.link_width(:).*signal_speed;
        node_signal=node_signal+accumarray(up,signal,[mesh.channel.n_nodes 1],@sum,0);
        node_signal=node_signal+accumarray(down,signal,[mesh.channel.n_nodes 1],@sum,0);
        channel_dt=min(channel_dt,config.courant*min(mesh.channel.link_length(:)./signal_speed));
    end
    if mesh.channel.n_transitions > 0
        node=mesh.channel.transition_node(:); cell_id=mesh.channel.transition_cell(:);
        channel_depth=max(channel_volume./mesh.channel.plan_area(:),0);
        node_wse=mesh.channel.bed(node)+channel_depth(node);
        cell_wse=mesh.surface_bed(cell_id)+depth(cell_id);
        transition_depth=max(max(node_wse,cell_wse)-max(mesh.channel.bed(node),mesh.surface_bed(cell_id)),config.dry_tolerance_m);
        area=mesh.channel.transition_width(:).*transition_depth;
        signal_speed=abs(transition_q(:))./max(area,eps)+sqrt(config.gravity.*transition_depth);
        signal=mesh.channel.transition_width(:).*signal_speed;
        node_signal=node_signal+accumarray(node,signal,[mesh.channel.n_nodes 1],@sum,0);
        cell_signal=cell_signal+accumarray(cell_id,signal,[mesh.n_cells 1],@sum,0);
        channel_dt=min(channel_dt,config.courant*min(mesh.channel.transition_length(:)./signal_speed));
    end
    active_node=node_signal>0;
    if any(active_node)
        channel_dt=min(channel_dt,config.courant*min(mesh.channel.plan_area(active_node)./node_signal(active_node)));
    end
end
active = cell_signal > 0;
if any(active)
    courant=config.courant;
    if strcmpi(config.routing_solver,'full_momentum'), courant=min(courant,0.45); end
    active_id=find(active);
    [surface_ratio,control_index]=min(mesh.surface_area(active) ./ cell_signal(active));
    control_cell=active_id(control_index);
    dt = courant * surface_ratio;
else
    dt = inf; control_cell=1;
end
diagnostic=struct('cell_id',control_cell,'x_m',mesh.cell_x(control_cell), ...
    'y_m',mesh.cell_y(control_cell),'depth_m',depth(control_cell), ...
    'surface_dt_s',dt,'channel_dt_s',channel_dt);
dt=min(dt,channel_dt);
end

function signal = kinematic_boundary_signal(mesh,depth,roughness,boundary,gravity,dry)
signal=zeros(mesh.n_cells,1);
[edge_id,types,values]=surface_boundary_parameters(mesh,boundary);
if isempty(edge_id), return; end
owner=mesh.edge_owner(edge_id); h=depth(owner); n=roughness(owner);
n(~isfinite(n) | n<=0)=1e-6;
normal=types=="normal_flow";
critical=types=="critical_flow";
if any(types=="stage")
    error('HydroPol2D:KinematicStageBoundaryUnavailable', ...
        'A kinematic-wave boundary cannot prescribe stage.');
end
speed=zeros(numel(edge_id),1);
q=h(normal).^(5/3)./n(normal).*sqrt(max(values(normal),0));
speed(normal)=(5/3).*abs(q)./max(h(normal),dry);
speed(critical)=(3/2).*sqrt(gravity.*h(critical));
signal=accumarray(owner,mesh.edge_length(edge_id).*speed,[mesh.n_cells 1],@sum,0);
end

function signal = diffusive_boundary_signal(mesh,depth,roughness,boundary,gravity,dry,slope_floor)
signal=zeros(mesh.n_cells,1);
[edge_id,types,values]=surface_boundary_parameters(mesh,boundary);
if isempty(edge_id), return; end
owner=mesh.edge_owner(edge_id); h=depth(owner); n=roughness(owner);
n(~isfinite(n) | n<=0)=1e-6;
stage=types=="stage"; normal=types=="normal_flow"; critical=types=="critical_flow";
conductance=zeros(numel(edge_id),1);
if any(stage)
    eta=mesh.surface_bed(owner(stage))+h(stage);
    hface=max(max(eta,values(stage))-mesh.surface_bed(owner(stage)),0);
    slope=(eta-values(stage))./mesh.edge_distance(edge_id(stage));
    coefficient=hface.^(5/3)./n(stage);
    conductance(stage)=mesh.edge_length(edge_id(stage)).*( ...
        0.5.*coefficient./sqrt(max(abs(slope),slope_floor))./mesh.edge_distance(edge_id(stage)) + ...
        (5/3).*hface.^(2/3)./n(stage).*sqrt(abs(slope)));
end
if any(normal)
    conductance(normal)=mesh.edge_length(edge_id(normal)).*(5/3).*h(normal).^(2/3)./n(normal).*sqrt(max(values(normal),0));
end
if any(critical)
    conductance(critical)=mesh.edge_length(edge_id(critical)).*(3/2).*sqrt(gravity.*h(critical));
end
conductance(~isfinite(conductance) | h<=dry)=0;
signal=accumarray(owner,conductance,[mesh.n_cells 1],@sum,0);
end

function [edge_id,types,values] = surface_boundary_parameters(mesh,boundary)
edge_id=zeros(0,1); types=strings(0,1); values=zeros(0,1);
if isempty(boundary) || ~isfield(boundary,'edge_id') || isempty(boundary.edge_id), return; end
edge_id = HydroPol2D_Voronoi_Boundary_Ids(boundary, mesh.n_edges);
assert(all(edge_id>=1 & edge_id<=mesh.n_edges & mesh.edge_neighbor(edge_id)==0), ...
    'HydroPol2D:InvalidBoundaryEdge', ...
    'Boundary edge ids must identify boundary faces of this mesh.');
types=string(boundary.type(:)); if isscalar(types), types=repmat(types,numel(edge_id),1); end
if isfield(boundary,'value')
    values=double(boundary.value(:)); if isscalar(values), values=repmat(values,numel(edge_id),1); end
else
    values=zeros(numel(edge_id),1);
end
assert(numel(types)==numel(edge_id) && numel(values)==numel(edge_id) && all(isfinite(values)));
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
