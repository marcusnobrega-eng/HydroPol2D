function [surface_volume, channel_volume, groundwater, diagnostics] = Voronoi_Groundwater_Advance( ...
    mesh, surface_volume, channel_volume, groundwater, duration_s, config)
%VORONOI_GROUNDWATER_ADVANCE Apply recharge then subcycle polygon Boussinesq flow.

if ~groundwater.enabled || duration_s <= 0
    diagnostics = empty_diagnostics();
    return
end
if ~isfield(config,'groundwater_relative_tolerance'), config.groundwater_relative_tolerance=1e-3; end
if ~isfield(config,'groundwater_head_tolerance_m'), config.groundwater_head_tolerance_m=1e-3; end
if ~isfield(config,'groundwater_courant'), config.groundwater_courant=0.25; end
if ~isfield(config,'min_dt_s'), config.min_dt_s=0.01; end

pending = groundwater.pending_exchange_m3;
assert(all(groundwater.storage + pending >= -1e-8), ...
    'HydroPol2D:GroundwaterMass', 'Vadose exchange requested more groundwater than is stored.');
groundwater.storage = max(groundwater.storage + pending,0);
groundwater.pending_exchange_m3(:) = 0;
groundwater.head = groundwater.bottom + groundwater.storage ./ ...
    (groundwater.specific_yield .* mesh.cell_area(:));

remaining = duration_s; substeps = 0; rejected_steps=0; maximum_flux = 0; river_exchange = 0;
river_exchange_by_cell=zeros(mesh.n_cells,1);
fixed_head_exchange=0;
candidate=duration_s;
if isfield(groundwater,'suggested_dt_s') && isfinite(groundwater.suggested_dt_s) && ...
        groundwater.suggested_dt_s>0
    candidate=min(candidate,groundwater.suggested_dt_s);
end
while remaining > 10*eps(max(duration_s,1))
    if isfield(config,'groundwater_lateral_flow') && ~config.groundwater_lateral_flow
        stable=remaining;
    else
        stable = groundwater_stable_timestep(mesh,groundwater,config.groundwater_courant);
    end
    dt = min([remaining,stable,candidate]);
    if isfield(config,'groundwater_lateral_flow') && ~config.groundwater_lateral_flow
        lateral=struct('max_flux_m3_s',0);
    else
        [coarse,~]=Voronoi_Boussinesq_Step(mesh,groundwater.head,groundwater.bottom, ...
            groundwater.hydraulic_conductivity,groundwater.specific_yield,dt);
        [half,first]=Voronoi_Boussinesq_Step(mesh,groundwater.head,groundwater.bottom, ...
            groundwater.hydraulic_conductivity,groundwater.specific_yield,0.5*dt);
        [fine,second]=Voronoi_Boussinesq_Step(mesh,half,groundwater.bottom, ...
            groundwater.hydraulic_conductivity,groundwater.specific_yield,0.5*dt);
        scale=config.groundwater_head_tolerance_m+config.groundwater_relative_tolerance.* ...
            max(abs(fine-groundwater.bottom),abs(coarse-groundwater.bottom));
        local_error=max(abs(fine-coarse)./max(scale,realmin),[],'omitnan');
        floor_s=min(max(min(config.min_dt_s,1),eps),remaining);
        if local_error>1
            if dt<=floor_s*(1+1e-12)
                error('HydroPol2D:GroundwaterMinimumStep', ...
                    'Voronoi groundwater reached its minimum step before meeting tolerance.');
            end
            rejected_steps=rejected_steps+1;
            candidate=max(floor_s,dt*max(0.2,0.9*local_error^(-1/3)));
            continue
        end
        groundwater.head=fine;
        lateral=struct('max_flux_m3_s',max(first.max_flux_m3_s,second.max_flux_m3_s));
        if local_error==0
            growth=2;
        else
            growth=min(2,max(0.5,0.9*local_error^(-1/3)));
        end
        candidate=max(floor_s,dt*growth);
    end
    groundwater.storage = groundwater.specific_yield .* max(groundwater.head-groundwater.bottom,0) .* mesh.cell_area(:);
    [surface_volume,channel_volume,groundwater,exchange] = river_exchange_step( ...
        mesh,surface_volume,channel_volume,groundwater,dt,config);
    if isfield(config,'groundwater_fixed_head_mask') && ~isempty(config.groundwater_fixed_head_mask)
        fixed=logical(config.groundwater_fixed_head_mask(:));
        target=groundwater.specific_yield(fixed).*(double(config.groundwater_fixed_head_m(fixed))-groundwater.bottom(fixed)).*mesh.cell_area(fixed);
        fixed_head_exchange=fixed_head_exchange+sum(max(target,0)-groundwater.storage(fixed));
        groundwater.storage(fixed)=max(target,0);
        groundwater.head(fixed)=double(config.groundwater_fixed_head_m(fixed));
    end
    maximum_flux=max([maximum_flux,lateral.max_flux_m3_s,exchange.maximum_flux_m3_s]);
    river_exchange=river_exchange+exchange.net_groundwater_to_river_m3;
    river_exchange_by_cell=river_exchange_by_cell+exchange.by_cell_m3;
    remaining=remaining-dt; substeps=substeps+1;
end

% Saturated storage above the polygon terrain becomes seepage/saturation
% excess. The vadose step will subsequently reduce its available capacity.
excess = groundwater.specific_yield .* max(groundwater.head-mesh.surface_bed(:),0) .* mesh.cell_area(:);
groundwater.storage = groundwater.storage-excess;
groundwater.head = min(groundwater.head,mesh.surface_bed(:));
surface_volume = surface_volume+excess;
groundwater.last_river_exchange_rate_m_s=river_exchange_by_cell./mesh.cell_area(:)./duration_s;
groundwater.last_seepage_rate_m_s=excess./mesh.cell_area(:)./duration_s;
groundwater.cumulative_river_exchange_m=groundwater.cumulative_river_exchange_m+ ...
    river_exchange_by_cell./mesh.cell_area(:);
groundwater.cumulative_seepage_m=groundwater.cumulative_seepage_m+ ...
    excess./mesh.cell_area(:);
diagnostics = struct('max_flux_m3_s',maximum_flux,'substep_count',substeps, ...
    'river_exchange_m3',river_exchange,'seepage_volume_m3',sum(excess), ...
    'fixed_head_exchange_m3',fixed_head_exchange,'mass_change_m3',0, ...
    'rejected_steps',rejected_steps,'suggested_dt_s',candidate);
end

function dt = groundwater_stable_timestep(mesh,groundwater,courant)
thickness=max(groundwater.head-groundwater.bottom,0);
internal=mesh.edge_neighbor>0;
o=mesh.edge_owner(internal); d=mesh.edge_neighbor(internal);
face_transmissivity=0.5.*(thickness(o)+thickness(d)).* ...
    0.5.*(groundwater.hydraulic_conductivity(o)+groundwater.hydraulic_conductivity(d));
conductance=face_transmissivity.*mesh.edge_length(internal)./mesh.edge_distance(internal);
total=accumarray(o,conductance,[mesh.n_cells 1],@sum,0)+accumarray(d,conductance,[mesh.n_cells 1],@sum,0);
active=total>0;
if any(active)
    dt=courant.*min(groundwater.specific_yield(active).*mesh.cell_area(active)./total(active));
else
    dt=inf;
end
end

function [surface_volume,channel_volume,groundwater,diagnostics] = river_exchange_step( ...
    mesh,surface_volume,channel_volume,groundwater,dt,config)
thickness=max(double(config.riverbed_thickness_m),eps);
net=0; maximum=0;
enabled=isfield(config,'river_groundwater_exchange') && logical(config.river_groundwater_exchange);
if ~enabled
    diagnostics=struct('net_groundwater_to_river_m3',0,'maximum_flux_m3_s',0, ...
        'by_cell_m3',zeros(mesh.n_cells,1));
    return
end
assert(isfield(mesh,'cell_riverbed_area') && all(isfinite(mesh.cell_riverbed_area)), ...
    'HydroPol2D:MissingRiverbedArea', ...
    'River exchange requires cell_riverbed_area_m2 in the mesh.');
assert(isfield(mesh,'cell_riverbed_hydraulic_conductivity') && ...
    all(isfinite(mesh.cell_riverbed_hydraulic_conductivity)) && ...
    all(mesh.cell_riverbed_hydraulic_conductivity>=0), ...
    'HydroPol2D:MissingRiverbedConductivity', ...
    ['River exchange requires nonnegative ' ...
     'cell_riverbed_hydraulic_conductivity_m_s in the mesh.']);
kbed=mesh.cell_riverbed_hydraulic_conductivity(:);
riverbed_area=mesh.cell_riverbed_area(:);
by_cell=zeros(mesh.n_cells,1);

resolved=riverbed_area>0;
if mesh.channel.n_nodes>0, resolved(mesh.channel.host_cell(:))=false; end
cells=find(resolved);
if ~isempty(cells)
    stage_all=mesh.surface_bed(:)+surface_volume./mesh.surface_area(:);
    if isfield(config,'subgrid_tables') && ~isempty(fieldnames(config.subgrid_tables))
        [subgrid_stage,~]=hp2d_voronoi_subgrid_cell_stage( ...
            config.subgrid_tables,surface_volume);
        selected=logical(config.subgrid_tables.cell_is_subgrid(:));
        stage_all(selected)=subgrid_stage(selected);
    end
    stage=stage_all(cells);
    [volume,q]=hp2d_river_groundwater_exchange_step( ...
        kbed(cells),thickness,riverbed_area(cells),groundwater.head(cells),stage, ...
        groundwater.storage(cells),surface_volume(cells), ...
        groundwater.specific_yield(cells).*mesh.cell_area(cells), ...
        mesh.surface_area(cells),dt,true(size(cells)));
    groundwater.storage(cells)=groundwater.storage(cells)-volume;
    surface_volume(cells)=surface_volume(cells)+volume;
    by_cell(cells)=by_cell(cells)+volume;
    net=net+sum(volume); maximum=max(maximum,max(abs(q),[],'omitnan'));
end

if mesh.channel.n_nodes>0
    cells=mesh.channel.host_cell(:);
    stage=mesh.channel.bed(:)+channel_volume./mesh.channel.plan_area(:);
    [volume,q]=hp2d_river_groundwater_exchange_step( ...
        kbed(cells),thickness,riverbed_area(cells),groundwater.head(cells),stage, ...
        groundwater.storage(cells),channel_volume, ...
        groundwater.specific_yield(cells).*mesh.cell_area(cells), ...
        mesh.channel.plan_area(:),dt,true(size(cells)));
    groundwater.storage(cells)=groundwater.storage(cells)-volume;
    channel_volume=channel_volume+volume;
    by_cell=by_cell+accumarray(cells,volume,[mesh.n_cells 1],@sum,0);
    net=net+sum(volume); maximum=max(maximum,max(abs(q),[],'omitnan'));
end
groundwater.head=groundwater.bottom+groundwater.storage./(groundwater.specific_yield.*mesh.cell_area(:));
diagnostics=struct('net_groundwater_to_river_m3',net,'maximum_flux_m3_s',maximum,'by_cell_m3',by_cell);
end

function diagnostics = empty_diagnostics()
diagnostics=struct('max_flux_m3_s',0,'substep_count',0,'river_exchange_m3',0, ...
    'seepage_volume_m3',0,'fixed_head_exchange_m3',0,'mass_change_m3',0, ...
    'rejected_steps',0,'suggested_dt_s',inf);
end
