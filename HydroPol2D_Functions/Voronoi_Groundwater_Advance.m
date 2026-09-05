function [surface_volume, channel_volume, groundwater, diagnostics] = Voronoi_Groundwater_Advance( ...
    mesh, surface_volume, channel_volume, groundwater, duration_s, config)
%VORONOI_GROUNDWATER_ADVANCE Apply recharge then subcycle polygon Boussinesq flow.

if ~groundwater.enabled || duration_s <= 0
    diagnostics = empty_diagnostics();
    return
end

pending = groundwater.pending_exchange_m3;
assert(all(groundwater.storage + pending >= -1e-8), ...
    'HydroPol2D:GroundwaterMass', 'Vadose exchange requested more groundwater than is stored.');
groundwater.storage = max(groundwater.storage + pending,0);
groundwater.pending_exchange_m3(:) = 0;
groundwater.head = groundwater.bottom + groundwater.storage ./ ...
    (groundwater.specific_yield .* mesh.cell_area(:));

remaining = duration_s; substeps = 0; maximum_flux = 0; river_exchange = 0;
river_exchange_by_cell=zeros(mesh.n_cells,1);
while remaining > 10*eps(max(duration_s,1))
    stable = groundwater_stable_timestep(mesh,groundwater,config.courant);
    dt = min(remaining,stable);
    [groundwater.head,lateral] = Voronoi_Boussinesq_Step(mesh,groundwater.head,groundwater.bottom, ...
        groundwater.hydraulic_conductivity,groundwater.specific_yield,dt);
    groundwater.storage = groundwater.specific_yield .* max(groundwater.head-groundwater.bottom,0) .* mesh.cell_area(:);
    [surface_volume,channel_volume,groundwater,exchange] = river_exchange_step( ...
        mesh,surface_volume,channel_volume,groundwater,dt,config);
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
diagnostics = struct('max_flux_m3_s',maximum_flux,'substep_count',substeps, ...
    'river_exchange_m3',river_exchange,'seepage_volume_m3',sum(excess), ...
    'mass_change_m3',0);
end

function dt = groundwater_stable_timestep(mesh,groundwater,courant)
thickness=max(groundwater.head-groundwater.bottom,0);
transmissivity=groundwater.hydraulic_conductivity.*thickness;
internal=mesh.edge_neighbor>0;
o=mesh.edge_owner(internal); d=mesh.edge_neighbor(internal);
conductance=2.*transmissivity(o).*transmissivity(d)./max(transmissivity(o)+transmissivity(d),eps) ...
    .*mesh.edge_length(internal)./mesh.edge_distance(internal);
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
kbed=double(config.riverbed_hydraulic_conductivity_m_s);
thickness=max(double(config.riverbed_thickness_m),eps);
net=0; maximum=0;
if kbed<=0
    diagnostics=struct('net_groundwater_to_river_m3',0,'maximum_flux_m3_s',0, ...
        'by_cell_m3',zeros(mesh.n_cells,1));
    return
end
by_cell=zeros(mesh.n_cells,1);

if isfield(config,'resolved_river_mask') && ~isempty(config.resolved_river_mask)
    resolved=logical(config.resolved_river_mask(:));
else
    resolved=false(mesh.n_cells,1);
end
cells=find(resolved);
if ~isempty(cells)
    stage=mesh.surface_bed(cells)+surface_volume(cells)./mesh.surface_area(cells);
    q=kbed.*mesh.cell_area(cells)./thickness.*(groundwater.head(cells)-stage);
    [q,volume]=limit_exchange(q,dt,groundwater.storage(cells),surface_volume(cells));
    groundwater.storage(cells)=groundwater.storage(cells)-volume;
    surface_volume(cells)=surface_volume(cells)+volume;
    by_cell(cells)=by_cell(cells)+volume;
    net=net+sum(volume); maximum=max(maximum,max(abs(q),[],'omitnan'));
end

if mesh.channel.n_nodes>0
    cells=mesh.channel.host_cell(:);
    stage=mesh.channel.bed(:)+channel_volume./mesh.channel.plan_area(:);
    q=kbed.*mesh.channel.plan_area(:)./thickness.*(groundwater.head(cells)-stage);
    [q,volume]=limit_exchange(q,dt,groundwater.storage(cells),channel_volume);
    groundwater.storage(cells)=groundwater.storage(cells)-volume;
    channel_volume=channel_volume+volume;
    by_cell=by_cell+accumarray(cells,volume,[mesh.n_cells 1],@sum,0);
    net=net+sum(volume); maximum=max(maximum,max(abs(q),[],'omitnan'));
end
groundwater.head=groundwater.bottom+groundwater.storage./(groundwater.specific_yield.*mesh.cell_area(:));
diagnostics=struct('net_groundwater_to_river_m3',net,'maximum_flux_m3_s',maximum,'by_cell_m3',by_cell);
end

function [q,volume] = limit_exchange(q,dt,groundwater_volume,surface_or_channel_volume)
volume=q.*dt;
to_river=volume>0; volume(to_river)=min(volume(to_river),groundwater_volume(to_river));
to_groundwater=volume<0; volume(to_groundwater)=-min(-volume(to_groundwater),surface_or_channel_volume(to_groundwater));
q=volume./dt;
end

function diagnostics = empty_diagnostics()
diagnostics=struct('max_flux_m3_s',0,'substep_count',0,'river_exchange_m3',0,'seepage_volume_m3',0,'mass_change_m3',0);
end
