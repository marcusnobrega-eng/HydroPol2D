function summary = run_voronoi_channel_regime_validation(mesh_file)
%RUN_VORONOI_CHANNEL_REGIME_VALIDATION Storage, overbank, reversal, wet/dry.

mesh = HydroPol2D_Read_UGRID(mesh_file); channel = mesh.channel;
assert(channel.n_nodes > 2 && channel.n_links > 1);
node = ceil(channel.n_nodes/2); host = channel.host_cell(node);
bank_volume = channel.plan_area(node) * (channel.bank(node) - channel.bed(node));

surface = zeros(mesh.n_cells,1); storage = zeros(channel.n_nodes,1);
surface(host) = 0.5 * bank_volume;
[surface, storage] = Voronoi_Equilibrate_Channel_Storage(mesh,surface,storage);
assert(surface(host) == 0 && abs(storage(node)-0.5*bank_volume) <= 1e-12*bank_volume);

surface(:) = 0; storage(:) = 0; surface(host) = 1.5 * bank_volume;
initial = sum(surface) + sum(storage);
[surface, storage] = Voronoi_Equilibrate_Channel_Storage(mesh,surface,storage);
channel_stage = channel.bed(node) + storage(node)/channel.plan_area(node);
surface_stage = mesh.surface_bed(host) + surface(host)/mesh.surface_area(host);
assert(surface(host) > 0 && abs(channel_stage-surface_stage) <= 1e-12);
assert(abs(sum(surface)+sum(storage)-initial) <= 1e-12*initial);

% A downstream water-level rise must reverse the signed graph discharge.
storage = channel.plan_area .* 0.5;
last = channel.link_down(end); storage(last) = channel.plan_area(last) * 5;
before = sum(storage);
[storage,q,~] = Voronoi_Neal_Channel_Step(channel,storage,zeros(channel.n_links,1),0.01);
assert(q(end) < 0 && all(storage >= 0));
reversal_mass_error = abs(sum(storage)-before)/before;

% The draining limiter may reduce a requested flux but may not create
% negative storage or lose volume in a nearly dry graph.
storage(:) = 0; storage(channel.link_up(1)) = 1;
before = sum(storage);
[storage,~,~] = Voronoi_Neal_Channel_Step(channel,storage,100*ones(channel.n_links,1),10);
drying_mass_error = abs(sum(storage)-before)/before;
assert(all(storage >= 0) && reversal_mass_error <= 1e-12 && drying_mass_error <= 1e-12);

summary = struct('overbank_common_stage_m',channel_stage, ...
    'reversal_mass_error',reversal_mass_error,'drying_mass_error',drying_mass_error);
end
