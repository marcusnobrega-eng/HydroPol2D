function Summary = run_voronoi_groundwater_river_exchange_validation(mesh_file)
%RUN_VORONOI_GROUNDWATER_RIVER_EXCHANGE_VALIDATION Gaining/losing Neal river.

mesh=HydroPol2D_Read_UGRID(mesh_file);
assert(mesh.channel.n_nodes>0,'Validation mesh must contain Neal channel nodes.');
config=struct('courant',0.2,'riverbed_hydraulic_conductivity_m_s',1e-6, ...
    'riverbed_thickness_m',0.5,'resolved_river_mask',false(mesh.n_cells,1));
duration=600;

% Gaining reach: groundwater stage is above an initially dry channel, but
% remains below the land surface so the transfer is river exchange, not seepage.
groundwater=make_groundwater(mesh,[11;10;9;6],-10);
surface=zeros(mesh.n_cells,1); channel=zeros(mesh.channel.n_nodes,1);
mass_before=sum(groundwater.storage)+sum(surface)+sum(channel);
[surface,channel,groundwater,gaining]=Voronoi_Groundwater_Advance( ...
    mesh,surface,channel,groundwater,duration,config);
gaining_residual=sum(groundwater.storage)+sum(surface)+sum(channel)-mass_before;

% Losing reach: a bankfull channel leaks into a lower aquifer.
groundwater=make_groundwater(mesh,-3,-10);
surface=zeros(mesh.n_cells,1);
channel=mesh.channel.plan_area.*max(mesh.channel.bank-mesh.channel.bed,0);
mass_before=sum(groundwater.storage)+sum(surface)+sum(channel);
[surface,channel,groundwater,losing]=Voronoi_Groundwater_Advance( ...
    mesh,surface,channel,groundwater,duration,config);
losing_residual=sum(groundwater.storage)+sum(surface)+sum(channel)-mass_before;

passed=gaining.river_exchange_m3>0 && losing.river_exchange_m3<0 && ...
    abs(gaining_residual)<1e-8 && abs(losing_residual)<1e-8;
Summary=table(gaining.river_exchange_m3,losing.river_exchange_m3, ...
    gaining_residual,losing_residual,passed, ...
    'VariableNames',{'gaining_exchange_m3','losing_exchange_m3', ...
    'gaining_mass_residual_m3','losing_mass_residual_m3','passed'});
disp(Summary);
assert(passed,'HydroPol2D:VoronoiGroundwaterRiverValidation', ...
    'Groundwater-river exchange validation failed.');
end

function groundwater = make_groundwater(mesh,head,bottom)
if isscalar(head), head=repmat(head,mesh.n_cells,1); else, head=head(:); end
groundwater=struct('enabled',true,'head',head, ...
    'bottom',repmat(bottom,mesh.n_cells,1),'hydraulic_conductivity',zeros(mesh.n_cells,1), ...
    'specific_yield',0.2*ones(mesh.n_cells,1),'pending_exchange_m3',zeros(mesh.n_cells,1));
groundwater.storage=groundwater.specific_yield.*(groundwater.head-groundwater.bottom).*mesh.cell_area;
end
