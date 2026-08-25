function summary = run_voronoi_full_momentum_validation(mesh_file, output_directory)
%RUN_VORONOI_FULL_MOMENTUM_VALIDATION CPU acceptance checks for HLL UGRID.

arguments
    mesh_file (1,:) char
    output_directory (1,:) char
end
root=fullfile(fileparts(mfilename('fullpath')),'..','..'); addpath(fullfile(root,'HydroPol2D_Functions'));
if exist(output_directory,'dir')~=7, mkdir(output_directory); end
mesh=HydroPol2D_Read_UGRID(mesh_file);
assert(mesh.channel.n_nodes==0,'Use a fully resolved mesh for full-momentum validation.');

% Lake at rest over the V-shaped, non-flat terrain: validates hydrostatic
% reconstruction and reflective-wall pressure fluxes.
eta0=min(mesh.surface_bed)+0.2; h0=max(eta0-mesh.surface_bed,0);
lake_config=base_config(h0,300,1,fullfile(output_directory,'full-momentum-lake.nc'));
lake=HydroPol2D_Voronoi_Run(mesh_file,lake_config,struct());
lake_error=max(abs(lake.final_surface_volume_m3./mesh.surface_area-h0));
lake_velocity=max(lake.surface_velocity_m_s,[],'all');
assert(lake_error<=1e-10 && lake_velocity<=1e-10, ...
    'HydroPol2D:FullMomentumLakeAtRest','Lake-at-rest balance failed.');

% Rainfall-runoff smoke test with the exact 60 m V-tilted outlet.
outlet=find(mesh.edge_neighbor==0 & abs(mesh.edge_midpoint_y)<1e-8 & ...
    mesh.edge_midpoint_x>=380 & mesh.edge_midpoint_x<=440);
assert(abs(sum(mesh.edge_length(outlet))-60)<=1e-8,'V-tilted outlet must be 60 m wide.');
run_config=base_config(0,3600,2,fullfile(output_directory,'full-momentum-rainfall.nc'));
forcing=struct('surface_source_m_s',@(time_s,~,~) 10.8e-3/3600*double(time_s<1800), ...
    'surface_boundary',struct('edge_id',outlet,'type',"normal_flow",'value',0.004));
run=HydroPol2D_Voronoi_Run(mesh_file,run_config,forcing);
d=run.diagnostics; input_volume=10.8e-3/3600*1800*sum(mesh.cell_area);
mass_error=max(abs([d.step_mass_residual_m3]))/input_volume;
assert(all(isfinite(run.final_surface_volume_m3)) && all(run.final_surface_volume_m3>=0));
assert(all(isfinite(run.surface_velocity_m_s),'all') && all(run.surface_velocity_m_s>=0,'all'));
assert(mass_error<=1e-8,'HydroPol2D:FullMomentumMass','Rainfall-runoff mass balance failed.');

% The same test at a 1 s maximum timestep is the temporal reference for
% this first explicit implementation.  Stability alone is not acceptance.
reference_config=base_config(0,3600,1,fullfile(output_directory,'full-momentum-rainfall-fixed1.nc'));
reference=HydroPol2D_Voronoi_Run(mesh_file,reference_config,forcing);
depth_convergence=norm(run.surface_depth_m(:,end)-reference.surface_depth_m(:,end)) / ...
    max(norm(reference.surface_depth_m(:,end)),eps);
outlet_volume=-sum([d.boundary_net_inflow_volume_m3]);
reference_outlet_volume=-sum([reference.diagnostics.boundary_net_inflow_volume_m3]);
outlet_convergence=abs(outlet_volume-reference_outlet_volume)/max(reference_outlet_volume,eps);
assert(depth_convergence<=0.02 && outlet_convergence<=0.02, ...
    'HydroPol2D:FullMomentumTemporalConvergence','2 s full-momentum solution differs excessively from 1 s reference.');

summary=table(lake_error,lake_velocity,mass_error,depth_convergence,outlet_convergence, ...
    numel(d),min([d.dt_s]),max([d.dt_s]), ...
    max([d.max_surface_depth_m]),max([d.max_surface_velocity_m_s]), ...
    outlet_volume, ...
    'VariableNames',{'lake_depth_error_m','lake_velocity_m_s','rainfall_mass_error_fraction', ...
    'final_depth_error_vs_1s','outlet_volume_error_vs_1s', ...
    'rainfall_timestep_count','rainfall_minimum_dt_s','rainfall_maximum_dt_s', ...
    'rainfall_maximum_depth_m','rainfall_maximum_velocity_m_s','rainfall_outlet_volume_m3'});
writetable(summary,fullfile(output_directory,'full-momentum-validation-summary.csv'));
end

function config=base_config(initial_depth,duration_s,maximum_dt_s,output_file)
config=struct('routing_solver','full_momentum','duration_s',duration_s, ...
    'initial_surface_depth_m',initial_depth,'min_dt_s',0.01,'max_dt_s',maximum_dt_s, ...
    'output_interval_s',60,'forcing_interval_s',inf,'courant',0.4, ...
    'surface_roughness',0.04,'compute_backend','cpu','output_netcdf',output_file, ...
    'overwrite_output',true,'full_momentum_maximum_velocity_m_s',10);
end
