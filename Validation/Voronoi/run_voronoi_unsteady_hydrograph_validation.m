function summary = run_voronoi_unsteady_hydrograph_validation(mesh_file)
%RUN_VORONOI_UNSTEADY_HYDROGRAPH_VALIDATION Time-varying channel inflow volume.

mesh = HydroPol2D_Read_UGRID(mesh_file); channel = mesh.channel;
inlet = setdiff((1:channel.n_nodes)',channel.link_down);
assert(numel(inlet) == 1);
config = struct('duration_s',100,'max_dt_s',2,'min_dt_s',1e-4, ...
    'forcing_interval_s',50,'output_interval_s',10);
forcing.channel_boundary = @(time_s,state,active_mesh) pulse(time_s,inlet); %#ok<NASGU,INUSD>
results = HydroPol2D_Voronoi_Run(mesh_file,config,forcing);
expected = 5 * 50;
actual = sum(results.final_surface_volume_m3) + sum(results.final_channel_volume_m3);
recorded = sum([results.diagnostics.boundary_net_inflow_volume_m3]);
volume_error = abs(actual-expected)/expected;
record_error = abs(recorded-expected)/expected;
assert(volume_error <= 1e-10 && record_error <= 1e-10);
summary = struct('prescribed_volume_m3',expected,'relative_volume_error',volume_error, ...
    'relative_record_error',record_error);
end

function boundary = pulse(time_s,inlet)
boundary = struct('node_id',inlet,'type',"inflow",'value',5*double(time_s < 50));
end
