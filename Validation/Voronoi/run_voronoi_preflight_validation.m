function summary = run_voronoi_preflight_validation(mesh_file)
%RUN_VORONOI_PREFLIGHT_VALIDATION Resource and unsupported-backend gates.

config = struct('routing_solver','local_inertial','compute_backend','cpu', ...
    'maximum_adjacent_size_ratio',2,'available_memory_bytes',1e9, ...
    'allow_legacy_mesh',true);
[~,report] = HydroPol2D_Voronoi_Preflight(mesh_file,config);

gpu_failed = false;
try
    bad = config; bad.compute_backend = 'gpu';
    HydroPol2D_Voronoi_Preflight(mesh_file,bad);
catch error
    gpu_failed = strcmp(error.identifier,'HydroPol2D:VoronoiGPUNotValidated');
end
assert(gpu_failed,'Unvalidated Voronoi GPU backend was accepted.');

memory_failed = false;
try
    bad = config; bad.available_memory_bytes = 1;
    HydroPol2D_Voronoi_Preflight(mesh_file,bad);
catch error
    memory_failed = strcmp(error.identifier,'HydroPol2D:VoronoiMemoryGate');
end
assert(memory_failed,'Memory gate did not reject an unaffordable mesh.');
summary = report;
end
