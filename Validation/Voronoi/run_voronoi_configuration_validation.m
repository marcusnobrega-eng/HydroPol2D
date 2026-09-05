function summary = run_voronoi_configuration_validation()
%RUN_VORONOI_CONFIGURATION_VALIDATION Unsafe D8 and mixed flags must fail.

flags = struct('flag_D8', 1, 'flag_inertial', 1);
d8_failed = false;
try
    HydroPol2D_validate_numerical_configuration(flags);
catch ME
    d8_failed = strcmp(ME.identifier, 'HydroPol2D:LegacyD8Retired');
end
assert(d8_failed, 'Legacy D8 configuration was not rejected.');

flags = struct('flag_D8', 0, 'flag_inertial', 1, 'flag_diffusive', 1);
mixed_failed = false;
try
    HydroPol2D_validate_numerical_configuration(flags);
catch ME
    mixed_failed = strcmp(ME.identifier, 'HydroPol2D:ConflictingRoutingSolvers');
end
assert(mixed_failed, 'Conflicting routing flags were not rejected.');
HydroPol2D_validate_numerical_configuration(struct('flag_D8', 0, 'flag_inertial', 1));
defaults = HydroPol2D_Voronoi_Options(1, struct());
assert(strcmp(defaults.unresolved_river_policy, 'neal_subgrid'));
none_policy = HydroPol2D_Voronoi_Options(1, struct('unresolved_river_policy','none'));
assert(strcmp(none_policy.unresolved_river_policy, 'none'));
neal = HydroPol2D_Normalize_Subgrid_Flags(struct( ...
    'flag_voronoi',0,'flag_neal_channel',1));
assert(neal.flag_subgrid == 1 && neal.flag_overbanks == 1);
lookup = HydroPol2D_Normalize_Subgrid_Flags(struct( ...
    'flag_voronoi',0,'flag_structured_lookup_subgrid',1));
assert(lookup.flag_subgrid == 1 && lookup.flag_overbanks == 0);
voronoi_subgrid = HydroPol2D_Normalize_Subgrid_Flags(struct( ...
    'flag_voronoi',1,'flag_voronoi_subgrid',1));
assert(voronoi_subgrid.flag_subgrid == 0 && voronoi_subgrid.flag_overbanks == 0);

conflict_failed=false;
try
    HydroPol2D_Normalize_Subgrid_Flags(struct( ...
        'flag_voronoi',0,'flag_neal_channel',1,'flag_structured_lookup_subgrid',1));
catch ME
    conflict_failed=strcmp(ME.identifier,'HydroPol2D:ConflictingSubgridMethods');
end
assert(conflict_failed,'Conflicting public subgrid flags were not rejected.');

topology_failed=false;
try
    HydroPol2D_Normalize_Subgrid_Flags(struct( ...
        'flag_voronoi',0,'flag_voronoi_subgrid',1));
catch ME
    topology_failed=strcmp(ME.identifier,'HydroPol2D:VoronoiSubgridRequiresVoronoi');
end
assert(topology_failed,'Voronoi subgrid without Voronoi topology was not rejected.');

legacy=HydroPol2D_Normalize_Subgrid_Flags(struct( ...
    'flag_voronoi',0,'flag_subgrid',1,'flag_overbanks',1));
assert(legacy.flag_neal_channel == 1 && legacy.flag_subgrid == 1);
inferred=HydroPol2D_Voronoi_Options(1,struct('subgrid_table_path','tables.nc'));
assert(inferred.voronoi_subgrid_enabled);
summary = struct('d8_rejected', d8_failed, 'mixed_flags_rejected', mixed_failed, ...
    'voronoi_options_validated', true,'subgrid_flags_validated',true, ...
    'subgrid_conflict_rejected',conflict_failed, ...
    'subgrid_topology_rejected',topology_failed);
end
