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
summary = struct('d8_rejected', d8_failed, 'mixed_flags_rejected', mixed_failed);
end
