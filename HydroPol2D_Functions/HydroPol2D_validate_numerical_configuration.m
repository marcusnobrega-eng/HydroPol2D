function HydroPol2D_validate_numerical_configuration(flags)
%HYDROPOL2D_VALIDATE_NUMERICAL_CONFIGURATION Reject unsafe solver mixtures.

if ~isstruct(flags)
    error('HydroPol2D:InvalidFlags', 'flags must be a structure.');
end
flags = HydroPol2D_Normalize_Subgrid_Flags(flags);
if isfield(flags, 'flag_D8') && gather_scalar(flags.flag_D8) == 1
    error('HydroPol2D:LegacyD8Retired', [ ...
        'The historical raster D8 implementation is incomplete and has been retired. ' ...
        'Use the validated raster D4 model or the separate Voronoi finite-volume runner.']);
end
voronoi = 0;
if isfield(flags, 'flag_voronoi')
    voronoi = gather_scalar(flags.flag_voronoi);
    if ~ismember(voronoi, [0 1])
        error('HydroPol2D:InvalidVoronoiFlag', 'flag_voronoi must be 0 or 1.');
    end
end
if voronoi == 1
    if ~isfield(flags, 'flag_inertial') || gather_scalar(flags.flag_inertial) ~= 1
        error('HydroPol2D:VoronoiRequiresLocalInertial', ...
            'Voronoi routing currently requires flag_inertial=1.');
    end
    if isfield(flags, 'flag_GPU') && gather_scalar(flags.flag_GPU) == 1
        error('HydroPol2D:VoronoiGPUNotValidated', ...
            'Voronoi routing is currently validated on CPU only; set flag_GPU=0.');
    end
    error('HydroPol2D:VoronoiStandaloneRunnerRequired', [ ...
        'flag_voronoi=1 cannot use the raster preprocessing loop. Generate the UGRID mesh ' ...
        'with HydroBathyDEM and run HydroPol2D_Voronoi_Run.']);
end

names = {'flag_CA','flag_inertial','flag_diffusive','flag_kinematic','flag_full_momentum'};
active = false(size(names));
for k = 1:numel(names)
    if isfield(flags, names{k})
        active(k) = gather_scalar(flags.(names{k})) == 1;
    end
end
if sum(active) > 1
    error('HydroPol2D:ConflictingRoutingSolvers', ...
        'Select exactly one routing solver; active flags: %s.', strjoin(names(active), ', '));
end

% D8 local-inertial routing (D8_Local_Inertial_Model) is a local-inertial
% option. Like HydroPol2D-Python it runs with rainfall and the outlets; the
% modules below have no D8 coupling yet and are refused rather than run
% without it.
if flag_value(flags, 'flag_d8_local_inertial') == 1
    if flag_value(flags, 'flag_inertial') ~= 1
        error('HydroPol2D:D8RequiresLocalInertial', ...
            'flag_d8_local_inertial = 1 requires flag_inertial = 1.');
    end
    if flag_value(flags, 'flag_numerical_scheme') > 1
        error('HydroPol2D:D8RequiresBates', ...
            'flag_d8_local_inertial = 1 requires flag_numerical_scheme = 1 (Bates).');
    end
    unsupported = {'flag_subgrid', 'flag_waterquality', 'flag_human_instability', ...
        'flag_snow_modeling', 'flag_inflow', 'flag_stage_hydrograph', 'flag_obs_gauges', ...
        'flag_river_groundwater_exchange', 'flag_reduce_DEM', 'flag_urban_drainage', ...
        'flag_detention_ponds', 'flag_green_infrastructure'};
    enabled = false(size(unsupported));
    for k = 1:numel(unsupported)
        enabled(k) = flag_value(flags, unsupported{k}) ~= 0;
    end
    if any(enabled)
        error('HydroPol2D:D8UnsupportedModule', ...
            'flag_d8_local_inertial = 1 does not support: %s.', strjoin(unsupported(enabled), ', '));
    end
end
end

function value = flag_value(flags, name)
if isfield(flags, name)
    value = gather_scalar(flags.(name));
else
    value = 0;
end
end

function value = gather_scalar(value)
if isa(value, 'gpuArray')
    value = gather(value);
end
if ~isscalar(value) || ~isnumeric(value)
    error('HydroPol2D:InvalidFlagValue', 'Numerical flags must be numeric scalars.');
end
value = double(value);
end
