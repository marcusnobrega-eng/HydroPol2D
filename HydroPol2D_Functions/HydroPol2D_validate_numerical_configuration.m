function HydroPol2D_validate_numerical_configuration(flags)
%HYDROPOL2D_VALIDATE_NUMERICAL_CONFIGURATION Reject unsafe solver mixtures.

if ~isstruct(flags)
    error('HydroPol2D:InvalidFlags', 'flags must be a structure.');
end
if isfield(flags, 'flag_D8') && gather_scalar(flags.flag_D8) == 1
    error('HydroPol2D:LegacyD8Retired', [ ...
        'The historical raster D8 implementation is incomplete and has been retired. ' ...
        'Use the validated raster D4 model or the separate Voronoi finite-volume runner.']);
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
