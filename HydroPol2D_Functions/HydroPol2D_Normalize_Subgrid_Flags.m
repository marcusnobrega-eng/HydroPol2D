function flags = HydroPol2D_Normalize_Subgrid_Flags(flags)
%HYDROPOL2D_NORMALIZE_SUBGRID_FLAGS Map public subgrid choices to legacy internals.

if ~isstruct(flags)
    error('HydroPol2D:InvalidFlags','flags must be a structure.');
end

public_names = {'flag_neal_channel','flag_structured_lookup_subgrid','flag_voronoi_subgrid'};
public_present = isfield(flags,public_names);
for k = 1:numel(public_names)
    name = public_names{k};
    if ~isfield(flags,name)
        flags.(name) = 0;
    else
        flags.(name) = binary_flag(flags.(name),name);
    end
end

if ~isfield(flags,'flag_voronoi'), flags.flag_voronoi = 0; end
flags.flag_voronoi = binary_flag(flags.flag_voronoi,'flag_voronoi');

legacy_subgrid_present = isfield(flags,'flag_subgrid');
legacy_overbanks_present = isfield(flags,'flag_overbanks');
legacy_subgrid = 0;
legacy_overbanks = 0;
if legacy_subgrid_present
    legacy_subgrid = binary_flag(flags.flag_subgrid,'flag_subgrid');
end
if legacy_overbanks_present
    legacy_overbanks = binary_flag(flags.flag_overbanks,'flag_overbanks');
end

if ~any(public_present) && legacy_subgrid == 1
    if legacy_overbanks == 1
        flags.flag_neal_channel = 1;
        migrated_name = 'flag_neal_channel';
    else
        flags.flag_structured_lookup_subgrid = 1;
        migrated_name = 'flag_structured_lookup_subgrid';
    end
    warning('HydroPol2D:LegacySubgridFlags', ...
        ['flag_subgrid/flag_overbanks are legacy inputs. This case was interpreted as %s=1. ' ...
         'Save the case with the explicit public flag.'], migrated_name);
end

active = [flags.flag_neal_channel, flags.flag_structured_lookup_subgrid, ...
    flags.flag_voronoi_subgrid];
if sum(active) > 1
    error('HydroPol2D:ConflictingSubgridMethods', ...
        'Select at most one subgrid method: Neal channel, structured lookup, or Voronoi subgrid.');
end
if flags.flag_voronoi == 1 && any(active(1:2))
    error('HydroPol2D:SubgridTopologyMismatch', ...
        'Neal channel and structured lookup subgrid require the structured raster model.');
end
if flags.flag_voronoi_subgrid == 1 && flags.flag_voronoi ~= 1
    error('HydroPol2D:VoronoiSubgridRequiresVoronoi', ...
        'flag_voronoi_subgrid=1 requires flag_voronoi=1.');
end

expected_subgrid = double(any(active(1:2)));
expected_overbanks = double(flags.flag_neal_channel == 1);
legacy_conflict = ...
    (legacy_subgrid_present && legacy_subgrid ~= expected_subgrid) || ...
    (legacy_overbanks_present && legacy_overbanks ~= expected_overbanks);
if any(public_present) && legacy_conflict
    error('HydroPol2D:ConflictingLegacySubgridFlags', ...
        'Legacy flag_subgrid/flag_overbanks conflict with the explicit public subgrid method.');
end

% The validated raster kernels retain these internal names.
flags.flag_subgrid = expected_subgrid;
flags.flag_overbanks = expected_overbanks;
end

function value = binary_flag(value,name)
if isa(value,'gpuArray'), value = gather(value); end
if islogical(value), value = double(value); end
if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ~ismember(double(value),[0 1])
    error('HydroPol2D:InvalidFlagValue','%s must be 0 or 1.',name);
end
value = double(value);
end
