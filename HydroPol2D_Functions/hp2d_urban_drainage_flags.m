function flags = hp2d_urban_drainage_flags(flags, bypass_mode)
%HP2D_URBAN_DRAINAGE_FLAGS Default flag_urban_drainage to 0 and refuse the
%   combinations HydroPol2D-Python refuses (core/preprocess.py): urban
%   drainage runs with the local-inertial solver on the structured grid,
%   without the channel subgrid and without water quality. The network is
%   read in bypass mode only: HydroPol2D-Python's workbook reader has no
%   drainage inputs.

if ~isfield(flags, 'flag_urban_drainage') || isempty(flags.flag_urban_drainage)
    flags.flag_urban_drainage = 0;
end
if flags.flag_urban_drainage ~= 1
    flags.flag_urban_drainage = 0;
    return
end
if ~bypass_mode
    error('HydroPol2D:UrbanDrainage', ['flag_urban_drainage = 1 needs bypass mode ' ...
        '(InputData_Bypass.Urban_Drainage); the workbook has no drainage inputs.']);
end
refused = {'flag_CA', 'flag_full_momentum', 'flag_diffusive', 'flag_kinematic', ...
    'flag_subgrid', 'flag_waterquality', 'flag_voronoi'};
reasons = {'the CA solver', 'the full-momentum solver', 'the diffusive solver', ...
    'the kinematic solver', 'the channel subgrid', 'water quality', 'Voronoi or quadtree meshes'};
for i = 1:numel(refused)
    if isfield(flags, refused{i}) && flags.(refused{i}) == 1
        error('HydroPol2D:UrbanDrainage', ...
            'Urban drainage runs with the local-inertial solver only; it is not coupled to %s (%s = 1).', ...
            reasons{i}, refused{i});
    end
end
end
