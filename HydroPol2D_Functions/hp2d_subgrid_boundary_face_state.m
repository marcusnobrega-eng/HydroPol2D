function face = hp2d_subgrid_boundary_face_state(SubgridTables, eta_face, side, Resolution)
%HP2D_SUBGRID_BOUNDARY_FACE_STATE Query boundary-face lookup hydraulics.
%
% side is north, south, west, or east. Returned arrays match the coarse
% cell grid, so callers can select outlet cells after querying.

side = lower(char(side));
valid = {'north','south','west','east'};
if ~any(strcmp(side, valid))
    error('hp2d_subgrid_boundary_face_state:badSide', ...
        'side must be north, south, west, or east.');
end

area_tab = SubgridTables.(['area_' side]);
width_tab = get_table_or_default(SubgridTables, ['width_' side], []);
wetfrac_tab = get_table_or_default(SubgridTables, ['wetfrac_' side], []);
hrep_tab = get_table_or_default(SubgridTables, ['hrep_' side], []);
phi_tab = get_table_or_default(SubgridTables, ['phi_' side], []);
perimeter_tab = get_table_or_default(SubgridTables, ['perimeter_' side], []);
Rh_tab = get_table_or_default(SubgridTables, ['Rh_' side], []);
n_tab = get_table_or_default(SubgridTables, ['nrep_' side], ...
    SubgridTables.(['n_' side]));
K_tab = get_table_or_default(SubgridTables, ['K_' side], []);
invert = SubgridTables.(['invert_' side]);

area = hp2d_subgrid_lookup_eta(area_tab, eta_face, invert, ...
    SubgridTables.dz, SubgridTables.maxDepth);
n_eff = hp2d_subgrid_lookup_eta(n_tab, eta_face, invert, ...
    SubgridTables.dz, SubgridTables.maxDepth);

if isempty(width_tab)
    width = Resolution .* double(area > 0);
else
    width = hp2d_subgrid_lookup_eta(width_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

if isempty(wetfrac_tab)
    wetfrac = width ./ max(Resolution, eps_like(width));
else
    wetfrac = hp2d_subgrid_lookup_eta(wetfrac_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

if isempty(phi_tab)
    phi = wetfrac;
else
    phi = hp2d_subgrid_lookup_eta(phi_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

if isempty(hrep_tab)
    hrep = area ./ max(width, eps_like(area));
else
    hrep = hp2d_subgrid_lookup_eta(hrep_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

if isempty(perimeter_tab)
    perimeter = area ./ max(hrep, eps_like(area));
    perimeter(~isfinite(perimeter)) = 0;
else
    perimeter = hp2d_subgrid_lookup_eta(perimeter_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

if isempty(Rh_tab)
    Rh = area ./ max(perimeter, eps_like(area));
else
    Rh = hp2d_subgrid_lookup_eta(Rh_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

if isempty(K_tab)
    K = (1 ./ n_eff) .* width .* max(hrep, 0).^(5/3);
else
    K = hp2d_subgrid_lookup_eta(K_tab, eta_face, invert, ...
        SubgridTables.dz, SubgridTables.maxDepth);
end

area(~isfinite(area)) = 0;
width(~isfinite(width)) = 0;
wetfrac(~isfinite(wetfrac)) = 0;
hrep(~isfinite(hrep)) = 0;
phi(~isfinite(phi)) = 0;
perimeter(~isfinite(perimeter)) = 0;
Rh(~isfinite(Rh)) = 0;
n_eff(~isfinite(n_eff) | n_eff <= 0) = NaN;
K(~isfinite(K)) = 0;

face = struct();
face.area = max(area, 0);
face.width = max(width, 0);
face.wetfrac = max(min(wetfrac, 1), 0);
face.phi = max(min(phi, 1), 0);
face.hrep = max(hrep, 0);
face.perimeter = max(perimeter, 0);
face.Rh = max(Rh, 0);
face.n = n_eff;
face.K = max(K, 0);
face.Hgrid = face.area ./ max(Resolution, eps_like(face.area));
face.H = face.hrep;
face.invert = invert;
end

function value = get_table_or_default(S, field, default_value)
if isfield(S, field) && ~isempty(S.(field))
    value = S.(field);
else
    value = default_value;
end
end

function e = eps_like(x)
if isa(x, 'gpuArray')
    e = eps(classUnderlying(x));
else
    e = eps(class(x));
end
end
