function face = hp2d_subgrid_face_state(SubgridTables, eta_face, direction, Resolution)
%HP2D_SUBGRID_FACE_STATE Query shared-face lookup hydraulics.
%
% direction is 'x' or 'y'. Returned arrays match the corresponding internal
% face table shape, not the padded HydroPol2D face-array shape.

direction = lower(char(direction));
switch direction
    case 'x'
        area_tab = SubgridTables.area_x;
        width_tab = get_table_or_default(SubgridTables, 'width_x', []);
        wetfrac_tab = get_table_or_default(SubgridTables, 'wetfrac_x', []);
        hrep_tab = get_table_or_default(SubgridTables, 'hrep_x', []);
        phi_tab = get_table_or_default(SubgridTables, 'phi_x', []);
        perimeter_tab = get_table_or_default(SubgridTables, 'perimeter_x', []);
        Rh_tab = get_table_or_default(SubgridTables, 'Rh_x', []);
        n_tab = get_table_or_default(SubgridTables, 'nrep_x', SubgridTables.n_x);
        K_tab = get_table_or_default(SubgridTables, 'K_x', []);
        invert = SubgridTables.invert_x;
    case 'y'
        area_tab = SubgridTables.area_y;
        width_tab = get_table_or_default(SubgridTables, 'width_y', []);
        wetfrac_tab = get_table_or_default(SubgridTables, 'wetfrac_y', []);
        hrep_tab = get_table_or_default(SubgridTables, 'hrep_y', []);
        phi_tab = get_table_or_default(SubgridTables, 'phi_y', []);
        perimeter_tab = get_table_or_default(SubgridTables, 'perimeter_y', []);
        Rh_tab = get_table_or_default(SubgridTables, 'Rh_y', []);
        n_tab = get_table_or_default(SubgridTables, 'nrep_y', SubgridTables.n_y);
        K_tab = get_table_or_default(SubgridTables, 'K_y', []);
        invert = SubgridTables.invert_y;
    otherwise
        error('hp2d_subgrid_face_state:badDirection', ...
            'direction must be x or y.');
end

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
    perimeter = area ./ max(area ./ max(width, eps_like(area)), eps_like(area));
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
    K = (1 ./ n_eff) .* area .* max(Rh, 0).^(2/3);
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
