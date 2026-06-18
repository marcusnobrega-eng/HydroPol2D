function face = hp2d_sfincs_boundary_velocity_state(SubgridTables, zu, side)
%HP2D_SFINCS_BOUNDARY_VELOCITY_STATE Query outlet/open-boundary face tables.

side = lower(char(side));
allowed = {'north', 'south', 'west', 'east'};
if ~ismember(side, allowed)
    error('hp2d_sfincs_boundary_velocity_state:badSide', ...
        'side must be north, south, west, or east.');
end

face = query_table_set(SubgridTables, side, zu);
end

function face = query_table_set(S, prefix, zu)
zmin = S.([prefix '_zmin']);
zmax = S.([prefix '_zmax']);
havg_tab = S.([prefix '_havg']);
nrep_tab = S.([prefix '_nrep']);
pwet_tab = S.([prefix '_pwet']);
navg = S.([prefix '_navg']);
ffit = S.([prefix '_ffit']);

havg = zeros(size(zu), 'like', zu);
nrep = Inf(size(zu), 'like', zu);
pwet = zeros(size(zu), 'like', zu);

valid = isfinite(zu) & isfinite(zmin) & isfinite(zmax);
inside = valid & zu >= zmin & zu < zmax;
above = valid & zu >= zmax;

if any(inside(:))
    tmp = lookup_levels(zmin, zmax, havg_tab, zu, inside);
    havg(inside) = tmp(inside);
    tmp = lookup_levels(zmin, zmax, nrep_tab, zu, inside);
    nrep(inside) = tmp(inside);
    tmp = lookup_levels(zmin, zmax, pwet_tab, zu, inside);
    pwet(inside) = tmp(inside);
end

if any(above(:))
    hM = havg_tab(:, :, end);
    nM = nrep_tab(:, :, end);
    dz_above = zu(above) - zmax(above);
    havg(above) = hM(above) + dz_above;
    nM_above = nM(above);
    navg_above = navg(above);
    bad_top_n = ~isfinite(nM_above) | nM_above <= 0;
    nM_above(bad_top_n) = navg_above(bad_top_n);
    denom = ffit(above) .* dz_above + 1;
    nrep_above = navg_above - (navg_above - nM_above) ./ denom;
    bad_n = ~isfinite(nrep_above) | nrep_above <= 0;
    nrep_above(bad_n) = navg_above(bad_n);
    nrep(above) = nrep_above;
    pwet(above) = 1;
end

dry = valid & zu < zmin;
havg(dry) = 0;
pwet(dry) = 0;
nrep(dry) = Inf;

havg(~isfinite(havg)) = 0;
pwet(~isfinite(pwet)) = 0;
nrep(~isfinite(nrep) | nrep <= 0) = Inf;

face.HG = max(havg, 0);
face.n = nrep;
face.phi = max(min(pwet, 1), 0);
face.zmin = zmin;
face.zmax = zmax;
end

function out = lookup_levels(zmin, zmax, table, zu, mask)
[ny, nx, nz] = size(table);
out = zeros(ny, nx, 'like', zu);
if ~any(mask(:))
    return;
end

flat = zmax <= zmin;
flat_mask = mask & flat;
if any(flat_mask(:))
    top = table(:, :, end);
    out(flat_mask) = top(flat_mask);
end

interp_mask = mask & ~flat;
if any(interp_mask(:))
    pos = (zu - zmin) ./ (zmax - zmin);
    pos = max(min(pos, 1), 0);
    fidx = pos .* (nz - 1) + 1;
    i1 = floor(fidx);
    i1 = max(min(i1, nz - 1), 1);
    i2 = i1 + 1;
    w = fidx - i1;
    [I, J] = ndgrid(1:ny, 1:nx);
    ind1 = sub2ind([ny, nx, nz], I(interp_mask), J(interp_mask), i1(interp_mask));
    ind2 = sub2ind([ny, nx, nz], I(interp_mask), J(interp_mask), i2(interp_mask));
    v1 = table(ind1);
    v2 = table(ind2);
    out(interp_mask) = v1 + w(interp_mask) .* (v2 - v1);
end
end
