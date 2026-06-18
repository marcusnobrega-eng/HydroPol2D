function zs = hp2d_sfincs_zs_from_cell_volume(SubgridTables, V)
%HP2D_SFINCS_ZS_FROM_CELL_VOLUME Recover water level from SFINCS volume table.

zmin = SubgridTables.z_zmin;
zmax = SubgridTables.z_zmax;
Vmax = SubgridTables.z_volmax;
cell_area = SubgridTables.cell_area;
z_levels = SubgridTables.z_level;
V_levels = SubgridTables.z_volume;

V = max(V, 0);
zs = zmin;
valid = isfinite(V) & isfinite(zmin) & isfinite(zmax);

above = valid & V >= Vmax;
zs(above) = zmax(above) + (V(above) - Vmax(above)) ./ cell_area;

inside = valid & ~above & V > 0;
if any(inside(:))
    tmp = lookup_cell_table_by_volume(V_levels, z_levels, V, inside);
    zs(inside) = tmp(inside);
end

zs(~valid) = NaN;
end

function out = lookup_cell_table_by_volume(V_levels, z_levels, V, mask)
[ny, nx, nz] = size(V_levels);
out = NaN(ny, nx, 'like', V);
if ~any(mask(:))
    return;
end

assigned = false(ny, nx);
for k = 1:(nz-1)
    v1 = V_levels(:, :, k);
    v2 = V_levels(:, :, k+1);
    z1 = z_levels(:, :, k);
    z2 = z_levels(:, :, k+1);
    interval = mask & ~assigned & isfinite(v1) & isfinite(v2) & ...
        isfinite(z1) & isfinite(z2) & v2 > v1 & V >= v1 & V <= v2;
    if any(interval(:))
        w = (V(interval) - v1(interval)) ./ (v2(interval) - v1(interval));
        out(interval) = z1(interval) + w .* (z2(interval) - z1(interval));
        assigned(interval) = true;
    end
end

remaining = mask & ~assigned;
if any(remaining(:))
    v_first = V_levels(:, :, 1);
    z_first = z_levels(:, :, 1);
    v_last = V_levels(:, :, end);
    z_last = z_levels(:, :, end);

    below = remaining & isfinite(v_first) & isfinite(z_first) & V <= v_first;
    out(below) = z_first(below);
    assigned(below) = true;

    above = remaining & ~assigned & isfinite(v_last) & isfinite(z_last) & V >= v_last;
    out(above) = z_last(above);
    assigned(above) = true;

    flat = remaining & ~assigned & isfinite(z_last);
    out(flat) = z_last(flat);
end
end
