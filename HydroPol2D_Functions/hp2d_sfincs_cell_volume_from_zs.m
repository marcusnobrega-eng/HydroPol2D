function V = hp2d_sfincs_cell_volume_from_zs(SubgridTables, zs)
%HP2D_SFINCS_CELL_VOLUME_FROM_ZS Exact cell storage from water level.

zmin = SubgridTables.z_zmin;
zmax = SubgridTables.z_zmax;
Vmax = SubgridTables.z_volmax;
cell_area = SubgridTables.cell_area;
z_levels = SubgridTables.z_level;
V_levels = SubgridTables.z_volume;

V = zeros(size(zs), 'like', zs);
valid = isfinite(zs) & isfinite(zmin) & isfinite(zmax);

below = valid & zs <= zmin;
V(below) = 0;

above = valid & zs >= zmax;
V(above) = Vmax(above) + (zs(above) - zmax(above)) .* cell_area;

inside = valid & ~below & ~above;
if any(inside(:))
    tmp = lookup_cell_table_by_z(z_levels, V_levels, zs, inside);
    V(inside) = tmp(inside);
end

V(~isfinite(V)) = 0;
V = max(V, 0);
end

function out = lookup_cell_table_by_z(z_levels, value_levels, zs, mask)
[ny, nx, nz] = size(z_levels);
out = zeros(ny, nx, 'like', zs);
if ~any(mask(:))
    return;
end

assigned = false(ny, nx);
for k = 1:(nz-1)
    z1 = z_levels(:, :, k);
    z2 = z_levels(:, :, k+1);
    v1 = value_levels(:, :, k);
    v2 = value_levels(:, :, k+1);
    interval = mask & ~assigned & isfinite(z1) & isfinite(z2) & ...
        isfinite(v1) & isfinite(v2) & z2 > z1 & zs >= z1 & zs <= z2;
    if any(interval(:))
        w = (zs(interval) - z1(interval)) ./ (z2(interval) - z1(interval));
        out(interval) = v1(interval) + w .* (v2(interval) - v1(interval));
        assigned(interval) = true;
    end
end

remaining = mask & ~assigned;
if any(remaining(:))
    z_first = z_levels(:, :, 1);
    v_first = value_levels(:, :, 1);
    z_last = z_levels(:, :, end);
    v_last = value_levels(:, :, end);

    below = remaining & isfinite(z_first) & isfinite(v_first) & zs <= z_first;
    out(below) = v_first(below);
    assigned(below) = true;

    above = remaining & ~assigned & isfinite(z_last) & isfinite(v_last) & zs >= z_last;
    out(above) = v_last(above);
    assigned(above) = true;

    flat = remaining & ~assigned & isfinite(v_last);
    out(flat) = v_last(flat);
end
end
