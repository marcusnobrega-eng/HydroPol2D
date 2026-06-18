function [drep_new,V_old,V_new,A_new,residual_m3] = hp2d_subgrid_apply_volume_change(drep_old, delta_depth_mm, SubgridTables, cell_area)
%HP2D_SUBGRID_APPLY_VOLUME_CHANGE Apply areal depth change through storage table.
%
% delta_depth_mm is interpreted as an areal water-depth change over the
% coarse model cell area, so the volume increment is
% delta_depth_mm/1000 * cell_area.

if isscalar(cell_area)
    area = cell_area .* ones(size(drep_old), 'like', drep_old);
else
    area = cell_area;
end

delta_volume = delta_depth_mm ./ 1000 .* area;

if isfield(SubgridTables, 'sfincs_exact') && SubgridTables.sfincs_exact
    zs_old = SubgridTables.z_zmin + max(drep_old, 0);
    V_old = hp2d_sfincs_cell_volume_from_zs(SubgridTables, zs_old);
    V_new = max(V_old + delta_volume, 0);
    zs_new = hp2d_sfincs_zs_from_cell_volume(SubgridTables, V_new);
    drep_new = max(zs_new - SubgridTables.z_zmin, 0);
    A_new = hp2d_sfincs_cell_wet_area_from_zs(SubgridTables, zs_new);
else
    V_old = hp2d_subgrid_lookup_depth( ...
        SubgridTables.volume_cell, drep_old, SubgridTables.dz, SubgridTables.maxDepth);

    V_new = max(V_old + delta_volume, 0);

    drep_new = hp2d_subgrid_inverse_volume( ...
        SubgridTables.volume_cell, V_new, SubgridTables.dz, SubgridTables.maxDepth, sqrt(max(area(:))));

    A_new = hp2d_subgrid_lookup_depth( ...
        SubgridTables.area_cell, drep_new, SubgridTables.dz, SubgridTables.maxDepth);
end

residual_m3 = V_new - V_old - delta_volume;
residual_m3(V_old + delta_volume < 0) = V_new(V_old + delta_volume < 0);
end
