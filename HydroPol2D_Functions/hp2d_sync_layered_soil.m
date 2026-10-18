function Soil_Properties = hp2d_sync_layered_soil(Soil_Properties, Subsurface, head_m, idx_nan)
%HP2D_SYNC_LAYERED_SOIL Mirror slots 1-3 into Soil_Properties.Layers and I_t (readers, maps).
ny = Subsurface.ny; nx = Subsurface.nx;
V = hp2d_soil_view(Subsurface.column, head_m(:), Subsurface.water_mm, zeros(ny*nx, 1));
L = Soil_Properties.Layers;
L.near_surface_storage_mm = reshape(Subsurface.water_mm(:, 1), ny, nx);
L.root_zone_storage_mm = reshape(Subsurface.water_mm(:, 2), ny, nx);
L.transmission_storage_mm = reshape(Subsurface.water_mm(:, 3), ny, nx);
L.near_surface_thickness_m = reshape(V.near_surface_thickness_m, ny, nx);
L.root_zone_thickness_m = reshape(V.root_zone_thickness_m, ny, nx);
L.transmission_thickness_m = reshape(V.transmission_thickness_m, ny, nx);
L.near_surface_capacity_mm = reshape(V.capacity_mm(:, 1), ny, nx);
L.root_zone_capacity_mm = reshape(V.capacity_mm(:, 2), ny, nx);
L.transmission_capacity_mm = reshape(V.capacity_mm(:, 3), ny, nx);
L.total_vadose_capacity_mm = reshape(V.total_capacity_mm, ny, nx);
L.water_table_depth_m = reshape(V.water_table_depth_m, ny, nx);
L.vadose_depth_m = reshape(V.vadose_depth_m, ny, nx);
names = fieldnames(L);
for i = 1:numel(names)
    value = L.(names{i});
    if isnumeric(value) && isequal(size(value), [ny, nx]) && isfloat(value)
        value(idx_nan) = nan;
        L.(names{i}) = value;
    end
end
Soil_Properties.Layers = L;
Soil_Properties.I_t = L.near_surface_storage_mm + L.root_zone_storage_mm + L.transmission_storage_mm;
Soil_Properties.I_t(idx_nan) = nan;
end
