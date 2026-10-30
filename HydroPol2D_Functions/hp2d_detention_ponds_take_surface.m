function [Detention_Ponds, d_t] = hp2d_detention_ponds_take_surface(Detention_Ponds, d_t)
%HP2D_DETENTION_PONDS_TAKE_SURFACE Move surface water on pond cells into the ponds [mm].
%   A pond cell is not surface: whatever the vertical step left on it (groundwater
%   exfiltration, saturation excess, an inflow cell) joins the pond's level pool.
%   For exfiltration this is add_groundwater_exfiltration_to_ponds in
%   HydroPol2D-Python: groundwater that reaches the surface under a pond fills the
%   pond first.

P = Detention_Ponds;
depth_mm = double(gather(d_t(P.cells)));
depth_mm(~isfinite(depth_mm)) = 0;
volume = accumarray(P.cell_pond, max(depth_mm, 0) / 1000 * P.cell_area_m2, [P.n_ponds 1]);
P.storage_m3 = P.storage_m3 + volume;
P.cumulative_surface_in_m3 = P.cumulative_surface_in_m3 + volume;
d_t(P.cells) = 0;
Detention_Ponds = P;
end
