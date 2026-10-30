function [area_m2, conveyance, depth_m] = hp2d_d8_link_hydraulics(Faces, Links, stage_m)
%HP2D_D8_LINK_HYDRAULICS Flow area, conveyance and depth of each D8 link.
%   The link water level is the higher of its two cell stages. Its depth is
%   that level minus the link datum (the higher channel invert); area and
%   conveyance K = A R^(2/3) / n interpolate the link table linearly between
%   breakpoints and hold the last value above the top one. A dry link
%   (depth <= 0) has zero area and conveyance. This is
%   _river_table_hydraulics_jax in HydroPol2D-Python.
%   Faces: datum, area, conveyance, count (one row per link) and either
%   zeta (one breakpoint row per link) or shared_zeta (one strictly
%   increasing breakpoint row used by every link; zeta is then []).
%   Links: source and target linear cell indices.

edge_stage = max(stage_m(Links.source), stage_m(Links.target));
zeta = edge_stage - Faces.datum;
count = Faces.count;
slot = hp2d_d8_row_slot(Faces.zeta, zeta, count, Faces.shared_zeta);
next_slot = min(slot + 1, count);
rows = (1:numel(zeta)).';
lin0 = sub2ind(size(Faces.area), rows, slot);
lin1 = sub2ind(size(Faces.area), rows, next_slot);
if isempty(Faces.shared_zeta)
    z0 = Faces.zeta(lin0);
    z1 = Faces.zeta(lin1);
else
    z0 = reshape(Faces.shared_zeta(slot), [], 1);
    z1 = reshape(Faces.shared_zeta(next_slot), [], 1);
end
weight = min(max((zeta - z0) ./ max(z1 - z0, 1e-12), 0), 1);
weight(slot >= count) = 0;

area_m2 = Faces.area(lin0) + weight .* (Faces.area(lin1) - Faces.area(lin0));
conveyance = Faces.conveyance(lin0) + weight .* (Faces.conveyance(lin1) - Faces.conveyance(lin0));
depth_m = max(zeta, 0);
wet = depth_m > 0;
area_m2 = max(area_m2, 0);
conveyance = max(conveyance, 0);
area_m2(~wet) = 0;
conveyance(~wet) = 0;

end
