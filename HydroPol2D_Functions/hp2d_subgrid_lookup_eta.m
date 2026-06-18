function Vq = hp2d_subgrid_lookup_eta(Vtab, eta, invert_el, dz, maxDepth)
%HP2D_SUBGRID_LOOKUP_ETA Interpolate a subgrid table by water-surface elevation.

Vq = hp2d_subgrid_lookup_depth(Vtab, eta - invert_el, dz, maxDepth);
end
