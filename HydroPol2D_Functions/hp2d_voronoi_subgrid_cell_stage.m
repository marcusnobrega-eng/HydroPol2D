function [stage, wet_area] = hp2d_voronoi_subgrid_cell_stage(tables, volume_m3)
%HP2D_VORONOI_SUBGRID_CELL_STAGE Water surface elevation from the cell volume curve.
%
%   Exact inverse of the piecewise-linear elevation-volume curve: the wetted area
%   is piecewise constant in elevation, so volume is piecewise linear and the
%   inversion is analytic, not iterative.
%
%   FULLY VECTORISED. An earlier version looped over cells, which is called every
%   timestep: on a 12,330-cell Pune mesh that is 12,330 interpreted iterations per
%   step, and MATLAB pays for every one. The bracketing index is found instead by
%   comparing the whole padded volume table against the volume vector at once.
%
%   Above the top breakpoint the cell is fully wet, so the plan area is the correct
%   slope; a zero slope in the first slot would divide by nothing, so the plan area
%   is used there too.

volume = max(double(volume_m3(:)), 0);
n = numel(volume);
assert(n == numel(tables.cell_datum_m), 'HydroPol2D:SubgridStageSize', ...
    'Volume vector length must match the sub-grid cell count.');

V = tables.cell_volume_m3;          % n x width, padded by repeating the last value
Z = tables.cell_zeta_m;
A = tables.cell_wet_area_m2;
width = size(V, 2);
count = max(tables.cell_point_count, 1);

% Mask the padding so it never wins the comparison, then count how many
% breakpoints the volume clears. That count IS the bracketing index.
valid = (1:width) <= count;                       % n x width logical
j = max(sum((volume >= V) & valid, 2), 1);
rows = (1:n).';
lin = rows + (j - 1) * n;                         % column-major linear index

base_v = V(lin);
base_z = Z(lin);
slope = A(lin);
top = j >= count;
plan = tables.cell_plan_area_m2;
slope(top | slope <= 0) = plan(top | slope <= 0);

stage = tables.cell_datum_m + base_z + (volume - base_v) ./ max(slope, realmin);

wet_area = A(lin);
dry = wet_area <= 0;
wet_area(dry) = plan(dry);
end
