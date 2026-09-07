function [area, perimeter, conveyance] = hp2d_voronoi_subgrid_face_state(tables, stage_m)
%HP2D_VORONOI_SUBGRID_FACE_STATE Flow area, wetted perimeter and conveyance at a stage.
%
%   Area and conveyance interpolate linearly inside the bracketing interval.
%   Perimeter does NOT: a sub-cell panel of the face is either wet or dry, so the
%   wetted perimeter is a STEP function of elevation. Interpolating it linearly
%   returned 12 m for a 20 m wide flat channel in the Python implementation, which
%   is why it is held piecewise constant here.
%
%   FULLY VECTORISED, for the same reason as the cell lookup: this runs every
%   timestep over every face -- 28,055 on the Pune mesh.

stage = double(stage_m(:));
n = numel(stage);
assert(n == numel(tables.face_datum_m), 'HydroPol2D:SubgridFaceSize', ...
    'Stage vector length must match the sub-grid face count.');

Z = tables.face_zeta_m;
Aa = tables.face_flow_area_m2;
Pp = tables.face_perimeter_m;
Cc = tables.face_conveyance;
width = size(Z, 2);
count = max(tables.face_point_count, 1);

zeta = stage - tables.face_datum_m;
valid = (1:width) <= count;
j = max(sum((zeta >= Z) & valid, 2), 1);
rows = (1:n).';
lin = rows + (j - 1) * n;
j1 = min(j + 1, count);
lin1 = rows + (j1 - 1) * n;

z0 = Z(lin);
z1 = Z(lin1);
span = max(z1 - z0, realmin);
w = min(max((zeta - z0) ./ span, 0), 1);
w(j >= count) = 0;                       % at or above the top, hold the last value

area = Aa(lin) + w .* (Aa(lin1) - Aa(lin));
conveyance = Cc(lin) + w .* (Cc(lin1) - Cc(lin));
perimeter = Pp(lin);                     % piecewise CONSTANT, see above

dry = zeta <= 0;
area(dry) = 0; perimeter(dry) = 0; conveyance(dry) = 0;
area = max(area, 0); perimeter = max(perimeter, 0); conveyance = max(conveyance, 0);
end
