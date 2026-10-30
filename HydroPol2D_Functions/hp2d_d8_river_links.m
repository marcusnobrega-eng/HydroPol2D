function [volume_m3, discharge_m3_s, cardinal_m3, diagonal_m3] = hp2d_d8_river_links( ...
    volume_m3, discharge_previous_m3_s, D8, dt_s)
%HP2D_D8_RIVER_LINKS One local-inertial update of the D8 links on cell storage.
%   Each link carries one discharge Q from its source cell to its downstream
%   (target) cell. With the link area A, conveyance K and water-surface
%   slope S = (stage_target - stage_source) / L (Bates et al. 2010, as in
%   Wflow.jl):
%       Q = (Q_prev - g A dt S) / (1 + g dt |Q_prev| A / K^2)
%   Q = 0 where the link depth is at or below D8.min_depth_m or A or K is 0;
%   then |Q| <= A sqrt(g h) (Froude number 1). Each donor cell scales its
%   outgoing links so one step cannot take more than it holds, and the
%   transfer moves Q dt from one cell to the other, so mass is conserved.
%   This is d8_river_local_inertial_step_jax in HydroPol2D-Python.
%   volume_m3 is the cell storage [m3] (column, D8.N cells); discharge in
%   m3/s, one value per link of D8.Links. cardinal_m3 and diagonal_m3 are the
%   volumes moved along cardinal and diagonal links in this step.

g = 9.81;
volume_m3 = max(volume_m3(:), 0);
source = D8.Links.source;
target = D8.Links.target;
stage = hp2d_d8_cell_stage(D8.Cells, volume_m3);
[area, conveyance, depth] = hp2d_d8_link_hydraulics(D8.Faces, D8.Links, stage);
slope = (stage(target) - stage(source)) ./ D8.Links.length;
q0 = discharge_previous_m3_s(:);
denominator = 1 + g * dt_s .* abs(q0) .* area ./ max(conveyance.^2, 1e-30);
discharge_m3_s = (q0 - g * area * dt_s .* slope) ./ denominator;
discharge_m3_s(~(depth > D8.min_depth_m & area > 0 & conveyance > 0)) = 0;
critical = area .* sqrt(g * depth);
discharge_m3_s = min(max(discharge_m3_s, -critical), critical);

n = D8.N;
outgoing = accumarray(source, max(discharge_m3_s, 0), [n 1]) + ...
    accumarray(target, max(-discharge_m3_s, 0), [n 1]);
scale = min(1, volume_m3 ./ max(outgoing * dt_s, 1e-30));
donor_scale = scale(target);
forward = discharge_m3_s >= 0;
donor_scale(forward) = scale(source(forward));
discharge_m3_s = discharge_m3_s .* donor_scale;

transfer = discharge_m3_s * dt_s;
volume_m3 = volume_m3 - accumarray(source, transfer, [n 1]) + accumarray(target, transfer, [n 1]);
volume_m3 = max(volume_m3, 0);
cardinal_m3 = sum(abs(transfer(~D8.Links.diagonal)));
diagonal_m3 = sum(abs(transfer(D8.Links.diagonal)));

end
