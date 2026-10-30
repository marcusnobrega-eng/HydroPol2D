function [area, radius, top_width] = hp2d_urban_drainage_section(depth_m, diameter_m)
%HP2D_URBAN_DRAINAGE_SECTION Wetted area [m2], hydraulic radius [m] and top
%   width [m] of a circular pipe at a depth (circular_section in
%   HydroPol2D-Python). The depth is clipped to [0, diameter].

D = diameter_m + zeros(size(depth_m));
depth = min(max(depth_m, 0), D);
theta = 2 * acos(min(max(1 - 2 * depth ./ D, -1), 1));
area = D.^2 / 8 .* (theta - sin(theta));
perimeter = 0.5 * D .* theta;
radius = zeros(size(area));
wet = perimeter > 0;
radius(wet) = area(wet) ./ perimeter(wet);
top_width = zeros(size(area));
open = depth > 0 & depth < D;
top_width(open) = 2 * sqrt(max(depth(open) .* (D(open) - depth(open)), 0));
end
