function [volume_m3, storage_area_m2] = hp2d_urban_drainage_storage(N, head_m)
%HP2D_URBAN_DRAINAGE_STORAGE Water stored at each node [m3] and its plan
%   storage area [m2] at a given head (node_storage_m3 and
%   _node_storage_area_m2 in HydroPol2D-Python).
%   Storage is the shaft area times the depth above the invert plus half of
%   each incident pipe's wetted volume at its end depth. Boundary (outfall)
%   nodes store nothing; their storage area is still returned.

n = numel(N.node_ids);
head = head_m(:);
[start_area, ~, start_width] = hp2d_urban_drainage_section( ...
    head(N.pipe_start) - N.pipe_start_invert_m, N.pipe_diameter_m);
[end_area, ~, end_width] = hp2d_urban_drainage_section( ...
    head(N.pipe_end) - N.pipe_end_invert_m, N.pipe_diameter_m);
volume_m3 = N.node_shaft_area_m2 .* max(head - N.node_invert_m, 0) ...
    + accumarray(N.pipe_start, N.pipe_half_length_m .* start_area, [n 1]) ...
    + accumarray(N.pipe_end, N.pipe_half_length_m .* end_area, [n 1]);
volume_m3(~N.dynamic) = 0;
if nargout > 1
    storage_area_m2 = N.node_shaft_area_m2 ...
        + accumarray(N.pipe_start, N.pipe_half_length_m .* start_width, [n 1]) ...
        + accumarray(N.pipe_end, N.pipe_half_length_m .* end_width, [n 1]);
end
end
