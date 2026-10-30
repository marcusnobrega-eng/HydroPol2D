function head_m = hp2d_urban_drainage_head(N, volume_m3, initial_head_m)
%HP2D_URBAN_DRAINAGE_HEAD Node head [m] that stores a volume [m3]
%   (node_head_m in HydroPol2D-Python). Without a first guess: 40 bisection
%   steps between the invert and the crown plus the volume over the shaft
%   area plus 1 m. With a guess: 6 safeguarded Newton steps. Boundary nodes
%   keep their boundary head.

tiny = 1e-12;
volume = max(volume_m3(:), 0);
lower = N.node_invert_m;
upper = N.node_crown_head_m + volume ./ max(N.node_shaft_area_m2, tiny) + 1;
if nargin < 3 || isempty(initial_head_m)
    low = lower;
    high = upper;
    for iteration = 1:40
        middle = 0.5 * (low + high);
        below = hp2d_urban_drainage_storage(N, middle) < volume;
        low(below) = middle(below);
        high(~below) = middle(~below);
    end
    solved = 0.5 * (low + high);
else
    low = lower;
    high = upper;
    current = min(max(initial_head_m(:), lower), upper);
    for iteration = 1:6
        [stored, area] = hp2d_urban_drainage_storage(N, current);
        below = stored < volume;
        low(below) = current(below);
        high(~below) = current(~below);
        candidate = current + (volume - stored) ./ max(area, tiny);
        inside = isfinite(candidate) & candidate >= low & candidate <= high;
        current(inside) = candidate(inside);
        current(~inside) = 0.5 * (low(~inside) + high(~inside));
    end
    solved = current;
end
head_m = solved;
head_m(~N.dynamic) = N.node_boundary_head_m(~N.dynamic);
end
