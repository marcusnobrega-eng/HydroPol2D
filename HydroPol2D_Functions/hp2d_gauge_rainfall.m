function out = hp2d_gauge_rainfall(action, varargin)
%HP2D_GAUGE_RAINFALL Rain-gauge series spread by inverse-distance weighting.
%   W = hp2d_gauge_rainfall('weights', gauges, x, y) sets up the weights of
%   points (x, y) once. gauges has x, y [m, same CRS as the points], time_s
%   (row times, each the start of its interval), intensity_mm_h (one row per
%   time, one column per gauge; NaN = not reporting) and an optional power
%   (default 2).
%   rate = hp2d_gauge_rainfall('rate', W, time_s) is the rain at the points
%   [mm/h] at time_s.
%   As in HydroPol2D-Python (forcing/gauge_rainfall.py): normalised
%   inverse-distance weights over the gauges that report, a point on a
%   reporting gauge takes that gauge's value, and the rain stops one interval
%   after the last row. Before the first row the rate is 0, as for every
%   start-stamped series in MATLAB (update_spatial_BC).

switch action
    case 'weights'
        [gauges, x, y] = varargin{:};
        power = 2;
        if isfield(gauges, 'power') && ~isempty(gauges.power), power = double(gauges.power); end
        gx = double(gauges.x(:))'; gy = double(gauges.y(:))';
        distance = hypot(double(x(:)) - gx, double(y(:)) - gy);
        exact = distance <= 100*eps(max(abs([gx gy x(:)' y(:)']), [], 'all') + 1);
        inverse = distance.^(-power);
        inverse(exact) = 0;
        times = double(gauges.time_s(:));
        if numel(times) >= 2
            until = times(end) + (times(end) - times(end-1));
        else
            until = inf;
        end
        out = struct('inverse', inverse, 'exact', double(exact), 'times', times, ...
            'until', until, 'intensity', double(gauges.intensity_mm_h));
    case 'rate'
        [W, time_s] = varargin{:};
        n = size(W.inverse, 1);
        index = find(W.times <= time_s, 1, 'last');
        if isempty(index) || time_s >= W.until
            out = zeros(n, 1);
            return
        end
        values = W.intensity(index, :)';
        reporting = isfinite(values);
        out = zeros(n, 1);
        if ~any(reporting), return; end
        v = values(reporting);
        exact = W.exact(:, reporting);
        on_gauge = any(exact > 0, 2);
        out(on_gauge) = (exact(on_gauge, :) * v) ./ sum(exact(on_gauge, :), 2);
        weight = W.inverse(~on_gauge, reporting);
        out(~on_gauge) = (weight * v) ./ max(sum(weight, 2), realmin);
        out = max(out, 0);
    otherwise
        error('hp2d_gauge_rainfall: unknown action %s.', action);
end
end
