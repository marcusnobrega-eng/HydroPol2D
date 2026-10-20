function sections = hp2d_gauge_cross_sections(rows, cols, widths, directions, elevation, valid)
%HP2D_GAUGE_CROSS_SECTIONS Line of cells across the flow, centred on each gauge cell.
%   rows, cols: gauge cell indices. widths: cells in the line (NaN or < 1 -> 5; an even
%   width becomes the next odd one so the gauge stays centred). directions: 'auto' or
%   the downstream direction 'N', 'S', 'E' or 'W'. Auto takes the downhill direction of
%   the terrain averaged over a width x width window around the gauge (rows grow
%   southward). An east- or west-flowing gauge gets a north-south line; a north- or
%   south-flowing one an east-west line. Cells outside the grid or the domain are dropped.
%   The discharge is the signed sum over the line's downstream faces
%   (hp2d_gauge_section_discharge), the same rule as HydroPol2D-Python routing gauges
%   (signed D4 faces).
[ny, nx] = size(elevation);
n = numel(rows);
if nargin < 3 || isempty(widths), widths = nan(n, 1); end
if nargin < 4 || isempty(directions), directions = repmat({'auto'}, n, 1); end
directions = cellstr(string(directions));
sections = repmat(struct('rows', [], 'cols', [], 'downstream', 'E', 'width', 5, ...
    'auto', true), n, 1);
for i = 1:n
    w = double(widths(min(i, numel(widths))));
    if ~isfinite(w) || w < 1
        w = 5;
    end
    w = round(w);
    if mod(w, 2) == 0
        warning('HydroPol2D:GaugeSection', ['Gauge %d: cross-section width %d is even; ', ...
            'using %d so the gauge cell stays centred.'], i, w, w + 1);
        w = w + 1;
    end
    k = (w - 1)/2;
    r = rows(i);
    c = cols(i);
    d = upper(strtrim(directions{min(i, numel(directions))}));
    is_auto = isempty(d) || strcmp(d, 'AUTO') || strcmp(d, 'NAN') || strcmp(d, 'MISSING');
    if is_auto
        d = downhill_direction(elevation, valid, r, c, max(k, 1));
        if isempty(d)
            warning('HydroPol2D:GaugeSection', ['Gauge %d: the terrain around it is flat; ', ...
                'set its downstream direction. Using east.'], i);
            d = 'E';
        end
    end
    assert(any(strcmp(d, {'N', 'S', 'E', 'W'})), 'HydroPol2D:GaugeSection', ...
        'Gauge %d: downstream direction must be auto, N, S, E or W (got %s).', i, d);
    if any(strcmp(d, {'E', 'W'}))
        line_rows = (r - k:r + k)';
        line_cols = c + zeros(w, 1);
    else
        line_rows = r + zeros(w, 1);
        line_cols = (c - k:c + k)';
    end
    inside = line_rows >= 1 & line_rows <= ny & line_cols >= 1 & line_cols <= nx;
    line_rows = line_rows(inside);
    line_cols = line_cols(inside);
    active = valid(sub2ind([ny, nx], line_rows, line_cols));
    sections(i).rows = line_rows(active);
    sections(i).cols = line_cols(active);
    sections(i).downstream = d;
    sections(i).width = w;
    sections(i).auto = is_auto;
end
end

function d = downhill_direction(z, valid, r, c, k)
% Mean central-difference gradient over the window; '' where it is flat.
[ny, nx] = size(z);
z(~valid) = nan;
gx = []; gy = [];
for rr = max(r - k, 1):min(r + k, ny)
    for cc = max(c - k, 1):min(c + k, nx)
        if cc > 1 && cc < nx && isfinite(z(rr, cc - 1)) && isfinite(z(rr, cc + 1))
            gx(end + 1) = (z(rr, cc + 1) - z(rr, cc - 1))/2; %#ok<AGROW>
        end
        if rr > 1 && rr < ny && isfinite(z(rr - 1, cc)) && isfinite(z(rr + 1, cc))
            gy(end + 1) = (z(rr + 1, cc) - z(rr - 1, cc))/2; %#ok<AGROW>
        end
    end
end
gx = mean([gx 0]); gy = mean([gy 0]);   % dz per cell, x east, y south
if abs(gx) < 1e-12 && abs(gy) < 1e-12
    d = '';
elseif abs(gx) >= abs(gy)
    if gx < 0, d = 'E'; else, d = 'W'; end
else
    if gy < 0, d = 'S'; else, d = 'N'; end
end
end
