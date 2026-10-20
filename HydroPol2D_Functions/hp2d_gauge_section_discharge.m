function Q = hp2d_gauge_section_discharge(sections, q, cell_area_m2)
%HP2D_GAUGE_SECTION_DISCHARGE Signed discharge across each gauge cross-section (m3/s).
%   sections: from hp2d_gauge_cross_sections. q: the per-direction cell outflows of the
%   routing step in mm/h over the cell area (fields left, right, up, down and, with D8
%   routing, ne, se, sw, nw; rows grow southward). The discharge crosses the line's
%   downstream faces. Positive is downstream.
%   Two conventions, set by q.signed_faces (default false):
%   - one-sided outflows >= 0 (CA routing, D8 inertial): what the line cells send
%     downstream minus what the cells beyond send back;
%   - signed face fluxes (D4 inertial, diffusive, kinematic, full momentum): each face is
%     stored twice (qout_right(r,c) = -qout_left(r,c+1)), so the line cells' own
%     downstream value is already the net flux and nothing is subtracted.
names = {'left', 'right', 'up', 'down', 'ne', 'se', 'sw', 'nw'};
for k = 1:numel(names)
    if ~isfield(q, names{k}) || isempty(q.(names{k}))
        q.(names{k}) = [];
    end
end
signed = isfield(q, 'signed_faces') && q.signed_faces;
[ny, nx] = size(q.right);
Q = zeros(numel(sections), 1);
to_m3_s = cell_area_m2/1000/3600;
for i = 1:numel(sections)
    s = sections(i);
    if isempty(s.rows)
        Q(i) = nan;
        continue
    end
    switch s.downstream
        case 'E'   % across the east faces
            out = pick(q.right, s.rows, s.cols) + pick(q.ne, s.rows, s.cols) + pick(q.se, s.rows, s.cols);
            back = pick(q.left, s.rows, s.cols + 1) + pick(q.nw, s.rows, s.cols + 1) + pick(q.sw, s.rows, s.cols + 1);
        case 'W'   % across the west faces
            out = pick(q.left, s.rows, s.cols) + pick(q.nw, s.rows, s.cols) + pick(q.sw, s.rows, s.cols);
            back = pick(q.right, s.rows, s.cols - 1) + pick(q.ne, s.rows, s.cols - 1) + pick(q.se, s.rows, s.cols - 1);
        case 'S'   % across the south faces
            out = pick(q.down, s.rows, s.cols) + pick(q.se, s.rows, s.cols) + pick(q.sw, s.rows, s.cols);
            back = pick(q.up, s.rows + 1, s.cols) + pick(q.ne, s.rows + 1, s.cols) + pick(q.nw, s.rows + 1, s.cols);
        otherwise  % 'N', across the north faces
            out = pick(q.up, s.rows, s.cols) + pick(q.ne, s.rows, s.cols) + pick(q.nw, s.rows, s.cols);
            back = pick(q.down, s.rows - 1, s.cols) + pick(q.se, s.rows - 1, s.cols) + pick(q.sw, s.rows - 1, s.cols);
    end
    if signed
        back = 0;
    end
    Q(i) = (out - back)*to_m3_s;
end

    function v = pick(field, r, c)
        % Sum of field over the given cells; 0 for missing fields or cells off the grid.
        if isempty(field)
            v = 0;
            return
        end
        keep = r >= 1 & r <= ny & c >= 1 & c <= nx;
        if ~any(keep)
            v = 0;
            return
        end
        vals = double(gather(field(sub2ind([ny, nx], r(keep), c(keep)))));
        vals(~isfinite(vals)) = 0;
        v = sum(vals);
    end
end
