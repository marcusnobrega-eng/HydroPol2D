function Outlet_Properties = define_outlet_faces_steepest_descent(z, idx_nan, Resolution, varargin)
%DEFINE_OUTLET_FACES_STEEPEST_DESCENT
% -------------------------------------------------------------------------
% Defines face-based outlet directions from DEM steepest descent.
%
% Purpose:
%   Convert outlet/perimeter candidate cells into outlet FACES.
%
% Main idea:
%   A cell is not simply an "outlet cell".
%   Instead, one of its exterior faces can be an outlet:
%
%       face_right
%       face_left
%       face_up
%       face_down
%
% Option A:
%   For each candidate boundary cell, activate ONLY ONE outlet face:
%   the exterior face with the strongest outward terrain descent.
%
% This avoids corner cells draining through multiple faces and prevents
% topographic divides from becoming artificial sinks.
%
% Inputs:
%   z          : DEM / bed elevation [m], ny x nx
%   idx_nan    : invalid-domain mask, ny x nx
%   Resolution : grid size [m]
%
% Optional name-value inputs:
%   'CandidateMask' : logical ny x nx mask of candidate outlet cells.
%                     If empty, all active perimeter cells are candidates.
%
%   'MinSlope'      : minimum outward terrain slope required to open a face.
%                     Default = 1e-6.
%
%   'TieTolerance'  : tolerance used when comparing slopes.
%                     Default = 1e-12.
%
%   'TiePriority'   : cell array controlling which face wins exact ties.
%                     Default = {'right','left','up','down'}.
%
% Outputs:
%   Outlet_Properties.face_right : true if right/east face is outlet
%   Outlet_Properties.face_left  : true if left/west face is outlet
%   Outlet_Properties.face_up    : true if upper/north face is outlet
%   Outlet_Properties.face_down  : true if lower/south face is outlet
%
%   Outlet_Properties.outlet_cell : true if cell has one active outlet face
%   Outlet_Properties.face_code   : 0 = no outlet
%                                   1 = right
%                                   2 = left
%                                   3 = up
%                                   4 = down
%
% Sign convention:
%   right/east face:  column + 1
%   left/west face:   column - 1
%   up/north face:    row - 1
%   down/south face:  row + 1
%
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Parse optional inputs
% -------------------------------------------------------------------------
CandidateMask = [];
MinSlope      = 1.0e-6;
TieTolerance  = 1.0e-12;
TiePriority   = {'right','left','up','down'};

if mod(numel(varargin),2) ~= 0
    error('Optional inputs must be provided as name-value pairs.');
end

for ii = 1:2:numel(varargin)
    name  = lower(varargin{ii});
    value = varargin{ii+1};

    switch name
        case {'candidatemask','candidate_mask'}
            CandidateMask = value;

        case {'minslope','min_slope'}
            MinSlope = value;

        case {'tietolerance','tie_tolerance'}
            TieTolerance = value;

        case {'tiepriority','tie_priority'}
            TiePriority = value;

        otherwise
            error('Unknown optional input: %s', varargin{ii});
    end
end

%% ------------------------------------------------------------------------
% Basic setup
% -------------------------------------------------------------------------
[ny,nx] = size(z);
dx = Resolution;

active = isfinite(z) & ~idx_nan;

if isempty(CandidateMask)
    candidate_cell = active;
else
    candidate_cell = active & logical(CandidateMask);
end

%% ------------------------------------------------------------------------
% Neighbor masks and neighbor elevations
% -------------------------------------------------------------------------
% Right neighbor: same row, column + 1
activeR = active; activeR(:) = false;
zR      = z;      zR(:)      = NaN;
if nx > 1
    activeR(:,1:end-1) = active(:,2:end);
    zR(:,1:end-1)      = z(:,2:end);
end

% Left neighbor: same row, column - 1
activeL = active; activeL(:) = false;
zL      = z;      zL(:)      = NaN;
if nx > 1
    activeL(:,2:end) = active(:,1:end-1);
    zL(:,2:end)      = z(:,1:end-1);
end

% Up neighbor: row - 1
activeU = active; activeU(:) = false;
zU      = z;      zU(:)      = NaN;
if ny > 1
    activeU(2:end,:) = active(1:end-1,:);
    zU(2:end,:)      = z(1:end-1,:);
end

% Down neighbor: row + 1
activeD = active; activeD(:) = false;
zD      = z;      zD(:)      = NaN;
if ny > 1
    activeD(1:end-1,:) = active(2:end,:);
    zD(1:end-1,:)      = z(2:end,:);
end

%% ------------------------------------------------------------------------
% Exterior candidate faces
% -------------------------------------------------------------------------
% A face is exterior if the neighbor on that face is inactive/outside.
% CandidateMask restricts which cells are allowed to become outlets.
% -------------------------------------------------------------------------
cand_right = candidate_cell & ~activeR;
cand_left  = candidate_cell & ~activeL;
cand_up    = candidate_cell & ~activeU;
cand_down  = candidate_cell & ~activeD;

candidate_exterior_cell = cand_right | cand_left | cand_up | cand_down;

%% ------------------------------------------------------------------------
% Internal descent slopes
% -------------------------------------------------------------------------
% These represent ordinary terrain descent toward active neighbor cells.
% Positive value means the neighbor is lower than the current cell.
% -------------------------------------------------------------------------
negInf = z;
negInf(:) = -Inf;

slope_int_right = negInf;
slope_int_left  = negInf;
slope_int_up    = negInf;
slope_int_down  = negInf;

idx = active & activeR;
slope_int_right(idx) = (z(idx) - zR(idx)) ./ dx;

idx = active & activeL;
slope_int_left(idx) = (z(idx) - zL(idx)) ./ dx;

idx = active & activeU;
slope_int_up(idx) = (z(idx) - zU(idx)) ./ dx;

idx = active & activeD;
slope_int_down(idx) = (z(idx) - zD(idx)) ./ dx;

max_internal_slope = max(max(slope_int_right, slope_int_left), ...
                         max(slope_int_up,    slope_int_down));

%% ------------------------------------------------------------------------
% Estimated outward descent slopes for exterior faces
% -------------------------------------------------------------------------
% Since the outside cell has no DEM value, estimate outward descent using
% the opposite interior neighbor.
%
% Example:
%   right boundary face:
%       if z(left neighbor) > z(current cell), terrain is descending toward
%       the right/outside direction.
%
% This avoids treating every perimeter cell as an outlet.
% -------------------------------------------------------------------------
slope_out_right = negInf;
slope_out_left  = negInf;
slope_out_up    = negInf;
slope_out_down  = negInf;

% Right outlet face: use left neighbor as opposite interior support.
idx = cand_right & activeL;
slope_out_right(idx) = (zL(idx) - z(idx)) ./ dx;

% Left outlet face: use right neighbor as opposite interior support.
idx = cand_left & activeR;
slope_out_left(idx) = (zR(idx) - z(idx)) ./ dx;

% Up outlet face: use down neighbor as opposite interior support.
idx = cand_up & activeD;
slope_out_up(idx) = (zD(idx) - z(idx)) ./ dx;

% Down outlet face: use up neighbor as opposite interior support.
idx = cand_down & activeU;
slope_out_down(idx) = (zU(idx) - z(idx)) ./ dx;

max_external_slope = max(max(slope_out_right, slope_out_left), ...
                         max(slope_out_up,    slope_out_down));

%% ------------------------------------------------------------------------
% Decide which cells are allowed to drain outward
% -------------------------------------------------------------------------
% A cell drains through an exterior face only if:
%
%   1. it has at least one exterior candidate face
%   2. the outward slope is larger than MinSlope
%   3. the outward slope is at least as steep as any internal descent
%
% This means topographic dividers remain walls.
% -------------------------------------------------------------------------
eligible = candidate_exterior_cell & ...
           isfinite(max_external_slope) & ...
           max_external_slope > MinSlope & ...
           max_external_slope >= (max_internal_slope - TieTolerance);

%% ------------------------------------------------------------------------
% Option A: choose only one outlet face per eligible cell
% -------------------------------------------------------------------------
face_right = active; face_right(:) = false;
face_left  = active; face_left(:)  = false;
face_up    = active; face_up(:)    = false;
face_down  = active; face_down(:)  = false;

selected = active; selected(:) = false;

for pp = 1:numel(TiePriority)

    switch lower(TiePriority{pp})

        case 'right'
            pick = eligible & ~selected & cand_right & ...
                   slope_out_right >= (max_external_slope - TieTolerance);
            face_right(pick) = true;
            selected(pick) = true;

        case 'left'
            pick = eligible & ~selected & cand_left & ...
                   slope_out_left >= (max_external_slope - TieTolerance);
            face_left(pick) = true;
            selected(pick) = true;

        case 'up'
            pick = eligible & ~selected & cand_up & ...
                   slope_out_up >= (max_external_slope - TieTolerance);
            face_up(pick) = true;
            selected(pick) = true;

        case 'down'
            pick = eligible & ~selected & cand_down & ...
                   slope_out_down >= (max_external_slope - TieTolerance);
            face_down(pick) = true;
            selected(pick) = true;

        otherwise
            error('Unknown TiePriority face: %s', TiePriority{pp});
    end
end

%% ------------------------------------------------------------------------
% Assemble output structure
% -------------------------------------------------------------------------
outlet_cell = face_right | face_left | face_up | face_down;

face_code = zeros(ny,nx,'uint8');
face_code(face_right) = uint8(1);
face_code(face_left)  = uint8(2);
face_code(face_up)    = uint8(3);
face_code(face_down)  = uint8(4);

Outlet_Properties = struct();

Outlet_Properties.method = 'perimeter_candidate_steepest_descent_option_A';

Outlet_Properties.face_right = face_right;
Outlet_Properties.face_left  = face_left;
Outlet_Properties.face_up    = face_up;
Outlet_Properties.face_down  = face_down;

Outlet_Properties.outlet_cell = outlet_cell;
Outlet_Properties.face_code   = face_code;

Outlet_Properties.active = active;
Outlet_Properties.candidate_cell = candidate_cell;
Outlet_Properties.candidate_exterior_cell = candidate_exterior_cell;

Outlet_Properties.candidate_face_right = cand_right;
Outlet_Properties.candidate_face_left  = cand_left;
Outlet_Properties.candidate_face_up    = cand_up;
Outlet_Properties.candidate_face_down  = cand_down;

Outlet_Properties.slope_out_right = slope_out_right;
Outlet_Properties.slope_out_left  = slope_out_left;
Outlet_Properties.slope_out_up    = slope_out_up;
Outlet_Properties.slope_out_down  = slope_out_down;

Outlet_Properties.max_internal_slope = max_internal_slope;
Outlet_Properties.max_external_slope = max_external_slope;

Outlet_Properties.MinSlope = MinSlope;
Outlet_Properties.TieTolerance = TieTolerance;
Outlet_Properties.TiePriority = TiePriority;

Outlet_Properties.n_outlet_faces = nnz(face_right) + nnz(face_left) + ...
                                   nnz(face_up)    + nnz(face_down);
Outlet_Properties.n_outlet_cells = nnz(outlet_cell);

[Outlet_Properties.row_outlet_face, Outlet_Properties.col_outlet_face] = find(outlet_cell);

end