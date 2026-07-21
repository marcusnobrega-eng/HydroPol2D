function [Q,Qc,Qf,Qci,Qfi,C_a] = subgrid_channel( ...
    h, w, zf, zc, dx, nc, nf, Qc_prev, Qf_prev, Qci_prev, Qfi_prev, g, dt, idx_rivers, outlet_index)
%SUBGRID_CHANNEL Neal et al. (2012) rectangular channel-floodplain routing.
%   This implementation follows the D4 channel + floodplain split described
%   in Neal et al. (2012). Channels are rectangular and sub-resolution, the
%   floodplain uses the coarse-cell width that remains outside the channel,
%   and the total face discharge is Q = Qc + Qf.

if nargin < 15
    error('subgrid_channel requires 15 inputs.');
end

[ny,nx] = size(h);

% D4-only state arrays. Diagonal memory is retained for interface
% compatibility but is not used in the Neal-mode implementation.
Qc = zeros(ny,nx,4,'like',h);
Qf = zeros(ny,nx,4,'like',h);
Q  = zeros(ny,nx,2,'like',h);
Qci = zeros(ny,nx,2,'like',h);
Qfi = zeros(ny,nx,2,'like',h);

nan_row = nan(1,nx,'like',h);
nan_col = nan(ny,1,'like',h);

w = max(w, 0);
h = max(h, 0);
idx_rivers = idx_rivers & (w > 0);

% Preserve non-channel cells as ordinary floodplain cells.
w_channel = w;
w_channel(~idx_rivers) = 0;

% Channel head and slopes [Eq. 6].
yc = h + zc;
Sc = slope_function(yc, dx, nan_col, nan_row);

% Channel flow depth [Eq. 7].
hfcflow = hf_function(yc, zc, nan_col, nan_row);
hfcflow = max(hfcflow, 0);

% Face channel width [Eq. 8].
wc_flow = width_function_min(w_channel, nan_col, nan_row);
wc_flow(~isfinite(wc_flow)) = 0;

Acflow = wc_flow .* hfcflow;                     % Eq. 9
Rc = Acflow ./ max(wc_flow + 2 .* hfcflow, eps); % Eq. 10

% A channel that fills the face is an explicitly resolved wide-flow cell.
% Using hydraulic depth here makes the Neal momentum equation reduce
% exactly to the ordinary local-inertial equation as w approaches dx.
full_width_face = wc_flow >= dx .* (1 - 10 .* eps(dx));
Rc(full_width_face) = hfcflow(full_width_face);

Qc_prev = ensure_d4_memory(Qc_prev, Qci_prev, ny, nx, 'like', h);
mask_channel = (wc_flow <= 0) | (hfcflow <= 0);

Qc = (Qc_prev - g .* Acflow .* dt .* Sc) ./ ...
    (1 + g .* dt .* nc.^2 .* abs(Qc_prev) ./ max((Rc .^ (4/3)) .* Acflow, eps)); % Eq. 11

Qc(mask_channel) = 0;
Qc(~isfinite(Qc)) = 0;

% Floodplain depth above the coarse DEM [Eq. 12].
hflood = max(0, h + zc - zf);
yf = hflood + zf;

Sf = slope_function(yf, dx, nan_col, nan_row);
hfflood = hf_function(yf, zf, nan_col, nan_row);
hfflood = max(hfflood, 0);

floodplain_width = max(dx - wc_flow, 0);
Qf_prev = ensure_d4_memory(Qf_prev, Qfi_prev, ny, nx, 'like', h);
qf_prev = zeros(size(Qf_prev), 'like', h);
wet_floodplain = floodplain_width > 0;
qf_prev(wet_floodplain) = Qf_prev(wet_floodplain) ./ floodplain_width(wet_floodplain);

mask_floodplain = (floodplain_width <= 0) | (hfflood <= 0);
Qf = ((qf_prev - g .* hfflood .* dt .* Sf) ./ ...
    (1 + g .* dt .* nf.^2 .* abs(qf_prev) ./ max(hfflood .^ (7/3), eps))) .* floodplain_width; % Eq. 13

Qf(mask_floodplain) = 0;
Qf(~isfinite(Qf)) = 0;

Q(:,:,1) = Qc(:,:,1) + Qf(:,:,1);
Q(:,:,2) = Qc(:,:,2) + Qf(:,:,2);

Q(:,end,1) = 0;
Q(1,:,2) = 0;
Qc(:,end,1) = 0;
Qc(1,:,2) = 0;
Qf(:,end,1) = 0;
Qf(1,:,2) = 0;

% Active water-surface area used by continuity [after Eqs. 14-16].
C_a = hp2d_neal_cell_area(h, w, zf - zc, dx);

% Outlet cells keep their internal storage area. The sink is applied later.
Q(repmat(idx_nan_or_outside(zf, outlet_index), 1, 1, 2)) = 0;
Qc(repmat(idx_nan_or_outside(zf, outlet_index), 1, 1, 4)) = 0;
Qf(repmat(idx_nan_or_outside(zf, outlet_index), 1, 1, 4)) = 0;

Qc = Qc(:,:,1:2);
Qf = Qf(:,:,1:2);

end

function mask = idx_nan_or_outside(zf, outlet_index)
mask = isnan(zf);
if nargin > 1 && ~isempty(outlet_index)
    mask = mask | false(size(outlet_index));
end
end

function Qd4 = ensure_d4_memory(Qprev, Qi_prev, ny, nx, varargin)
Qd4 = zeros(ny,nx,4,varargin{:});

if isempty(Qprev)
    return;
end

if isscalar(Qprev)
    return;
end

if ndims(Qprev) == 3
    ncopy = min(size(Qprev,3), 2);
    Qd4(:,:,1:ncopy) = Qprev(:,:,1:ncopy);
end

if nargin >= 2 && ~isempty(Qi_prev) && ~isscalar(Qi_prev) && ndims(Qi_prev) == 3
    ncopy = min(size(Qi_prev,3), 2);
    Qd4(:,:,3:(2+ncopy)) = Qi_prev(:,:,1:ncopy);
end
end

function S = slope_function(M, dx, nan_col, nan_row)
S = zeros(size(M,1), size(M,2), 4, 'like', M);
S(:,:,1) = [(M(:,2:end) - M(:,1:end-1)) ./ dx, nan_col]; % East
S(:,:,2) = [nan_row; (M(1:end-1,:) - M(2:end,:)) ./ dx]; % North
end

function hf = hf_function(y, z, nan_col, nan_row)
hf = zeros(size(y,1), size(y,2), 4, 'like', y);
hf(:,:,1) = [max(y(:,2:end), y(:,1:end-1)) - max(z(:,2:end), z(:,1:end-1)), nan_col];
hf(:,:,2) = [nan_row; max(y(1:end-1,:), y(2:end,:)) - max(z(1:end-1,:), z(2:end,:))];
end

function wc = width_function_min(w, nan_col, nan_row)
wc = zeros(size(w,1), size(w,2), 4, 'like', w);
wc(:,:,1) = [min(w(:,1:end-1), w(:,2:end)), nan_col];
wc(:,:,2) = [nan_row; min(w(1:end-1,:), w(2:end,:))];
end
