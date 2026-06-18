function [Q_face,Hf_x,Hf_y,Wf_x,Wf_y,face_x,face_y] = hp2d_subgrid_local_inertial_flux( ...
    eta_n,Q_prev,nc,Resolution,dt,h_min,idx_nan,SubgridTables,flag_numerical_scheme,V_available)
%HP2D_SUBGRID_LOCAL_INERTIAL_FLUX Local-inertial update in discharge form.
%
% Q_prev and Q_face are m3/s on padded HydroPol2D faces:
%   (:,:,1) x-face, positive from j to j+1
%   (:,:,2) y-face, positive from i to i+1

if nargin < 9 || isempty(flag_numerical_scheme)
    flag_numerical_scheme = 1;
end
if nargin < 10
    V_available = [];
end

[ny,nx] = size(eta_n);
g = 9.81;
small = eps_like(eta_n);

if isempty(idx_nan)
    inactive = false(ny,nx);
else
    inactive = logical(idx_nan);
end

if isempty(Q_prev)
    Q_prev = zeros(ny,nx,2,'like',eta_n);
end

Q_face = zeros(ny,nx,2,'like',eta_n);
Hf_x = zeros(ny,nx,'like',eta_n);
Hf_y = zeros(ny,nx,'like',eta_n);
Wf_x = zeros(ny,nx,'like',eta_n);
Wf_y = zeros(ny,nx,'like',eta_n);
face_x = struct();
face_y = struct();

if nx > 1
    etaL = eta_n(:,1:nx-1);
    etaR = eta_n(:,2:nx);
    eta_face_x = max(etaL, etaR);
    face_x = hp2d_subgrid_face_state(SubgridTables, eta_face_x, 'x', Resolution);
    Sx = (etaR - etaL) ./ Resolution;
    active_x = ~(inactive(:,1:nx-1) | inactive(:,2:nx)) & ...
        face_x.area > 0 & face_x.hrep > 0 & face_x.width > 0 & face_x.phi > 0 & ...
        isfinite(Sx) & isfinite(face_x.n) & face_x.n > 0;
    if h_min > 0
        active_x = active_x & face_x.hrep > h_min;
    end
    Qold = Q_prev(:,1:nx-1,1);
    Qnew = local_inertial_update(Qold, Resolution, face_x.hrep, face_x.n, Sx, dt, g, active_x);
    if flag_numerical_scheme == 2 || flag_numerical_scheme == 3
        Qnew = apply_adaptive_weighting_x(Qnew, Qold, face_x.area, face_x.width, dt, Resolution, g, flag_numerical_scheme);
        Qnew(~active_x) = 0;
    end
    Q_face(:,1:nx-1,1) = Qnew;
    Hf_x(:,1:nx-1) = face_x.hrep;
    Wf_x(:,1:nx-1) = face_x.width;
end

if ny > 1
    etaS = eta_n(1:ny-1,:);
    etaN = eta_n(2:ny,:);
    eta_face_y = max(etaS, etaN);
    face_y = hp2d_subgrid_face_state(SubgridTables, eta_face_y, 'y', Resolution);
    Sy = (etaN - etaS) ./ Resolution;
    active_y = ~(inactive(1:ny-1,:) | inactive(2:ny,:)) & ...
        face_y.area > 0 & face_y.hrep > 0 & face_y.width > 0 & face_y.phi > 0 & ...
        isfinite(Sy) & isfinite(face_y.n) & face_y.n > 0;
    if h_min > 0
        active_y = active_y & face_y.hrep > h_min;
    end
    Qold = Q_prev(1:ny-1,:,2);
    Qnew = local_inertial_update(Qold, Resolution, face_y.hrep, face_y.n, Sy, dt, g, active_y);
    if flag_numerical_scheme == 2 || flag_numerical_scheme == 3
        Qnew = apply_adaptive_weighting_y(Qnew, Qold, face_y.area, face_y.width, dt, Resolution, g, flag_numerical_scheme);
        Qnew(~active_y) = 0;
    end
    Q_face(1:ny-1,:,2) = Qnew;
    Hf_y(1:ny-1,:) = face_y.hrep;
    Wf_y(1:ny-1,:) = face_y.width;
end

if ~isempty(V_available)
    Q_face = cap_face_discharge_by_donor_volume(Q_face, V_available, dt);
end

Q_face(~isfinite(Q_face)) = 0;
Hf_x(~isfinite(Hf_x) | Hf_x <= small) = 0;
Hf_y(~isfinite(Hf_y) | Hf_y <= small) = 0;
Wf_x(~isfinite(Wf_x)) = 0;
Wf_y(~isfinite(Wf_y)) = 0;

end

function Qnew = local_inertial_update(Qold, full_width, H, n, S, dt, g, active)
% SFINCS-style subgrid local-inertial closure.
%
% The update is evaluated in grid-average unit-discharge form, matching the
% SFINCS subgrid LIE. HydroPol2D stores total face discharge Q [m3/s], so
% Qold is divided by the full coarse face width before the update and
% multiplied by that same width afterward.
Qnew = zeros(size(Qold), 'like', Qold);
qold = Qold ./ max(full_width, eps_like(Qold));
den = 1 + g .* dt .* n.^2 .* abs(qold) ./ ...
    max(max(H, 0).^(7/3), eps_like(Qold));
qnew = zeros(size(Qold), 'like', Qold);
qnew(active) = (qold(active) - g .* dt .* H(active) .* S(active)) ./ den(active);
Qnew(active) = qnew(active) .* full_width;
Qnew(~active) = 0;
end

function Q_face = cap_face_discharge_by_donor_volume(Q_face, V_available, dt)
%CAP_FACE_DISCHARGE_BY_DONOR_VOLUME Conservatively prevent donor overdraft.
[ny,nx,~] = size(Q_face);
out_rate = zeros(ny,nx,'like',V_available);

if nx > 1
    Qx = Q_face(:,1:nx-1,1);
    out_rate(:,1:nx-1) = out_rate(:,1:nx-1) + max(Qx, 0);
    out_rate(:,2:nx) = out_rate(:,2:nx) + max(-Qx, 0);
end

if ny > 1
    Qy = Q_face(1:ny-1,:,2);
    out_rate(1:ny-1,:) = out_rate(1:ny-1,:) + max(Qy, 0);
    out_rate(2:ny,:) = out_rate(2:ny,:) + max(-Qy, 0);
end

scale = ones(ny,nx,'like',V_available);
needs_cap = out_rate .* dt > V_available & out_rate > 0;
scale(needs_cap) = max(V_available(needs_cap), 0) ./ ...
    max(out_rate(needs_cap) .* dt, eps_like(V_available));
scale(~isfinite(scale)) = 0;
scale = max(min(scale, 1), 0);

if nx > 1
    Qx = Q_face(:,1:nx-1,1);
    scale_left = scale(:,1:nx-1);
    scale_right = scale(:,2:nx);
    donor_scale = ones(size(Qx), 'like', Qx);
    pos = Qx >= 0;
    donor_scale(pos) = scale_left(pos);
    donor_scale(~pos) = scale_right(~pos);
    Q_face(:,1:nx-1,1) = Qx .* donor_scale;
end

if ny > 1
    Qy = Q_face(1:ny-1,:,2);
    scale_up = scale(1:ny-1,:);
    scale_down = scale(2:ny,:);
    donor_scale = ones(size(Qy), 'like', Qy);
    pos = Qy >= 0;
    donor_scale(pos) = scale_up(pos);
    donor_scale(~pos) = scale_down(~pos);
    Q_face(1:ny-1,:,2) = Qy .* donor_scale;
end
end

function Qw = apply_adaptive_weighting_x(Qnew, Qold, A, W, dt, dx, g, scheme)
Qw = Qnew;
H = A ./ max(W, eps_like(A));
celerity = sqrt(g .* max(H, 0));
velocity = abs(Qold) ./ max(A, eps_like(A));
theta = 1 - dt ./ dx .* min(velocity, celerity);
theta(~isfinite(theta)) = 1;
theta = max(min(theta, 1), 0.7);
neighbor = zeros(size(Qold), 'like', Qold);
if size(Qold,2) > 2
    if scheme == 2
        neighbor(:,2:end-1) = Qold(:,1:end-2);
        neg = Qold < 0;
        right = zeros(size(Qold), 'like', Qold);
        right(:,2:end-1) = Qold(:,3:end);
        neighbor(neg) = right(neg);
    else
        neighbor(:,2:end-1) = 0.5 .* (Qold(:,1:end-2) + Qold(:,3:end));
    end
end
Qw = theta .* Qnew + (1 - theta) .* neighbor;
end

function Qw = apply_adaptive_weighting_y(Qnew, Qold, A, W, dt, dx, g, scheme)
Qw = Qnew;
H = A ./ max(W, eps_like(A));
celerity = sqrt(g .* max(H, 0));
velocity = abs(Qold) ./ max(A, eps_like(A));
theta = 1 - dt ./ dx .* min(velocity, celerity);
theta(~isfinite(theta)) = 1;
theta = max(min(theta, 1), 0.7);
neighbor = zeros(size(Qold), 'like', Qold);
if size(Qold,1) > 2
    if scheme == 2
        neighbor(2:end-1,:) = Qold(1:end-2,:);
        neg = Qold < 0;
        down = zeros(size(Qold), 'like', Qold);
        down(2:end-1,:) = Qold(3:end,:);
        neighbor(neg) = down(neg);
    else
        neighbor(2:end-1,:) = 0.5 .* (Qold(1:end-2,:) + Qold(3:end,:));
    end
end
Qw = theta .* Qnew + (1 - theta) .* neighbor;
end

function e = eps_like(x)
if isa(x, 'gpuArray')
    e = eps(classUnderlying(x));
else
    e = eps(class(x));
end
end
