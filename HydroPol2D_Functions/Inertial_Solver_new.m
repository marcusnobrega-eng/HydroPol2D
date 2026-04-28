function q = Inertial_Solver(flag_numerical_scheme,q_p,dt,Hf,S,n_face_sq,dx,idx_nan)
% Calculates solution of Local Inertial Model
%
% Input:
% flag_numerical_scheme: chooses which numerical scheme is used
% 1: original Bates
% 2: s-upwind
% 3: s-centered
% q_p: discharge per face in [m2/s]
% dt: time-step in seconds
% Hf: effective flow depth in the interface [m]
% S: water surface elevation slope [m/m]
% n_face_sq: Manning's roughness coefficient squared at faces [s^2 m^(-2/3)]
% dx: space discretization [m]
% idx_nan: cells outside of the domain
%
% Output:
% q: discharge per face in [m2/s] at the end of the time-step

% -------------------------------------------------------------------------
% 1) Keep control/scalar values on CPU
% -------------------------------------------------------------------------
if isa(flag_numerical_scheme,'gpuArray')
    flag_numerical_scheme = gather(flag_numerical_scheme);
end
if isa(dt,'gpuArray')
    dt = gather(dt);
end
if isa(dx,'gpuArray')
    dx = gather(dx);
end

% -------------------------------------------------------------------------
% 2) Typed constants matching q_p
% -------------------------------------------------------------------------
g         = cast(9.81,'like',q_p);
dt_t      = cast(dt,'like',q_p);
dx_t      = cast(dx,'like',q_p);
alpha     = dt_t ./ dx_t;
gdt       = g .* dt_t;
eps_depth = cast(1e-12,'like',q_p);
one_t     = cast(1,'like',q_p);
zero_t    = cast(0,'like',q_p);
half_t    = cast(0.5,'like',q_p);

% -------------------------------------------------------------------------
% 3) Safe effective depth
% -------------------------------------------------------------------------
% Hf_safe = max(Hf, eps_depth);
Hf_safe = Hf;

% -------------------------------------------------------------------------
% 4) Outside-domain face mask (keep as logical)
% -------------------------------------------------------------------------
face_invalid = false(size(q_p),'like',idx_nan);
face_invalid(:,1:end-1,1) = idx_nan(:,1:end-1) | idx_nan(:,2:end);
face_invalid(:,end,1)     = idx_nan(:,end);
face_invalid(2:end,:,2)   = idx_nan(1:end-1,:) | idx_nan(2:end,:);
face_invalid(1,:,2)       = idx_nan(1,:);

% -------------------------------------------------------------------------
% 5) Split directions once (scheme 2/3 benefit most from this)
% -------------------------------------------------------------------------
qx = q_p(:,:,1);
qy = q_p(:,:,2);

Hfx = Hf_safe(:,:,1);
Hfy = Hf_safe(:,:,2);

Sx = S(:,:,1);
Sy = S(:,:,2);

nfx = n_face_sq(:,:,1);
nfy = n_face_sq(:,:,2);

invalid_x = face_invalid(:,:,1);
invalid_y = face_invalid(:,:,2);

qx(invalid_x) = 0;
qy(invalid_y) = 0;

[ny,nx] = size(qx);

% -------------------------------------------------------------------------
% Scheme 1: Original Bates
% -------------------------------------------------------------------------
if flag_numerical_scheme == 1

    abs_q = abs(q_p);
    invH73 = Hf_safe.^(-7/3);

    q = (q_p - gdt .* Hf .* S) ./ ...
        (one_t + gdt .* n_face_sq .* abs_q .* invH73);

% -------------------------------------------------------------------------
% Scheme 2: s-upwind
% -------------------------------------------------------------------------
elseif flag_numerical_scheme == 2

    % ----- X direction -----
    abs_qx  = abs(qx);
    celer_x = sqrt(g .* Hfx);
    theta_x = one_t - alpha .* min(abs_qx ./ Hfx, celer_x);
    theta_x(~isfinite(theta_x)) = one_t;

    qx_up = zeros(size(qx),'like',qx);

    % positive flow -> take left/upwind cell
    if nx > 2
        qx_up(:,2:nx-1) = qx(:,1:nx-2);
    end

    % negative flow -> take right/upwind cell
    neg_x = (qx < zero_t);
    if nx > 2
        qx_right = qx(:,3:nx);
        qx_mid   = qx_up(:,2:nx-1);
        neg_x_mid = neg_x(:,2:nx-1);
        qx_mid(neg_x_mid) = qx_right(neg_x_mid);
        qx_up(:,2:nx-1) = qx_mid;
    end

    invH73_x = Hfx.^(-7/3);
    qx_new = (theta_x .* qx + (one_t - theta_x) .* qx_up - gdt .* Hfx .* Sx) ./ ...
             (one_t + gdt .* nfx .* abs_qx .* invH73_x);

    % ----- Y direction -----
    abs_qy  = abs(qy);
    celer_y = sqrt(g .* Hfy);
    theta_y = one_t - alpha .* min(abs_qy ./ Hfy, celer_y);
    theta_y(~isfinite(theta_y)) = one_t;

    qy_up = zeros(size(qy),'like',qy);

    % positive flow -> take lower/upwind cell
    if ny > 2
        qy_up(2:ny-1,:) = qy(1:ny-2,:);
    end

    % negative flow -> take upper/upwind cell
    neg_y = (qy < zero_t);
    if ny > 2
        qy_upper = qy(3:ny,:);
        qy_mid   = qy_up(2:ny-1,:);
        neg_y_mid = neg_y(2:ny-1,:);
        qy_mid(neg_y_mid) = qy_upper(neg_y_mid);
        qy_up(2:ny-1,:) = qy_mid;
    end

    invH73_y = Hfy.^(-7/3);
    qy_new = (theta_y .* qy + (one_t - theta_y) .* qy_up - gdt .* Hfy .* Sy) ./ ...
             (one_t + gdt .* nfy .* abs_qy .* invH73_y);

    q = cat(3,qx_new,qy_new);

% -------------------------------------------------------------------------
% Scheme 3: s-centered
% -------------------------------------------------------------------------
elseif flag_numerical_scheme == 3

    % ----- X direction -----
    abs_qx  = abs(qx);
    celer_x = sqrt(g .* Hfx);
    theta_x = one_t - alpha .* min(abs_qx ./ Hfx, celer_x);
    theta_x(~isfinite(theta_x)) = one_t;
    theta_x(isnan(theta_x)) = 1;
    qx_avg = zeros(size(qx),'like',qx);
    if nx > 2
        qx_avg(:,2:nx-1) = half_t .* (qx(:,1:nx-2) + qx(:,3:nx));
    end

    invH73_x = Hfx.^(-7/3);
    qx_new = (theta_x .* qx + (one_t - theta_x) .* qx_avg - gdt .* Hfx .* Sx) ./ ...
             (one_t + gdt .* nfx .* abs_qx .* invH73_x);

    % ----- Y direction -----
    abs_qy  = abs(qy);
    celer_y = sqrt(g .* Hfy);
    theta_y = one_t - alpha .* min(abs_qy ./ Hfy, celer_y);
    theta_y(~isfinite(theta_y)) = one_t;
    theta_y(isnan(theta_y)) = 1;

    qy_avg = zeros(size(qy),'like',qy);
    if ny > 2
        qy_avg(2:ny-1,:) = half_t .* (qy(1:ny-2,:) + qy(3:ny,:));
    end

    invH73_y = Hfy.^(-7/3);
    qy_new = (theta_y .* qy + (one_t - theta_y) .* qy_avg - gdt .* Hfy .* Sy) ./ ...
             (one_t + gdt .* nfy .* abs_qy .* invH73_y);

    q = cat(3,qx_new,qy_new);

else
    error('Invalid flag_numerical_scheme. Use 1 (Bates), 2 (s-upwind), or 3 (s-centered).');
end

% -------------------------------------------------------------------------
% Final hygiene
% -------------------------------------------------------------------------
q(face_invalid) = 0;
q(~isfinite(q)) = 0;

end