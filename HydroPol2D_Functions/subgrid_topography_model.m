function [q_face,Hf_x,Hf_y,Wf_x,Wf_y] = subgrid_topography_model( ...
    flag_numerical_scheme,eta_n,q_prev,nc,Resolution,dt,h_min,idx_nan,Subgrid_Properties,SubgridTables)
%SUBGRID_TOPOGRAPHY_MODEL Compatibility wrapper for lookup-table subgrid flow.
%
% This function preserves the older HydroPol2D subgrid helper contract:
% callers pass previous unit discharge q [m2/s] and receive updated q plus
% effective face depth/width. Internally, the repaired subgrid method solves
% the local-inertial equation in discharge form Q [m3/s] using shared-face
% area, wetted width, hydraulic radius, and roughness lookup tables.

% Compatibility input retained for older callers.
unused_Subgrid_Properties = Subgrid_Properties; %#ok<NASGU>

[ny,nx] = size(eta_n);

Q_prev = zeros(ny,nx,2,'like',eta_n);
if ~isempty(q_prev)
    Q_prev(:,:,1) = q_prev(:,:,1) .* Resolution;
    Q_prev(:,:,2) = q_prev(:,:,2) .* Resolution;
end

[Q_face,Hf_x,Hf_y,Wf_x,Wf_y] = hp2d_subgrid_local_inertial_flux( ...
    eta_n,Q_prev,nc,Resolution,dt,h_min,idx_nan,SubgridTables,flag_numerical_scheme);

q_face = zeros(ny,nx,2,'like',eta_n);
q_face(:,:,1) = Q_face(:,:,1) ./ max(Wf_x, eps_like(Wf_x));
q_face(:,:,2) = Q_face(:,:,2) ./ max(Wf_y, eps_like(Wf_y));
q_face(:,:,1) = zero_dry_width(q_face(:,:,1), Wf_x);
q_face(:,:,2) = zero_dry_width(q_face(:,:,2), Wf_y);
q_face(~isfinite(q_face)) = 0;
end

function q = zero_dry_width(q, W)
q(W <= 0 | ~isfinite(W)) = 0;
end

function e = eps_like(x)
if isa(x, 'gpuArray')
    e = eps(classUnderlying(x));
else
    e = eps(class(x));
end
end
