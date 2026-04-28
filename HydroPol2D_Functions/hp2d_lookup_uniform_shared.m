function Vq = hp2d_lookup_uniform_shared(Vtab, q, dz, maxDepth)
%HP2D_LOOKUP_UNIFORM_SHARED
% Fast uniform-depth lookup for CPU/GPU
%
% INPUTS
%   Vtab     : [ny x nx x nz] lookup table sampled on uniform depth bins
%   q        : [ny x nx] query depth
%   dz       : scalar depth interval
%   maxDepth : scalar maximum tabulated depth
%
% OUTPUT
%   Vq       : [ny x nx] interpolated value
%
% ASSUMPTION
%   depth grid is:
%       0, dz, 2*dz, ..., maxDepth
%   and nz = size(Vtab,3)

    [ny, nx, nz] = size(Vtab);
    N = ny * nx;

    qv = reshape(q, N, 1);

    qv(~isfinite(qv)) = 0;
    qv = max(qv, 0);
    qv = min(qv, maxDepth);

    idxL = floor(qv ./ dz) + 1;
    idxL = min(max(idxL, 1), nz-1);

    q0 = (idxL - 1) .* dz;
    a  = (qv - q0) ./ dz;

    atTop = (qv >= maxDepth);
    idxL(atTop) = nz-1;
    a(atTop)    = 1;

    V = reshape(Vtab, N, nz);

    base = (0:N-1)' * nz;
    ind1 = base + idxL;
    ind2 = ind1 + 1;

    v1 = V(ind1);
    v2 = V(ind2);

    Vqi = v1 + a .* (v2 - v1);

    Vqi(~isfinite(Vqi)) = 0;
    Vqi = max(Vqi, 0);

    Vq = reshape(Vqi, ny, nx);
end