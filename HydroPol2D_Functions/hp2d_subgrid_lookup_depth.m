function Vq = hp2d_subgrid_lookup_depth(Vtab, depth, dz, maxDepth)
%HP2D_SUBGRID_LOOKUP_DEPTH Interpolate a uniform-depth subgrid table.

[ny,nx,nz] = size(Vtab);
Vq = NaN(size(depth), 'like', Vtab);

if ~isequal(size(depth), [ny,nx])
    error('hp2d_subgrid_lookup_depth:sizeMismatch', ...
        'Query depth size must match the first two dimensions of the table.');
end

d = max(depth, 0);
d = min(d, maxDepth);

idxL = floor(d ./ dz) + 1;
idxL = max(idxL, 1);
idxL = min(idxL, nz - 1);
idxU = idxL + 1;

dL = (idxL - 1) .* dz;
w = (d - dL) ./ dz;

atTop = d >= maxDepth;
idxL(atTop) = nz - 1;
idxU(atTop) = nz;
w(atTop) = 1;

[I,J] = ndgrid(1:ny, 1:nx);
indL = sub2ind([ny,nx,nz], I, J, idxL);
indU = sub2ind([ny,nx,nz], I, J, idxU);

vL = Vtab(indL);
vU = Vtab(indU);
ok = isfinite(depth) & isfinite(vL) & isfinite(vU) & isfinite(w);

Vq(ok) = vL(ok) + w(ok) .* (vU(ok) - vL(ok));
Vq(isfinite(depth) & depth <= 0) = 0;
end
