function drep = hp2d_subgrid_inverse_volume(Vtab, Vq, dz, maxDepth, Resolution)
%HP2D_SUBGRID_INVERSE_VOLUME Invert a monotone volume-depth lookup table.

[ny,nx,nz] = size(Vtab);
if ~isequal(size(Vq), [ny,nx])
    error('hp2d_subgrid_inverse_volume:sizeMismatch', ...
        'Volume query size must match the first two dimensions of the table.');
end

drep = NaN(ny,nx, 'like', Vq);
depth_axis = reshape((0:nz-1) * dz, 1, 1, nz);

finite_count = sum(isfinite(Vtab), 3);
Vtop = Vtab(:,:,end);
valid_cell = finite_count >= 2 & isfinite(Vtop) & isfinite(Vq);

if ~any(valid_cell(:))
    return;
end

Vwork = Vq;
Vwork(valid_cell) = max(Vwork(valid_cell), 0);
Vwork(~valid_cell) = NaN;

Vclamp = Vwork;
Vclamp(valid_cell) = min(Vwork(valid_cell), Vtop(valid_cell));

cmp = Vtab <= Vclamp;
cmp(~isfinite(Vtab)) = false;
idxL = sum(cmp, 3);
idxL = max(idxL, 1);
idxL = min(idxL, nz - 1);
idxU = idxL + 1;

[I,J] = ndgrid(1:ny, 1:nx);
indL = sub2ind([ny,nx,nz], I, J, idxL);
indU = sub2ind([ny,nx,nz], I, J, idxU);

V1 = reshape(Vtab(indL), ny, nx);
V2 = reshape(Vtab(indU), ny, nx);

D = repmat(depth_axis, ny, nx, 1);
d1 = reshape(D(indL), ny, nx);
d2 = reshape(D(indU), ny, nx);

ok = valid_cell & isfinite(V1) & isfinite(V2) & isfinite(d1) & isfinite(d2);
denom = V2 - V1;
bad = ok & abs(denom) <= eps(class_underlying_like(Vtab));
denom(bad) = 1;

a = zeros(ny,nx, 'like', Vq);
a(ok) = (Vclamp(ok) - V1(ok)) ./ denom(ok);
a(ok) = max(min(a(ok), 1), 0);
a(bad) = 0;

drep(ok) = d1(ok) + a(ok) .* (d2(ok) - d1(ok));

above_top = valid_cell & Vwork > Vtop;
drep(above_top) = maxDepth + (Vwork(above_top) - Vtop(above_top)) ./ (Resolution^2);
drep(valid_cell & Vq <= 0) = 0;
drep(valid_cell & ~isfinite(drep)) = 0;
end

function c = class_underlying_like(x)
if isa(x, 'gpuArray')
    c = classUnderlying(x);
else
    c = class(x);
end
end
