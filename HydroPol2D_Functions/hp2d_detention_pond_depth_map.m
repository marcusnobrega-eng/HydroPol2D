function depth_m = hp2d_detention_pond_depth_map(Detention_Ponds, storage_m3)
%HP2D_DETENTION_POND_DEPTH_MAP Level-pool depth on the pond cells [m], NaN elsewhere.
%   max(stage - cell elevation, 0) on every pond cell (pond_depth_map_m in
%   HydroPol2D-Python). storage_m3 defaults to the current pond storage.

P = Detention_Ponds;
if nargin < 2
    storage_m3 = P.storage_m3;
end
stage = hp2d_detention_pond_stage(storage_m3, P);
depth_m = nan(size(P.mask));
depth_m(P.cells) = max(stage(P.cell_pond) - P.cell_elevation_m, 0);
end
