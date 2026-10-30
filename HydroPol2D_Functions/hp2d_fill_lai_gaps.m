function [lai, filled] = hp2d_fill_lai_gaps(lai, idx_nan, fallback)
%HP2D_FILL_LAI_GAPS Give active cells without an LAI value the fallback LAI.
%   A NaN LAI in an active cell made the canopy store NaN there, so that
%   cell's rain never reached the ground. HydroPol2D-Python fills such cells
%   with the configured fallback (static.lai, default 0) and warns
%   (core/preprocess.py, _prepare_leaf_area_index). Cells outside the domain
%   stay NaN. filled is the number of cells filled.
if nargin < 3 || isempty(fallback), fallback = 0; end
missing = ~idx_nan & ~isfinite(lai);
filled = nnz(missing);
if filled > 0
    warning('HydroPol2D:LAIFallback', ...
        'The LAI raster has %d missing active cells; those cells use the LAI fallback (%g).', ...
        filled, fallback);
    lai(missing) = fallback;
end
end
