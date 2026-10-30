function [x, y] = hp2d_overlap_raster_centres(mapping)
%HP2D_OVERLAP_RASTER_CENTRES Centres of the overlap raster cells, in the
%   overlap's own order (south-up, row-major: index (row - 1)*cols + col with
%   row 1 the southernmost).
x_edges = mapping.x_edges(:); y_edges = mapping.y_edges(:);
x_centre = (x_edges(1:end-1) + x_edges(2:end))/2;
y_centre = (y_edges(1:end-1) + y_edges(2:end))/2;
if y_centre(1) > y_centre(end), y_centre = flipud(y_centre); end
[X, Y] = meshgrid(x_centre, y_centre);   % rows south to north
x = reshape(X.', [], 1); y = reshape(Y.', [], 1);
end
