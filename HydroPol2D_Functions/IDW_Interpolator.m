function [Var_Interpolated,X,Y,cache] = IDW_Interpolator( ...
    x_coordinate, y_coordinate, var_obs, x_grid, y_grid, method, cache)
% Variable Spatial Interpolator with IDW or Nearest Neighbor method
%
% Input:
% x_coordinate : column vector of x coordinates
% y_coordinate : column vector of y coordinates
% var_obs      : column vector of observed values
% x_grid       : x grid vector
% y_grid       : y grid vector
% method       : 'idw' or 'nearest'
% cache        : optional struct for precomputed data
%
% Output:
% Var_Interpolated : interpolated field on grid
% X, Y             : meshgrid coordinates
% cache            : precomputed cache for reuse
%
% Notes:
% - 'nearest' is optimized for repeated calls.
% - If station coordinates and grid do not change, the nearest-neighbor
%   search is done only once.
%
% Example:
% cache = [];
% [V,X,Y,cache] = IDW_Interpolator(xc,yc,f,xg,yg,'nearest',cache);
% [V,X,Y,cache] = IDW_Interpolator(xc,yc,f2,xg,yg,'nearest',cache);

%% Defaults
if nargin < 6 || isempty(method)
    method = 'idw';
end

if nargin < 7 || isempty(cache)
    cache = struct();
end

%% Force column vectors and convert to double
x_coordinate = double(x_coordinate(:));
y_coordinate = double(y_coordinate(:));
var_obs      = double(var_obs(:));
x_grid       = double(x_grid(:));
y_grid       = double(y_grid(:));

%% Input checks
if ~isvector(x_coordinate) || ~isvector(y_coordinate) || ~isvector(var_obs)
    error('x_coordinate, y_coordinate, and var_obs must be vectors.');
end

if numel(x_coordinate) ~= numel(y_coordinate) || numel(x_coordinate) ~= numel(var_obs)
    error('x_coordinate, y_coordinate, and var_obs must have the same length.');
end

%% Reuse or build meshgrid
need_new_grid = true;

if isfield(cache,'x_grid') && isfield(cache,'y_grid') && ...
   isfield(cache,'X') && isfield(cache,'Y')
    if isequal(cache.x_grid, x_grid) && isequal(cache.y_grid, y_grid)
        X = cache.X;
        Y = cache.Y;
        need_new_grid = false;
    end
end

if need_new_grid
    [X,Y] = meshgrid(x_grid,y_grid);
    cache.X = X;
    cache.Y = Y;
    cache.x_grid = x_grid;
    cache.y_grid = y_grid;
else
    X = cache.X;
    Y = cache.Y;
end

%% Interpolation
switch lower(method)

    case 'idw'
        % Keep original IDW behavior
        Var_Interpolated = idw_function([x_coordinate, y_coordinate], var_obs, X, Y);

    case 'nearest'
        % Fast repeated nearest-neighbor:
        % precompute nearest observation index for each grid cell once

        need_new_nearest = true;

        if isfield(cache,'method') && strcmpi(cache.method,'nearest') && ...
           isfield(cache,'x_coordinate') && isfield(cache,'y_coordinate') && ...
           isfield(cache,'nearest_idx') && isfield(cache,'grid_size')

            if isequal(cache.x_coordinate, x_coordinate) && isequal(cache.y_coordinate, y_coordinate)
                need_new_nearest = false;
            end
        end

        if need_new_nearest
            obs_points  = [x_coordinate, y_coordinate];
            grid_points = [X(:), Y(:)];

            % knnsearch expects double; inputs already cast above
            cache.nearest_idx = knnsearch(obs_points, grid_points);

            cache.x_coordinate = x_coordinate;
            cache.y_coordinate = y_coordinate;
            cache.grid_size    = size(X);
            cache.method       = 'nearest';
        end

        % Fast lookup only
        Var_Interpolated = reshape(var_obs(cache.nearest_idx), cache.grid_size);

    otherwise
        error('Unknown interpolation method. Use ''idw'' or ''nearest''.');
end

end