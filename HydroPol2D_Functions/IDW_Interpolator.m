% Variable Spatial Interpolator with IDW method
% Developer: Marcus Nobrega
% Goal - Spatially interpolate observed values in a 2-D space with known
% coordinates of the observations
%
% Input:
% x_coordinate: col vector of x easting coordinates [L]
% y_coordinate: col vector of y coordinates [L]
% var_obs: col vector of observed values for each point 
% x_grid: vector with x_grid discretization, oriented to
% rightward direction [L]
% y_grid: vector with y_grid discretization, not coordinates, oriented from
% upward direction
%
% Output:
% Var_Interpolated: meshgrid with interpolated values in each pixel
% X: Matrix with coordinates in x
% Y: Matrixd with coordinates in y
%
% Example
%
% Interpolate the raifall data in a watershed with the following data:
% x_coordinate = [2 20 40]'; 
% y_coordinate = [5 10 80]';
% x_grid = [1:1:50]'; % 
% y_grid = [1:1:100]'; %
% var_obs = [10 50 0]'; % mm/h
% [Var_Interpolated,X,Y] = Rainfall_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid);
% contourf(X, Y, Var_Interpolated, 20)
% hold on
% plot(x_coordinate, y_coordinate, 'r.', 'MarkerSize', 30)
% colormap('jet')
% cb = colorbar;
% cb.Label.String = 'Rainfall [mm/h]';
% xlabel('x coordinates','interpreter','latex');
% ylabel('y coordinates','interpreter','latex');


function [Var_Interpolated,X,Y] = IDW_Interpolator(x_coordinate,y_coordinate,var_obs,x_grid,y_grid)
%% Meshgrid
[X,Y] = meshgrid(x_grid,y_grid);

%% Observations
X1 = x_coordinate; % Col Vector of easting (m)
Y1 = y_coordinate; % Col Vector of northing (m)
F = var_obs; % Col Vector with observed value

Var_Interpolated = idw_function([X1,Y1],F,X,Y); % Call idw_function
end