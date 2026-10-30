function [ETP, Ep, avg_temp, wind_speed, min_temp] = ETP_model(~, ~, ~, ~, x_grid, ~, varargin)
%ETP_MODEL Test double: a cold, humid day with negative Penman-Monteith values.
shape = size(x_grid');
ETP = -0.4 * ones(shape); ETP(1) = 2.0;
Ep = -0.6 * ones(shape); Ep(1) = 3.0;
avg_temp = 1 * ones(shape); wind_speed = 2 * ones(shape); min_temp = -2 * ones(shape);
end
