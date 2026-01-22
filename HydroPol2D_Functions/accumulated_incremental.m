function [accum_variable, delta_variable, variable_intensity] = accumulated_incremental(steps, variable_discretized, time_step_variable)
% Optimized version of accumulated_incremental
% Author: Marcus Nobrega Gomes Junior
% Last Updated: Optimized by ChatGPT - July 2025

% Inputs:
%   steps                - Number of time steps
%   variable_discretized - Matrix of variable values (rows = series, cols = time)
%   time_step_variable   - Time step in minutes

% Outputs:
%   accum_variable       - Accumulated variable over time
%   delta_variable       - Incremental change at each time step
%   variable_intensity   - Intensity (e.g., mm/h)

% Precompute constants
dt_hours = time_step_variable / 60;

% Preallocate
[n_series, n_times] = size(variable_discretized);
delta_variable = zeros(n_series, steps);

% Only compute up to min(steps, n_times)
valid_steps = min(steps, n_times);

% Cumulative sum for accumulation
accum_variable = cumsum(variable_discretized(:,1:valid_steps) * dt_hours, 2);

% Compute delta values efficiently
delta_variable(:,1:valid_steps) = [accum_variable(:,1), diff(accum_variable, 1, 2)];

% Fill remaining steps (if any) with zeros
if steps > n_times
    delta_variable(:, n_times+1:steps) = 0;
end

% Intensity in mm/h
variable_intensity = delta_variable / dt_hours;

end
