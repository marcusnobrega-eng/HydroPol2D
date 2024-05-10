function [accum_variable,delta_variable,variable_intensity] = accumulated_incremental(steps,variable_discretized,time_step_variable)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 11 September, 2021                %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

% Calculates accumulated and incremental variable, assuming the model's
% time-step and the timeseries time-step
% 09/30/2020
%
% delta_p is the variable discretized into model's time-step
delta_variable = zeros(size(variable_discretized,1),steps);
for z = 1:size(variable_discretized,1)
    accum_variable = cumsum(variable_discretized(z,:)*time_step_variable/60);    
    for t = 1:(steps)
        perc_rainfall = t/steps*100
        if t == 1
            delta_variable(z,t) = accum_variable(t);
        elseif t >length(variable_discretized)
            delta_variable(z,t) = 0;
        else
            delta_variable(z,t) = accum_variable(t)-accum_variable(t-1); % mm
        end
    end

    variable_intensity = delta_variable/(time_step_variable/60); % mm/h  through all time-steps
    delta_variable = round(delta_variable,6); % ROUNDED VALUE!!! be careful
end

