function [recharge_rate, updated_soil_moisture] = simulate_groundwater_recharge(infiltration_rate, initial_soil_moisture, alpha, dt, min_soil_moisture, max_soil_moisture, idx_imp)
    % Simulate groundwater recharge using a linear reservoir approach.
    %
    % INPUTS:
    % infiltration_rate      - The infiltration rate at the surface (mm/h)
    % initial_soil_moisture  - The initial soil moisture in the unsaturated zone (mm)
    % alpha                  - The linear coefficient representing the recharge rate (1/h)
    % dt                     - The time step (hours)
    % min_soil_moisture      - The minimum soil moisture (m)
    % max_soil_moisture      - The maximum soil moisture (m)
    % idx_impervious         - Logical mask showing surface impervious
    % areas
    %
    % OUTPUTS:
    % recharge_rate          - The computed recharge rate (m/s)
    % updated_soil_moisture  - The updated soil moisture (m/s)
    
    % Calculate the recharge rate based on the linear reservoir approach
    recharge_rate = alpha .* initial_soil_moisture;
    % recharge_rate(idx_imp) = 0;
    
    % Update soil moisture state
    updated_soil_moisture = initial_soil_moisture + dt .* (infiltration_rate - recharge_rate);
    
    % Ensure that soil moisture does not go below zero
    updated_soil_moisture = max(updated_soil_moisture, min_soil_moisture);
    updated_soil_moisture = min(updated_soil_moisture, max_soil_moisture);

    % Recalculating the recharge rate to ensure no physical violations
    % dS/dt = f - R
    % R = f - dS/dt
    recharge_rate =  infiltration_rate - (updated_soil_moisture -initial_soil_moisture)/dt;

    % Mass balance check
    error = nansum(nansum((updated_soil_moisture - initial_soil_moisture) - dt*(infiltration_rate - recharge_rate)));

    if error > 1/100*(dt*infiltration_rate)
        error('Mass balance error in recharge too large')
    end
    
end
