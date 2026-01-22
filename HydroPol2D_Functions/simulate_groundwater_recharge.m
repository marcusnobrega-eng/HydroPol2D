function [recharge_rate, updated_soil_moisture, cumulative_recharge] = simulate_groundwater_recharge(infiltration_rate, initial_soil_moisture, alpha, dt, min_soil_moisture, max_soil_moisture, idx_imp, current_recharge)
%% ═══════════════════════════════════════════════════════════════════════
%  Function: simulate_groundwater_recharge
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 03/06/2025
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Simulate groundwater recharge using a linear reservoir approach.
%      This function computes the recharge rate based on infiltration,
%      updates the soil moisture state, and tracks cumulative recharge,
%      while ensuring mass balance and enforcing physical soil moisture limits.
%
%  ➤ Inputs:
%      • infiltration_rate      - Infiltration rate at the surface (mm/h)
%      • initial_soil_moisture  - Initial soil moisture in the unsaturated zone (mm)
%      • alpha                  - Linear coefficient representing the recharge rate (1/h)
%      • dt                     - Time step (hours)
%      • min_soil_moisture      - Minimum allowable soil moisture (m)
%      • max_soil_moisture      - Maximum allowable soil moisture (m)
%      • idx_imp                - Logical mask for impervious areas (if applicable)
%      • current_recharge       - Current cumulative recharge (mm)
%
%  ➤ Outputs:
%      • recharge_rate          - Computed recharge rate (m/s)
%      • updated_soil_moisture  - Updated soil moisture state (m)
%      • cumulative_recharge    - Total cumulative recharge (mm)
%
%  ➤ Notes:
%      • The recharge rate is initially computed using a linear reservoir
%        approach, then recalculated to maintain mass balance.
%      • Soil moisture is bounded between minimum and maximum limits to
%        prevent non-physical values.
%      • A mass balance check is performed, and an error is raised if the 
%        discrepancy exceeds 1% of the product (dt*infiltration_rate).
% ═══════════════════════════════════════════════════════════════════════
    
    % Calculate the recharge rate based on the linear reservoir approach
    effective_moisture = initial_soil_moisture; effective_moisture(effective_moisture < 0.005) = 0; % No recharge for inf depth < 5 mm
    recharge_rate = alpha .* effective_moisture;
    % recharge_rate(idx_imp) = 0;
    
    % Update soil moisture state
    updated_soil_moisture = initial_soil_moisture + dt .* (infiltration_rate - recharge_rate);
    
    % Ensure that soil moisture does not go below zero
    updated_soil_moisture = max(updated_soil_moisture, min_soil_moisture);
    updated_soil_moisture = min(updated_soil_moisture, max_soil_moisture);

    % Recalculating the recharge rate to ensure no physical violations
    % dS/dt = f - R
    % R = f - dS/dt
    recharge_rate =  infiltration_rate - (updated_soil_moisture -initial_soil_moisture)/dt; % m/s

    % Mass balance check
    error = nansum(nansum((updated_soil_moisture - initial_soil_moisture) - dt*(infiltration_rate - recharge_rate)));

    % Cumulative recharge
    cumulative_recharge = current_recharge + recharge_rate*1000*dt; % Cumulative recharge in mm
    
end
