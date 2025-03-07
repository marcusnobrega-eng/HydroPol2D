% Infiltration Module

% Current UZ + depth Storage
S_UZ_inf_0 = nansum(nansum(C_a.*Soil_Properties.I_t/1000));
S_p_inf_0 = nansum(nansum(C_a.*depths.d_t/1000));
S_inf_0 = S_UZ_inf_0 + S_p_inf_0;

if flags.flag_infiltration == 1
    % Inflow Rate
    Hydro_States.i_a = (depths.d_t)./(time_step/60);  % Inflow rate [mm/h]

    % Infiltration Capacity
    C = Soil_Properties.ksat.*(1 + ((depths.d_t + ...
        Soil_Properties.psi).*(Soil_Properties.teta_sat - Soil_Properties.teta_i))./Soil_Properties.I_t); % matrix form of Infiltration Capacity [mm/h]
    
    % if flags.flag_baseflow ~= 1 % Deactivating this method due to the linear reservoir approach
    %     % Cells with top layer exceeding root zone
    %     Non_C_idx = Soil_Properties.I_t > Soil_Properties.Lu*1000; % Cells that exceed the top layer infiltrated depth
    % 
    %     % Limiting infiltration in these areas
    %     C(Non_C_idx) = Soil_Properties.ksat(Non_C_idx); % If that condition ocurrs, reduce the capacity to the saturation
    % 
    %     % Cells where infiltration capacity is higher than inflow
    %     Hydro_States.idx_C = (Hydro_States.i_a <= C); % Cells where the inflow rate is smaller than the infiltration capacity
    % 
    %     % Cells where the inflow rate is larger than the infiltration capacity
    %     idx_i_a = (Hydro_States.i_a > C); % Cells where the inflow rate is larger than the infiltration capacity
    % end
    % % Recovery time
    % if k == 1
    %     Soil_Properties.T = zeros(size(elevation,1),size(elevation,2)); Soil_Properties.T(idx_nan) = nan;
    % else
    %     Soil_Properties.T = Soil_Properties.T - time_step/60; % Recoverying time (hours
    % end

    % if flags.flag_baseflow ~= 1
    %     Soil_Properties.T(idx_i_a) = Soil_Properties.Tr(idx_i_a); % Saturated Areas we begin the process again
    %     Soil_Properties.idx_T = Soil_Properties.T < 0; % Cells where the replenishing time is reached
    % 
    %     % Refreshing Soil_Properties.I_t to I_0 for cases where idx_T > 0
    %     Soil_Properties.I_t(Soil_Properties.idx_T) = min_soil_moisture(Soil_Properties.idx_T);
    % 
    % end

    % Infiltration rate
    if time_step*60 < 1/60 % If we are using high resolution time-steps
        % We can approximate the soil infiltration capacity curve
        Hydro_States.f = min(C,Hydro_States.i_a); % Infiltration rate (mm/hr)
        Hydro_States.f(LULC_Properties.idx_imp) = 0;

        % % Soil matrix mass balance
        % if flags.flag_baseflow ~= 1
        % 
        %     Soil_Properties.I_t = max(Soil_Properties.I_t + Hydro_States.f*time_step/60 - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
        % 
        % else
        %     Soil_Properties.I_t = Soil_Properties.I_t + Hydro_States.f*time_step/60; % Assuming all infiltrated depth into the soil matrix (recharge will be done with no f)
        % end

        % Infiltrated Volume
        inf_volume_cell = Hydro_States.f*(time_step/60); % Infiltrated volume in mm per cell
        inf_volume = 1/1000*nansum(nansum((inf_volume_cell).*C_a)); % m3 total per domain

    else
        % We need to solve the implicit GA equation
        [Soil_Properties.I_t,Hydro_States.f] = GA_Newton_Raphson(Soil_Properties.I_t,time_step/60 ...
            ,Soil_Properties.ksat,Soil_Properties.psi,Soil_Properties.teta_sat - Soil_Properties.teta_i,depths.d_t,Hydro_States.i_a,12,LULC_Properties.idx_imp);

        % if flags.flag_baseflow ~= 1 % Deep percolation
        %     Soil_Properties.I_t = max(Soil_Properties.I_t - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
        % end

        % Soil matrix balance already considered in the GA model, so
        % recharge occurs with f = 0
        
        % Infiltrated Volume
        inf_volume_cell = Hydro_States.f*(time_step/60); % Infiltrated volume in mm per cell
        inf_volume = 1/1000*nansum(nansum((inf_volume_cell).*C_a)); % m3 total per domain

    end

    depths.d_t = depths.d_t - inf_volume_cell; % Taking care of the infiltration depth
end

% Final UZ + depth Storage
S_UZ_inf_t = nansum(nansum(C_a.*Soil_Properties.I_t/1000));
S_p_inf_t = nansum(nansum(C_a.*depths.d_t/1000));
S_inf_t = S_UZ_inf_t + S_p_inf_t;

% Reservoir 1: surface
% dh1/dt = -f

% Reservoir 2: soil matrix
% dh2/dt = +f

% dh1/dt + dh2/dt = 0

% Storage(t + dt) - Storage(t) = 0

% Error 
dS_inf = (S_inf_t - S_inf_0);
error = dS_inf;
