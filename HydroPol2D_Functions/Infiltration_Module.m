%% ═══════════════════════════════════════════════════════════════════════
%  Function: infiltrationModule
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 03/06/2025
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Compute the infiltration process for the hydrological model. This 
%      module calculates the effective infiltration rate, updates the soil
%      moisture storage, and computes the error in the infiltration volume.
%      It handles both standard and subgrid (inbank/overbank) processes.
%
%  ➤ Inputs:
%      • Wshed_Properties:
%           - Resolution: Spatial resolution (used to compute cell area)
%           - River_Width, River_Depth: Parameters for subgrid river modeling.
%
%      • Soil_Properties:
%           - I_t: Current UZ soil moisture (mm)
%           - ksat: Saturated hydraulic conductivity (mm/h)
%           - psi: Soil suction head (mm)
%           - teta_sat, teta_i: Saturated and initial soil moisture contents
%
%      • depths:
%           - d_t: Current total water depth (mm)
%           - d_p: Ponded water depth from previous time-step (mm) [for subgrid corrections]
%
%      • LULC_Properties:
%           - idx_imp: Logical mask for impervious areas (if applicable)
%
%      • flags:
%           - flag_infiltration: Activates the infiltration process
%           - flag_subgrid: Indicates subgrid modeling for finer channel detail
%           - flag_overbanks: Activates inbank/overbank corrections
%
%      • time_step: Duration of the current time-step (minutes)
%      • k: Time-step counter (used for initialization)
%      • DEM_raster: Digital Elevation Model (for domain dimensions and NaN mask)
%      • C_a: Cell area (may be modified during subgrid corrections)
%      • cumulative_infiltration: Matrix storing the cumulative infiltration
%      • errors: Array to store infiltration error (mass balance error)
%
%  ➤ Outputs:
%      • Updates Hydro_States with:
%           - i_a: Computed inflow rate (mm/h)
%           - f: Infiltration rate (mm/h)
%
%      • Updates Soil_Properties.I_t (soil moisture storage) based on 
%        infiltration computations.
%
%      • Updates depths.d_t (water depth) by subtracting the infiltrated 
%        volume.
%
%      • Updates cumulative_infiltration with the new infiltration volume.
%
%      • Updates errors(4) with the infiltration mass balance error (m³)
%
%  ➤ Local Functions:
%      • inbank_to_overbank(B,H,h_p,R)
%           - Converts inbank depths to overbank conditions.
%
%      • overbank_to_inbank(B,H,h_p,R)
%           - Converts overbank depths back to inbank conditions.
%
%  ➤ Notes:
%      • The module distinguishes between standard infiltration and scenarios
%        requiring subgrid corrections (i.e., inbank/overbank transitions).
%      • For high-resolution time-steps, the infiltration rate is approximated
%        using a simple min() function; otherwise, the GA Newton-Raphson 
%        method is used to solve the implicit infiltration equation.
%      • The code applies unit conversions as necessary (e.g., converting mm 
%        to m³, mm/h adjustments).
%      • Some sections (e.g., baseflow-related code) are commented out but 
%        indicate potential extensions.
%      • Ensure that all input structures (e.g., Wshed_Properties, Soil_Properties)
%        are properly initialized in the workspace.
% ═══════════════════════════════════════════════════════════════════════


% Infiltration Module

Coarse_Area = Wshed_Properties.Resolution^2;

% Current UZ + depth Storage
S_UZ_inf_0 = nansum(nansum(Coarse_Area.*Soil_Properties.I_t/1000));
if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    Vol = ((Wshed_Properties.Resolution - Wshed_Properties.River_Width).*Wshed_Properties.Resolution.*max((depths.d_t/1000 - Wshed_Properties.River_Depth),0) + ...
                      (Wshed_Properties.Resolution.*Wshed_Properties.River_Width.*depths.d_t/1000)); % m3 per cell
    S_p_inf_0 = nansum(nansum(Vol));
    eff_depth = 1000*(Vol / Coarse_Area); % Depth in the coarse model that matches the subgrid volume (mm)
else
    eff_depth = depths.d_t;
    S_p_inf_0 = nansum(nansum(Coarse_Area.*depths.d_t/1000));
end

S_inf_0 = S_UZ_inf_0 + S_p_inf_0;


if flags.flag_infiltration == 1
    % Inflow Rate
    Hydro_States.i_a = (eff_depth)./(time_step/60);  % Inflow rate [mm/h]

    % Effective Soil Moisture for Green-Ampt Calculation
    I_t_GA = max(Soil_Properties.I_t,5); % Limiting to 5 mm as minimum

    % Infiltration Capacity
    C = Soil_Properties.ksat.*(1 + ((eff_depth + ...
        Soil_Properties.psi).*(Soil_Properties.teta_sat - Soil_Properties.teta_i))./I_t_GA); % matrix form of Infiltration Capacity [mm/h]
    
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
        inf_volume = 1/1000*nansum(nansum((inf_volume_cell.*Coarse_Area))); % m3 total per domain

    else
        % We need to solve the implicit GA equation
        [Soil_Properties.I_t,Hydro_States.f] = GA_Newton_Raphson(Soil_Properties.I_t,time_step/60 ...
            ,Soil_Properties.ksat,Soil_Properties.psi,Soil_Properties.teta_sat - Soil_Properties.teta_i,eff_depth,Hydro_States.i_a,12,LULC_Properties.idx_imp);

        % if flags.flag_baseflow ~= 1 % Deep percolation
        %     Soil_Properties.I_t = max(Soil_Properties.I_t - Soil_Properties.k_out.*double(Hydro_States.idx_C)*time_step/60,min_soil_moisture);
        % end

        % Soil matrix balance already considered in the GA model, so
        % recharge occurs with f = 0
        
        % Infiltrated Volume
        inf_volume_cell = Hydro_States.f*(time_step/60); % Infiltrated volume in mm per coarse cell
        inf_volume = 1/1000*nansum(nansum((inf_volume_cell/1000 .* Coarse_Area))); % m3 total per domain

    end

    depths.d_t = depths.d_t - inf_volume_cell .* Coarse_Area ./ C_a; % Taking care of the infiltration depth
end

if flags.flag_subgrid == 1 && flags.flag_overbanks % Maybe we have a change from inbank <-> overbank
    % Eq. 15 and 16 of A subgrid channel model for simulating river hydraulics andfloodplain inundation over large and data sparse areas
    % Inbank - Overbank
    idx = depths.d_t/1000 > Wshed_Properties.River_Depth & depths.d_p/1000 <= Wshed_Properties.River_Depth & (Wshed_Properties.River_Width>0); % Cells in which there is a change from inbank to overbank
    if sum(sum(idx)) > 0
        [depths.d_t(idx)] = inbank_to_overbank(Wshed_Properties.River_Width(idx),Wshed_Properties.River_Depth(idx),depths.d_t(idx)/1000,Wshed_Properties.Resolution);
        C_a(idx) = Wshed_Properties.Resolution^2;
    end
    % Overbank - Inbank
    idx = depths.d_t/1000 <= Wshed_Properties.River_Depth & depths.d_p/1000 >= Wshed_Properties.River_Depth & (Wshed_Properties.River_Width>0); % Cells in which there is a change from overbank to inbank
    if sum(sum(idx)) > 0
        [depths.d_t(idx)] = overbank_to_inbank(Wshed_Properties.River_Width(idx),Wshed_Properties.River_Depth(idx),depths.d_t(idx)/1000,Wshed_Properties.Resolution);
        C_a(idx) = Wshed_Properties.River_Width(idx)*Wshed_Properties.Resolution;
    end 
end

% Final UZ + depth Storage
S_UZ_inf_t = nansum(nansum(Coarse_Area.*Soil_Properties.I_t/1000));

if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    S_p_inf_t = nansum(nansum((Wshed_Properties.Resolution - Wshed_Properties.River_Width).*Wshed_Properties.Resolution.*max((depths.d_t/1000 - Wshed_Properties.River_Depth),0))) + ...
                      nansum(nansum(Wshed_Properties.Resolution.*Wshed_Properties.River_Width.*depths.d_t/1000));
else
    S_p_inf_t = nansum(nansum(Coarse_Area.*depths.d_t/1000));
end
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

errors(4) = error; % m3 

% Cumulative Infiltration
if k == 1
    cumulative_infiltration = ones(size(DEM_raster.Z)); cumulative_infiltration(isnan(cumulative_infiltration)) = nan;
end
cumulative_infiltration = cumulative_infiltration + Hydro_States.f/1000/3600 * (time_step/60); % [mm]


%%%% Local Functions %%%%

function [h] = inbank_to_overbank(B,H,h_p,R)
    h = 1000*(B.* H + (R - B).* h_p) ./ R;
end

function [h] = overbank_to_inbank(B,H,h_p,R)
    h = 1000*(B .* h_p + (R - B) .* (H - h_p)) ./ B;
end