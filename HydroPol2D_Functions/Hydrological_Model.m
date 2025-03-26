% Infiltration Routine
% Developer: Marcus Nobrega, Ph.D.
% Goal: Estimate effective precipitation and infiltration at cells domain
% Date: 2/22/2025
% Added recharge calculation. Added 2D boussinesq solution to the
% groundwater model


% Interception Module
interception_module

% Evaporation / Evapotranspiration Module
Evaporation_Evapotranspiration_Module

% Infiltration Module 
Infiltration_Module

% Recharge and Groundwater Module
Groundwater_Module

% Inflow
BC_States.inflow(:,150) = BC_States.inflow(:,151);
BC_States.inflow(:,151) = 0;

depths.d_t = depths.d_t + BC_States.inflow;

% Depths
if min(min(depths.d_t)) < -1e-8
    catch_index = catch_index + 1;
    warning('Negative depths. Please reduce the time-step.')
else
    depths.d_t(depths.d_t < 1e-6) = 0;
end

depths.d_t(depths.d_t < 0) = 0;

% Effective Precipitation (Available depth at the beginning of the time-step
depths.d_tot = depths.d_t;

% Correcting depths in case it reach overbanks
if flags.flag_subgrid == 1 % Maybe we have a change from inbank <-> overbank
    % Eq. 15 and 16 of A subgrid channel model for simulating river hydraulics andfloodplain inundation over large and data sparse areas
    % Inbank - Overbank
    idx = depths.d_tot/1000 > Wshed_Properties.River_Depth & depths.d_p/1000 <= Wshed_Properties.River_Depth & idx_rivers; % Cells in which there is a change from inbank to overbank
    if sum(sum(idx)) > 0
        depths.d_tot(idx) = 1000*(depths.d_tot(idx)/1000 - (depths.d_tot(idx)/1000 - Elevation_Properties.elevation_cell(idx) + (Elevation_Properties.elevation_cell(idx) - Wshed_Properties.River_Depth(idx))).*(1 - Wshed_Properties.River_Width(idx)/Wshed_Properties.Resolution)); % mm
        C_a(idx) = Wshed_Properties.Resolution^2;
    end
    % Overbank - Inbank
    idx = depths.d_tot/1000 <= Wshed_Properties.River_Depth & depths.d_p/1000 >= Wshed_Properties.River_Depth & idx_rivers; % Cells in which there is a change from overbank to inbank
    if sum(sum(idx)) > 0
        factor = depths.d_tot(idx)/1000 - Elevation_Properties.elevation_cell(idx) + (Elevation_Properties.elevation_cell(idx) - Wshed_Properties.River_Depth(idx));
        depths.d_tot(idx) = 1000*(depths.d_tot(idx)/1000 + (Wshed_Properties.Resolution*(factor))./Wshed_Properties.River_Width(idx) ...
            - (depths.d_tot(idx)/1000 - Elevation_Properties.elevation_cell(idx) + (Elevation_Properties.elevation_cell(idx) - Wshed_Properties.River_Depth(idx)))); % mm
        C_a(idx) = Wshed_Properties.River_Width(idx)*Wshed_Properties.Resolution;
    end
end

