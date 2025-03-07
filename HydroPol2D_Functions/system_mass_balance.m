% Storage Computation
% Developer: Marcus Nobrega, Ph.D.,

function [dS, fluxes, Storage, error] = system_mass_balance(P, Qin, E_int, ETR, E_ow, Qout, S, S_p, S_UZ, S_GW,S_prev)
% Input:
% P: Precipitation [m3]
% Qin: Inflow Hydrograph [m3]
% E_int: Evaporation in canopy [m3]
% ETR: Real evapotranspiration [m3]
% E_ow: Evaporation in open water [m3]
% Qout: Outlet flow [m3]
% S: Plant canopy storage in [m3]
% S_p: Ponded depth storage [m3]
% S_UZ: UZ storage in [m3]
% S_GW: GW Storage in [m3]
% S_prev: Previous system storage [m3]
% 
% Output: 
% dS: Variation in storage [m3]
% fluxes: Total flux in [m3]
% Storage: Current storage in [m3]
% error: Mass balance error in [m3]

% Flux Computation [m3]
fluxes = P + Qin - E_int - ETR - E_ow - Qout;

% Variation in storage dS [m3]
dS = (S + S_p + S_UZ + S_GW) - S_prev;

% Storge [m3]
Storage = S_prev + dS;

% Error Computation
error = dS - fluxes;

if abs(error) > 100
    ttt = 1;
end
% if error/P > 1/100
%     error('Mass balance error to large')
% end
end