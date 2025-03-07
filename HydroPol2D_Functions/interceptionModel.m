function [S, T, E, F] = interceptionModel(P, Ep, LAI, S_prev, C)
    % Estimates plant canopy interception following the Rutter Model
    %
    % Inputs:
    % P - Precipitation (mm)
    % Ep - Potential evaporation (mm)
    % LAI - Leaf Area Index (unitless)
    % S_prev - Interception storage from previous timestep (mm)
    % C - Canopy storage coefficient (typically 0.2-0.5 mm per LAI unit)
    % 
    % Outputs:
    % S - Updated canopy interception storage (mm)
    % T - Throughfall (mm)
    % E - Evaporation from interception storage (mm)
    % F: stemflow (mm)

    % Step 1: Compute Maximum Canopy Storage
    LAI(isnan(LAI)) = 0;
    LAI(isnan(P)) = nan;
    S_max = C * LAI; % Maximum interception storage (mm)    

    % Compute stemflow
    alpha = 0.3; % alpha = Decay parameter (depends on bark roughness, e.g., 0.3)
    maximum_delay = 20; % maximum delay in minutes
    F = computeStemflow(P, maximum_delay, alpha, S_prev, S_max);
    F = 0*F; % Not activated !

    % Compute Evaporation
    beta = S_prev ./ S_max;  beta(isnan(beta)) = 0;
    E = beta .* Ep;

    % Compute canopy storage change
    % dS = P - F - E - 0 (Assuming initally no throughfall)
    dS = P - F - E;

    % Compute S^{t + 1/2}
    S  = S_prev + dS; % No throughfall considered

    % dS = P - F - E - T, T = P - dS - F - E

    % Compute Throughfall
    T = max(S - S_max,0); % Throughfall is the difference between storage - maximum storage

    S = min(S,S_max); % Final canopy storage [mm]   

    % Mass Balance Check
    error = nansum(nansum(((S - S_prev) - (P - E - F - T))))/nansum(nansum(P))*100;

    if error > 0.1
        warning('Mass balance errors in the interception model of %.d (%)',error)
    end
    

    % Nested function to calculate stemflow
    function stemflowRate = computeStemflow(P, T0, alpha, canopyStorage, maxCanopyStorage)
        % P = Rainfall intensity (mm/hr)
        % T0 = Maximum possible delay (minutes) (e.g., 20 min)
        % alpha = Decay parameter (depends on bark roughness, e.g., 0.3)
        % canopy_storage = canopy water storage (mm)
        % maximum maxCanopyStorage (mm)
        % Computes the stemflow delay and rate
        stemflowDelay = T0 .* exp(-alpha .* P); % Time delay in minutes
        stemflowFraction = 0.05; % Fraction of intercepted water that becomes stemflow
        stemflowRate = stemflowFraction .* canopyStorage ./ (stemflowDelay / 60); % Convert to mm/hr
        stemflowRate = min(stemflowRate,maxCanopyStorage);
        
    end
end

