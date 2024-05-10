function [DEM_raster,DEM,S] = DEM_smoothening(DEM_raster,min_area,flag_trunk,tau,K_value)
% Developer: Marcus Nobrega
% Smooth DEM in main flow accumulation paths
% Input: DEM_raster (GridOBj), min_area (km2), flag_trunk = [0,1], tau, K_value
% min_area is the flow accumulation threshold
% If flag_trunk == 1, we smooth only the main channel
% tau = [0, ... 1], and K_value = [0 - 10] is the degree of smoothing

%% DEM Smoothening
% Basic Calculations (Flow Direction and Accumulation)
FD = FLOWobj(DEM_raster,'preprocess','fill'); % Fill sinks and calculate flow direction

% Stream Network
% We then extract the stream network for a threshold upstream area 
area_km2 = min_area; % km2
area_cells = area_km2./((DEM_raster.cellsize/1000)^2); % pixels
S = STREAMobj(FD,'minarea',area_cells); % Streams with faccum > area_threshold

% If you want to smooth only the largest river, use flag_trunk == 1
DEMs = DEM_raster; % Smoothed DEM = Original DEM (we still will do the smoothing)
if flag_trunk == 1
    S = klargestconncomps(S,1); % Largest Stream
    S = trunk(S); % Only largest river
end

% plotdz(S,DEM_raster,'color','k'); % Trunked River
% pause(2)

% Smoothing Calculation
% To hydrologically correct and smooth the river profile, we use the function crs:
% min_gradient = nan; % Minimum Gradient
% zs = crs(S,DEM,'K',K_value,'tau',tau_value,'mingradient',min_gradient); % Hydrologically corrected DEM
zs = crs(S,DEM_raster,'K',K_value,'tau',tau); % Hydrologically corrected DEM

% Rebuilding the DEM using zs (Smoothed DEM)
DEMs.Z(S.IXgrid) = zs; % Not quite sure if 100% correct, but its supposed to get all streams
DEM_raster = DEMs;

% Raster Value
DEM = DEM_raster.Z;
end