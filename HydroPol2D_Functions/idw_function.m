% Inverse Distance Weighting Method
% Developer: Marcus Nobrega
% Description:  Interpolate spatial variables
%
% Inputs:
% X0 = vector with x and y projected coordinates [x,y] where observations
% are available
% F0 = value of the observations at points X0
% Xint = Coordinates [xint,yint] where interpolated points are calculated
% p = exponent of the inverse weight function. Default = 2.
% rad = radius of interest [L]. Default = infinity
% L = norm used for distance calculation. Default is the Euclidean norm = 2
%
% Outputs:
% Fint = values at interpolated point Xint.


function Fint = idw_function(X0,F0,X,Y,p,rad,L)

%
% Default input parameters
if nargin < 7 % Less than 6 parameters
    L = 2; % 2nd norm
    if nargin < 6 % Less than 6 parameters
        rad = inf; % Radius of influence is infinity
        if nargin < 5 % Less than 4 parameters
            p = 2; % Normal weight
        end
    end
end

if isgpuarray(X0) % Checking if we have gpuarrays
    L = gpuArray(L);
    p = gpuArray(p);
    rad = gpuArray(rad);
end

% Inverse distance weight output
if isgpuarray(X0) % Checking if we have gpuarrays
    DabsL = gpuArray(zeros(size(X,1),size(X,2),size(X0,1))); % Distance from each observed point
else
    DabsL = zeros(size(X,1),size(X,2),size(X0,1)); % Distance from each observed point
end

for i = 1:size(X0,1) % Number of Observations
    % P-Norm from Observations
    DeltaX = abs(X - (X0(i,1)));
    DeltaY = abs(Y - (X0(i,2)));
    DabsL(:,:,i) = (DeltaX.^L + DeltaY.^L).^(1/L);
end

% Filters
DabsL(DabsL == 0) = eps; % Small Number
DabsL(DabsL > rad) = inf; % Outside of radius influence

% Weights
W = 1./(DabsL.^p); % Weights for each observed point for each point in the mesh

% Interpolation
values = reshape(F0,1,1,length(F0)); % Observed Values in a 3D array

Fint = sum(W.*values,3)./sum(W,3); % Interpolated Values

end