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


%%%%%%%%%%%%%%%%%%% Another Source Code %%%%%%%%%%%%%%%%%%%%%%%%

% function Fint = idw_function(X0, F0, X, Y, p, rad, L)
%     % Default input parameters
%     if nargin < 7
%         L = 2; % 2nd norm (Euclidean)
%         if nargin < 6
%             rad = inf; % Radius of influence is infinity
%             if nargin < 5
%                 p = 2; % Normal weight
%             end
%         end
%     end
% 
%     % Check if we have GPU arrays
%     if isgpuarray(X0)
%         L = gpuArray(L);
%         p = gpuArray(p);
%         rad = gpuArray(rad);
%     end
% 
%     % Get the number of points and the number of observations
%     [nPts, nObs] = size(X0, 1), size(F0, 1);
%     nMesh = size(X, 1) * size(X, 2); % Total number of mesh points
% 
%     % Precompute the absolute distances (avoid looping)
%     % DeltaX and DeltaY are the differences between the mesh points and observation points
%     DeltaX = bsxfun(@minus, X, reshape(X0(:, 1), [], 1)); % X mesh points - X0
%     DeltaY = bsxfun(@minus, Y, reshape(X0(:, 2), [], 1)); % Y mesh points - Y0
% 
%     % Compute the distance matrix (vectorized)
%     DabsL = (abs(DeltaX).^L + abs(DeltaY).^L).^(1/L);
% 
%     % Filters: Set values greater than the radius to inf and zero distances to eps
%     DabsL(DabsL == 0) = eps;
%     DabsL(DabsL > rad) = inf;
% 
%     % Compute the weights (inverse distance to the power p)
%     W = 1 ./ (DabsL.^p);
% 
%     % Reshape the observed values for broadcasting
%     F0_reshaped = reshape(F0, [1, 1, nObs]);
% 
%     % Perform the weighted sum and normalization to get the interpolated values
%     W_sum = sum(W, 3); % Sum of weights for normalization
%     W_values = sum(W .* F0_reshaped, 3); % Weighted sum of values
% 
%     % Compute the interpolated values (avoid divide by zero errors)
%     Fint = W_values ./ W_sum;
%     Fint(W_sum == 0) = NaN; % If no valid weights, assign NaN
% end
