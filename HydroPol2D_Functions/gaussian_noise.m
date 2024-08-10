function [noise] = gaussian_noise(mean,variance,rows,cols)
% Generates a gaussian noise with bias given by the mean
%
% Input:
% mean: bias in the noise
% variance
% rows: integer indicating the number of rows
% cols: integer indicatin the number of cols
% 
% Output:
% noise
% 
% Equation:
% N~(mu,var) = mu + sqrt(var).random number of known size

noise = mean + sqrt(variance).*randn(rows,cols);
end