%%% Gaussian Noise %%%
% Goal - Create a single number representing a gaussian noise
function [g_noise] = gaussian_noise_generator(variance,average)
std_deviation = sqrt(variance);
% g_noise = (average + std_deviation*randn)
g_noise = average + std_deviation*randn;
end