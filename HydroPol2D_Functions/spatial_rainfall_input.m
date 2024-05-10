% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2022                        %
%                                                                 %
%                 Last Updated: 1 July, 2022                      %
%   Goal: Set the model to solve spatial-temporal rainfall         %
function [spatial_rainfall] = spatial_rainfall_input(rainfall_spatial_duration,t)
%% 1.0 - Read Rainfall Maps
% We are assuming that the input maps are following the
% "rainfall_timestep_maps" duration, in minutes.
% Please, enter the name of the rasters and concatenate them into a
% multi-dimensional array
% --------- Convention ----------- %
% Let's assume you have n maps of rainfall. You should enter them here as:
% spatial_rainfall_1, spatial_rainfall_2 ... spatial_rainfall_n"
% Therefore, you create you spatial_rainfall 3D array as:
% spatial_rainfall = cat(spatial_rainfall_1, spatial_rainfall_2, ...
% spatial_rainfall_n)

% ---------- We need a code to read these files better -------------- %
i = find(rainfall_spatial_duration <= t,1,'last');
files = dir('*/*.asc');
z1 = files(i).folder;
z2 = '/';
z3 = files(i).name;
file_address = strcat(z1,z2,z3);
rain = load(file_address);
rain(rain<0) = nan;
spatial_rainfall = rain;
end

