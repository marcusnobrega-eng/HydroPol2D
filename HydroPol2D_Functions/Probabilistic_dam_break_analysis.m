%% Probabilistic Dam-Break Analysis
clear all

load maximum_depths_18.08.mat
load outputs_18.08.mat

% Mask with maximum flood extent
threshold = 0.3; % m
mask = maximum_depths > 0.3;

% Maximum
total_mask = sum(mask,3);
total_mask = logical(total_mask);
% Minimum
% total_mask = sum(mask,3);
% total_mask = logical(total_mask);

for i = 1:10
    subplot(5,2,i)
    surf(double(mask(:,:,i)))
    view(0,90)
    shading interp
    colormap jet
end

%% Analysis per pixel (mask)
close all
stdev = zeros(size(mask,1),size(mask,2));
mean_cell = zeros(size(mask,1),size(mask,2));
stdev_mask = zeros(size(mask,1),size(mask,2));
mean_cell_mask = zeros(size(mask,1),size(mask,2));
for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        stdev_mask(i,j) = std(squeeze(mask(i,j,:)));
        mean_cell_mask(i,j) = mean(squeeze(mask(i,j,:)));
    end
end
subplot(2,1,1)
surf(mean_cell_mask);
view(0,90);
shading interp
colormap(turbo(5))
colorbar

subplot(2,1,2)
surf(stdev_mask);
view(0,90);
shading interp
colormap turbo
colorbar

%% Analysis per pixel (depth)
stdev = zeros(size(mask,1),size(mask,2));
mean_cell = zeros(size(mask,1),size(mask,2));
for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        stdev(i,j) = std(squeeze(maximum_depths(i,j,:)));
        mean_cell(i,j) = mean(squeeze(maximum_depths(i,j,:)));
    end
end
subplot(2,1,1)
surf(mean_cell);
view(0,90);
shading interp
colormap('parula')
colorbar

subplot(2,1,2)
surf(stdev);
view(0,90);
shading interp
colormap turbo
colorbar

%% Analise por celula
points = [190 317; 225 252];
depth_points = zeros(size(points,1),size(outputs,1));

for j = 1:size(points,1)
    subplot(2,1,j)
    row = points(j,1);
    col = points(j,2);
    depth_points(j,:) = sort(squeeze(maximum_depths(row,col,:)),"descend");
    pd = makedist('Normal',mean(depth_points(j,:)),std(depth_points(j,:)));
    % cdf_point(j,:) = cdf(pd,depth_points(j,:));
    x_plot = [0:0.01:max(depth_points(j,:))];
    p_ = cdf(pd,depth_points(j,:));
    p = cdf(pd,x_plot);
    cdfplot(depth_points(j,:))
    hold on
    plot(depth_points(j,:),p_);
end


% % Percentis
% pecentile = 0.01:0.01:1;
% 
% Y = cdf(depth_points);


