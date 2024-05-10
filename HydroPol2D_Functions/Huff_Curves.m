function [t, i, idf] = Huff_Curves(td, dt, K, a, b, c, RP, flag_plot)
% Developer: Marcus Nobrega, Ph.D
% Goal: Create a hyetograph from an IDF curve
close all
% Inputs:
% td: Rainfall duration (min)
% dt: time-step (min)
% K, a, b, c, RP: IDF parameters such that i(mm/hr) = K*(RP^a)/((b+td)^c)
% RP: retun periods in years
% flag_plot: boolean [0,1] where 1 we plot the results and 0 we don't.
%
% Outputs:
% t: vector with the time from 0 to td, with dt intervals
% i: rainfall intensity for each value of t in mm/h
% idf: IDF rainfall intensity for each t
%
% Example: Determine the alternated blocks hyetograph of a rainfall with 60
% min duration,
% and plot the hyetograph with 5-min intervals
%
% [t, i, idf] = huff_curves(60, 5, 819, 0.138, 10, 0.75, 10, 1)
% Problems with the function? Contact me at marcusnobrega.engcivil@gmail.com

%% 1. Set duration
t = dt:dt:td;
nsteps = length(t);

%% 2. Calculate intensity
idf = K.*(RP^a)./((b + t).^c); % Rainfall Intensity for all t durations
rainfal_volume = idf(end)*td/60; % mm
delta_t = zeros(1,length(idf));
for k = 1:(length(idf))
    if k == 1
        delta_t(k) = t(k);
    else
        delta_t(k) = t(k) - t(k-1);
    end
end
p = idf .* delta_t * 60; % Rainfall Volume for each duration

%% 3. Huff Curves Method
% Huff Curves can be written as a 6-th order polynomial for each quartile
% Therefore, we enter a matrix with each of the coefficients of a
% polynomial such that:
% p/pt = c6(t/d)^6 + c5(t/td)^5 + c4(t/td)^4 + c3(t/td)^3 + c2(t/td)^2 +
% c1(t/td)1 + c0
% Rows of this matrix are 1st, 2nd, 3rd, and 4th quartiles
huff_coefficients = [-0.9633	3.8869	-7.895	10.089	-8.0108	3.8936	-0.0032
-39.436	125.18	-146.04	73.604	-13.936	1.6243	-0.0068
46.542	-131.55	132.63	-57.315	10.796	-0.1107	0.005
-25.289	67.54	-64.926	28.031	-5.2061	0.8535	-0.0042];

t_prime = t / td;
p_prime_quartiles = zeros(4, nsteps);
for k = 1:nsteps
    for j = 1:size(huff_coefficients,2)
        exp = size(huff_coefficients,2) - j;
        p_prime_quartiles(:,k) = huff_coefficients(:,j)*(t_prime(k))^exp + p_prime_quartiles(:,k);
    end
end

% Calculating Rainfall Intensities
rainfall_intensity = zeros(size(p_prime_quartiles));
for k = 1:nsteps
    if k == 1
        rainfall_intensity(:,k) = rainfal_volume*p_prime_quartiles(:,k)/(dt/60); % mm/h
    else
        rainfall_intensity(:,k) = rainfal_volume*(p_prime_quartiles(:,k) - p_prime_quartiles(:,k-1))/(dt/60); % mm/h
    end
end
rainfall_intensity(rainfall_intensity < 0) = 0; % Boundary Condition

% i = zeros(1, nsteps);
% for k = 1:nsteps
% % Sort quartiles in descending order
% sorted_quartiles = sort(p_prime_quartiles(:,k), 'descend');
% % Calculate hyetograph based on Huff Curves method
% i(k) = max(sorted_quartiles(1) + sorted_quartiles(2), sorted_quartiles(3) + sorted_quartiles(4));
% end

%% 4. Plotting Hyetograph
colors = linspecer(4);
if flag_plot == 1
    % Creating Modeling Results Folder
        % Create the folder name
        folderName = 'Modeling_Results';
        
        % Check if the folder already exists
        if ~exist(folderName, 'dir')
            % If it doesn't exist, create the folder
            mkdir(folderName);
            disp('Folder "Modeling_Results" created successfully!');
        else
            disp('Data sucessfully exported in Modeling_Results Folder');
        end  


    figure % Plotting Intensities
    set(gcf,'units','inches','position',[3,3,6.5,4])  
    plot(t_prime,p_prime_quartiles(1,:),'color',colors(1,:),'LineStyle','-','LineWidth',2,'Marker','*'); hold on;
    plot(t_prime,p_prime_quartiles(2,:),'color',colors(2,:),'LineStyle','--','LineWidth',2,'Marker','.'); hold on;
    plot(t_prime,p_prime_quartiles(3,:),'color',colors(3,:),'LineStyle','-.','LineWidth',2,'Marker','<'); hold on;
    plot(t_prime,p_prime_quartiles(4,:),'color',colors(4,:),'LineStyle',':','LineWidth',2,'Marker','>');
    hold off
    xlabel('$t/t_d$','Interpreter','latex','FontSize',12)
    ylabel('$P/P_t$','Interpreter','latex','FontSize',12)
    title('Huff Cumulative Rainfall','Interpreter','latex','FontSize',12)
    legend('IDF Curve','Quartile 1','Quartile 2','Quartile 3','Quartile 4','Location','Best','Interpreter','latex','FontSize',12)
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)    
    grid on
    box on
    exportgraphics(gcf,fullfile(folderName,'Quartiles_Huff.pdf'),'ContentType','vector') 

    figure % Plotting Cumulative Volumes
    set(gcf,'units','inches','position',[3,3,6.5,4])      
    hold on
    plot(t_prime*td,rainfall_intensity(1,:),'LineStyle','-','LineWidth',2,'Marker','*');
    plot(t_prime*td,rainfall_intensity(2,:),'LineStyle','--','LineWidth',2,'Marker','.');
    plot(t_prime*td,rainfall_intensity(3,:),'LineStyle','-.','LineWidth',2,'Marker','<');
    plot(t_prime*td,rainfall_intensity(4,:),'LineStyle',':','LineWidth',2,'Marker','>');
    hold off
    xlabel('Elapsed Time [min]','Interpreter','latex','FontSize',12)
    ylabel('Rainfall Intensity [$\mathrm{mm.h^{-1}}$]','Interpreter','latex','FontSize',12);
    title('Huff Hyetograph','Interpreter','latex')
    legend('Quartile 1','Quartile 2','Quartile 3','Quartile 4','Location','Best','Interpreter','latex')
    grid on
    axis tight
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)    
    grid on
    box on
    title(strcat('RP = ',num2str(RP,'%.2d'),' years'),'interpreter','latex')   
    exportgraphics(gcf,fullfile(folderName,'Intensity_Huff_All_Quartiles.pdf'),'ContentType','vector') 


    figure % Plotting Cumulative Volumes
    set(gcf,'units','inches','position',[3,3,6.5,4])      
    hold on
    if td <= 120
        rain = rainfall_intensity(1,:);
    elseif td <= 12*60
        rain = rainfall_intensity(2,:);
    elseif td <= 24*60
        rain = rainfall_intensity(3,:);
    else
        rain = rainfall_intensity(4,:);
    end
    bar(t_prime*td,rain,'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
    hold off
    xlabel('Elapsed Time [min]','Interpreter','latex','FontSize',12)
    ylabel('Rainfall Intensity [$\mathrm{mm.h^{-1}}$]','Interpreter','latex','FontSize',12);
    title('Huff Hyetograph','Interpreter','latex')
    legend('Huff Hyetograph','Location','Best','Interpreter','latex')
    grid on
    axis tight
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', 12)    
    grid on
    box on
    title(strcat('RP = ',num2str(RP,'%.2d'),' years'),'interpreter','latex')   
    exportgraphics(gcf,fullfile(folderName,'Intensity_Huff.pdf'),'ContentType','vector') 
% Outputs 
i = rain;
end
end
