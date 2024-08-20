function [t,i,P,idf] = alternated_blocks(td,dt,K,a,b,c,RP,flag_plot)
% %------- Alternated Blocks Method ------- %
% Developer: Marcus Nobrega, Ph.D
% Goal: Create a hyetograph from an IDF curve
%
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
% P: cumulative rainfall intensity
% idf: IDF rainfall intensity for each t
%
% Example: Determine the alternated blocks hyetograph of a rainfall with 60
% min duration, RP = 10 years, K = 819, a = 0.138 b = 10, and c = 0.75 and plot the
% hyetograph with 5-min intervals
%
% [t,i,P,idf] = alternated_blocks(60,5,819,0.138,10,0.75,10,1)
% Problems with the function? Contact me at marcusnobrega.engcivil@gmail.com

%% 1.0 - Set duration
t = dt:dt:td; % duration in min
nsteps = length(t); % number of time-steps
%% 2.0 - Calculate intensity
idf = K.*(RP^a)./((b + t).^c); % mm/hr
P_cum = (idf.*t/60); % mm
delta_p = zeros(1,length(P_cum));
for k = 1:(length(P_cum))
    if k == 1
        delta_p(k) = P_cum(k);
    else
        delta_p(k) = P_cum(k) - P_cum(k-1); % Delta precipitation in mm
    end
end
P = cumsum(delta_p); % Cumulative volume of rainfall
intensity_alternated = delta_p/(dt/60); % mm/h

%% 3.0 - Alternated Blocks Method
i = zeros(1,nsteps);
if mod(nsteps,2) == 1 % ODD
    for k = 1:(floor(nsteps/2))
        z1 = nsteps - k*2 + 1;
        z2 = nsteps - k*2;
        i(nsteps - k) = intensity_alternated(z2);
        i(k) = intensity_alternated(z1);
    end
    i(nsteps) = min(intensity_alternated);
else % EVEN number of steps
    for k = 1:(floor(nsteps/2)-1)
        z1 = nsteps - k*2 + 1;
        z2 = nsteps - k*2;
        i(nsteps - k) = intensity_alternated(z2);
        i(k) = intensity_alternated(z1);
    end
    i(nsteps-1) = intensity_alternated(nsteps-2);
    i(nsteps) = min(intensity_alternated);
    mid = floor((td/2)/dt); % middle of the duration
    i(1,mid) = max(idf);
end
i = i';
t = t';

colors = linspecer(2);


if flag_plot == 1
%% Creating Modeling Results Folder
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

    % Plots
    close all
    set(gcf,'units','inches','position',[3,3,6.5,4])    
    plot(t,i,'color',colors(1,:),'linewidth',2,'Marker','*','MarkerSize',5,'MarkerEdgeColor','black');
    hold on
    plot(t,idf,'color',colors(2,:),'linewidth',2,'Marker','^','MarkerSize',5,'MarkerEdgeColor','black','LineStyle','--');
    xlabel('Elapsed Time [min]','Interpreter','latex','FontSize',12)
    ylabel('Rainfall Intensity [$\mathrm{mm.h^{-1}}$]','Interpreter','latex','FontSize',12);
    grid on
    axis tight
    font_size = 12;
    set(gca, 'TickLength', [0.015 0.01]);
    set(gca,'Tickdir','out')
    set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
    legend('Alternated Blocks','IDF curve','interpteter','latex')
    title(strcat('RP = ',num2str(RP,'%.2d'),' years'),'interpreter','latex')
    exportgraphics(gcf,fullfile(folderName,'Alternated_Blocks_Hyetograph.pdf'),'ContentType','vector')    
end
end
