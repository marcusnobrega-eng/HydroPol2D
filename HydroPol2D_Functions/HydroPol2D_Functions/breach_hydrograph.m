% Instantaneous Dam-Breach Hydrograph
% Developer: Marcus Nobrega
% Goal - Derive a function to determine the inflow hydrograph as a function
% of the Dam dimensions

function [time_min,flow_m3_s] = breach_hydrograph(h,W,B,Breach,dt,tf,flag_plot,ID,Name)

% Calculating tp
tp = 15/8*W/sqrt(9.81*h); % sec

% Calculating Initial Volume
S0 = h*W*B; % m3

% Calculating Peak Flow
% Qp = 8/27*B*sqrt(9.81)*h^(3/2); % m3/s
Qp = 8/27*Breach*sqrt(9.81)*h^(3/2); % m3/s

% Calculating storage at the end of peak time
Sf = S0 - Qp*tp;

% Discretizing time-domain
time = [0:dt:tf]; % seconds

% Hs
hs = h - Qp*tp/(W*B);

h0 = h;

for i = 1:length(time)
    current_time = (i-1)*dt; % sec
    if current_time <= tp
        flow_m3_s(i,1) = Qp; % m3/s
        depth_m(i,1) = h0 - Qp*dt/(W*B);
    else
        flow_m3_s(i,1) = (-1/2*(B*sqrt(9.81))^(2/3)*h*(tp - current_time)/(Sf) + Qp^(-1/3))^(-3); % m3/s
        depth_m(i,1) = h0 - mean(flow_m3_s(i,1) + flow_m3_s(i-1,1))*dt/(W*B);
    end
    h0 = depth_m(i,1);
end

% Time Minutes
time_min = time'/60; % minutes

% Plotting

if flag_plot == 1
    close all
    set(gcf,'units','inches','position',[3,3,2.5,2.5])
    plot(time_min,flow_m3_s,'Color','black','LineStyle','--','LineWidth',2);
    xlabel('Elapsed time [min]','Interpreter','latex','FontSize',10);
    ylabel('Flow discharge [$\mathrm{m^3s^{-1}}]$','Interpreter','latex','FontSize',10);
    hold on
    set(gca,'TickDir','out')
    yyaxis right
    set(gca,'ycolor','black')
    ylabel('Water Depth [m]','Interpreter','latex','FontSize',10);
    plot(time_min,depth_m,'Color','black','LineStyle','-','LineWidth',2);
    title(strcat(Name,{' '},'-',{' '},ID),'Interpreter','latex','FontSize',12,'Color','black')
    set(gca,'FontName','Garamond')
    set(gca,'TickDir','out')


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

    label_plot = strcat(ID,'_',Name,'_','Breach_Hydrograph.pdf');

    % Export graphics
    exportgraphics(gcf,fullfile(folderName,label_plot),'ContentType','vector')

    label_plot = strcat(ID,'_',Name,'_','Breach_Hydrograph.m');

    saveas(gcf,fullfile(folderName,label_plot));

    close all

end
end