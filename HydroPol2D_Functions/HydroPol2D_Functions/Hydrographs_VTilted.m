% Plotting GSSHA x HydroPol2D

hydrographs_data = readtable('HydroPol2D_GSSHA_Comparison_VTILTED.xlsx','Sheet','New_Data');

time_steps_data = readtable('HydroPol2D_GSSHA_Comparison_VTILTED.xlsx','Sheet','Time-steps');

Labels = {'$\Delta t = 60~$sec','$\Delta t = 50~$sec', '$\Delta t = 40~$sec','$\Delta t = 30~$sec','$\Delta t = 20~$sec','$\Delta t = 5~$sec','$\Delta t = 2~$sec',' $\Delta t = 1~$sec', '$\Delta t = 0.5~$sec', '$\Delta t = 0.1~$sec'};


%% Plots 

% Data Observed GSSHA
time_GSSHA = table2array(hydrographs_data(:,1));
flow_GSSHA = table2array(hydrographs_data(:,2));

time_Adaptive_Diffusive = table2array(hydrographs_data(:,3));
flow_Adaptive_Diffusive = table2array(hydrographs_data(:,4));

time_Adaptive_CA = table2array(hydrographs_data(:,5));
flow_Adaptive_CA = table2array(hydrographs_data(:,6));

set(gcf,'units','inches','position',[3,0,6.5,6])
font_size = 12; 
lw = 1.5;
mk_type = {'*','square','o','v','d','h'};
ln_type = {'--','-.','-',':','-','--'};
mk_size = 4;
subplot(2,2,1)
colors_plot = linspecer(5);
marker_colors = linspecer(5);
% Unstable Simulations - Diffusive
for i = 1:5
    begin = 7;
    offset = (i-1)*4;
    time_plot = table2array(hydrographs_data(:,begin  + offset));
    data_plot = table2array(hydrographs_data(:,begin  + offset + 1));
    plot(time_plot,data_plot,'linewidth',lw,'marker',mk_type(i),'MarkerSize',mk_size,'LineStyle',ln_type(i),'MarkerEdgeColor',colors_plot(i,:),'MarkerFaceColor',colors_plot(i,:),'color',colors_plot(i,:));
    hold on    
end
% Plotting GSSHA
s1 = scatter(time_GSSHA,flow_GSSHA,25,'d','MarkerEdgeColor',[0,0,128]/255,'MarkerFaceColor',[0,0,128]/255,'LineWidth',1.5);
s1.MarkerFaceAlpha = 'flat';
xlabel('Elapsed Time (min)','Interpreter','latex','FontSize',font_size)
ylabel('Flow Discharge ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
set(gca, 'TickLength', [0.01 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
% title('Flood Depths','Interpreter','Latex');

subplot(2,2,2)
% Unstabe CA
for i = 1:5
    begin = 9;
    offset = (i-1)*4;
    time_plot = table2array(hydrographs_data(:,begin  + offset));
    data_plot = table2array(hydrographs_data(:,begin  + offset + 1));
    plot(time_plot,data_plot,'linewidth',lw,'marker',mk_type(i),'MarkerSize',mk_size,'LineStyle',ln_type(i),'MarkerEdgeColor',colors_plot(i,:),'MarkerFaceColor',colors_plot(i,:),'color',colors_plot(i,:));
    hold on
end
s2 = scatter(time_GSSHA,flow_GSSHA,25,'d','MarkerEdgeColor',[0,0,128]/255,'MarkerFaceColor',[0,0,128]/255,'LineWidth',1.5);xlabel('Elapsed Time (min)','Interpreter','latex','FontSize',font_size)
s2.MarkerFaceAlpha = 'flat';
ylabel('Flow Discharge ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
set(gca, 'TickLength', [0.01 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
legend(Labels{1:5},'interpreter','latex')

% title('Flood Depths','Interpreter','Latex');

color_plot = linspecer(7);
color_plot = color_plot(6:7,:);

subplot(2,2,[3])
plot(time_Adaptive_Diffusive,flow_Adaptive_Diffusive,'linewidth',lw,'marker',mk_type(1),'MarkerSize',mk_size,'LineStyle',ln_type(1),'MarkerEdgeColor',colors_plot(1,:),'MarkerFaceColor',colors_plot(1,:),'color',colors_plot(1,:));
hold on
plot(time_Adaptive_CA,flow_Adaptive_CA,'linewidth',lw,'marker',mk_type(2),'MarkerSize',mk_size,'LineStyle',ln_type(2),'MarkerEdgeColor',colors_plot(2,:),'MarkerFaceColor',colors_plot(2,:),'color',colors_plot(2,:));
hold on;
s3 = scatter(time_GSSHA,flow_GSSHA,25,'d','MarkerEdgeColor',[0,0,128]/255,'MarkerFaceColor',[0,0,128]/255,'LineWidth',1.5);
s3.MarkerFaceAlpha = 'flat';
xlabel('Elapsed Time (min)','Interpreter','latex','FontSize',font_size)
ylabel('Flow Discharge ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
set(gca, 'TickLength', [0.01 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
% title('Flood Depths','Interpreter','Latex');


subplot(2,2,[4])
scatter(flow_Adaptive_Diffusive,flow_Adaptive_CA,25,'o','filled','black');

for i = 1:5
    scatter(flow_Adaptive_Diffusive,flow_Adaptive_CA,25,'o','filled','black');
end

xlabel('HydroPol2D (a) ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
ylabel('HydroPol2D (b) Discharge ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
% 45 deg
x_max = max(max(flow_Adaptive_Diffusive,flow_Adaptive_CA));
x_45 = [0,x_max]; y_45 = [0,x_max];
hold on
plot(x_45,y_45,"Color",'black','LineStyle','--','LineWidth',2);
set(gca, 'TickLength', [0.01 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
% title('Flood Depths','Interpreter','Latex');
set(gca, 'TickLength', [0.01 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
box on

% Insert Adaptive Time-step Graphs% 
xstart=.5;
xend=.7;
ystart=.1;
yend=.4;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
% Read data
elapsed_time_diffusive = table2array(time_steps_data(:,1));
time_diffusive = table2array(time_steps_data(:,2));

elapsed_time_CA = table2array(time_steps_data(:,3));
time__CA = table2array(time_steps_data(:,4));

colorplots = linspecer(12);
colorplots = colorplots(10:12,:);

plot(elapsed_time_diffusive,time_diffusive,'LineWidth',1.5,'Linestyle','--','Color',colorplots(1,:))% here i am plotting sub part of same figure. you plot another figure
grid on
hold on
plot(elapsed_time_CA,time__CA,'LineWidth',1.5,'Color',colorplots(2,:))% here i am plotting sub part of same figure. you plot another figure
xlabel('Elapsed Time (min)','FontSize',12,'Interpreter','latex')
ylabel('$\Delta t(t)$','FontSize',12,'Interpreter','latex')
set(gca, 'TickLength', [0.01 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
legend('HydroPol2D (a)','HydroPol2D (b)','interpreter','latex');

%% Calculating Functions
% Nash-Suctclife-Efficiency and Objective function
NSE = 1 - sum((obs_discharge - Qmod).^2)/(sum((obs_discharge - mean(obs_discharge)).^2)); % We want to maximize it

n_elements = length(obs_discharge);
RMSE = sqrt(sum((obs_discharge - Qmod).^2/n_elements));

r2 = corrcoef(obs_discharge,Qmod);


% subplot(2,2,3)
% % Stable Diffusive
% color_plot = linspecer(5 + 4);
% color_plot = color_plot(6:9,:);
% for i = 1:4
%     begin = 27;
%     offset = (i-1)*4;
%     time_plot = table2array(hydrographs_data(:,begin  + offset));
%     data_plot = table2array(hydrographs_data(:,begin  + offset + 1));
%     plot(time_plot,data_plot,'linewidth',lw,'marker',mk_type(i),'MarkerSize',mk_size,'LineStyle',ln_type(i),'MarkerEdgeColor','red','MarkerFaceColor','red','color',colors_plot(i,:));
%     hold on
% end
% scatter(time_GSSHA,flow_GSSHA,25,'d','filled','color',[0,0,128]/255);
% xlabel('Elapsed Time (min)','Interpreter','latex','FontSize',font_size)
% ylabel('Flow Discharge ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
% set(gca, 'TickLength', [0.01 0.01]);
% set(gca,'Tickdir','out')
% set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
% % title('Flood Depths','Interpreter','Latex');
% 
% 
% subplot(2,2,4)
% % Stable CA
% color_plot = linspecer(5 + 4);
% color_plot = color_plot(6:9,:);
% for i = 1:4
%     begin = 29;
%     offset = (i-1)*4;
%     time_plot = table2array(hydrographs_data(:,begin  + offset));
%     data_plot = table2array(hydrographs_data(:,begin  + offset + 1));
%     plot(time_plot,data_plot,'linewidth',lw,'marker',mk_type(i),'MarkerSize',mk_size,'LineStyle',ln_type(i),'MarkerEdgeColor','red','MarkerFaceColor','red','color',colors_plot(i,:));
%     hold on
% end
% scatter(time_GSSHA,flow_GSSHA,25,'d','filled','color',[0,0,128]/255);
% xlabel('Elapsed Time (min)','Interpreter','latex','FontSize',font_size)
% ylabel('Flow Discharge ($\mathrm{m^3/s}$)','Interpreter','latex','FontSize',font_size)
% set(gca, 'TickLength', [0.01 0.01]);
% set(gca,'Tickdir','out')
% set(gca, 'FontName', 'Garamond', 'FontSize', font_size)
% subplot(2,2,5)