%% Plot Results of Calibration and Validation

input_data = readtable("Resultados_Calibracao_Validacao.xlsx");
input_data = table2array(input_data);
n_events = 8;
for i = 1:n_events
    indicators(i,1:4) = input_data(1,((7*(i-1)+1):(7*(i-1)+1)+3));
end

for i = 1:n_events
    % time, modeled, observed
    results(:,(i-1)*3+1:(i-1)*3+3) = input_data(:,((i-1)*7 + 5):((i-1)*7 + 7));
end

%% ScatterPlot
sz = 3;
for i = 1:n_events
    modeled = results(:,(i-1)*3 + 2);
    observed = results(:,(i-1)*3 + 3);
%     scatter(modeled,observed,'.k','MarkerSize',8);
    plot(modeled,observed,'.k','MarkerSize',8);
    hold on
end
axis equal
axis tight
box on
max_value = 4*10^4;
xlim([0 max_value])
ylim([0 max_value])
% 45 degree line
x_45 = [0:max_value];
y_45 = x_45;
plot(x_45,y_45,'Color','black','LineStyle','--','LineWidth',1.5)
% axis properties
h = gca;
h.TickLength = [0.025 0.010];
h.XMinorTick = 'on' ; h.YMinorTick = 'on';
xlabel('Modeled TSS Concentration (mg/L)','interpreter','latex')
ylabel('Observed TSS Concentration (mg/L)','interpreter','latex')
legend('Event 1','Event 2','Event 3','Event 4','Event 4','Event 5','Event 6', 'Event 7','Event 8','interpreter','latex')

%% Pollutographs
figure(2)
set(gcf,'units','inches','position',[2,0,6.5,6])

c_1 = [0,0,204]/256;
c_2 = [0,255,255]/256;
c_3 = [250,6,6]/256;
c_4 = [55,155,213]/256;
c_5 = [0,0,0]/256;

for i = 1:n_events
    h = subplot(4,2,i);
    modeled = results(:,(i-1)*3 + 2);
    observed = results(:,(i-1)*3 + 3);
    time = results(:,(i-1)*3 + 1);
    plot(time,modeled,'.r','MarkerSize',8);
    hold on
    plot(time,observed,'.black','MarkerSize',8); 
    h.FontName = 'Garamond';
    h.FontSize = 12;
    h.TickDir = 'out';
    h.TickLength = [0.025 0.015];
    h.XTick = [0 2 4 6 8 10];
end
ylabel('TSS Concentration (mg/L)','Interpreter','Latex')
xlabel('Elapsed Time (min)','Interpreter','Latex')
legend('Modeled','Observed');