%% plotter between dem fixed and raw data

% Specify the filename and path of the CSV file
results = fullfile(strcat(cd,"\",folderName),'Rating_Curve_Gauges.csv');
general_data = fullfile(strcat(cd,"\"),'Input_Data_sheets/General_Data_HydroPol2D.xlsx');

% Read the CSV files
simulated = readtable(results);
general_data = readtable(general_data);
date_begin = table2array(general_data(10,2));
date_end = table2array(general_data(11,2));
date_begin = datetime(datestr(date_begin+datenum('30-Dec-1899')));
date_end = datetime(datestr(date_end+datenum('30-Dec-1899')));
% % data_raw.(1)=datetime(data_raw.(1),'InputFormat','dd-MMM-yyyy HH:mm:ss');
% data_fix.(1)=datetime(data_fix.(1),'InputFormat','dd-MMM-yyyy HH:mm:ss');

% Reading the gauges codes to dowbload data
gauges_codes = table2array(general_data(2:end,36));
gauges_codes = gauges_codes(~isnan(gauges_codes));

T = [];
for i =1:length(gauges_codes)
   % request data (html)
   html = webread(strcat('http://telemetriaws1.ana.gov.br/ServiceANA.asmx/DadosHidrometeorologicos?CodEstacao=',string(gauges_codes(i)),'&DataInicio=', ...
       string(day(date_begin)),'/',string(month(date_begin)),'/',string(year(date_begin)),'&DataFim=', ...
       string(day(date_end)),'/',string(month(date_end)),'/',string(year(date_end))));
   % Gather each row of codes
   code = regexp(html,'<CodEstacao>(.*?)</CodEstacao>','tokens');
   code = cellfun(@(x) x{1}, code, 'UniformOutput', false);
   code = vertcat(code');
   % Gather each row of dates
   date = regexp(html,'<DataHora>(.*?)</DataHora>','tokens');
   date = cellfun(@(x) x{1}, date, 'UniformOutput', false);
   date = datetime(vertcat(date'));
   % Filling flow missing data, and gather each row
   html = regexprep(html,'<Vazao />','<Vazao> </Vazao>');
   flow = regexp(html,'<Vazao>(.*?)</Vazao>','tokens');
   flow = cellfun(@(x) x{1}, flow, 'UniformOutput', false);
   flow = cellfun(@str2double,flow');
   % Filling stage missing data, and gather each row
   html = regexprep(html,'<Nivel />','<Nivel> </Nivel>');
   stage = regexp(html,'<Nivel>(.*?)</Nivel>','tokens');
   stage = cellfun(@(x) x{1}, stage, 'UniformOutput', false);
   stage = cellfun(@str2double,stage');
   % Filling chuva missing data, and gather each row
   html = regexprep(html,'<Chuva />','<Chuva> </Chuva>');
   rainfall = regexp(html,'<Chuva>(.*?)</Chuva>','tokens');
   rainfall = cellfun(@(x) x{1}, rainfall, 'UniformOutput', false);
   rainfall = cellfun(@str2double,rainfall');
   % Create table
   T.(strcat('G',string(gauges_codes(i)))) = table(code, date, flow, stage, rainfall);
   strcat('Gauge data No.-', string(i), '-of-', string(length(gauges_codes)))
end

% Create the plot
figure;
set(gcf, 'Units', 'inches', 'Position', [1 1 16 8]);  

for i =1:length(gauges_codes)
    subplot(1,2,i)
    % plot observed data
    plot(T.(strcat('G',string(gauges_codes(i)))).date, T.(strcat('G',string(gauges_codes(i)))).flow, LineWidth=1,color='black');hold on
    % plot simulated data
    plot(simulated.(1), simulated.(i+1), LineWidth=1,color='#FF6666'); hold on
    % labels
    ylabel('Streamflow m$$^3$$/s','interpreter','latex', 'FontSize',10);
    title(string(gauges_codes(i)),'interpreter','latex', 'FontSize',10)
    % yyaxis right; set(gca,'ydir','reverse','ycolor','black');
    % bar(data.Var1, data.Var3,'FaceColor','#99CCFF','EdgeColor','#99CCFF','LineWidth',1.5)
    
    % Add labels and title
    xlabel('','interpreter','latex', 'FontSize',10);
    % ylabel({'Rainfall','Intensity (mm/h)'},'interpreter','latex','FontSize',14); ylim([0,max(data.Var3)*5]);
    % Adjust the plot as needed
    box on;
end
legend('Observed data','Simulated data','interpreter','latex', 'FontSize',10);

% Exporting results into PDF

exportgraphics(gcf,fullfile(strcat(cd,"\",folderName),'Plot_simulated_vs_Observed.pdf'),'ContentType', 'vector'); 
close all
