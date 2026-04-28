function [t,i,P,idf] = alternated_blocks(td,dt,K,a,b,c,RP,flag_plot,path)
% %------- Alternated Blocks Method ------- %
% Developer: Marcus Nobrega, Ph.D
% Updated to save outputs in the user-defined "path"
%
% Goal:
%   Create a hyetograph from an IDF curve using the Alternated Blocks Method
%
% Inputs:
%   td        : Rainfall duration (min)
%   dt        : Time-step (min)
%   K,a,b,c   : IDF parameters such that
%               i(mm/hr) = K*(RP^a)/((b+td)^c)
%   RP        : Return period (years)
%   flag_plot : Boolean [0,1]
%               1 = plot and export results
%               0 = do not plot
%   path      : Folder path where the output file will be saved
%
% Outputs:
%   t   : Vector with time from dt to td, with dt intervals (min)
%   i   : Alternated-block rainfall intensity for each t (mm/h)
%   P   : Cumulative rainfall depth (mm)
%   idf : IDF rainfall intensity for each t (mm/h)
%
% Example:
%   [t,i,P,idf] = alternated_blocks(60,5,819,0.138,10,0.75,10,1,'C:\Results')
%
% Notes:
%   - The exported figure will be saved as:
%       <path>\Alternated_Blocks.pdf
%   - If the folder in "path" does not exist, it will be created.
%   - The function linspecer must be available in your MATLAB path if used.

%% 0.0 - Check input path
if nargin < 9 || isempty(path)
    path = pwd; % Use current folder if no path is provided
end

% Create the output folder if it does not exist
if ~exist(path, 'dir')
    mkdir(path);
    disp(['Folder created successfully: ', path]);
end

%% 1.0 - Set duration
t = dt:dt:td;              % duration in min
nsteps = length(t);        % number of time-steps

%% 2.0 - Calculate intensity
idf = K.*(RP^a)./((b + t).^c);   % mm/hr
P_cum = (idf .* t / 60);         % cumulative rainfall depth in mm

delta_p = zeros(1,length(P_cum));
for k = 1:length(P_cum)
    if k == 1
        delta_p(k) = P_cum(k);
    else
        delta_p(k) = P_cum(k) - P_cum(k-1); % Delta precipitation in mm
    end
end

P = cumsum(delta_p);                  % cumulative rainfall depth (mm)
intensity_alternated = delta_p/(dt/60); % mm/h

%% 3.0 - Alternated Blocks Method
i = zeros(1,nsteps);

if mod(nsteps,2) == 1   % ODD number of steps
    for k = 1:floor(nsteps/2)
        z1 = nsteps - k*2 + 1;
        z2 = nsteps - k*2;
        i(nsteps - k) = intensity_alternated(z2);
        i(k) = intensity_alternated(z1);
    end
    i(nsteps) = min(intensity_alternated);

else                    % EVEN number of steps
    for k = 1:(floor(nsteps/2)-1)
        z1 = nsteps - k*2 + 1;
        z2 = nsteps - k*2;
        i(nsteps - k) = intensity_alternated(z2);
        i(k) = intensity_alternated(z1);
    end

    i(nsteps-1) = intensity_alternated(nsteps-2);
    i(nsteps)   = min(intensity_alternated);

    mid = floor((td/2)/dt); % middle of the duration
    i(mid) = max(idf);
end

i = i(:);
t = t(:);
idf = idf(:);
P = P(:);

%% Colors
colors = linspecer(2);

%% 4.0 - Plot and export
if flag_plot == 1
    close all
    figure('Units','inches','Position',[3,3,6.5,4]);

    plot(t,i,'Color',colors(1,:),'LineWidth',2, ...
        'Marker','*','MarkerSize',5,'MarkerEdgeColor','black');
    hold on

    plot(t,idf,'Color',colors(2,:),'LineWidth',2, ...
        'Marker','^','MarkerSize',5,'MarkerEdgeColor','black', ...
        'LineStyle','--');

    xlabel('Elapsed Time [min]','Interpreter','latex','FontSize',12)
    ylabel('Rainfall Intensity [$\mathrm{mm.h^{-1}}$]', ...
        'Interpreter','latex','FontSize',12);

    grid on
    axis tight

    font_size = 12;
    set(gca,'TickLength',[0.015 0.01]);
    set(gca,'TickDir','out');
    set(gca,'FontName','Garamond','FontSize',font_size);

    legend('Alternated Blocks','IDF curve','Interpreter','latex','Location','best');
    title(['RP = ', num2str(RP,'%.2f'), ' years'],'Interpreter','latex');

    % Export to the user-defined path
    exportFile = fullfile(path,'Alternated_Blocks.pdf');
    exportgraphics(gcf, exportFile, 'ContentType', 'vector');

    disp(['Figure exported successfully to: ', exportFile]);
end

end