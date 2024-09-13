%% Plotting

load Workspace_Optimization.mat

% Number of effective generations
generations = 10;

% Plotting Optimization Results

for i = 1:generations
    filename = sprintf('gen_%04d.mat',i);
    data_generation(:,:,i) = load(filename);
end

generation_array = 1:1:generations;

index = 1;
for i = 1:2
    for j = 1:n_LULC
        if i == 1
            x_labels{index,1} = sprintf('$n_%d$',j);
        elseif i == 2
            x_labels{index,1} = sprintf('$h_{0,%d}$',j);
%         elseif i == 3
%             x_labels{index,1} = sprintf('$C_{1,%d}$',j);
%         elseif i == 4
%             x_labels{index,1} = sprintf('$C_{2,%d}$',j);
%         elseif i == 5
%             x_labels{index,1} = sprintf('$C_{3,%d}$',j);
%         elseif i == 6
%             x_labels{index,1} = sprintf('$C_{4,%d}$',j);
        end
        index = index + 1;
    end
end

for i = 1:3
    for j = 1:n_SOIL
        if i == 1
            x_labels{index,1} = sprintf('$k_{sat,%d}$',j);
        elseif i == 2
            x_labels{index,1} = sprintf('Soil Moisture %d',j);
        else
            x_labels{index,1} = sprintf('$psi_{%d}$',j);
        end
        index = index + 1;
    end
end
for i = 1:generations
    score_gen(i,:) = data_generation(:,:,i).Score_gen;
    idx_min_score = find(score_gen(i,:) == min(score_gen(i,:)),1,'first');
    idx_max_score = find(score_gen(i,:) == max(score_gen(i,:)),1,'first');    
    population_gen  = data_generation(:,:,i).Population_gen;
    best_individual_gen(i,:) = population_gen(idx_min_score,:);  
    worst_individual_gen(i,:) = population_gen(idx_max_score,:);        
    percentage(i,:) = (best_individual_gen(i,:) - x_known)./x_known*100;
end

% Saving Output Functions for best individuals
% Run HydroPol2D Routing Solver
for i = 1:generations 
    x_best = best_individual_gen(i,:); %  Decision Vector
    [~,Qmod_best(i,:),Cmod_best(i,:)] = HydroPol2D_Solver_with_x(x_best,flag_calibrate_wq,idx_imp,idx_lulc,idx_soil,n_LULC,n_SOIL,ADD_events,delta_p_obs,warmup_depths,warmup_I_0,n_events,time_observed,observed_flow,pollutant_concentration , alfa_albedo_input	, alfa_max	, alfa_min	, alfa_save	, avgtemp_stations	, B_t , Bmax	, Bmin	, C	 , C_3	, C_4	, Cd	, cell_area	, climatologic_spatial_duration	, col_outlet	, coordinate_x	, coordinate_y	, coordinates_stations	, d_0	, d_p	, date_begin	, delta_p_agg	, dem	, DEM_etp	, DEM_raster	, depth_tolerance	, elevation	, elevation_down_t	, elevation_left_t	, elevation_right_t	, elevation_up_t	, ETP	, ETP_save	, factor_cells	, flag_critical	, flag_D8	, flag_ETP	, flag_infiltration	, flag_inflow	, flag_rainfall	, flag_spatial_rainfall	, flag_timestep		, flag_waterquality	, flag_wq_model	, flow_acceleration	, flow_tolerance	, flows_cells	, G_stations	, gravity	, h_0	, h_0_fulldomain		, I_tot_end	, I_tot_end_cell, idx_nan	, idx_nan_5	, inflow	, inflow_cells	, k	, k_out	, Krs	, ksat_fulldomain	, last_record_maps	, lat	, mass_lost	, mass_outlet	, max_time_step	, maxtemp_stations	, min_Bt	, min_time_step	, mintemp_stations	, mu	, n_stream_gauges	, nx_max	, ny_max		, Out_Conc	, outlet_index	, outlet_index_fulldomain	, outlet_type	, P_conc	, psi_fulldomain	, rainfall_matrix	, rainfall_matrix_full_domain	, Resolution	, ro_water	, roughness	, roughness_fulldomain		, row_outlet	, slope_alfa	, slope_outlet		, spatial_domain	, t	, t_previous		, teta_i_fulldomain	, teta_sat	, teta_sat_fulldomain	, time_calculation_routing	, time_change_matrices	, time_change_records	, time_deltap	, time_ETP	, time_records	, time_save_previous	, time_step	, time_step_change	, time_step_increments	, time_step_model	, time_step_save	, tmin_wq	, Tot_Washed	, Tr	, u2_stations	, ur_stations	, v_threshold	, vel_down	, vel_left	, vel_right	, vel_up	, vol_outlet	, weight_person, width1_person, width2_person);
end

% Run HydroPol2D Routing Solver
for i = 1:generations 
    x_worst = worst_individual_gen(i,:); %  Decision Vector
    [~,Qmod_worst(i,:),Cmod_worst(i,:)] = HydroPol2D_Solver_with_x(x_worst,flag_calibrate_wq,idx_imp,idx_lulc,idx_soil,n_LULC,n_SOIL,ADD_events,delta_p_obs,warmup_depths,warmup_I_0,n_events,time_observed,observed_flow,pollutant_concentration , alfa_albedo_input	, alfa_max	, alfa_min	, alfa_save	, avgtemp_stations	, B_t , Bmax	, Bmin	, C	 , C_3	, C_4	, Cd	, cell_area	, climatologic_spatial_duration	, col_outlet	, coordinate_x	, coordinate_y	, coordinates_stations	, d_0	, d_p	, date_begin	, delta_p_agg	, dem	, DEM_etp	, DEM_raster	, depth_tolerance	, elevation	, elevation_down_t	, elevation_left_t	, elevation_right_t	, elevation_up_t	, ETP	, ETP_save	, factor_cells	, flag_critical	, flag_D8	, flag_ETP	, flag_infiltration	, flag_inflow	, flag_rainfall	, flag_spatial_rainfall	, flag_timestep		, flag_waterquality	, flag_wq_model	, flow_acceleration	, flow_tolerance	, flows_cells	, G_stations	, gravity	, h_0	, h_0_fulldomain		, I_tot_end	, I_tot_end_cell, idx_nan	, idx_nan_5	, inflow	, inflow_cells	, k	, k_out	, Krs	, ksat_fulldomain	, last_record_maps	, lat	, mass_lost	, mass_outlet	, max_time_step	, maxtemp_stations	, min_Bt	, min_time_step	, mintemp_stations	, mu	, n_stream_gauges	, nx_max	, ny_max		, Out_Conc	, outlet_index	, outlet_index_fulldomain	, outlet_type	, P_conc	, psi_fulldomain	, rainfall_matrix	, rainfall_matrix_full_domain	, Resolution	, ro_water	, roughness	, roughness_fulldomain		, row_outlet	, slope_alfa	, slope_outlet		, spatial_domain	, t	, t_previous		, teta_i_fulldomain	, teta_sat	, teta_sat_fulldomain	, time_calculation_routing	, time_change_matrices	, time_change_records	, time_deltap	, time_ETP	, time_records	, time_save_previous	, time_step	, time_step_change	, time_step_increments	, time_step_model	, time_step_save	, tmin_wq	, Tot_Washed	, Tr	, u2_stations	, ur_stations	, v_threshold	, vel_down	, vel_left	, vel_right	, vel_up	, vol_outlet	, weight_person, width1_person, width2_person);
end

%% Plotting
% Define the number of colors
numColors = generations;


colors = linspecer(generations);
% Plotting example (scatter plot)
figure;
hold on;


set(gcf,'units','inches','position',[3,3,6,4])
subplot(2,1,1)
for i = 1:numColors
    if i == numColors
        h1 = scatter(time_observed, Qmod_best(i,:),25,'*', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :));        
    else
        scatter(time_observed, Qmod_best(i,:),10,'*', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :));
    end
        hold on
    if i == numColors
        plot(time_observed, Qmod_best(i,:),'linewidth',.5,'color',colors(i, :),'LineStyle','-'); 
    end
end

for i = 1:numColors
    hold on
    if i == numColors
        h3 = scatter(time_observed, Qmod_worst(i,:),10,'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :));
    else
        scatter(time_observed, Qmod_worst(i,:),10,'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :));
    end
end

hold on
h2 = scatter(time_observed, observed_flow,25,'diamond', 'MarkerFaceColor','red', 'MarkerEdgeColor', 'red');

colormap(linspecer(40))
cbar = colorbar;
cbar.Ticks = [0 1];
cbar.TickLabels = [0 generations];
cbar.TickDirection = 'out';
cbar.TickLength = [0.02];
ylabel(cbar,'Generations','Interpreter','Latex','FontSize',12)
hold off

xlabel('Elapsed Time (min)','Interpreter','latex')
ylabel('Concentration (mg/L)','Interpreter','latex');
legend([h1 h3 h2],{'Best Individual','Worst Individual','Observed'},'interpreter','latex');
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
grid off
box on ;


% Insert Optimization Function Chart
xstart=.3;
xend=.7;
ystart=.7;
yend=.9;
axes('position',[xstart ystart xend-xstart yend-ystart ])
box on
% define range of subpart
best_scores = min(score_gen');
worst_scores = max(score_gen');
color_plots = linspecer(2);
plot(generation_array,best_scores,'LineWidth',1.5,'Marker','*','MarkerEdgeColor','black','Color',color_plots(1,:))% here i am plotting sub part of same figure. you plot another figure
grid on
hold on
plot(generation_array,worst_scores,'LineWidth',1.5,'Marker','+','MarkerEdgeColor','black','Color',color_plots(2,:))% here i am plotting sub part of same figure. you plot another figure

xlabel('Generations','FontSize',12,'Interpreter','latex')
ylabel('$O_f~[\mathrm{mg/L}]$','FontSize',12,'Interpreter','latex')
legend('Best Individual','Worst Individual','interpreter','latex');


subplot(2,1,2)
color_plot = linspecer(12);
for i = 1:12
    plot(generation_array,percentage(:,(i)),'color',color_plot(i,:),'marker','*');
    hold on
end
plot(generation_array,0*generation_array,'linewidth',1.5,'LineStyle','--','color','black')
legend(x_labels,'interpreter','latex')
% legend('$C_1$','$C_2$','$C_3$','$C_4$','interpreter','latex')
xlabel('Generations','Interpreter','latex')
ylabel('Relative Error (\%)','Interpreter','latex');
set(gca, 'TickLength', [0.02 0.01]);
set(gca,'Tickdir','out')
set(gca, 'FontName', 'Garamond', 'FontSize', 12)
grid off
box on ;