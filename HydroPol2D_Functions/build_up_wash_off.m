function [B_t,P_conc,Out_Conc,tmin_wq,tot_W_out,mass_lost,Tot_Washed] = build_up_wash_off(C_3,C_4,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow,B_t,time_step,nx_max,ny_max,cell_area,outlet_index,idx_nan_5,flag_wq_model,mass_lost,Tot_Washed,Bmin,Bmax,min_Bt)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 30 October, 2026                  %
%          Goal: Solve the Build-up and Wash-off Model            %
%
% The wash-off W = C_3 Q^C_4 B [kg/h] (or the rating curve) leaves each cell
% through its four faces and, on outlet cells, the outlet. The step is split
% into n = ceil(max(W_tot/B) dt) equal sub-steps that sum to dt, and in each
% sub-step a cell's outflow is scaled so it never exceeds the mass it holds,
% as HydroPol2D-Python's mass_based_washoff_step does. The pollutant
% therefore sub-cycles inside the hydraulic step: tmin_wq is the full step
% and no longer limits the model time step.

%% Pollutant Rate Flux
% 3 dim array with pages 1, 2, 3, 4, and 5 being left, right, up, down, and
% outlet flow rates in mm/hr
q_out_t = cat(3,qout_left_t,qout_right_t,qout_up_t,qout_down_t,outlet_flow); % mm/h
dt_h = time_step/60; % hours
Br_kg = Bmin/1000*cell_area; % kg
Bm_kg = Bmax/1000*cell_area; % kg
valid = ~idx_nan_5(:,:,1);
not_a_number = isnan(B_t);
B_t(~isfinite(B_t)) = 0;
B_t = max(B_t,0);
outlet_mask = logical(outlet_index) & (outlet_flow > 0);
if any(outlet_mask(:))
    outlet_water_volume_L = max(nansum(outlet_flow(outlet_mask)) * cell_area * time_step/60, 0);
else
    outlet_water_volume_L = 0;
end

% Sub-steps from the gross outflow rate at the start of the step
W_out_t = washoff_rates(B_t);
rate_ratio = zeros(size(B_t),'like',B_t);
has_mass = B_t > 0 & valid;
tot_rate = sum(W_out_t,3);
rate_ratio(has_mass) = tot_rate(has_mass)./B_t(has_mass); % 1/h
steps = max(1, ceil(max(rate_ratio(:))*dt_h));
dt_sub = dt_h/steps; % hours

outflow_mass = zeros(size(B_t),'like',B_t); % kg that left each cell
outlet_mass = 0; % kg captured by the outlet sink during the step
for i = 1:steps
    W_out_t = washoff_rates(B_t);
    tot_W_sub = sum(W_out_t,3); % kg/hr
    % A cell cannot release more than it holds in one sub-step
    scale = ones(size(B_t),'like',B_t);
    over = tot_W_sub*dt_sub > B_t;
    scale(over) = B_t(over)./max(tot_W_sub(over)*dt_sub,1e-30);
    W_out_t = W_out_t.*scale;
    tot_W_sub = sum(W_out_t,3);
    % Inflows from the neighbours' faces
    W_in_t = [zeros(ny_max,1),W_out_t(:,1:(nx_max-1),2)] + ... % from the left cell's right face
        [W_out_t(:,2:nx_max,1), zeros(ny_max,1)] + ...           % from the right cell's left face
        [zeros(1,nx_max); W_out_t(1:(ny_max-1),:,4)] + ...       % from the upper cell's down face
        [W_out_t(2:ny_max,:,3); zeros(1,nx_max)];                % from the lower cell's up face
    W_in_t(~valid) = 0;
    outflow_mass = outflow_mass + tot_W_sub*dt_sub;
    B_t = max(B_t - (tot_W_sub - W_in_t)*dt_sub, 0);
    if any(outlet_mask(:))
        outlet_capture = W_out_t(:,:,5);
        outlet_capture = outlet_capture(outlet_mask) * dt_sub;
        outlet_mass = outlet_mass + sum(outlet_capture(isfinite(outlet_capture) & outlet_capture > 0));
    end
end
% Average Pollutant Flux in the Time-Step
tot_W_out = outflow_mass/dt_h; % kg/hr
tot_W_out(~valid) = 0;
B_t(not_a_number) = nan;
tmin_wq = time_step*60; % s: the pollutant sub-cycles inside the step

% Concentrations
B_t(abs(B_t) < 1e-12) = 0;
tot_q_out = sum(q_out_t,3); % mm/hr
tot_q_out(~valid) = 0;
tot_q_out = max(tot_q_out,10); % 10 mm/hr
P_conc = max(10^6*(tot_W_out)./(tot_q_out*cell_area),0); % mg/L
if outlet_water_volume_L > 0
    Out_Conc = max(1e6*outlet_mass/outlet_water_volume_L,0); % mg/L
else
    Out_Conc = 0;
end
Tot_Washed = tot_W_out*time_step/60 + Tot_Washed; % (kg) Actually, total received

    function W = washoff_rates(B)
        % Directional wash-off [kg/h]; no flux from cells below min_Bt.
        if flag_wq_model == 1 % Rating Curve
            B_extra = B;
            B_extra(B*1000/cell_area <= Bmin) = 0;
            B_extra(B*1000/cell_area >= Bmax) = Bm_kg;
            f_bt = (1 + (max(B_extra - Br_kg,0)));
            W = C_3.*(q_out_t/1000/3600*cell_area).^(C_4).*f_bt; % kg/hr (rating curve)
        else
            W = C_3.*(q_out_t/1000/3600*cell_area).^(C_4).*B; % kg/hr (mass based)
        end
        below = (B*1000/cell_area < min_Bt) | ~isfinite(B) | ~valid;
        W(repmat(below,1,1,5)) = 0;
        W(idx_nan_5) = 0;
        W(~isfinite(W)) = 0;
    end
end
