% Human Instability Module

% ---- Calculation of Human Stability --- %
if flags.flag_human_instability == 1
    if k == 1
        Human_Instability.max_risk = 0*DEM_raster.Z;
    end
    h = depths.d_t/1000;
    v = velocities.total_velocity;

    Human_Instability.F_person = Human_Instability.weight_person*Human_Instability.gravity; % N
    % Buyoance
    Human_Instability.F_buoy = (Human_Instability.width1_person*Human_Instability.width2_person)*h*Human_Instability.ro_water*Human_Instability.gravity; % N
    % Available Friction
    Human_Instability.available_friction = Human_Instability.mu*(max(Human_Instability.F_person - Human_Instability.F_buoy,0));
    % Hydrodynamic Force
    Human_Instability.hydro_force = max(1/2*(Human_Instability.Cd.*Human_Instability.ro_water.*Human_Instability.width1_person*h.*v.^2),0);
    % Risk Factor
    Human_Instability.risk_t = min(Human_Instability.hydro_force./Human_Instability.available_friction,1);
    Human_Instability.risk_t(idx_nan) = nan;
    Human_Instability.max_risk = max(Human_Instability.risk_t,Human_Instability.max_risk); % Max value

elseif flags.flag_human_instability == 2
    % Something has to be here
elseif flags.flag_human_instability == 3
    % Please add some meaning of what are these maps
    Human_Instability.risk_t_cm = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_c_m,Human_Instability.y_c_m,Human_Instability.w_c_m,Human_Instability.d_c_m);
    Human_Instability.risk_t_tm = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_t_m,Human_Instability.y_t_m,Human_Instability.w_t_m,Human_Instability.d_t_m);
    Human_Instability.risk_t_am = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_a_m,Human_Instability.y_a_m,Human_Instability.w_a_m,Human_Instability.d_a_m);
    Human_Instability.risk_t_om = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_o_m,Human_Instability.y_o_m,Human_Instability.w_o_m,Human_Instability.d_o_m);
    Human_Instability.risk_t_cf = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_c_f,Human_Instability.y_c_f,Human_Instability.w_c_f,Human_Instability.d_c_f);
    Human_Instability.risk_t_tf = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_t_f,Human_Instability.y_t_f,Human_Instability.w_t_f,Human_Instability.d_t_f);
    Human_Instability.risk_t_af = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_a_f,Human_Instability.y_a_f,Human_Instability.w_a_f,Human_Instability.d_a_f);
    Human_Instability.risk_t_of = Human_risk(flags.flag_human_instability,velocities.velocity_raster,depths.d_t./1000,Human_Instability.ro_water,Human_Instability.gravity,Human_Instability.mu,Human_Instability.Cd,Human_Instability.slope,Human_Instability.m_o_f,Human_Instability.y_o_f,Human_Instability.w_o_f,Human_Instability.d_o_f);
end