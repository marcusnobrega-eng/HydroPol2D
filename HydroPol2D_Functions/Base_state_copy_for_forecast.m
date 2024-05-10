function [bs_k,bs_t,bs_t2,bs_t_previous,bs_t_save,bs_t_save_2,bs_t_store,bs_time_calculation_routing,bs_time_step,bs_d_t,bs_I_t,flag_backup_state] = Base_state_copy_for_forecast(k,t,t2,t_previous,t_save,t_save_2,t_store,time_calculation_routing,time_step,d_t,I_t)
    bs_k = k;
    bs_t = t;
    bs_t2 = t2;
    bs_t_previous = t_previous;
    bs_t_save = t_save;
    bs_t_save_2 = t_save_2;
    bs_t_store = t_store;
    bs_time_calculation_routing = time_calculation_routing;
    bs_time_step = time_step;
    bs_d_t = d_t;
    bs_I_t = I_t;

    flag_backup_state = 1;
end