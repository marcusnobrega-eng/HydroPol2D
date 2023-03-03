###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                    September 2021                               #
#                                                                 #
#       Last update : 7 July, 2021                                #
###################################################################

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# -------- DESCRIPTION ----------------
# This function estimates the transferred volume from a cell to the
# neighbour cells using a weighting approach in terms of the available
# volume in the neighbour cells
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np

def CA_Routing(elevation_cell,  elevation_left_t,  elevation_right_t,  elevation_up_t,  elevation_down_t,  d_t_cell,  d_left_cell,  d_right_cell,  d_up_cell,  d_down_cell,  roughness_cell,  cell_area,  time_step,  h_0_cell,  Resolution,  I_tot_end_cell,  outlet_index,  outlet_type,  slope_outlet,  row_outlet,  col_outlet,  idx_nan,  flag_critical):
    # ----------- Finding Nans or Values Below a Threshold ------------ #
    d_t_min = 1e-6  # m
    idx_2 = d_t_cell < d_t_min*1000
    idx_3 = idx_nan
    idx2 = np.logical_or(idx_2, idx_3)

    # ----------- Adding minimum slope to do calculation ------------ #
    h_min = d_t_min  # m
    I_tot_end_previous = I_tot_end_cell
    # ----------- Notation ------------ #
    #  <-I-> (left right up down) = [left ; right ; up ; down]
    # ----------- Cell Depth and Water surface Elevation ------------ #
    depth_cell = np.array(d_t_cell)/1000  # meters
    wse_cell = elevation_cell + depth_cell

    # ----------- Depth Differences for All Cells ------------ #
    delta_h = np.zeros((5, np.size(elevation_cell[:, 0]), np.size(elevation_cell[0, :])))
    delta_h[0, :, :] = wse_cell - elevation_left_t - np.array(d_left_cell)/1000  # left
    delta_h[1, :, :] = wse_cell - elevation_right_t - np.array(d_right_cell)/1000  # right
    delta_h[2, :, :] = wse_cell - elevation_up_t - np.array(d_up_cell)/1000  # up
    delta_h[3, :, :] = wse_cell - elevation_down_t - np.array(d_down_cell)/1000  # down

    # %%% CHECKED 9/28/2021 %%%% #
    # ----------- Outlet Calculations ------------ #
    if outlet_type == 1:
        S_0 = slope_outlet  # Normal slope
        # V = h * area
        delta_h[4, :, :] = S_0*Resolution*outlet_index  # normal slope depth difference
        # %%% CHECKED 9/28/2021 %%%% #
    else:
        S_0 = np.zeros((np.shape(elevation_cell)))
        for i in range(len(row_outlet)):
            # Checking left, right, up, down
            row = int(row_outlet[i])  # Row of outlet
            col = int(col_outlet[i])  # Col of outlet
            S_0[row, col] = np.array(np.power(np.array(depth_cell[row, col]),-1/6))*np.array(np.sqrt(9.81))*np.array(roughness_cell[row, col])  # Critical slope
        delta_h[4, :, :] = S_0*Resolution

    # ----------- Correcting delta_h ------------ #
    delta_h[np.isnan(delta_h)] = 0  # TESTING
    delta_h[np.isinf(delta_h)] = 0  # TESTING

    # ----------- Boundary Conditions ------------ #
    # These only applied when we are modeling a rectangular grid watershed
    delta_h[0, :, 0] = 0  # left
    delta_h[1, :, -1] = 0  # Right
    delta_h[2, 0, :] = 0  # Up
    delta_h[3, -1, :] = 0  # Down
    idx = delta_h < h_min
    delta_h[idx] = 0

    # ----------- Available volumes ------------ #
    Vol = cell_area*delta_h  # 3-D Array with pages following direction notation
    # Max h
    maxh = np.max(delta_h, axis=0)
    # %%% CHECKED 9/28/2021 %%%% #
    # ----------- Minimum and Ma volumes ------------ #
    # Minimum
    c = 1000000  # Large number
    Vol_nan = Vol.copy()
    Vol_nan[idx] = c
    Vol_min_matrix = np.min(Vol_nan, axis=0)
    # Total Volume
    Vol_tot = np.sum(Vol, axis=0)

    # ----------- Weights ------------ #
    weight = np.array(Vol)/np.array(Vol_tot + Vol_min_matrix)
    weight_max = np.max(weight, axis=0)

    # ----------- Velocity Calculation ------------ #
    # Velocity to the cell with the highest gradient
    if flag_critical == 1:
        v_m = np.minimum(np.sqrt(9.81*depth_cell), np.array(1/roughness_cell) * np.array(np.power(np.array(np.maximum(depth_cell - h_0_cell/1000, 0)), 2/3))*np.array(np.sqrt(maxh/Resolution)))
    else:
        v_m = np.array(1/roughness_cell) * np.array(np.power(np.array(np.maximum(depth_cell - h_0_cell/1000, 0)), 2/3))*np.array(np.sqrt(maxh/Resolution))

    # ----------- Inter-cell Volume Calculation ------------ #
    # Estimated Volume that goes to the cell with the highest gradient
    # I_m = np.array(v_m) * depth-cell * Resolution * time_step * 60  # Velocity * area * Time_step (m3)

    # ----------- Total outflow volume ------------ #
    # I_tot_end = np.minimum(depth_cell*cell_area, I_m/np.max(weight), Vol_min_matrix + I_tot_begin)
    I_tot_end_cell = np.minimum(np.array(depth_cell) * cell_area, np.array(v_m)*np.divide(np.array(depth_cell), np.array(weight_max), where=weight_max != 0)*(Resolution*time_step*60))
    I_tot_end_cell = np.minimum(I_tot_end_cell, I_tot_end_previous + Vol_min_matrix)
    I_tot_end_cell[np.isnan(I_tot_end_cell)] = 0  # Numerical constraint at I_tot
    I_tot_end_cell[idx2] = 0  # if the depth is smaller than dmin ATTENTION HERE

    # ----------- Outflow Volume to each direction ------------ #
    I = np.array(weight)*np.array(I_tot_end_cell)  # m3
    I_tot_end_cell = np.sum(I, axis=0)  # Volume that leaves the cell in m^3
    qout = np.array(I)/(time_step*60)/(cell_area)*1000*3600  # mm/hr
    qout_left = qout[0, :, :]
    qout_right = qout[1, :, :]
    qout_up = qout[2, :, :]
    qout_down = qout[3, :, :]
    outlet_flow = qout[4, :, :]

    # ----------- Final depth at the cell ------------ #
    d_t_cell = d_t_cell - (I_tot_end_cell/cell_area)*1000  # final depth in mm

    return qout_left, qout_right, qout_up, qout_down, outlet_flow, d_t_cell, I_tot_end_cell, I

def CA_Routing_GPU(elevation_cell,  elevation_left_t,  elevation_right_t,
               elevation_up_t,  elevation_down_t,  d_t_cell,  d_left_cell,
               d_right_cell,  d_up_cell,  d_down_cell,  roughness_cell,
               cell_area,  time_step,  h_0_cell,  Resolution,  I_tot_end_cell,
               outlet_index,  outlet_type,  slope_outlet,  row_outlet,
               col_outlet,  idx_nan,  flag_critical):
    import cupy as cp
    # ----------- Finding Nans or Values Below a Threshold ------------ #
    d_t_min = 0.000001  # m
    d_t_min = cp.asarray(d_t_min)
    idx_2 = d_t_cell < d_t_min*1000
    idx_3 = idx_nan
    idx2 = idx_2 + idx_3
    idx2 = idx2 > 0
    # ----------- Adding minimum slope to do calculation ------------ #
    h_min = d_t_min  # m
    I_tot_end_previous = I_tot_end_cell
    # ----------- Notation ------------ #
    #  <-I-> (left right up down) = [left ; right ; up ; down]
    # ----------- Cell Depth and Water surface Elevation ------------ #
    depth_cell = cp.divide(d_t_cell, 1000)  # meters
    wse_cell = elevation_cell + depth_cell

    # ----------- Depth Differences for All Cells ------------ #
    delta_h = cp.zeros((5, cp.size(elevation_cell[:, 0]), cp.size(elevation_cell[0, :])))
    delta_h[0, :, :] = wse_cell - elevation_left_t - cp.divide(d_left_cell, 1000)  # left
    delta_h[1, :, :] = wse_cell - elevation_right_t - cp.divide(d_right_cell, 1000)  # right
    delta_h[2, :, :] = wse_cell - elevation_up_t - cp.divide(d_up_cell, 1000)  # up
    delta_h[3, :, :] = wse_cell - elevation_down_t - cp.divide(d_down_cell, 1000)  # down

    # %%% CHECKED 9/28/2021 %%%% #
    # ----------- Outlet Calculations ------------ #
    if outlet_type == 1:
        S_0 = slope_outlet  # Normal slope
        # V = h * area
        delta_h[4, :, :] = cp.multiply(S_0, cp.multiply(Resolution, outlet_index))  # normal slope depth difference
        # %%% CHECKED 9/28/2021 %%%% #
    else:
        S_0 = cp.zeros((cp.shape(elevation_cell)))
        for i in range(len(row_outlet)):
            # Checking left, right, up, down
            row = int(row_outlet[i])  # Row of outlet
            col = int(col_outlet[i])  # Col of outlet
            S_0[row, col] = cp.multiply(cp.power(depth_cell[row, col],-1/6), cp.multiply((cp.sqrt(9.81)),(roughness_cell[row, col])))  # Critical slope
        delta_h[4, :, :] = S_0*Resolution

    # ----------- Correcting delta_h ------------ #
    delta_h[cp.isnan(delta_h)] = 0  # TESTING
    delta_h[cp.isinf(delta_h)] = 0  # TESTING

    # ----------- Boundary Conditions ------------ #
    # These only applied when we are modeling a rectangular grid watershed
    delta_h[0, :, 0] = 0  # left
    delta_h[1, :, -1] = 0  # Right
    delta_h[2, 0, :] = 0  # Up
    delta_h[3, -1, :] = 0  # Down
    idx = delta_h < h_min
    delta_h[idx] = 0

    # ----------- Available volumes ------------ #
    Vol = cp.multiply(cell_area,delta_h)  # 3-D Array with pages following direction notation
    # Max h
    maxh = cp.max(delta_h, axis=0)
    # %%% CHECKED 9/28/2021 %%%% #
    # ----------- Minimum and Ma volumes ------------ #
    # Minimum
    c = 10000  # Large number
    c = cp.asarray(c)
    Vol_nan = Vol.copy()
    Vol_nan[idx] = c
    Vol_min_matrix = cp.min(Vol_nan, axis=0)
    # Total Volume
    Vol_tot = cp.sum(Vol, axis=0)

    # ----------- Weights ------------ #
    weight = cp.divide(Vol, (Vol_tot + Vol_min_matrix))
    weight_max = cp.max(weight, axis=0)

    # ----------- Velocity Calculation ------------ #
    # Velocity to the cell with the highest gradient
    if flag_critical == 1:
        v_m = cp.multiply(cp.minimum(cp.sqrt(9.81*depth_cell), cp.array(1/roughness_cell)), cp.multiply((cp.power((cp.maximum(depth_cell - h_0_cell/1000, 0)), 2/3)), cp.sqrt(maxh/Resolution)))
    else:
        v_m = cp.multiply((1/roughness_cell), cp.multiply(cp.power(cp.maximum(depth_cell - h_0_cell/1000, 0), 2/3), (cp.sqrt(maxh/Resolution))))

    # ----------- Inter-cell Volume Calculation ------------ #
    # Estimated Volume that goes to the cell with the highest gradient
    # I_m = np.array(v_m) * depth-cell * Resolution * time_step * 60  # Velocity * area * Time_step (m3)

    # ----------- Total outflow volume ------------ #
    # I_tot_end = np.minimum(depth_cell*cell_area, I_m/np.max(weight), Vol_min_matrix + I_tot_begin)
    I_tot_end_cell = cp.minimum(cp.multiply(depth_cell, cell_area), cp.multiply(v_m, cp.multiply(cp.where(weight_max == 0, 0, cp.divide(depth_cell,weight_max)), (Resolution*time_step*60))))
    I_tot_end_cell = cp.minimum(I_tot_end_cell, I_tot_end_previous + Vol_min_matrix)
    I_tot_end_cell[idx2] = 0  # if the depth is smaller than dmin ATTENTION HERE

    # ----------- Outflow Volume to each direction ------------ #
    I = cp.multiply(weight, I_tot_end_cell)  # m3
    I_tot_end_cell = cp.sum(I, axis=0)  # Volume that leaves the cell in m^3
    qout = (I)/(time_step*60)/(cell_area)*1000*3600  # mm/hr
    qout_left = qout[0, :, :]
    qout_right = qout[1, :, :]
    qout_up = qout[2, :, :]
    qout_down = qout[3, :, :]
    outlet_flow = qout[4, :, :]

    # ----------- Final depth at the cell ------------ #
    d_t_cell = d_t_cell - (I_tot_end_cell/cell_area)*1000  # final depth in mm

    return qout_left, qout_right, qout_up, qout_down, outlet_flow, d_t_cell, I_tot_end_cell, I