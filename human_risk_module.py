###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                        March 2023                               #
#                                                                 #
#                 Last update : 1 March, 2023                     #
###################################################################

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# -------- DESCRIPTION ----------------
# This function calculates human risk of being dragged by water forces
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import cupy as cp
import math

def human_risk(vel_t_zzz,h, slope, m, y, D, d, density, gravity, Cc, Fc):
    m_c, m_t, m_a, m_o = 22.4, 45.5, 60.0, 71.0
    y_c, y_t, y_a, y_o = 1.21, 1.50, 1.70, 1.65
    D_c, D_t, D_a, D_o = 0.17, 0.20, 0.26, 0.24
    d_c, d_t, d_a, d_o = 0.085, 0.1, 0.13, 0.13

    density, gravity = 1000, 9.81
    Cc, Fc= 1, 0.46

    # Person's Geometry for h<y
    sub_temp_1 = cp.logical_and(cp.divide(h,cp.cos(slope)) < y/2)
    V_s1 = cp.multiply(sub_temp_1,2*cp.divide(h,cp.cos(slope))*math.pi*((d**2)/4))
    A_s1 = cp.multiply(sub_temp_1,2*cp.divide(h,cp.cos(slope))*d)
    X_gs1 = cp.multiply(sub_temp_1,d/2)
    Y_gs1 = cp.multiply(sub_temp_1,cp.divide(h,(2*cp.cos(slope))))

    # Person's Geometry for h>y
    sub_temp_2 = cp.logical_and(cp.divide(h, cp.cos(slope)) > y_c/2)
    V_s2 = cp.multiply(sub_temp_2, 2*(y/2)*math.pi*((d**2)/4) + (cp.divide(h, cp.cos(slope))-(y/2))*math.pi*((D**2)/4))
    A_s2 = cp.multiply(sub_temp_2, 2*(y/2)*d+((cp.divide(h, cp.cos(slope)))-(y/2))*D)
    X_gs2 = cp.multiply(sub_temp_2, cp.divide(2*((math.pi*d**2)/4)*((y/2)*(d/2)) + ((math.pi*D**2)/4)*(cp.divide(h, cp.cos(slope))-(y/2))*(D/2), 2*((math.pi*d**2)/4)*(y/2) + ((math.pi*D**2)/4)*(cp.divide(h, cp.cos(slope))-(y/2))))
    Y_gs2 = cp.multiply(sub_temp_2, cp.divide(2*((math.pi*d**2)/4)*((y**2)/8)+((math.pi*D**2)/4)*cp.multiply((cp.divide(h, cp.cos(slope))-(y/2)), ((0.5*(cp.divide(h, cp.cos(slope))-(y/2)))+(y/2))), 2*((math.pi*d**2)/4)*(y/2) + ((math.pi*D**2)/4)*(cp.divide(h, cp.cos(slope))-(y/2))))

    # Drag forces by the flow
    L1 = 0.5*density*Cc*cp.multiply(cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 2), cp.cos((math.pi/2)-slope)), cp.multiply(cp.power(vel_t_zzz, 2), A_s1))
    L2 = 0.5*density*Cc*cp.multiply(cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 2), cp.cos((math.pi/2)-slope)), cp.multiply(cp.power(vel_t_zzz, 2), A_s2))
    L = L1+L2
    D1 = 0.5*density*Cc*cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 3), cp.multiply(cp.power(vel_t_zzz, 2), A_s1))
    D2 = 0.5*density*Cc*cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 3), cp.multiply(cp.power(vel_t_zzz, 2), A_s2))
    Drag = D1+D2
    # Person weight forces acting
    Wp = m*gravity*cp.sin(slope)
    Wn = m*gravity*cp.cos(slope)
    # Buoyancy effect
    Bn1 = density*gravity*cp.multiply(V_s1, cp.cos(slope))
    Bn2 = density*gravity*cp.multiply(V_s2, cp.cos(slope))
    Bn = Bn1+Bn2
    T = Fc*(Wn-Bn-L) #  Friccion available
    M = cp.multiply(Drag, cp.divide(h,2)) + cp.multiply(Wp, (7/12)*y*cp.cos(slope)) +  # Momemtum acting
    # Slipping stability
    slide = cp.multiply(cp.logical(Drag + Wp > T), Drag + Wp)
    # Toppling stability
    Topple =
    # Drawing risk