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

def human_risk(v, h, slope, m, y, D, d, density, gravity, Cc, Fc):

    # Person's Geometry for h<y
    sub_temp_1 = cp.less(cp.divide(h,cp.cos(slope)), y/2).astype(int)*1
    V_s1 = cp.multiply(sub_temp_1,2*cp.divide(h,cp.cos(slope))*math.pi*((d**2)/4))
    A_s1 = cp.multiply(sub_temp_1,2*cp.divide(h,cp.cos(slope))*d)
    X_gs1 = cp.multiply(sub_temp_1,d/2)
    Y_gs1 = cp.multiply(sub_temp_1,cp.divide(h,(2*cp.cos(slope))))

    # Person's Geometry for h>y
    sub_temp_2 = cp.greater(cp.divide(h, cp.cos(slope)), y/2).astype(int)*1
    V_s2 = cp.multiply(sub_temp_2, 2*(y/2)*math.pi*((d**2)/4) + (cp.divide(h, cp.cos(slope))-(y/2))*math.pi*((D**2)/4))
    A_s2 = cp.multiply(sub_temp_2, 2*(y/2)*d+((cp.divide(h, cp.cos(slope)))-(y/2))*D)
    X_gs2 = cp.multiply(sub_temp_2, cp.divide(2*((math.pi*d**2)/4)*((y/2)*(d/2)) + ((math.pi*D**2)/4)*(cp.divide(h, cp.cos(slope))-(y/2))*(D/2), 2*((math.pi*d**2)/4)*(y/2) + ((math.pi*D**2)/4)*(cp.divide(h, cp.cos(slope))-(y/2))))
    Y_gs2 = cp.multiply(sub_temp_2, cp.divide(2*((math.pi*d**2)/4)*((y**2)/8)+((math.pi*D**2)/4)*cp.multiply((cp.divide(h, cp.cos(slope))-(y/2)), ((0.5*(cp.divide(h, cp.cos(slope))-(y/2)))+(y/2))), 2*((math.pi*d**2)/4)*(y/2) + ((math.pi*D**2)/4)*(cp.divide(h, cp.cos(slope))-(y/2))))

    # Drag forces by the flow
    L1 = 0.5*density*Cc*cp.multiply(cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 2), cp.cos((math.pi/2)-slope)), cp.multiply(cp.power(v, 2), A_s1))
    L2 = 0.5*density*Cc*cp.multiply(cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 2), cp.cos((math.pi/2)-slope)), cp.multiply(cp.power(v, 2), A_s2))
    L = L1+L2
    D1 = 0.5*density*Cc*cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 3), cp.multiply(cp.power(v, 2), A_s1))
    D2 = 0.5*density*Cc*cp.multiply(cp.power(cp.sin((math.pi/2)-slope), 3), cp.multiply(cp.power(v, 2), A_s2))
    Drag = D1+D2
    # Person weight forces acting
    Wp = m*gravity*cp.sin(slope)
    Wn = m*gravity*cp.cos(slope)
    # Buoyancy effect
    Bn1 = density*gravity*cp.multiply(V_s1, cp.cos(slope))
    Bn2 = density*gravity*cp.multiply(V_s2, cp.cos(slope))
    Bn = Bn1+Bn2
    X_gs = X_gs1+X_gs2
    Y_gs = Y_gs1+Y_gs2
    #  Friccion available
    T = Fc*(Wn-Bn-L)
    # Momemtum acting
    M = cp.multiply(Drag, cp.divide(h,2)) + cp.multiply(Wp, (7/12)*y*cp.cos(slope)) + (cp.multiply(Bn,cp.divide(X_gs,cp.cos(slope)))+cp.multiply(Y_gs,cp.sin(slope))) \
        + cp.multiply(L, cp.divide(X_gs, cp.cos(slope)) + cp.multiply(cp.divide(h, 2), cp.tan(slope)))
    # Slipping stability
    slide = cp.greater(Drag + Wp, T).astype(int)*1
    # Toppling stability
    topple = cp.greater(M, cp.multiply(Wn, cp.divide((5/6)*d, cp.cos(slope)) + cp.divide((7/12)*y, cp.sin(slope)))).astype(int)*2
    # Drawing risk
    drawing = cp.greater(h, (13/16)*y).astype(int)*3
    return cp.maximum(slide, topple, drawing)