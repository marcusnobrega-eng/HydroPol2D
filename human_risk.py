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
import math
import os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc
import matplotlib as mpl

# state some variables #
m_child = 22.4  # kg
y_child = 1.21  # m
D_child = 0.17  # m
d_child = 0.085  # m

m_adult = 71  # kg
y_adult = 1.71  # m
D_adult = 0.26  # m
d_adult = 0.13  # m

slope = 0  # degrees
density = 1000  # kg/m3
gravity = 9.81  # m/s2
Cc = 1  # Drag coefficient
Fc = 0.46  # Friction coefficient

# define some funcitons #
def submerged (h, Y, d, D, slope):
    if (h/math.cos(slope) < Y/2):
        V_s = 2 * (h / math.cos(slope)) * math.pi * ((d**2) / 4)
        A_s = 2*(h/math.cos(slope))*d
        X_gs = d/2
        Y_gs = h/(2*math.cos(slope))
        return V_s, A_s, X_gs, Y_gs
    elif (h/math.cos(slope) > Y/2):
        V_s = 2*(Y/2)*math.pi*((d**2)/4)+((h/math.cos(slope))-(Y/2))*math.pi*((D**2)/4)
        A_s = 2*(Y/2)*d+((h/math.cos(slope))-(Y/2))*D
        X_gs = (2*((math.pi*d**2)/4)*((Y/2)*(d/2))+((math.pi*D**2)/4)*((h/math.cos(slope))-(Y/2))*(D/2)) / (2*((math.pi*d**2)/4)*(Y/2) + ((math.pi*D**2)/4)*((h/math.cos(slope))-(Y/2)))
        Y_gs = (2*((math.pi*d**2)/4)*((Y**2)/8)+((math.pi*D**2)/4)*((h/math.cos(slope))-(Y/2))*((0.5*((h/math.cos(slope))-(Y/2)))+(Y/2))) / (2*((math.pi*d**2)/4)*(Y/2) + ((math.pi*D**2)/4)*((h/math.cos(slope))-(Y/2)))
        return V_s,A_s,X_gs,Y_gs

main = np.zeros((len(np.arange(0.1, y_child, 0.01)), 1+ len(np.arange(0.1, 5, 0.1))))
row = 0
for i in np.arange(0.1, y_child, 0.01):
    main[row, 0] = i
    col = 1
    for j in np.arange(0.1, 5, 0.1):
        Buoy_child, area_child, xc_child, yc_child = submerged(i, y_child, d_child, D_child, slope)
        L = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**2)*math.cos((math.pi/2)-slope)*(j**2)*area_child
        D = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**3)*(j**2)*area_child
        Wp = m_child*gravity*math.sin(slope)
        Wn = m_child*gravity*math.cos(slope)
        Bn = density*gravity*Buoy_child*math.cos(slope)
        T = Fc*(Wn-Bn-L)
        final = T - D - Wp
        main[row, col] = final
        col += 1
    row += 1














