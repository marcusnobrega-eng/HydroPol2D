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
from labellines import labelLines
from matplotlib import rc
import matplotlib as mpl
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Garamond']}, size = '13')
rc('text', usetex=False)

# state some variables #
for k in np.arange(0, 22, 2):
    m_child = 22.4  # kg
    y_child = 1.21  # m
    D_child = 0.17  # m
    d_child = 0.085  # m

    m_adult = 71  # kg
    y_adult = 1.71  # m
    D_adult = 0.26  # m
    d_adult = 0.13  # m

    slope_ = k  # %
    slope = math.atan(slope_/100)
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

    child = np.empty((len(np.arange(0.01, y_child+0.3, 0.01)), 5))*np.nan
    adult = np.empty((len(np.arange(0.01, y_adult+0.3, 0.01)), 5))*np.nan
    row = 0
    slide = 1
    toppling = 2
    drowning = 3

    child[:, drowning] = (13/16)*y_child
    adult[:, drowning] = (13/16)*y_adult
    for i in np.arange(0.01, y_child+0.3, 0.01):
        child[row, 0] = i
        flag_s = 0
        flag_t = 0
        for j in np.arange(0.0, 5, 0.001):
            Buoy_child, area_child, xgs_child, ygs_child = submerged(i, y_child, d_child, D_child, slope)
            L = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**2)*math.cos((math.pi/2)-slope)*(j**2)*area_child
            D = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**3)*(j**2)*area_child
            Wp = m_child*gravity*math.sin(slope)
            Wn = m_child*gravity*math.cos(slope)
            Bn = density*gravity*Buoy_child*math.cos(slope)
            T = Fc*(Wn-Bn-L)
            if (D + Wp > T) and (flag_s==0):  # no longer stability by slipping
                child[row, slide] = j
                flag_s = 1
            if (D*(i/2)+(Wp*((7/12)*y_child*math.cos(slope)))+(Bn*((xgs_child/math.cos(slope))+(ygs_child*math.sin(slope))))+(L*((xgs_child/math.cos(slope))+((i/2)*(math.tan(slope)))))>Wn*(((5/6)*d_child)/math.cos(slope)+((7/12)*y_child)*math.sin(slope))) and (flag_t==0): # no longer stability by toppling
                child[row, toppling] = j
                flag_t = 1
            if flag_s==1 and flag_t==1:
                break
        row += 1
    row = 0
    slide = 1
    toppling = 2
    for i in np.arange(0.1, y_adult+0.3, 0.01):
        adult[row, 0] = i
        flag_s = 0
        flag_t = 0
        for j in np.arange(0.0, 5, 0.001):
            Buoy_adult, area_adult, xgs_adult, ygs_adult = submerged(i, y_adult, d_adult, D_adult, slope)
            L = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**2)*math.cos((math.pi/2)-slope)*(j**2)*area_adult
            D = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**3)*(j**2)*area_adult
            Wp = m_adult*gravity*math.sin(slope)
            Wn = m_adult*gravity*math.cos(slope)
            Bn = density*gravity*Buoy_adult*math.cos(slope)
            T = Fc*(Wn-Bn-L)
            if (D + Wp > T) and (flag_s==0):  # no longer stability by slipping
                adult[row, slide] = j
                flag_s = 1
            if (D*(i/2)+(Wp*((7/12)*y_adult*math.cos(slope)))+(Bn*((xgs_adult/math.cos(slope))+(ygs_adult*math.sin(slope))))+(L*((xgs_adult/math.cos(slope))+((i/2)*(math.tan(slope)))))>Wn*(((5/6)*d_adult)/math.cos(slope)+((7/12)*y_adult)*math.sin(slope))) and (flag_t==0): # no longer stability by toppling
                adult[row, toppling] = j
                flag_t = 1
            if flag_s==1 and flag_t==1:
                break
        row += 1
    child[:,4] = np.minimum(np.nan_to_num(child[:,slide], nan=np.inf),np.nan_to_num(child[:,toppling],nan=np.inf))
    adult[:,4] = np.minimum(np.nan_to_num(adult[:,slide], nan=np.inf),np.nan_to_num(adult[:,toppling],nan=np.inf))
    child= np.delete(child, np.argwhere(child[:,4]==0)[1:],axis=0)
    child= np.delete(child, np.argwhere(child[:,4]==np.inf),axis=0)
    adult= np.delete(adult, np.argwhere(adult[:,4]==0)[1:],axis=0)
    adult= np.delete(adult, np.argwhere(adult[:,4]==np.inf),axis=0)
    #
   # plt.plot(child[:,slide], child[:,0])
    #plt.plot(child[:,toppling], child[:,0])
    #plt.plot(child[:,0],child[:,drowning])
    plt.plot(child[:,4],np.minimum(child[:,0],child[:,drowning]),label=r'$\theta$= ' + str(k) + ' %', c='Black')
    # plt.plot(adult[:, slide], adult[:, 0])
    # plt.plot(adult[:, toppling], adult[:, 0])
    # plt.plot(adult[:, 0], adult[:, drowning])
    plt.plot(adult[:,4],np.minimum(adult[:,0],adult[:,drowning]),label=r'$\theta$= ' + str(k) + ' %', c='Black',linestyle='dashed')
    plt.xlim(0,4)
    plt.xlabel('Flow velocity (m/s)')
    plt.ylabel('Flow depth (m)')
lines = plt.gca().get_lines()
labelLines(lines, xvals=np.linspace(0.5,3.5,len(lines)), fontsize=7)
plt.legend(['Child','Adult'])










