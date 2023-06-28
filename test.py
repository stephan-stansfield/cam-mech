"""
Author: Stephan Stansfield, adapted from code originally by Yuxiang Ma
Date:   11/14/2022 - 6/13/2022
Python version 3.10
"""
#
from cam_class import CamGeneration
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import sys

#gear_ratios = np.array([[0.25, 0.15, 0.15, 0.25,    1,    1, 0.75, 0.3],
#                        [0,    0.05,  0.25, 0.3, 0.45, 0.55, 0.6, 0.7]])

# gear_ratios = np.array([[0.2, 0.18, 0.18, 0.32,    1.7,    2, 1.8, 0.8],
#                        [0,    0.08,  0.3, 0.35, 0.45, 0.55, 0.6, 0.65]])

# gear_ratios = np.array([[0.3, 0.22, 0.3, 0.8, 1.7 , 2 , 1.8 , 0.8,  0.38, ],
#                         [   0, 0.1, 0.25, 0.38, 0.45, 0.5, 0.58, 0.66, 0.75]])
# ninterp = 360
# angles = np.arange(0, ninterp+1, 1)/ninterp*360
# grs = interpolate.interp1d(gear_ratios[1]*360, gear_ratios[0], kind='cubic')
# plt.figure()
# plt.scatter(angles, grs(angles), label='gear ratio')
# # plt.legend(loc="upper left")
# plt.xlabel('angle (deg)')
# plt.ylabel('gear ratio')
# # plt.xlim([0, 300])
# # plt.ylim([0, 2])

# column 0 of gear_ratios holds desired gear ratios
# column 1 of gear_ratios holds percentage of sitting motion associated with the gear ratios in col 0
gear_ratios = np.array([[0.18, 0.18, 0.3, 0.8, 1.7 , 2 , 1.8 , 0.8 ],
                        [0  , 0.08, 0.28, 0.4, 0.45, 0.5, 0.58, 0.66]])

# gear_ratios = np.array([[]
#                         []])

gear_ratios[1, :] = 2*np.pi*gear_ratios[1, :] # convert column 1 of gear_ratios to radians

initial_gr = 0.5

# generate cams using selected gear ratios, initial gear ratio, and sit angle.
# Sit angle is input in radians.
Cam = CamGeneration(gear_ratios.T, initial_gr, sit_angle=np.pi*220/180)

Cam.calculate_cam_radii()
# Cam.plot_cams(stroke=10)
# Cam.plot_forces(torque=False, stroke=0.1) #set to False if you don't want to see the torque plot
# Cam.normalize_angles(knee_csv)

# calculate cable force using knee angle data
filename = 'data/Knee-angle_Chugo_2006.csv'
Cam.plot_forces_percentage(filename, torque=False, stroke=0.1)


# Cam.alter_cam_profiles(routing=np.array([[np.pi/2, 0.2, ]]))
inner_pts, outer_pts = Cam.cam_pts(stroke=0.1) # output arrays of point coordinates (X, Y, Z) on perimeter of cams

np.savetxt('inner.txt', np.append(inner_pts, np.zeros((inner_pts.shape[0],1))*1000, axis=1), delimiter=',')  # mm
np.savetxt('outer.txt', np.append(outer_pts, np.zeros((outer_pts.shape[0],1))*1000, axis=1), delimiter=',')  # mm