"""
Author: Stephan Stansfield, adapted from code originally by Yuxiang Ma
Date:   11/14/2022 - 7/9/2023
Python version 3.10
"""
#
from cam_class_optimization import CamGeneration
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import sys

# row 0 of gear_ratios holds desired gear ratios
# row 1 of gear_ratios holds percentage of sitting motion associated with the gear ratios in row 0
# profiles as of committee meeting 1 (peak at 32.5%; peak force 175 N, sitting force 27.5 N):
gear_ratios = np.array([[0.18, 0.2, .23, 0.5, 1.4,  1.7,  1.5, 0.9], # is there a min and max value for these? - ratios are sensitive, convex check, 
                        [0  , 0.08, 0.15, 0.4, 0.45, 0.5, 0.60, 0.66]])

gear_ratios[1, :] = 2 * np.pi * gear_ratios[1, :] # convert row 1 of gear_ratios to radians
print(gear_ratios.shape)
initial_gr = 0.5 # initial overall gear ratio, inner cam/outer cam

# generate cams using selected gear ratios, initial gear ratio, and sit angle.
# sit angle is input in radians.
Cam = CamGeneration(gear_ratios.T, initial_gr, sit_angle = np.pi * 220 / 180)

# calculate the convex hull of the cams and get its max radius
max_radius = Cam.calculate_cam_radii() 

# calculate cable force using knee angle data
filename = 'data/Knee-angle_Chugo_2006.csv'
forces,percentages = Cam.plot_forces_percentage(filename, torque = False, stroke = 0.1,plot=True)

#print(forces,np.round(percentages,1))
ideal_perc_ind = np.where(np.round(percentages,1) == 32.3)[0][0]
ideal_perc_force = forces[ideal_perc_ind]
max_force = max(forces)
max_force_perc = percentages[np.where(forces == max_force)[0][0]]
sitting_perc_ind = np.where(np.round(percentages,1) == 0)[0][0]
sitting_force = forces[sitting_perc_ind]

print("All index value of 32.3 is: ", ideal_perc_force)
print('force at 32.3% = ',ideal_perc_force, 'max force of', max_force,' at ',max_force_perc, ' %')
print('envelope = ', 2 * max_radius)
print(f'sitting force: {sitting_force}')

# inner_pts, outer_pts = Cam.cam_pts(stroke=0.1) # output arrays of point coordinates (X, Y, Z) on perimeter of cams
# np.savetxt('inner.txt', np.append(inner_pts, np.zeros((inner_pts.shape[0], 1)) * 1000, axis = 1), delimiter = ',')  # mm
# np.savetxt('outer.txt', np.append(outer_pts, np.zeros((outer_pts.shape[0], 1)) * 1000, axis = 1), delimiter = ',')  # mm