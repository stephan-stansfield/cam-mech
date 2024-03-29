"""Using python 3.10"""
import math

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import make_interp_spline, BSpline
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from collections import defaultdict

font = {'size': 14}
matplotlib.rc('font', **font)


class CamGeneration:
    """Class that generates points on the cams for
    every degree given torque/gear ratios. 
    Can give energy output given stiffness as well.
    """
    
    def __init__(self, gear_ratios, input_angles, scaling, sit_angle=np.pi):
        """points: number of points to define on cams
        gear_ratios: (points*2)-dim ndarray [gear_ratios; angles]
        scaling: scalar multiple, e.g. 0.5, to determine overall size of cams
        """
        # Sort gear ratios and input angles in ascending order of input angles.
        sorted_inds = np.argsort(input_angles)
        self.input_angles = input_angles[sorted_inds]
        self.input_angles = input_angles * sit_angle / np.pi # SS: don't understand this scaling still...the implicit assumption is that the input angle is pi?
        self.gear_ratios = gear_ratios[sorted_inds]

        # Initialize class variables.
        self.scaling = scaling
        self.cam_radii = np.ones((self.gear_ratios.shape[0], 2))
        self.n_interp = 360 
        self.sit_ind = round(sit_angle * self.n_interp / (2*np.pi))
        self.angles = np.arange(0, self.n_interp+1, 1) * 2*np.pi / self.n_interp 
        self.pts_inner = None
        self.pts_outer = None
        # self.cam_radii = np.ones((self.points, 2)) @ np.diag([1, 1/scaling])  # initial guess for inner and outer cam radii
        # self.input_angles = self.gear_ratios[:, 1] # take the second row of gear ratios to be the input angles

    def calculate_cam_radii(self, stroke=np.pi, plot=False, index=0):
        """Calculates the cam radii for each gear ratio and input angle.
        Determines the convex hull of the given gear ratios and angles, then
        interpolates between these points to generate the cam radii.
        
        Parameters:
        stroke: scalar, full stroke of the pulling cable

        Returns:
        pts_inner: ndarray, inner cam points in Cartesian space
        pts_outer: ndarray, outer cam points in Cartesian space
        radius_max: scalar maximum radius of the cam envelopes
        """
        
        # Calculate initial cam sizes using gear ratios and scaling factor.
        r_min = 0.5
        for ind, ratio in enumerate(self.gear_ratios):
            # r = self.scaling * np.sqrt(ratio) * self.cam_radii[ind, 0]
            # R = self.scaling * self.cam_radii[ind, 1] / np.sqrt(ratio)
            [r, R] = self.scaling * np.array([np.sqrt(ratio), 1/np.sqrt(ratio)])
            if r < r_min:
                r = r_min
                R = r_min / ratio
            self.cam_radii[ind, :] = np.array([r, R])
            # self.cam_radii[ind, 1] = R

        # Pad arrays to ensure continuity and then interpolate between points.
        self.cam_radii = np.append(self.cam_radii, self.cam_radii[: 2, :],
                                   axis=0)
        self.input_angles = np.append(self.input_angles,
                                      self.input_angles[0 : 2] + 2*np.pi)
        self.cam_radii = np.insert(self.cam_radii, 0, self.cam_radii[-4 : -2, :],
                                   axis = 0)
        self.input_angles = np.insert(self.input_angles, 0,
                                       self.input_angles[-4 : -2] - 2*np.pi)
        cam_fcn_inner = interpolate.interp1d(self.input_angles,
                                             self.cam_radii[:, 0],
                                             kind='cubic',
                                             fill_value='extrapolate')
        cam_fcn_outer = interpolate.interp1d(self.input_angles,
                                             self.cam_radii[:, 1],
                                             kind='cubic',
                                             fill_value='extrapolate')
        self.cam_radii = np.vstack((cam_fcn_inner(self.angles),
                                    cam_fcn_outer(self.angles))).T
        
        if np.any(np.isnan(self.cam_radii)):
            print("Warning: NaN in cam radii")

        # Define cam points in Cartesian space.
        # SS: do we need the phase change here?
        # self.pts_inner = np.array([self.cam_radii[:, 0] * np.cos(self.angles),
        #                            self.cam_radii[:, 0] * np.sin(self.angles)]).T
        # self.pts_outer = np.array([R * np.cos(self.angles - np.pi / 2), R * np.sin(self.angles - np.pi / 2)]).T
        self.pts_inner = (self.cam_radii[:, 0] * [np.cos(self.angles),
                                                  np.sin(self.angles)
                                                  ]).T
        self.pts_outer = (self.cam_radii[:, 1] * [np.cos(self.angles - np.pi/2),
                                                  np.sin(self.angles - np.pi/2)
                                                  ]).T
        if plot:
            plt.scatter(self.pts_inner[:,0], self.pts_inner[:,1], lw = 2)
            plt.scatter(self.pts_outer[:,0], self.pts_outer[:,1], lw = 2)
            plt.legend(['inner cam','outer cam'])
            plt.axis('equal')
            plt.title('before convex function plot')
            plt.show()
        
        # Calculate points of convex hull for each cam.
        self.cam_radii[:, 0] = self.convex_cam_pts(self.pts_inner, plot)
        self.cam_radii[:, 1] = self.convex_cam_pts(self.pts_outer, plot)

        if np.any(np.isnan(self.cam_radii)):
            print("Warning: NaN in cam radii")

        if plot:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.plot(self.angles - np.pi, self.cam_radii[:, 0])
            ax.plot(self.angles - (3*np.pi/2), self.cam_radii[:, 1] )  # apply phase change to the outer cam
            ax.set_rticks([2, 4, 6, 8, 10])
            ax.set_rlabel_position(-22.5)  # move radial labels away from plotted line
            ax.grid(True)
            ax.set_title('interpolated cam')

        # Define cam points in Cartesian space after interpolating in polar
        # coordinates. Apply phase changes to counteract that made above
        # SS: do we need both...?
        self.pts_inner = (self.cam_radii[:, 0] * [np.cos(self.angles - np.pi),
                                                 np.sin(self.angles - np.pi)
                                                 ]).T
        self.pts_outer = (self.cam_radii[:, 1] * [np.cos(self.angles - 3*np.pi/2),
                                                  np.sin(self.angles - 3*np.pi/2)
                                                  ]).T

        # Check that Cartesian points match polar points.
        if plot:
            plt.figure()
            plt.plot(self.pts_inner[:,0], self.pts_inner[:,1], lw = 2)
            plt.plot(self.pts_outer[:,0], self.pts_outer[:,1], lw = 2)
            plt.legend(['inner cam','outer cam'])
            plt.axis('equal')
            plt.show()

        # Scale cam size to desired stroke length
        self.x_cable = np.cumsum(self.cam_radii[:, 0] * 2*np.pi / self.n_interp)
        ratio = stroke / self.x_cable[self.sit_ind]
        self.cam_radii *= ratio
        radius_max = np.max(self.cam_radii)
        if plot:
            self.plot_cams(self.cam_radii, stroke, index)
        return self.pts_inner, self.pts_outer, radius_max

    def convex_cam_pts(self, points, plot=False):
        hull = ConvexHull(points, incremental=True)
        points = np.vstack((points[hull.vertices, 0], points[hull.vertices, 1])).T

        if plot:
            plt.figure()
            plt.scatter(hull.points[hull.vertices][:,0],hull.points[hull.vertices][:,1])
            plt.title('hull points (in convex_cam_pts function)')
            plt.show()
        """ plt.plot(self.pts_inner[:,0], self.pts_inner[:,1], lw = 2)
        plt.plot(self.pts_outer[:,0], self.pts_outer[:,1], lw = 2)
        plt.legend(['inner cam','outer cam'])
        plt.axis('equal')
        plt.title('Inner & outer points (in convex_cam_pts function)')
        plt.show() """

        # Split cam points into 2 or 3 sub-arrays that each monotonically
        # increase or decrease in x. Flip any segments that are monotonically
        # decreasing, which is required by numpy.interp(), and interpolate.
        break_ind = []
        for ind in range(len(points) - 2):
            if (math.copysign(1, np.diff(points[ind : ind+2, 0]))
                != math.copysign(1, np.diff(points[ind+1 : ind+3, 0]))):
                    break_ind.append(ind+2)
        points_array = []
        points_array.append(points[: break_ind[0], :])
        if len(break_ind) == 2:
            points_array.append(points[break_ind[0] : break_ind[1], :])
            points_array.append(points[break_ind[1] :, :])
        elif len(break_ind) == 1:
            points_array.append(points[break_ind[0] :, :])
        for ind, sub_array in enumerate(points_array):
            if not np.all(np.diff(sub_array[:, 0]) > 0):
                points_array[ind] = np.flip(points_array[ind], 0)
        for ind, sub_array in enumerate(points_array):
            # x = np.linspace(sub_array[0, 0], sub_array[-1, 0], len(sub_array))
            x = np.linspace(sub_array[0, 0], sub_array[-1, 0], 100)
            y = np.interp(x, sub_array[:, 0], sub_array[:, 1])
            if ind == 0:
                interp_array = np.array([x, y]).T
            else:
                interp_array = np.concatenate((interp_array, np.array([x, y]).T))

        if plot:
            plt.scatter(interp_array[:,0], interp_array[:,1])
            plt.title('Interpolated convex cam points')
            plt.axis('equal')
            plt.show()

        # Convert interpolated convex cam points into polar coordinates.
        # Linearly interpolate between points in polar coordinates and return.
        polar_radii = np.empty(interp_array.shape[0])
        for ind, point in enumerate(interp_array):
            dist = np.linalg.norm(point)
            polar_radii[ind] = dist
        polar_angles = np.arctan2(interp_array[:, 1], interp_array[:, 0]) + np.pi
        polar_angles = np.array([polar_angles]).flatten() ## SS: check if needed

        if plot:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.scatter(polar_angles, polar_radii)
            ax.set_rticks([2, 4, 6, 8, 10])
            ax.set_rlabel_position(-22.5)  # move radial labels away from plotted line
            ax.grid(True)
            ax.set_title('radii and angles before interpolation')

        polar_cam_fcn = interpolate.interp1d(polar_angles,
                                             polar_radii,
                                             kind='linear',
                                             fill_value='extrapolate')
        
        if np.any(np.isnan(polar_cam_fcn(self.angles))):
            print("Warning: NaN in cam radii")

        return polar_cam_fcn(self.angles)

    def calc_forces_percentages(self, angle_data, torque=False, stroke=2*np.pi,
                               plot=False, index=0):
        """Calculate the force profiles vs stance percentage given the solved
        cam radii.
        stroke: scalar, full stroke of the pulling cable
        """
        # Calculate elastic band tension from displacement and experimental 
        # characterization assuming linear stiffness. Calculate cable tension
        # from elastic band tension and gear ratios.
        x_elastic = 2 * np.cumsum(self.cam_radii[:, 1] * 2*np.pi / self.n_interp)
        E_elastic = 14
        k_elastic = 2 * E_elastic / x_elastic[self.sit_ind]**2
        f_elastic = k_elastic * x_elastic
        f_cable =  f_elastic * self.cam_radii[:, 1] / self.cam_radii[:, 0]
        # print(f"Elastic band stiffness {k} N/m; Pulling distance {x_e[self.sit_ind]} m")

        # Account for nonlinear change in knee angle during standing.
        # Percentages go from ~0% to ~100% going from sitting to standing.
        knee_array = np.loadtxt(angle_data, delimiter = ',', ndmin = 2)
        percentages = knee_array[:, 0]
        knee_angles = knee_array[:, 1]
        knee_fcn = interpolate.interp1d(percentages, knee_angles,
                                        kind='cubic', fill_value='extrapolate')
        percentages = np.linspace(0, 100, self.sit_ind + 1)
        knee_angles = knee_fcn(percentages)

        # Find linear relationship between stand percentage and cable
        # displacement. Flip arrays to go from standing to sitting. Scale cable
        # displacement to knee angle percentage.
        knee_angles = np.flip(knee_angles)
        percentages = np.flip(percentages)
        x_cable_knee = (self.x_cable[self.sit_ind]
                        * (1 - (knee_angles-np.min(knee_angles))
                           / np.ptp(knee_angles))) 

        if plot:
            plt.figure()
            plt.plot(knee_angles, x_cable_knee)
            plt.xlabel('knee angle (degree)')
            plt.ylabel('cable displacement, xc (m)')
            plt.title('knee angle and cable displacement linear relationship')

            # plot stance percentage vs. cable displacement
            plt.figure()
            plt.plot(percentages, x_cable_knee)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cable Displacement (m)')

        # Add "negative" stance phase to account for forces that occur when cam 
        # is rotated past the sitting position. For this phase, assume a linear
        # relationship between cable displacement and stance percentage.
        perc_remaining = (self.sit_ind - self.n_interp) / self.n_interp * 100
        ind_perc_rem = (np.where(np.round(percentages)
                                 == np.round(-perc_remaining))[0][0])
        perc_neg = np.linspace(-1/self.sit_ind * 100,
                               perc_remaining,
                               self.n_interp - self.sit_ind)
        x_c_knee_neg = np.linspace(x_cable_knee[-1],
                                 2 * x_cable_knee[-1] - x_cable_knee[ind_perc_rem],
                                 self.n_interp - self.sit_ind)
        percentages = np.concatenate((percentages, perc_neg))
        x_cable_knee = np.concatenate((x_cable_knee, x_c_knee_neg))

        if np.any(np.isnan(self.x_cable)):
            print("Warning: NaN in cam radii")
        
        # Remove duplicate points to avoid errors in interpolation.
        repeat_dict = defaultdict(list)
        for ind, point in enumerate(self.x_cable):
            repeat_dict[point].append(ind)
        repeat_dict = {k:v for k,v in repeat_dict.items() if len(v)>1}
        repeat_ind = np.array(list(repeat_dict.values()))
        if len(repeat_dict) > 0:
            print("Repeated values and indices: ", repeat_dict)
            self.x_cable = np.delete(self.x_cable, repeat_ind)
            self.angles = np.delete(self.angles, repeat_ind)
            x_elastic = np.delete(x_elastic, repeat_ind)

        # Relate elastic force to stance percentage into account.
        cable_fcn = interpolate.interp1d(self.x_cable, self.angles,
                                         kind='cubic', fill_value='extrapolate')
        elastic_fcn = interpolate.interp1d(self.angles, x_elastic,
                                           kind='cubic', fill_value='extrapolate')
        angle_scaled = cable_fcn(x_cable_knee)
        x_elastic_scaled = elastic_fcn(angle_scaled)
        f_cable_scaled = (self.cam_radii[:, 1] / self.cam_radii[:, 0]
                          * k_elastic * x_elastic_scaled)
        
        # Plot and save values.
        if plot:
            plt.figure()
            plt.plot(percentages, angle_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cam Angle (rad)')

            # plot elastic band displacement vs. stance percentage
            plt.figure()
            plt.plot(percentages, x_elastic_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Elastic displacement (m)')
            
            # plot cable force vs. stance percentage
            plt.figure()
            plt.plot(percentages, f_cable_scaled,
                     label='Cable Force', linewidth=3)
            plt.legend(loc='upper right')
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Force (N)')
            plt.xlim([-60, 100])
            plt.ylim([0, 250])
            filename = 'force_plot_' + str(index) + '.png'
            plt.savefig(filename, dpi=300)
            
            np.savetxt('force_output.csv',
                       np.stack((self.angles, self.x_cable, f_elastic, f_cable),
                                axis=1),
                                header='angles(rad), pulling distance(m), \
                                elastic band force(n), cable force(n)')
            if torque:
                plt.figure()
                plt.scatter(self.angles * 360 / (2*np.pi),
                            f_elastic * self.cam_radii[:, 0])
                plt.xlabel('Angle (deg)')
                plt.ylabel('Torque (N-m)')
            plt.show()

        return f_cable_scaled[:-10], percentages[:-10] ## SS: why [:-10]?

    def plot_cams(self, cam_radii=0, stroke=np.pi, index=0):
        """
        Plots the cam points as a sanity check feature.
        stroke: scalar, full stroke of the pulling cable
        """
        # print("cam_radii: ", cam_radii)
        # if cam_radii.all() == 0:
        #     cam_radii = self.cam_radii
        r = self.cam_radii[:, 0]
        R = self.cam_radii[:, 1]
        # x_c = np.cumsum(r * 2 * np.pi / self.n_interp) 
        # ratio = stroke / x_c[self.sit_ind] # full stroke divided by distance along cam perimeter of sitting angle
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        # ax.plot(self.angles - np.pi / 2, R * ratio)  # apply phase change to the outer cam
        # ax.plot(self.angles, r*ratio)
        ax.plot(self.angles, r)
        ax.plot(self.angles, R)
        ax.set_rticks([2, 4, 6, 8, 10])
        ax.set_rlabel_position(-22.5)  # move radial labels away from plotted line
        ax.grid(True)
        ax.set_title(f"Cam Demonstration, min radius={np.min(self.cam_radii)}, max radius={np.max(self.cam_radii)}", va='bottom')
        filename = 'cam_plot_' + str(index) + '.png'
        plt.savefig(filename, dpi=300)