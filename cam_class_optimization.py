"""Using python 3.10"""
import math
from datetime import date
from os import path, makedirs

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
    
    def __init__(self, gear_ratios, input_angles, scaling, sit_angle, offset_angle):
        """
        Parameters:
        gear_ratios: desired gear ratio at each angle along the compound cam
        input_angles: angles at which to define the desired gear ratios
        scaling: scalar multiple, e.g. 0.5, to determine overall size of cams
        sit_angle: angle in radians along the compound cam at which force cable
        is acting while the user is seated
        offset_angle: angle in radians along the compound cam at which the 
        energy storage band acts relative to the force cable
        """
        # Sort gear ratios and input angles in ascending order of input angles.
        self.input_angles = input_angles
        self.gear_ratios = gear_ratios
        sorted_inds = np.argsort(self.input_angles)
        self.input_angles = self.input_angles[sorted_inds]
        self.gear_ratios = gear_ratios[sorted_inds]

        # Initialize class variables.
        self.scaling = scaling
        self.cam_radii = np.ones((self.gear_ratios.shape[0], 2))
        self.n_interp = 360 
        self.sit_ind = round(sit_angle * self.n_interp / (2*np.pi))
        self.offset_ind = round(offset_angle * self.n_interp / (2*np.pi))
        self.angles = np.arange(0, self.n_interp+1, 1) * 2*np.pi / self.n_interp 
        self.pts_inner = None
        self.pts_outer = None
        self.dateStr = date.today().strftime("%Y-%m-%d")

    def calculate_cam_radii(self, user_height=1.67, plot=False, index=0):
        """Calculates the cam radii for each gear ratio and input angle.
        Determines the convex hull of the given gear ratios and angles, then
        interpolates between these points to generate the cam radii.
        
        Parameters:
        user_height: height, in meters, of the person the device is designed for
        user_height: height, in meters, of the person the device is designed for

        Returns:
        pts_inner: ndarray, inner cam points in Cartesian space
        pts_outer: ndarray, outer cam points in Cartesian space
        radius_max: scalar maximum radius of the cam envelopes
        """
        
        # Calculate initial cam keypoints using gear ratios and scaling factor.
        # Check if any of the radii are smaller than the minimum allowed radius
        # and adjust as needed.
        radius_min = 0.5
        for ind, ratio in enumerate(self.gear_ratios):
            [r, R] = self.scaling * np.array([np.sqrt(ratio), 1/np.sqrt(ratio)])
            if ratio < 1:
                if r < radius_min:
                    r = radius_min
                    R = r / ratio
            if ratio > 1:
                if R < radius_min:
                    R = radius_min
                    r = R * ratio
            self.cam_radii[ind, :] = np.array([r, R])

        if plot:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.scatter(self.input_angles, self.cam_radii[:, 0])
            ax.scatter(self.input_angles, self.cam_radii[:, 1] )
            ax.grid(True)
            ax.set_title('original, scaled cam points')
            plt.show()
            

        # Convert cam points to Cartesian space.
        # SS: do we need the phase change here?
        # self.pts_inner = np.array([self.cam_radii[:, 0] * np.cos(self.angles),
        #                            self.cam_radii[:, 0] * np.sin(self.angles)]).T
        # self.pts_outer = np.array([R * np.cos(self.angles - np.pi / 2), R * np.sin(self.angles - np.pi / 2)]).T
        self.pts_inner = (self.cam_radii[:, 0] * [np.cos(self.input_angles),
                                                  np.sin(self.input_angles)
                                                  ]).T
        self.pts_outer = (self.cam_radii[:, 1] * [np.cos(self.input_angles), # SS: this had a -np.pi/2 here before (4/30/24)
                                                  np.sin(self.input_angles)
                                                  ]).T
        """ if plot:
            plt.figure()
            plt.scatter(self.pts_inner[:,0], self.pts_inner[:,1], lw = 2)
            plt.scatter(self.pts_outer[:,0], self.pts_outer[:,1], lw = 2)
            plt.legend(['inner cam','outer cam'])
            plt.axis('equal')
            plt.title('before convex function plot')
            plt.show() """
        
        # Calculate points of convex hull for each cam. Results are cam radii
        # (polar coordinates).
        cam_radii_inner = self.convex_cam_pts(self.pts_inner, plot)
        cam_radii_outer = self.convex_cam_pts(self.pts_outer, plot)
        self.cam_radii = np.vstack((cam_radii_inner, cam_radii_outer)).T

        if plot:
            """ fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.scatter(self.angles, self.cam_radii[:, 0])
            ax.scatter(self.angles, self.cam_radii[:, 1] )
            ax.grid(True)
            ax.set_title('Polar interpolated convex cam points')
            plt.show() """

        # Determine stroke length achieved by flexing knee from upright to 90
        # degrees. This calculation is based on measurements taken in Song et
        # al., 2022(https://doi.org/10.1177/15589250221138546), linearly scaled
        # to the user's height.
        unstretch_len = 0.05
        stretch_pct = 1.016
        user_height_ref = 1.67
        stroke = user_height / user_height_ref * unstretch_len * stretch_pct

        # Scale cam size to calculated stroke length. Iteratively check that the
        # minimum radius is not violated.
        ratio = 0
        threshold = 0.01
        r_min = 0.25*25.4
        while abs(ratio - 1) > threshold:
            self.x_cable = np.cumsum(self.cam_radii[:, 0] * 2*np.pi / self.n_interp)
            ratio = stroke / self.x_cable[self.sit_ind]
            self.cam_radii *= ratio
            self.x_cable *= ratio
            for cam in self.cam_radii:
                for point in cam:
                    if point < r_min:
                        point = r_min
        radius_max = np.max(self.cam_radii)

        # Rotate values in outer cam to account for where elastic band leaves
        # surface relative to where cable leaves surface.
        outer_1 = self.cam_radii[:self.offset_ind, 1]
        outer_2 = self.cam_radii[self.offset_ind:, 1]
        outer_rotated = np.concatenate((outer_2, outer_1))
        self.cam_radii[:, 1] = outer_rotated

        if plot:
            self.plot_cams(self.cam_radii, index)

        # Convert cam points to Cartesian space and return the final cam shapes.
        self.pts_inner = (self.cam_radii[:, 0] * [np.cos(self.angles),
                                                 np.sin(self.angles)
                                                 ]).T
        self.pts_outer = (self.cam_radii[:, 1] * [np.cos(self.angles),
                                                  np.sin(self.angles)
                                                  ]).T
    
        """ if plot:
            self.plot_cams_cartesian(self.pts_inner, self.pts_outer) """
        return self.pts_inner, self.pts_outer, radius_max

    def convex_cam_pts(self, points, plot=False):
        """
        Compute the convex cam points given a set of input points.

        Parameters:
        - points: numpy array
            The input points used to compute the convex cam points.
        - plot: bool, optional
            If True, plot the intermediate steps of the computation. Default is False.

        Returns:
        - numpy array
            The computed convex cam points in polar coordinates.

        This function computes the convex cam points by performing the following steps:
        1. Compute the convex hull of the input points.
        2. Convert the convex points into polar coordinates and interpolate between them.
        3. Convert the interpolated points back into Cartesian coordinates.
        4. Split the cam points into sub-arrays that monotonically increase or decrease in x.
        5. Flip any segments that are monotonically decreasing.
        5. Interpolate between the points in each sub-array.
        6. Convert the interpolated points into polar coordinates.
        7. Linearly interpolate between the points in polar coordinates.
        8. Return the computed convex cam points in polar coordinates.

        If the `plot` parameter is set to True, the function will also generate plots
        to visualize the intermediate steps of the computation.
        """
        # Compute the convex hull of the input points.
        hull = ConvexHull(points, incremental=True)
        points = np.vstack((points[hull.vertices, 0], points[hull.vertices, 1])).T

        # Convert convex cam points to polar coordinates in range [0, 2*pi].
        radii_convex = np.empty(points.shape[0])
        for ind, point in enumerate(points):
            dist = np.linalg.norm(point)
            radii_convex[ind] = dist
        angles_convex = np.arctan2(points[:, 1], points[:, 0])
        for ind, angle in enumerate(angles_convex):
            if angle < 0:
                angles_convex[ind] += 2*np.pi

        # Sort points in ascending order of angles.
        sorted_inds = np.argsort(angles_convex)
        angles_convex = angles_convex[sorted_inds]
        radii_convex = radii_convex[sorted_inds]

        # Pad arrays to ensure continuity at end points.
        radii_convex = np.append(radii_convex, radii_convex[0])
        angles_convex = np.append(angles_convex, angles_convex[0] + 2*np.pi)
        radii_convex = np.insert(radii_convex, 0, radii_convex[-2])
        angles_convex = np.insert(angles_convex, 0, angles_convex[-2] - 2*np.pi)

        # Interpolate between points in polar coordinates.
        polar_cam_fcn = interpolate.interp1d(angles_convex,
                                             radii_convex,
                                             kind='quadratic',
                                             fill_value='extrapolate')
        radii_interp = polar_cam_fcn(self.angles)

        # Convert to Cartesian coordinates and calculate convex hull again.
        points = (radii_interp * [np.cos(self.angles), np.sin(self.angles)]).T
        hull = ConvexHull(points, incremental=True)
        points = np.vstack((points[hull.vertices, 0], points[hull.vertices, 1])).T

        if plot:
            """ plt.figure()
            plt.scatter(points[:, 0], points[:, 1])
            plt.title('Cam points after second convex hull function')
            plt.axis('equal')
            plt.show() """

        # Split Cartesian cam points into 2 or 3 sub-arrays that each
        # monotonically increase or decrease in x.
        break_ind = []
        for ind in range(len(points) - 2):
            if (math.copysign(1, np.diff(points[ind : ind+2, 0]))
                != math.copysign(1, np.diff(points[ind+1 : ind+3, 0]))):
                    break_ind.append(ind+2)
        points_array = []
        points_array.append(points[:break_ind[0], :])
        if len(break_ind) == 2:
            points_array.append(points[break_ind[0]:break_ind[1], :])
            points_array.append(points[break_ind[1]:, :])
        elif len(break_ind) == 1:
            points_array.append(points[break_ind[0]:, :])

        # Flip any segments that are monotonically decreasing, which is required
        # by numpy.interp(). Then, linearly interpolate between 100 points over
        # the range of each sub-array and between the sub-arrays?. This fills
        # any gaps that resulted from the convex hull operation.
        max_gap = 0.1
        # gap_inds = []
        interp_array = np.empty((0, 2))
        for ind, sub_array in enumerate(points_array):
            # Identify any large gaps between sub-arrays and inerpolate over them.
            if ind != len(points_array) - 1:
                if math.dist(sub_array[-1], points_array[ind+1][0]) > max_gap:
                    # gap_inds.append(ind)
                    x1 = sub_array[-1, 0]
                    x2 = points_array[ind+1][0, 0]
                    y1 = sub_array[-1,1]
                    y2 = points_array[ind+1][0,1]
                    if x1 > x2:
                        x = np.linspace(x2, x1, 100)
                        y = np.interp(x, [x2, x1], [y2, y1])
                    else:
                        x = np.linspace(x1, x2, 100)
                        y = np.interp(x, [x1, x2], [y1, y2])
                    interp_array = np.concatenate((interp_array, np.array([x, y]).T))
            else:
                if math.dist(sub_array[-1], points_array[0][0]) > max_gap:
                    # gap_inds.append(ind)
                    x1 = sub_array[-1, 0]
                    x2 = points_array[0][0, 0]
                    y1 = sub_array[-1,1]
                    y2 = points_array[0][0,1]
                    if x1 > x2:
                        x = np.linspace(x2, x1, 100)
                        y = np.interp(x, [x2, x1], [y2, y1])
                    else:
                        x = np.linspace(x1, x2, 100)
                        y = np.interp(x, [x1, x2], [y1, y2])
                    interp_array = np.concatenate((interp_array, np.array([x, y]).T))
        for ind, sub_array in enumerate(points_array):
            # Flip segments that are monotonically decreasing.
            if not np.all(np.diff(sub_array[:, 0]) > 0):
                points_array[ind] = np.flip(points_array[ind], 0)
        for ind, sub_array in enumerate(points_array):
            x = np.linspace(sub_array[0, 0], sub_array[-1, 0], 100)
            y = np.interp(x, sub_array[:, 0], sub_array[:, 1])
            """ if ind == 0:
                interp_array = np.array([x, y]).T
            else:
                interp_array = np.concatenate((interp_array, np.array([x, y]).T)) """
            interp_array = np.concatenate((interp_array, np.array([x, y]).T))
            # if np.abs(points_array[ind+1][0,0] - sub_array[-1,0]) > max_gap:
            # if ind != len(points_array) - 1:
            #     if math.dist(sub_array[-1], points_array[ind+1][0]) > max_gap:
            #         x = np.linspace(sub_array[-1, 0], points_array[ind+1][0, 0], 100)
            #         y = np.interp(x, [sub_array[-1, 0], points_array[ind+1][0, 0]],
            #                     [sub_array[-1, 1], points_array[ind+1][0,1]])
            #         interp_array = np.concatenate((interp_array, np.array([x, y]).T))
            # else:
            #     if math.dist(sub_array[-1], points_array[0][0]) > max_gap:
            #         x = np.linspace(sub_array[-1, 0], points_array[0][0, 0], 100)
            #         y = np.interp(x, [sub_array[-1, 0], points_array[0][0, 0]],
            #                     [sub_array[-1, 1], points_array[0][0,1]])
            #         interp_array = np.concatenate((interp_array, np.array([x, y]).T))

        if plot:
            """ plt.figure()
            plt.scatter(interp_array[:, 0], interp_array[:, 1])
            plt.title('Cartesian interpolated convex cam points')
            plt.axis('equal')
            plt.show() """

        # Convert interpolated convex cam points into polar coordinates in range
        # [0, 2*pi].
        radii_interp = np.empty(interp_array.shape[0])
        for ind, point in enumerate(interp_array):
            dist = np.linalg.norm(point)
            radii_interp[ind] = dist
        angles_interp = np.arctan2(interp_array[:, 1], interp_array[:, 0])
        for ind, angle in enumerate(angles_interp):
            if angle < 0:
                angles_interp[ind] += 2*np.pi

        # # Sort points in ascending order of angles.
        # sorted_inds = np.argsort(angles_interp)
        # angles_interp = angles_interp[sorted_inds]
        # radii_interp = radii_interp[sorted_inds]

        # Remove duplicates to avoid errors in interpolation.
        dup_removed = self.remove_duplicates(angles_interp, radii_interp)
        angles_interp = dup_removed[0]
        radii_interp = dup_removed[1]

        # Pad arrays to ensure continuity at end points. Linearly interpolate
        # between points in polar coordinates and return.
        radii_interp = np.append(radii_interp, radii_interp[0])
        angles_interp = np.append(angles_interp, angles_interp[0] + 2*np.pi)
        radii_interp = np.insert(radii_interp, 0, radii_interp[-2])
        angles_interp = np.insert(angles_interp, 0, angles_interp[-2] - 2*np.pi)
        polar_cam_fcn = interpolate.interp1d(angles_interp,
                                             radii_interp,
                                             kind='cubic',
                                             fill_value='extrapolate')
        
        if np.any(np.isnan(polar_cam_fcn(self.angles))):
            print("Warning: NaN in cam radii")

        return polar_cam_fcn(self.angles)

    def calc_forces_percentages(self, angle_data, torque=False, plot=False,
                                index=0):
        """Calculate the force profiles vs stance percentage given the solved
        cam radii. The order of causality is:
        stance percentage -> knee angle -> cable displacement -> cam angle -> 
        elastic band displacement -> elastic band force -> cable force
        """
        # Calculate elastic band tension from displacement and experimental 
        # characterization assuming linear stiffness. Calculate cable tension
        # from elastic band tension and gear ratios. Values are in N, m, & N/m.
        x_elastic = 2 * np.cumsum(self.cam_radii[:, 1] * 2*np.pi / self.n_interp)
        k_elastic = 110
        f_elastic = k_elastic * x_elastic
        f_cable =  f_elastic * self.cam_radii[:, 1] / self.cam_radii[:, 0]
        # print(f"Elastic band stiffness {k_elastic} N/m; Pulling distance {x_elastic[self.sit_ind]} m")

        # Create a function that relates knee angle to stance percentage using
        # angle_data, accounting for nonlinear change in knee angle while 
        # standing. Percentages go from ~0% to ~100% going from sit to stand.
        knee_array = np.loadtxt(angle_data, delimiter = ',', ndmin = 2)
        percentages = knee_array[:, 0]
        knee_angles = knee_array[:, 1]
        knee_fcn = interpolate.interp1d(percentages, knee_angles,
                                        kind='cubic', fill_value='extrapolate')
        percentages = np.linspace(0, 100, self.sit_ind + 1)
        knee_angles = knee_fcn(percentages)

        # Flip arrays to go from standing (100%) to sitting (0%), since cable
        # displacement is 0 at 100% stance. Then scale cable displacement by the
        # normalized knee angle.
        knee_angles = np.flip(knee_angles)
        percentages = np.flip(percentages)
        x_cable_scaled = (self.x_cable[self.sit_ind]
                          * (1 - (knee_angles-np.min(knee_angles))
                                 / np.ptp(knee_angles))) 

        if plot:
            """
            plt.figure()
            plt.plot(knee_angles, x_cable_scaled)
            plt.xlabel('knee angle (degree)')
            plt.ylabel('cable displacement, xc (m)')
            plt.title('knee angle and cable displacement linear relationship')
            

            # plot stance percentage vs. cable displacement
            plt.figure()
            plt.plot(percentages, x_cable_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cable Displacement (m)')
            """

        if np.any(np.isnan(self.x_cable)):
            print("Warning: NaN in cam radii")
        
        # Remove duplicate points to avoid errors in interpolation.
        dup_removed = self.remove_duplicates(self.x_cable, self.angles, x_elastic)
        self.x_cable = dup_removed[0]
        self.angles = dup_removed[1]
        x_elastic = dup_removed[2]
        # repeat_dict = defaultdict(list)
        # for ind, point in enumerate(self.x_cable):
        #     repeat_dict[point].append(ind)
        # repeat_dict = {k:v for k,v in repeat_dict.items() if len(v)>1}
        # repeat_list = list(repeat_dict.values())
        # repeat_ind = np.array(sum(repeat_list, []))
        # if len(repeat_dict) > 0:
        #     print("Repeated values and indices: ", repeat_dict)
        #     self.x_cable = np.delete(self.x_cable, repeat_ind)
        #     self.angles = np.delete(self.angles, repeat_ind)
        #     x_elastic = np.delete(x_elastic, repeat_ind)

        # Find force in elastic band from stance percentage as input. Use cable
        # displacement scaled to knee angle to find cam angle. Then use cam
        # angle to find elastic band displacement and force. Finally, find cable
        # force based on the elastic tension and gear ratio at each point.
        cable_fcn = interpolate.interp1d(self.x_cable, self.angles,
                                         kind='cubic', fill_value='extrapolate')
        elastic_fcn = interpolate.interp1d(self.angles, x_elastic,
                                           kind='cubic', fill_value='extrapolate')
        angle_scaled = cable_fcn(x_cable_scaled)
        x_elastic_scaled = elastic_fcn(angle_scaled)
        f_cable_scaled = (self.cam_radii[: self.sit_ind+1, 1] /
                          self.cam_radii[: self.sit_ind+1, 0]
                          * k_elastic * x_elastic_scaled)
        
        # Plot and save values.
        if plot:
            """
            # plot cam angle vs. stance percentage
            plt.figure()
            plt.plot(percentages, angle_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cam Angle (rad)')

            # plot elastic band displacement vs. stance percentage
            plt.figure()
            plt.plot(percentages, x_elastic_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Elastic displacement (m)')
            """
            
            # plot cable force vs. stance percentage
            plt.figure()
            plt.plot(percentages, f_cable_scaled,
                     label='Cable Force', linewidth=3)
            plt.legend(loc='upper right')
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Force (N)')
            plt.xlim([0, 100])
            plt.ylim([0, 250])

            newpath = 'results/force_plots/force_plots_' + self.dateStr
            if not path.exists(newpath):
                makedirs(newpath)
            filename = newpath + '/force_plot_' + str(index) + '.png'
            plt.savefig(filename, dpi=300)
            if torque:
                plt.figure()
                plt.scatter(self.angles * 360 / (2*np.pi),
                            f_elastic * self.cam_radii[:, 0])
                plt.xlabel('Angle (deg)')
                plt.ylabel('Torque (N-m)')
            plt.show()

            np.savetxt('force_output.csv',
                       np.stack((self.angles, self.x_cable, f_elastic, f_cable),
                                axis=1),
                                header='angles(rad), pulling distance(m), \
                                    elastic band force(n), cable force(n)')

        return f_cable_scaled, percentages
    
    def remove_duplicates(self, x, y=None, z=None):
        """
        Remove duplicate values from the given array `x` and corresponding values from `y`.

        Args:
            x (ndarray): Array of values.
            y (ndarray or None): Array of corresponding values or None.

        Returns:
            ndarray: Updated array `x` with duplicate values removed.
            ndarray or None: Updated array `y` with corresponding values removed, or None if `y` is None.

        """
        # Store indices of repeated values in a dictionary.
        repeat_dict = defaultdict(list)
        for ind, point in enumerate(x):
            repeat_dict[point].append(ind)
        repeat_dict = {k:v for k,v in repeat_dict.items() if len(v)>1}
        repeat_list = list(repeat_dict.values())
        repeat_ind = np.array(sum(repeat_list, []))
        
        # If there are repeated values, print and remove them from the arrays.
        if len(repeat_dict) > 0:
            # print("Repeated values and indices: ", repeat_dict)
            x = np.delete(x, repeat_ind)
            if y is not None:
                y = np.delete(y, repeat_ind)
            if z is not None:
                z = np.delete(z, repeat_ind)
        
        return x, y, z

    def plot_cams(self, cam_radii=0, index=0):
        """
        Plots the cam points in polar coordinates.
        Plots the cam points in polar coordinates.
        """
        r = self.cam_radii[:, 0]
        R = self.cam_radii[:, 1]
        plt.figure()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(self.angles, r)
        ax.plot(self.angles, R)
        ax.set_rlabel_position(-22.5)
        ax.grid(True)
        ax.set_title(f"Cam Demonstration, min radius={np.min(self.cam_radii)}, max radius={np.max(self.cam_radii)}", va='bottom')

        newpath = 'results/cam_plots/cam_plots_' + self.dateStr
        if not path.exists(newpath):
            makedirs(newpath)
        filename = newpath + '/cam_plot_' + str(index) + '.png'
        plt.savefig(filename, dpi=300)
        plt.show()

    def plot_cams_cartesian(self, pts_inner, pts_outer):
        """
        Plots the cam points in Cartesian coordinates.
        Plots the cam points in Cartesian coordinates.
        """
        plt.figure()
        plt.plot(pts_inner[:, 0], pts_inner[:, 1], lw = 2)
        plt.plot(pts_outer[:, 0], pts_outer[:, 1], lw = 2)
        plt.axis('equal')
        plt.legend(['inner cam','outer cam'])
        plt.title('Cam Points in Cartesian Space')
        plt.show()