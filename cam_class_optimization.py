"""Using python 3.10"""
import math
from datetime import date
from os import path, makedirs

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from scipy.interpolate import make_interp_spline, BSpline
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.optimize import Bounds
from collections import defaultdict

font = {'size': 14}
matplotlib.rc('font', **font)


class CamGeneration:
    """Class that generates points on the cams for every degree given
    torque/gear ratios. Can also return force and energy output given stiffness.
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
        energy storage band acts relative to the force transmission cable
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
        n_eval = 0

    def calculate_cam_radii(self, user_height=1.67, k_elastic=0,
                            plot=False, index=0):
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
        # Check if any of the radii are smaller than the minimum allowed radius,
        # in meters, and adjust as needed.
        radius_min = 0.005
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
        self.pts_inner = (self.cam_radii[:, 0] * [np.cos(self.input_angles),
                                                  np.sin(self.input_angles)
                                                  ]).T
        self.pts_outer = (self.cam_radii[:, 1] * [np.cos(self.input_angles),
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
        # al., 2022 (https://doi.org/10.1177/15589250221138546), linearly scaled
        # to the user's height.
        unstretch_len = 0.05
        stretch_pct = 1.016
        user_height_ref = 1.67
        stroke = user_height / user_height_ref * unstretch_len * stretch_pct

        # Scale cam size to calculated stroke length. Iteratively check that the
        # minimum radius, in meters, is not violated.
        ratio = 0
        threshold = 0.01
        r_min = 0.25*.0254
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
        self.cam_offset = self.cam_radii[:, 1].copy()
        outer_1 = self.cam_radii[:self.offset_ind, 1]
        outer_2 = self.cam_radii[self.offset_ind:, 1]
        cam_rotated = np.concatenate((outer_2, outer_1))
        self.cam_radii[:, 1] = cam_rotated

        if plot:
            self.plot_cams(self.cam_radii, k_elastic, index)

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

        If the 'plot' parameter is set to True, the function will also generate plots
        to visualize the intermediate steps of the computation.
        """
        # Compute the convex hull of the input points.
        hull = ConvexHull(points, incremental=True)
        points = np.vstack((points[hull.vertices, 0], points[hull.vertices, 1])).T

        # Convert convex cam points to polar coordinates in range [0, 2*pi].
        # radii_convex = np.empty(points.shape[0])
        # for ind, point in enumerate(points):
        #     dist = np.linalg.norm(point)
        #     radii_convex[ind] = dist
        # angles_convex = np.arctan2(points[:, 1], points[:, 0])
        # for ind, angle in enumerate(angles_convex):
        #     if angle < 0:
        #         angles_convex[ind] += 2*np.pi
        radii_convex, angles_convex = self.to_polar(points)

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
        interp_array = np.empty((0, 2))
        for ind, sub_array in enumerate(points_array):
            # Identify any large gaps between sub-arrays and inerpolate over them.
            if ind != len(points_array) - 1:
                if math.dist(sub_array[-1], points_array[ind+1][0]) > max_gap:
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
            interp_array = np.concatenate((interp_array, np.array([x, y]).T))

        if plot:
            """ plt.figure()
            plt.scatter(interp_array[:, 0], interp_array[:, 1])
            plt.title('Cartesian interpolated convex cam points')
            plt.axis('equal')
            plt.show() """

        # Convert interpolated convex cam points into polar coordinates in range
        # [0, 2*pi].
        # radii_interp = np.empty(interp_array.shape[0])
        # for ind, point in enumerate(interp_array):
        #     dist = np.linalg.norm(point)
        #     radii_interp[ind] = dist
        # angles_interp = np.arctan2(interp_array[:, 1], interp_array[:, 0])
        # for ind, angle in enumerate(angles_interp):
        #     if angle < 0:
        #         angles_interp[ind] += 2*np.pi
        radii_interp, angles_interp = self.to_polar(interp_array)

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

    def calc_forces(self, angle_data, k_elastic, torque=False,
                                plot=False, index=0):
        """Calculate the force profiles vs stance percentage given the solved
        cam radii. The order of causality is:
        stance percentage -> knee angle -> cable displacement -> cam angle -> 
        elastic band displacement -> elastic band force -> cable force
        """
        # Calculate elastic band tension from displacement around storage cam 
        # and using experimental stiffness characterization, simplified to be
        # linear. Note that storage cam angle is offset. Calculate cable tension
        # from elastic band tension and gear ratios. Values are in N, m, & N/m.
        x_elastic = np.cumsum(self.cam_offset * 2*np.pi / self.n_interp)
        f_elastic = k_elastic * x_elastic
        f_cable =  f_elastic * self.cam_offset / self.cam_radii[:, 0]
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

        # Use cable displacement scaled to knee angle to find cam angle.
        # Then use cam angle to find elastic band displacement and force.
        # Finally, find cable force based on the elastic tension and gear ratio
        # at each point.
        cable_fcn = interpolate.interp1d(self.x_cable, self.angles,
                                         kind='cubic', fill_value='extrapolate')
        elastic_fcn = interpolate.interp1d(self.angles, x_elastic,
                                           kind='cubic', fill_value='extrapolate')
        angle_scaled = cable_fcn(x_cable_scaled)
        x_elastic_scaled = elastic_fcn(angle_scaled)
        f_cable_scaled = (self.cam_offset[: self.sit_ind+1] /
                          self.cam_radii[: self.sit_ind+1, 0]
                          * k_elastic * x_elastic_scaled)
        # f_cable_scaled = (self.cam_radii[: self.sit_ind+1, 1] /
        #                   self.cam_radii[: self.sit_ind+1, 0]
        #                   * k_elastic * x_elastic_scaled)
        
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

            filepath = 'results/force_plots/force_plots_' + self.dateStr
            if not path.exists(filepath):
                makedirs(filepath)
            filename = filepath + '/force_plot_' + str(index) + '.png'
            plt.savefig(filename, dpi=300)
            if torque:
                plt.figure()
                plt.scatter(self.angles * 360 / (2*np.pi),
                            f_elastic * self.cam_radii[:, 0])
                plt.xlabel('Angle (deg)')
                plt.ylabel('Torque (N-m)')
            plt.show()

            np.savetxt(filepath + '/_force_output.csv',
                       np.stack((self.angles, self.x_cable, f_elastic, f_cable),
                                axis=1),
                                header='angles(rad), pulling distance(m), \
                                    elastic band force(n), cable force(n)')
            
        return f_cable_scaled, percentages
    
    def generate_sit_cam(self, folder_path, cam_idx, n_params=6):
        # Load sit-to-stand cam point cloud in Cartesian coordinates.
        points_stand = np.array(pd.read_csv(f'{folder_path}inner_{cam_idx}.txt',
                                         header=None, usecols=[0,1])) / 1000

        # Convert point cloud to polar coordinates.
        radii_stand, angles_stand = self.to_polar(points_stand)

        def interp_cam(x):
            # Interpolate between radii keypoints
            key_angles = x[:n_params]
            key_radii = x[:n_params]
            spline_fcn = interpolate.interp1d(key_angles, key_radii,
                                              kind='cubic', fill_value='extrapolate')
            radii_interp = spline_fcn(angles_stand[:self.sit_ind+1])
            return radii_interp
        
        # Bounds:
        # Constrain the radii at sit and stand angles to be equal to those of
        # the sit-to-stand cam, within a threshold. Set lower and upper bounds
        # at all other angles to reasonable values.
        thresh_end = 0.001
        lb_ang = np.zeros(n_params)
        lb_rad = np.ones(n_params) * 0.00635
        lb_rad[0] = radii_stand[0] - thresh_end
        lb_rad[-1] = radii_stand[self.sit_ind] - thresh_end
        lb = np.concatenate((lb_ang, lb_rad))

        ub_ang = np.ones(n_params) * self.angles[self.sit_ind]
        ub_rad = np.ones(n_params) * 0.10
        ub_rad[0] = radii_stand[0] + thresh_end
        ub_rad[-1] = radii_stand[self.sit_ind] + thresh_end
        ub = np.concatenate((ub_ang, ub_rad))
        bounds = Bounds(lb, ub, keep_feasible=False)

        # Nonlinear contraint:
        # Constrain the path lengths to be equal within a threshold.
        x_cable_stand = np.cumsum(radii_stand * 2*np.pi / self.n_interp)
        path_length_stand = x_cable_stand[self.sit_ind]
        thresh_path = 0.01
        ub_path = path_length_stand + thresh_path
        lb_path = path_length_stand - thresh_path
        def cons_path(x):
            radii_interp = interp_cam(x)
            x_cable_sit = np.cumsum(radii_interp * 2*np.pi / self.n_interp)
            path_length_sit = x_cable_sit[-1]
            return path_length_sit
        path_constraint = NonlinearConstraint(cons_path, lb_path, ub_path,
                                              keep_feasible=False)

        # Objective:
        # Minimize sum of radii first derivative.
        def objective(x):
            radii_interp = interp_cam(x)
            return np.sum(np.abs(np.diff(radii_interp))) / 1000

        # Initial guess: 
        # x0 = np.linspace(radii_stand[0], radii_stand[self.sit_ind],
        #                  num=self.sit_ind+1) # a straight line between the sit-to-stand cam radii end points.
        rmax = (path_length_stand * (self.n_interp/(np.pi * self.sit_ind)) 
                - 0.5 * (radii_stand[0] + radii_stand[self.sit_ind]))
        rad1 = np.linspace(radii_stand[0], rmax, num=int(n_params/2))
        rad2 = np.linspace(rmax, radii_stand[self.sit_ind], num=int(n_params/2))
        rad0 = np.concatenate((rad1, rad2))
        ang0 = np.linspace(angles_stand[0], angles_stand[self.sit_ind], num=n_params)
        x0 = np.concatenate((ang0, rad0))
        path_length_init = np.cumsum(rad0 * 2*np.pi / n_params)[-1]
        print(f"INTIAL GUESS PATH LENGTH: {path_length_init}")
        print(f"INITIAL GUESS OBJECTIVE VALUE: {objective(x0)}")

        plt.figure()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(ang0, rad0, label='Initial guess')
        plt.show()

        # Callback: check the optimization progress intermittently.
        global n_eval
        n_eval = 1
        def callback_fun(x):
            global n_eval
            if n_eval % 10 == 0:
                print(f"Optimization iteration {n_eval}. Objective = {objective(x)}")
                print(f"Standing path length: {path_length_stand}")
                radii_interp = interp_cam(x)
                print(f"Current path length: {np.cumsum(radii_interp * 2*np.pi / self.n_interp)[-1]}")
                plt.figure()
                fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
                ax.plot(angles_stand[:self.sit_ind+1], radii_interp, label=f"Optimization iteration {n_eval}")
                plt.title(f"Optimization iteration {n_eval}")
                plt.show()
            n_eval += 1

        # Minimize: perform the optimization.
        result = minimize(objective, x0, method='COBYLA', 
                          options={'maxiter': 10000},
                          constraints=path_constraint,
                          bounds=bounds,
                          callback=callback_fun)
        
        # Report and plot the results.
        radii_interp = interp_cam(result.x)
        opt_path = np.cumsum(radii_interp * 2*np.pi / self.n_interp)
        print("Sit cam path length: ", opt_path[-1])
        print("Sit cam start radius: ", radii_interp[0])
        print("Sit cam end radius: ", radii_interp[-1])

        plt.figure()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(angles_stand[:self.sit_ind+1], radii_interp, label='Sitting cam')
        plt.show()

        print("check it out")
        return result
    
    def remove_duplicates(self, x, y=None, z=None):
        """
        Remove duplicate values from the given array 'x' and corresponding
        values from 'y'.

        Args:
            x (ndarray): Array of values.
            y (ndarray or None): Array of corresponding values or None.

        Returns:
            ndarray: Updated array 'x' with duplicate values removed.
            ndarray or None: Updated array 'y' with corresponding values removed,
            or None if 'y' is None.

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
    
    def to_polar(self, points):
        """
        Convert the given Cartesian points to polar coordinates.

        Args:
            points (ndarray): Array of Cartesian points.

        Returns:
            ndarray: Array of radii.
            ndarray: Array of angles.
        """
        # Calculate radii and angles of the points.
        radii = np.linalg.norm(points, axis=1)
        # radii = np.empty(points.shape[0])
        # for ind, point in enumerate(points):
            # dist = np.linalg.norm(point)
            # radii[ind] = dist
        angles = np.arctan2(points[:, 1], points[:, 0])
        for ind, angle in enumerate(angles):
            if angle < 0:
                angles[ind] += 2*np.pi
        return radii, angles

    def plot_cams(self, cam_radii=0, k_elastic=0, index=0):
        """
        Plots the cam points in polar coordinates.
        """
        r = 100*self.cam_radii[:, 0]
        R = 100*self.cam_radii[:, 1]
        plt.figure()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(self.angles, r, label='Transmission cam')
        ax.plot(self.angles, R, label='Storage cam')
        ax.legend(loc='lower right')
        ax.set_xticklabels([])
        ax.grid(True)
        ax.set_title(f"""Cam shapes\nmin radius={100*np.min(self.cam_radii):.2f} cm, max radius={100*np.max(self.cam_radii):.2f} cm\nK={k_elastic} N/m""")

        filepath = 'results/cam_plots/cam_plots_' + self.dateStr
        if not path.exists(filepath):
            makedirs(filepath)
        filename = filepath + '/cam_plot_' + str(index) + '.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.show()

    def plot_cams_cartesian(self, pts_inner, pts_outer):
        """
        Plots the cam points in Cartesian coordinates.
        """
        plt.figure()
        plt.plot(pts_inner[:, 0], pts_inner[:, 1], lw = 2)
        plt.plot(pts_outer[:, 0], pts_outer[:, 1], lw = 2)
        plt.axis('equal')
        plt.legend(['inner cam','outer cam'])
        plt.title('Cam Points in Cartesian Space')
        plt.show()