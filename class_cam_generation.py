"""
cam-mech: cam generation utilities

Provides CamGeneration and helper functions for reading CSV data,
computing curve lengths, and converting between Cartesian and polar
coordinates. This module is intended to be used as a library and
avoids doing heavy computation on import.

The public API includes:
- read_xy_from_csv
- curve_length_from_csv
- CamGeneration

"""
from __future__ import annotations

import math
import csv
from datetime import date
from pathlib import Path
from typing import Optional, Tuple
import logging

import numpy as np
import pandas as pd

# Module logger
logger = logging.getLogger(__name__)


def _get_plt():
    """Lazily import matplotlib.pyplot for plotting helpers.

    Raises a RuntimeError if matplotlib is not available.
    """
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        logger.error("matplotlib is required for plotting but is not available: %s", e)
        raise RuntimeError("matplotlib is required for plotting") from e
    return plt

from scipy import interpolate
from scipy.interpolate import make_interp_spline, BSpline
from scipy.spatial import ConvexHull, convex_hull_plot_2d
# from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.optimize import Bounds
# 'minimize' from cobyqa is imported lazily inside optimize routines to avoid importing heavy optional deps at module import time
from scipy.integrate import trapezoid
from collections import defaultdict


def read_xy_from_csv(csv_filepath, x_col=0, y_col=1, delimiter=','):
    """Read ordered x,y points from a CSV file.

    Parameters:
    csv_filepath: str
        Path to a CSV file with at least two columns of numeric data.
    x_col: int, optional
        Column index for x coordinates. Default is 0.
    y_col: int, optional
        Column index for y coordinates. Default is 1.
    delimiter: str, optional
        Field delimiter for the CSV file. Default is ','.

    Returns:
    x: numpy.ndarray
        Array of x coordinate values.
    y: numpy.ndarray
        Array of y coordinate values.
    """
    df = pd.read_csv(csv_filepath, header=None, comment='#', delimiter=delimiter)
    if df.shape[1] < 2:
        raise ValueError("CSV file must contain at least two columns of numeric x and y data.")
    df = df.iloc[:, [x_col, y_col]].apply(pd.to_numeric, errors='coerce')
    if df.isna().any().any():
        raise ValueError("CSV file contains non-numeric values in the x or y columns.")
    x = df.iloc[:, 0].to_numpy(dtype=float)
    y = df.iloc[:, 1].to_numpy(dtype=float)
    return x, y


def curve_length_from_csv(csv_filepath, n_eval=1000, spline_order=3,
                          x_col=0, y_col=1, delimiter=','):
    """Compute cumulative spline lengths for ordered CSV points.

    The function reads x/y points from the provided CSV file, constructs a
    parametric spline curve that passes through the points in order, and
    returns the cumulative arc length at each original point.

    Parameters:
    csv_filepath: str
        Path to the CSV file containing ordered x and y coordinates.
    n_eval: int, optional
        Number of points to evaluate along the spline for length estimation.
    spline_order: int, optional
        Order of the spline. Default is 3 (cubic) and is automatically reduced
        when there are too few data points.
    x_col: int, optional
        Index of the column containing x coordinates.
    y_col: int, optional
        Index of the column containing y coordinates.
    delimiter: str, optional
        Field delimiter used by the CSV file.

    Returns:
    numpy.ndarray
        Array of cumulative arc lengths at each data point, in the same order
        as the input points.
    """
    x, y = read_xy_from_csv(csv_filepath, x_col=x_col, y_col=y_col,
                            delimiter=delimiter)
    if x.size < 2:
        return np.zeros(x.shape, dtype=float)
    t = np.linspace(0.0, 1.0, x.size)
    k = min(spline_order, x.size - 1)
    spline = make_interp_spline(t, np.column_stack((x, y)), k=k, axis=0)

    n_eval = max(n_eval, x.size)
    t_eval = np.linspace(t[0], t[-1], n_eval)
    derivative = spline.derivative()
    dx_dy = derivative(t_eval)
    speeds = np.hypot(dx_dy[:, 0], dx_dy[:, 1])

    dt = t_eval[1] - t_eval[0]
    cumulative = np.concatenate(([0.0], np.cumsum((speeds[:-1] + speeds[1:]) * 0.5 * dt)))
    return np.interp(t, t_eval, cumulative)


font = {'size': 14}
try:
    import matplotlib
    matplotlib.rc('font', **font)
except Exception:
    logger.debug("matplotlib not available; skipping font rc setup")



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

    def calculate_cam_radii(self, user_height: float = 1.67, k_elastic: float = 0,
                            plot: bool = False, save: bool = False, save_dir: Optional[str] = None, index: int = 0):
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
            plt = _get_plt()
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
        cam_radii_inner = self.convex_cam_pts(self.pts_inner, self.angles, plot)
        cam_radii_outer = self.convex_cam_pts(self.pts_outer, self.angles, plot)
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
        r_min = 0.25 * .0254 # convert from inches to meters
        while abs(ratio - 1) > threshold:
            x_cable = np.cumsum(self.cam_radii[:, 0] * 2*np.pi / self.n_interp)
            ratio = stroke / x_cable[self.sit_ind]
            self.cam_radii *= ratio
            x_cable *= ratio
            for cam in self.cam_radii:
                for point in cam:
                    if point < r_min:
                        point = r_min
        radius_max = np.max(self.cam_radii)

        # Rotate values in outer cam to account for where elastic band leaves
        # surface relative to where cable leaves surface.
        cam_rotated, cam_original = self.rotate_cam(self.cam_radii[:, 1])        # self.cam_offset = self.cam_radii[:, 1].copy()
        # outer_1 = self.cam_radii[:self.offset_ind, 1]
        # outer_2 = self.cam_radii[self.offset_ind:, 1]
        # cam_rotated = np.concatenate((outer_2, outer_1))
        self.cam_radii[:, 1] = cam_rotated

        if plot:
            # Plot; saving is opt-in via save flag and save_dir
            self.plot_cams(self.cam_radii, k_elastic, index, save=save, save_dir=save_dir)

        # Convert cam points to Cartesian space and return the final cam shapes.
        self.pts_inner = (self.cam_radii[:, 0] * [np.cos(self.angles),
                                                  np.sin(self.angles)
                                                  ]).T
        self.pts_outer = (self.cam_radii[:, 1] * [np.cos(self.angles),
                                                  np.sin(self.angles)
                                                  ]).T
    
        """ if plot:
            self.plot_cams_cartesian(self.pts_inner, self.pts_outer) """
        return self.cam_radii, self.pts_inner, self.pts_outer, radius_max, x_cable

    def convex_cam_pts(self, points, angles, plot=False):
        """
        Compute the convex cam points given a set of input points.

        Parameters:
        - points: numpy array
            The input points used to compute the convex cam points.
        - angles: numpy array
            The angles used to generate the radii of the convex cam.
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
        radii_interp = polar_cam_fcn(angles)

        # Convert to Cartesian coordinates and calculate convex hull again.
        points = (radii_interp * [np.cos(angles), np.sin(angles)]).T
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
        radii_interp, angles_interp = self.to_polar(interp_array)

        # Remove duplicates to avoid errors in interpolation.
        dup_removed = self.remove_duplicates(angles_interp, radii_interp)
        angles_interp = dup_removed[0]
        radii_interp = dup_removed[1]

        # Pad arrays to ensure continuity at end points. Cubically interpolate
        # between points in polar coordinates and return.
        radii_interp = np.append(radii_interp, radii_interp[0])
        angles_interp = np.append(angles_interp, angles_interp[0] + 2*np.pi)
        radii_interp = np.insert(radii_interp, 0, radii_interp[-2])
        angles_interp = np.insert(angles_interp, 0, angles_interp[-2] - 2*np.pi)
        polar_cam_fcn = interpolate.interp1d(angles_interp,
                                             radii_interp,
                                             kind='cubic',
                                             fill_value='extrapolate')
        
        if np.any(np.isnan(polar_cam_fcn(angles))):
            logger.warning("NaN detected in polar interpolation of cam radii")

        return polar_cam_fcn(angles)
    
    def rotate_cam(self, radii):
        # Rotate values in outer cam to account for where elastic band leaves
        # surface relative to where cable leaves surface.
        cam_original = radii.copy()
        outer_1 = radii[: self.offset_ind]
        outer_2 = radii[self.offset_ind :]
        cam_rotated = np.concatenate((outer_2, outer_1))
        return cam_rotated, cam_original
    
    def derotate_cam(self, radii):
        outer_1 = radii[-self.offset_ind :]
        outer_2 = radii[: -self.offset_ind]
        cam_derotated = np.concatenate((outer_1, outer_2))
        return cam_derotated
    
    def calc_forces(self, cam_radii, x_cable, k_elastic,
                    torque: bool = False, plot: bool = False, save: bool = False, save_dir: Optional[str] = None, index: int = 0):
        """Calculate the force profiles vs stance percentage given the solved
        cam radii. The order of causality is:
        stance percentage -> knee angle -> scaled cable displacement -> cam angle -> 
        elastic band displacement -> elastic band force -> cable force
        """
        # Obtain original (unrotated) storage cam radii
        cam_original = self.derotate_cam(cam_radii[:, 1])

        # Calculate elastic band tension from displacement around storage cam 
        # and using experimental stiffness characterization, simplified to be
        # linear. Note that storage cam angle is offset. Calculate cable tension
        # from elastic band tension and gear ratios. Values are in N, m, & N/m.
        x_elastic = np.cumsum(cam_original * 2*np.pi / self.n_interp)
        f_elastic = k_elastic * x_elastic
        f_cable =  f_elastic * cam_original / cam_radii[:, 0]
        logger.info("Elastic band stiffness %s N/m; Pulling distance %s m", k_elastic, x_elastic[self.sit_ind])

        # Stored energy vs. elastic band displacement
        # E = 0.5 * k_elastic * x_elastic**2
        E = trapezoid(f_elastic, x_elastic)
        logger.info("Total energy stored in elastic band: %s J", E)

        if plot:
            plt = _get_plt()
            plt.figure()
            plt.plot(100*x_cable[:self.sit_ind], f_cable[:self.sit_ind],
                     label='Transmission cable force vs. displacement; unscaled', linewidth=3)
            plt.legend(loc='upper right')
            plt.xlabel('Transmission cable displacement (cm)')
            plt.ylabel('Force (N)')
            plt.show()

            plt.figure()
            plt.plot(100*x_elastic[:self.sit_ind], f_elastic[:self.sit_ind],
                     label='Elastic force vs. displacement; unscaled', linewidth=3)
            plt.legend(loc='upper right')
            plt.xlabel('Elastic band displacement (cm)')
            plt.ylabel('Force (N)')
            plt.show()

        # Generate knee angles corresponding to stance percentage
        percentages = np.linspace(0, 100, self.sit_ind + 1)
        knee_angles = self.knee_angle_to_stance(percentages)

        # Flip knee angles and percentages to correspond with cam angle order,
        # which increases from stand to sit.
        percentages = percentages[::-1]
        knee_angles = knee_angles[::-1]

        # Scale transmission cable displacement by the normalized knee angle.
        # At this point, x_cable_scaled should decrease with index.
        x_cable_scaled = (x_cable[self.sit_ind]
                          * (knee_angles-np.min(knee_angles))
                          / np.ptp(knee_angles)) 

        # Isolate "effective" portions of x_cable (from 0 to 220
        # degrees) and flip order to match the order of knee angles,
        #  which go from sit to stand. This means that x_cable_eff 
        # and x_cable_scaled_eff are now decreasing in value.
        x_cable = x_cable[:self.sit_ind+1]
        x_cable = x_cable[::-1]
        # x_cable_scaled_eff = x_cable_scaled[:self.sit_ind+1]
        # x_cable_scaled_eff = x_cable_scaled_eff[::-1] 

        if plot:
            plt = _get_plt()

            # Plot knee angle vs. cable displacement to verify
            plt.figure()
            plt.plot(knee_angles, 100 * x_cable_scaled)
            plt.xlabel('Knee Angle (degrees)')
            plt.ylabel('Transmission Cable Displacement, xc (cm)')

            # Plot stance percentage vs. cable displacement to verify
            plt.figure()
            plt.plot(percentages, 100 * x_cable_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Transmission Cable Displacement (cm)')

        if np.any(np.isnan(x_cable)):
            logger.warning("NaN detected in x_cable passed to calc_forces")
        
        # Remove duplicate points from cable path information (and corresponding
        # points from the angles and elastic band path) to avoid errors in
        # interpolation.
        dup_removed = self.remove_duplicates(x_cable, self.angles[:x_cable.size], x_elastic)
        x_cable = dup_removed[0]
        self.angles = dup_removed[1]
        x_elastic = dup_removed[2]

        # Create functions relating transmission cable displacement to cam 
        # angle and relating cam angle to storage cable displacement.
        cable_fcn = interpolate.interp1d(x_cable, self.angles[:x_cable.size],
                                         kind='cubic', fill_value='extrapolate')
        elastic_fcn = interpolate.interp1d(self.angles[:x_cable.size], x_elastic[:x_cable.size],
                                           kind='cubic', fill_value='extrapolate')
        
        # Use cable displacement scaled to knee angle to find cam angle as
        # driven by knee rotation.
        angle_scaled = cable_fcn(x_cable_scaled)

        # Use knee-driven cam angle to find storage cable displacement.
        x_elastic_scaled = elastic_fcn(angle_scaled)

        # Determine gear ratio at each point based on the scaled cam angle and the original cam radii.
        gear_ratios = cam_original[: self.sit_ind+1] / cam_radii[: self.sit_ind+1, 0]
        gear_ratio_fcn = interpolate.interp1d(self.angles[:gear_ratios.size], gear_ratios, kind='cubic', fill_value='extrapolate')
        gear_ratios_scaled = gear_ratio_fcn(angle_scaled)

        # Find transmission cable force based on the storage cable displacement,
        # spring stiffness, and gear ratio at each point.
        f_elastic_scaled = k_elastic * x_elastic_scaled
        f_cable_scaled = gear_ratios_scaled * f_elastic_scaled

        # f_cable_scaled = (cam_original[: self.sit_ind+1] /
        #                   cam_radii[: self.sit_ind+1, 0]
        #                   * k_elastic * x_elastic_scaled)
        
        # Plot and save values.
        if plot:
            plt = _get_plt()

            # plot cam angle vs. stance percentage
            plt.figure()
            plt.plot(percentages, angle_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cam Angle (rad)')

            # Plot storage cable displacement vs. stance percentage. Flip 
            # percentages, since storage cable displacement is opposite to
            # transmission cable.
            plt.figure()
            plt.plot(percentages, x_elastic_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Storage cable displacement (m)')

            # # Plot storage cable force vs. displacement.
            # plt.figure()
            # plt.plot(x_elastic_scaled, f_cable_scaled) # Note: this plot doesn't agree with the labels. Believe it was a typo.
            # plt.xlabel('Elastic displacement (m)')
            # plt.ylabel('Elastic force (N)')
            # plt.title('Elastic force vs. displacement)')
            
            # Plot & save transmission cable force vs. stance percentage.
            plt.figure()
            plt.plot(percentages, f_cable_scaled,
                     label='Cable Force', linewidth=3)
            plt.legend(loc='upper right')
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Force (N)')
            plt.xlim([0, 100])
            # plt.ylim([0, 250])
            if save:
                filepath = save_dir or ('results/force_plots/force_plots_' + self.dateStr)
                filename = str(Path(filepath) / f'force_plot_{index}.png')
                # Save the plotted figure
                logging.debug("Calling _save_plot with filename: %s", filename)
                self._save_plot(plt, filename)

                # Plot & save transmission cable force vs. scaled cable displacement
                plt.figure()
                plt.plot(100 * x_cable_scaled, f_cable_scaled)
                plt.xlabel('Transmission cable displacement (cm)')
                plt.ylabel('Transmission cable force (N)')
                plt.title('Transmission cable force vs. displacement; scaled to knee angle')
                filename = str(Path(filepath) / f'force_plot_unscaled_{index}.png')
                logging.debug("Calling _save_plot with filename: %s", filename)
                self._save_plot(plt, filename)

                out_file = str(Path(filepath) / f'_force_output{index}.csv')
                logging.debug("Calling _save_array with output file: %s", out_file)
                self._save_array(out_file,
                                 np.stack((self.angles,
                                           percentages,
                                           x_cable,
                                           f_cable[:len(self.angles)],
                                           f_elastic[:len(self.angles)],
                                           x_cable_scaled,
                                           f_cable_scaled,
                                           f_elastic_scaled),
                                           axis=1),
                                 header='Angles (rad),' \
                                 'Stance Percentage (%),' \
                                 'Transmission Cable Displacement (m),' \
                                 'Transmission Cable Force (N),' \
                                 'Storage Cable Force (N),' \
                                 'Scaled Transmission Cable Displacement (cm),' \
                                 'Scaled Transmission Cable Force (N),' \
                                 'Scaled Storage Cable Force (N)')
                logging.debug("Sucessfully saved force output to: %s", out_file)
            if torque:
                plt.figure()
                plt.scatter(self.angles * 360 / (2*np.pi),
                            f_elastic * cam_radii[:, 0])
                plt.xlabel('Angle (deg)')
                plt.ylabel('Torque (N-m)')
            plt.show()

            # fields = ['Angles(rad)', 'Transmission cable displacement (m)', \
            #     'Transmission cable force (N)', 'Storage cable displacment (m)', \
            #     'Storage cable force (N)',]
            # rows = np.stack((self.angles, x_cable, f_cable, x_elastic, f_elastic), axis=1)
            
        return f_cable_scaled, percentages
    
    def knee_angle_to_stance(self, percentages_in):
        '''
        Create a function that relates knee angle to stance percentage using
        angle_data, accounting for nonlinear change in knee angle while 
        standing. Percentages go from ~0% to ~100% going from sit to stand.
        If no input percentages are provided, produces knee angles using evenly
        spaced percentages with same number of indices as sit data.
        '''
        angle_data = 'data/Knee-angle_Chugo_2006.csv'
        knee_array = np.loadtxt(angle_data, delimiter = ',', ndmin = 2)
        percentages = knee_array[:, 0]
        knee_angles = knee_array[:, 1]
        knee_fcn = interpolate.interp1d(percentages, knee_angles,
                                        kind='cubic', fill_value='extrapolate')
        knee_angles = knee_fcn(percentages_in)
        return knee_angles
    
    def percentage_to_x(self, percentages_in, stroke):
        # Generate knee angles corresponding to stance percentage.
        knee_angles = self.knee_angle_to_stance(percentages_in)

        # Map transmission cable stroke to stance percentage.
        x_cable_original = (stroke * (1 - (knee_angles - np.min(knee_angles))
                                 / np.ptp(knee_angles))) 
        return x_cable_original
    
    def generate_sit_cam(self, radii_stand, angles_stand, radii_outer, k_elastic,
                         n_params=6):

        def interp_cam(x):
            # # Interpolate between radii keypoints
            # key_angles = x[:n_params]
            # key_radii = x[n_params:]
            # spline_fcn = interpolate.interp1d(key_angles, key_radii,
            #                                   kind='cubic', fill_value='extrapolate')
            # radii_interp = spline_fcn(angles_stand[: self.sit_ind+1])
            # return radii_interp
            # Convert radii to points in Cartesian space.

            # Take radii keypoints and angles and return the convex hull of the cam radii.
            key_angles = x[:n_params]
            key_radii = x[n_params:]
            points = (key_radii * [np.cos(key_angles),
                                      np.sin(key_angles)
                                      ]).T
            # Take convex hull of cam radii.
            radii_convex = self.convex_cam_pts(points, angles_stand[: self.sit_ind+1], False)
        
            return radii_convex
        
        # Bounds:
        # Constrain the radii at sit and stand angles to be equal to those of
        # the sit-to-stand inner cam, within a threshold. Set lower and upper
        # bounds at all other angles to reasonable values.
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

        # Nonlinear contraints:
        # 1. Constrain the path lengths of the stand-to-sit and sit-to-stand
        # cams to be equal within a threshold.
        x_cable_stand = np.cumsum(radii_stand * 2*np.pi / self.n_interp)
        path_length_stand = x_cable_stand[self.sit_ind]
        thresh_path = 0.001
        ub_path = path_length_stand + thresh_path
        lb_path = path_length_stand - thresh_path
        def constr_path(x):
            radii_interp = interp_cam(x)
            x_cable_sit = np.cumsum(radii_interp * 2*np.pi / self.n_interp)
            path_length_sit = x_cable_sit[-1]
            return path_length_sit
        path_constraint = NonlinearConstraint(constr_path, lb_path, ub_path,
                                              keep_feasible=False)
        
        # 2. Constrain the radii at the end points to be equal to the
        # sit-to-stand cam radii.
        thresh_end = 0.0001
        ub_end = thresh_end * np.ones(2)
        lb_end = -1 * ub_end
        def constr_ends(x):
            radii_interp = interp_cam(x)
            return np.array([radii_interp[0] - radii_stand[0],
                    radii_interp[-1] - radii_stand[self.sit_ind]])
        end_constraint = NonlinearConstraint(constr_ends, lb_end, ub_end,
                                             keep_feasible=False)
        
        cons = [path_constraint, end_constraint]

        # Objective:
        # Minimize maximum force generated by stand-to-sit cam in combination
        # with outer cam.
        def objective(x):
            radii_interp = interp_cam(x)

            # # Convert radii to points in Cartesian space.
            # angles = angles_stand[: self.sit_ind+1]
            # points = (radii_interp * [np.cos(angles),
            #                           np.sin(angles)
            #                           ]).T
            # # Take convex hull of cam radii.
            # radii_convex = self.convex_cam_pts(points, angles, False)

            # Prepare cam radii and cable path for force calculation.
            cam_radii_sit = np.vstack((radii_interp, radii_outer[:self.sit_ind+1])).T
            x_cable_sit = np.cumsum(radii_interp * 2*np.pi / self.n_interp)

            forces, percentages = self.calc_forces(cam_radii_sit, x_cable_sit,
                                                k_elastic, torque=False)

            return np.max(forces)

        # Initial guess: a triangular profile (in terms of radii) between the
        # end points with the correct path length
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
        logger.info("Initial guess path length: %s", path_length_init)
        logger.info("Initial guess objective value: %s", objective(x0))

        plt = _get_plt()
        plt.figure()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(ang0, rad0, label='Initial guess')
        plt.show()

        # Callback: intermittently check the optimization progress.
        n_eval = 1
        def callback_fun(x):
            nonlocal n_eval
            if n_eval % 10 == 0:
                logger.info("Optimization iteration %s. Objective = %s", n_eval, objective(x))
                logger.info("Standing path length: %s", path_length_stand)
                radii_interp = interp_cam(x)
                logger.info("Current path length: %s", np.cumsum(radii_interp * 2*np.pi / self.n_interp)[-1])
                plt.figure()
                fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
                ax.plot(angles_stand[: self.sit_ind+1], radii_interp, label=f"Optimization iteration {n_eval}")
                plt.title(f"Optimization iteration {n_eval}")
                plt.show()
            n_eval += 1

        # Minimize: perform the optimization.
        # Import minimize lazily so the module can be imported without optional deps
        try:
            from cobyqa import minimize
        except Exception:
            logger.error("cobyqa.minimize is required for generate_sit_cam optimization but is not available")
            raise
        result = minimize(objective, x0,
                          options={'maxiter': 10000},
                          constraints=cons)
                        #   bounds=bounds)
                        #   callback=callback_fun)
        
        # Add the "non-functional" cam radii from the sit-to-stand cam (from the
        # sit angle and above) to the stand-to-sit cam radii. Report and plot 
        # the results.
        radii_interp = interp_cam(result.x)
        radii_interp = np.concatenate((radii_interp, radii_stand[self.sit_ind+1 :]))
        x_cable_sit = np.cumsum(radii_interp * 2*np.pi / self.n_interp)
        logger.info("Sit cam path length (360 deg): %s", x_cable_sit[-1])
        logger.info("Sit cam path length (220 deg): %s", x_cable_sit[self.sit_ind])
        logger.info("Stand cam path length: %s", path_length_stand)
        logger.info("Sit cam start radius: %s", radii_interp[0])
        logger.info("Sit cam end radius: %s", radii_interp[-1])

        plt.figure()
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(angles_stand, 100*radii_interp, label='Sitting cam')
        ax.set_xticklabels([])
        ax.set_title("Stand-to-Sit Cam")
        plt.show()

        # Convert radii to Cartesian points for saving.
        result_points = (radii_interp * [np.cos(self.angles), np.sin(self.angles)]).T

        return result, radii_interp, x_cable_sit, result_points
    
    def remove_duplicates(self, x: np.ndarray, y: Optional[np.ndarray] = None, z: Optional[np.ndarray] = None, *, strategy: str = 'drop_all') -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Remove duplicates from x and corresponding entries in y and z.

        Parameters:
        - x: numpy array of keys
        - y, z: optional arrays of the same length as x
        - strategy: how to handle repeated values in x:
            * 'drop_all' (default): remove all occurrences of any value that
              appears more than once (preserves only values that are unique in x)
            * 'keep_first': keep the first occurrence of each duplicate value
            * 'keep_last': keep the last occurrence of each duplicate value

        Returns:
        Tuple of (x_new, y_new, z_new) where y_new or z_new may be None if
        corresponding inputs were None.
        """
        if strategy not in {'drop_all', 'keep_first', 'keep_last'}:
            raise ValueError("strategy must be one of 'drop_all', 'keep_first', 'keep_last'")

        # Build mapping from value to list of indices
        idx_map = defaultdict(list)
        for i, v in enumerate(x):
            idx_map[v].append(i)

        keep_indices = []
        if strategy == 'drop_all':
            # Keep indices whose value occurs exactly once
            for v, inds in idx_map.items():
                if len(inds) == 1:
                    keep_indices.append(inds[0])
        elif strategy == 'keep_first':
            for v, inds in idx_map.items():
                keep_indices.append(inds[0])
        else:  # keep_last
            for v, inds in idx_map.items():
                keep_indices.append(inds[-1])

        keep_indices = np.array(sorted(keep_indices), dtype=int)
        x_new = x[keep_indices]
        y_new = y[keep_indices] if y is not None else None
        z_new = z[keep_indices] if z is not None else None

        return x_new, y_new, z_new
    
    def to_polar(self, points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Convert Cartesian points to polar coordinates in [0, 2*pi).

        Parameters:
        - points: array of shape (N, 2)

        Returns:
        - radii: shape (N,)
        - angles: shape (N,) in radians in range [0, 2*pi)
        """
        radii = np.linalg.norm(points, axis=1)
        angles = np.mod(np.arctan2(points[:, 1], points[:, 0]), 2 * np.pi)
        return radii, angles

    # --- File I/O helpers -------------------------------------------------
    def _ensure_dir(self, dirpath: str | Path) -> None:
        """Ensure a directory exists (create if needed). Accepts a string or Path."""
        if not dirpath:
            return
        p = Path(dirpath)
        if not p.exists():
            p.mkdir(parents=True, exist_ok=True)

    def _save_plot(self, plt_module, filepath: str | Path, dpi: int = 300, bbox_inches: str | None = 'tight') -> None:
        """Save a plot using the provided matplotlib.pyplot module.

        Ensures directory exists before calling plt.savefig.
        """
        logging.debug("In _save_plot!")
        p = Path(filepath)
        self._ensure_dir(p.parent)
        # matplotlib accepts Path-like objects in recent versions but convert to str for compatibility
        outpath = str(p)
        if bbox_inches is not None:
            plt_module.savefig(outpath, dpi=dpi, bbox_inches=bbox_inches)
            logging.debug(f"Saved plot to {outpath} with bbox_inches={bbox_inches}")
        else:
            plt_module.savefig(outpath, dpi=dpi)
            logging.debug(f"Saved plot to {outpath} with tight bounding box")

    def _save_array(self, filepath: str | Path, arr: np.ndarray, header: str | None = None, delimiter: str = ',') -> None:
        """Save an array to a text file, ensuring the containing directory exists."""
        p = Path(filepath)
        self._ensure_dir(p.parent)
        np.savetxt(str(p), arr, delimiter=delimiter, header=header)

    def plot_cams(self, cam_radii=0, k_elastic=0, index=0, save: bool = False, save_dir: Optional[str] = None):
        """
        Plots the cam points in polar coordinates.
        """
        plt = _get_plt()
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

        if save:
            logging.debug("Saving cam plot for index %s", index)
            target_dir = save_dir or ('results/cams/cams_' + self.dateStr + '/plots')
            filename = str(Path(target_dir) / f'cam_plot_{index}.png')
            self._save_plot(plt, filename, bbox_inches='tight')
        plt.show()

    def plot_cams_cartesian(self, pts_inner, pts_outer):
        """
        Plots the cam points in Cartesian coordinates.
        """
        plt = _get_plt()
        plt.figure()
        plt.plot(pts_inner[:, 0], pts_inner[:, 1], lw = 2)
        plt.plot(pts_outer[:, 0], pts_outer[:, 1], lw = 2)
        plt.axis('equal')
        plt.legend(['inner cam','outer cam'])
        plt.title('Cam Points in Cartesian Space')
        plt.show()