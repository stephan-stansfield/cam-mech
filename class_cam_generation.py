"""
cam-mech: cam generation utilities

Provides CamGeneration and helper functions for reading CSV data,
computing curve lengths, and converting between Cartesian and polar
coordinates. This module is intended to be used as a library and
avoids doing heavy computation on import.

The public API includes:
- read_xy_from_csv
- curve_length_from_csv
- knee_angle_to_stance
- percentage_to_x
- CamGeneration

"""
from __future__ import annotations

import math
from datetime import date
from pathlib import Path
from typing import Optional, Tuple
import logging

import numpy as np
import pandas as pd

from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.spatial import ConvexHull
from scipy.optimize import NonlinearConstraint
from scipy.integrate import trapezoid
from collections import defaultdict

logger = logging.getLogger(__name__) # Initiate module logger

font = {'size': 14}
try:
    import matplotlib
    matplotlib.rc('font', **font)
except Exception:
    logger.debug("matplotlib not available; skipping font rc setup")

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

    The function reads x/y points from the provided CSV file, constructs
    a parametric spline curve that passes through the points in order,
    and returns the cumulative arc length at each original point.

    Parameters:
    csv_filepath: str
        Path to the CSV file containing ordered x and y coordinates.
    n_eval: int, optional
        Number of points to evaluate along the spline for length
        estimation.
    spline_order: int, optional
        Order of the spline. Default is 3 (cubic) and is automatically
        reduced when there are too few data points.
    x_col: int, optional
        Index of the column containing x coordinates.
    y_col: int, optional
        Index of the column containing y coordinates.
    delimiter: str, optional
        Field delimiter used by the CSV file.

    Returns:
    numpy.ndarray
        Array of cumulative arc lengths at each data point, in the same
        order as the input points.
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


def knee_angle_to_stance(percentages_in):
    """Map stance percentage to knee angle using reference gait data.

    Interpolates the knee angle vs. stance percentage relationship from
    Chugo et al., 2006, accounting for the nonlinear change in knee
    angle while standing. Percentages go from ~0% to ~100% going from
    sit to stand.

    Parameters:
    percentages_in: stance percentages (0 to 100) to convert

    Returns:
    numpy.ndarray
        Knee angle, in degrees, at each input percentage.
    """
    angle_data = 'data/Knee-angle_Chugo_2006.csv'
    knee_array = np.loadtxt(angle_data, delimiter=',', ndmin=2)
    percentages = knee_array[:, 0]
    knee_angles = knee_array[:, 1]
    knee_fcn = interpolate.interp1d(percentages, knee_angles,
                                    kind='cubic', fill_value='extrapolate')
    return knee_fcn(percentages_in)


def percentage_to_x(percentages_in, stroke):
    """Map stance percentage to transmission cable displacement.

    Parameters:
    percentages_in: stance percentages (0 to 100) to convert
    stroke: total transmission cable stroke length corresponding to
        the full knee angle range

    Returns:
    numpy.ndarray
        Transmission cable displacement at each input percentage.
    """
    # Generate knee angles corresponding to stance percentage.
    knee_angles = knee_angle_to_stance(percentages_in)

    # Map transmission cable stroke to stance percentage.
    return (stroke * (1 - (knee_angles - np.min(knee_angles))
                     / np.ptp(knee_angles)))

class CamGeneration:
    """Class that generates points on the cams for every degree given
    torque/gear ratios. Can also return force and energy output given 
    stiffness.
    """
    
    def __init__(self, gear_ratios, input_angles, scaling, sit_angle, offset_angle):
        """
        Parameters:
        gear_ratios: desired gear ratio at each keypoint along the
            compound cam
        input_angles: angles at which to define the desired gear ratios
        scaling: scalar multiple to determine overall size of cams
        sit_angle: angle in radians along the compound cam at which
            force cable is acting while the user is seated
        offset_angle: angle in radians along the compound cam at which 
            the energy storage band acts relative to the force
            transmission cable
        """
        # Sort gear ratios and input angles in ascending order of input
        # angles.
        self.input_angles = input_angles
        self.gear_ratios = gear_ratios
        sorted_inds = np.argsort(self.input_angles)
        self.input_angles = self.input_angles[sorted_inds]
        self.gear_ratios = gear_ratios[sorted_inds]

        # Initialize class variables.
        self.scaling = scaling
        self.cam_radii = np.ones((self.gear_ratios.shape[0], 2))
        self.n_interp = 360 
        self.sit_ind = round(np.rad2deg(sit_angle))
        self.offset_ind = round(np.rad2deg(offset_angle))
        self.angles = np.deg2rad(np.arange(0, self.n_interp+1, 1))
        self.pts_trans = None
        self.pts_stor = None
        self.dateStr = date.today().strftime("%Y-%m-%d")

    def calculate_cam_radii(self, user_height: float = 1.67, k: float = 0,
                            plot: bool = False, save: bool = False, 
                            save_dir: Optional[str] = None, index: int = 0):
        """Calculates the cam radii for each gear ratio and input angle.
        Determines the convex hull of the given gear ratios and angles,
        then interpolates between these points to generate the cam
        radii.

        Parameters:
        user_height: height, in meters, of the person the device is
            designed for
        k: spring stiffness, in N/m, used only to label saved plots
        plot: if True, display intermediate and final cam plots
        save: if True, save the final cam plot to disk (requires
            plot=True)
        save_dir: directory to save the plot to; defaults to a dated
            subfolder of results/cams if not provided
        index: identifier appended to the saved plot's filename

        Returns:
        pts_trans: ndarray, transmission cam points in Cartesian space
        pts_stor: ndarray, storage cam points in Cartesian space
        radius_max: scalar maximum radius of the cam envelopes
        """
        
        # Calculate initial cam keypoints using gear ratios and scaling
        # factor. Check if any of the radii are smaller than the
        # minimum allowed radius, in meters, and adjust as needed.
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
            ax.set_title('Optimization input parameters:\nscaled gear ratios and angles')
            plt.show()
            
        # Convert cam points to Cartesian space.
        self.pts_trans = (self.cam_radii[:, 0] * [np.cos(self.input_angles),
                                                  np.sin(self.input_angles)
                                                  ]).T
        self.pts_stor = (self.cam_radii[:, 1] * [np.cos(self.input_angles),
                                                  np.sin(self.input_angles)
                                                  ]).T

        # Calculate points of convex hull for each cam. Ouputs cam radii
        # in polar coordinates.
        cam_radii_trans = self.convex_cam_pts(self.pts_trans, self.angles)
        cam_radii_stor = self.convex_cam_pts(self.pts_stor, self.angles)
        self.cam_radii = np.vstack((cam_radii_trans, cam_radii_stor)).T

        # Determine stroke length achieved by flexing knee from upright
        # to 90 degrees. This calculation is based on measurements
        # taken in Song et al., 2022
        # (https://doi.org/10.1177/15589250221138546), linearly scaled
        # to the user's height.
        unstretch_len = 0.05
        stretch_pct = 1.016
        user_height_ref = 1.67
        stroke = user_height / user_height_ref * unstretch_len * stretch_pct

        # Scale cam size to calculated stroke length. Iteratively check
        # that the minimum radius, in meters, is not violated.
        ratio = 0
        threshold = 0.01
        r_min = 0.25 * .0254 # convert from inches to meters
        while abs(ratio - 1) > threshold:
            x_trans = np.cumsum(np.deg2rad(self.cam_radii[:, 0]))
            ratio = stroke / x_trans[self.sit_ind]
            self.cam_radii *= ratio
            x_trans *= ratio
            for cam in self.cam_radii:
                for point in cam:
                    if point < r_min:
                        point = r_min
        radius_max = np.max(self.cam_radii)

        # Rotate values in storage cam to account for where cable leaves
        # surface relative to where transmission cable leaves surface.
        cam_rotated, _ = self.rotate_cam(self.cam_radii[:, 1])
        self.cam_radii[:, 1] = cam_rotated

        if plot:
            # Plot; saving is opt-in via save flag and save_dir
            self.plot_cams(self.cam_radii, k, index, save=save)

        # Convert cam points to Cartesian space and return the final
        # cam shapes.
        self.pts_trans = (self.cam_radii[:, 0] * [np.cos(self.angles),
                                                  np.sin(self.angles)
                                                  ]).T
        self.pts_stor = (self.cam_radii[:, 1] * [np.cos(self.angles),
                                                  np.sin(self.angles)
                                                  ]).T

        return self.cam_radii, self.pts_trans, self.pts_stor, radius_max, x_trans

    def convex_cam_pts(self, points, angles):
        """
        Compute the convex cam points given a set of input points.

        Parameters:
        - points: numpy array
            The input points used to compute the convex cam points.
        - angles: numpy array
            The angles used to generate the radii of the convex cam.

        Returns:
        - numpy array
            The computed convex cam points in polar coordinates.

        This function computes the convex cam points by performing the
        following steps:
        1. Compute the convex hull of the input points.
        2. Convert the convex points into polar coordinates and
           interpolate between them.
        3. Convert the interpolated points back into Cartesian
           coordinates.
        4. Split the cam points into sub-arrays that monotonically
           increase or decrease in x.
        5. Flip any segments that are monotonically decreasing.
        6. Interpolate between the points in each sub-array.
        7. Convert the interpolated points into polar coordinates.
        8. Linearly interpolate between the points in polar coordinates.
        9. Return the computed convex cam points in polar coordinates.
        """
        # Compute the convex hull of the input points.
        hull = ConvexHull(points, incremental=True)
        points = np.vstack((points[hull.vertices, 0],
                            points[hull.vertices, 1])).T

        # Convert convex cam points to polar coordinates in range [0,
        # 2*pi].
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

        # Convert to Cartesian coordinates and calculate convex hull
        # again.
        points = (radii_interp * [np.cos(angles), np.sin(angles)]).T
        hull = ConvexHull(points, incremental=True)
        points = np.vstack((points[hull.vertices, 0],
                            points[hull.vertices, 1])).T

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

        # Flip any segments that are monotonically decreasing, which is
        # required by numpy.interp(). Then, linearly interpolate between
        # 100 points over the range of each sub-array and between the 
        # sub-arrays. This fills any gaps that resulted from the convex
        # hull operation.
        max_gap = 0.1
        interp_array = np.empty((0, 2))
        for ind, sub_array in enumerate(points_array):
            # Identify any large gaps between sub-arrays and inerpolate
            # over them.
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
                    interp_array = np.concatenate((interp_array, 
                                                   np.array([x, y]).T))
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
                    interp_array = np.concatenate((interp_array, 
                                                   np.array([x, y]).T))
        for ind, sub_array in enumerate(points_array):
            # Flip segments that are monotonically decreasing.
            if not np.all(np.diff(sub_array[:, 0]) > 0):
                points_array[ind] = np.flip(points_array[ind], 0)

        for ind, sub_array in enumerate(points_array):
            x = np.linspace(sub_array[0, 0], sub_array[-1, 0], 100)
            y = np.interp(x, sub_array[:, 0], sub_array[:, 1])
            interp_array = np.concatenate((interp_array, 
                                           np.array([x, y]).T))

        # Convert interpolated convex cam points into polar coordinates
        # in range [0, 2*pi].
        radii_interp, angles_interp = self.to_polar(interp_array)

        # Remove duplicates to avoid errors in interpolation.
        dup_removed = self.remove_duplicates(angles_interp, radii_interp)
        angles_interp = dup_removed[0]
        radii_interp = dup_removed[1]

        # Pad arrays to ensure continuity at end points. Cubically
        # interpolate between points in polar coordinates and return.
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
        """Rotate an array of cam radii by offset_ind indices (leftward
        cyclic shift) to account for where the storage cable leaves the
        cam surface relative to where the transmission cable leaves it.

        Parameters:
        radii: cam radii indexed by angle, before rotation

        Returns:
        cam_rotated: radii cyclically shifted by offset_ind
        cam_original: unmodified copy of the input radii
        """
        cam_original = radii.copy()
        radii_1 = radii[: self.offset_ind]
        radii_2 = radii[self.offset_ind :]
        cam_rotated = np.concatenate((radii_2, radii_1))
        return cam_rotated, cam_original

    def derotate_cam(self, radii):
        """Undo the cyclic shift applied by rotate_cam, restoring radii
        to their original angular indexing.
        """
        radii_1 = radii[-self.offset_ind :]
        radii_2 = radii[: -self.offset_ind]
        cam_derotated = np.concatenate((radii_1, radii_2))
        return cam_derotated
    
    def calc_forces(self, cam_radii, x_trans, k, torque: bool = False, 
                    plot: bool = False, save: bool = False,
                    save_dir: Optional[str] = None, index: int = 0):
        """
        Calculate the cable force profiles vs. stance percentage given
        the solved cam radii.

        The order of causality is:
        stance percentages -> 
        knee angles ->
        scaled transmission cable displacement ->
        cam angles -> 
        storage cable displacement (with spring stiffness) ->
        storage cable force (with gear ratios) -> 
        transmission cable force
        """

        # Define transmission and storage cam radii corresponding to 
        # each angle of rotation.
        cam_radii_trans = cam_radii[:, 0]
        cam_radii_stor = self.derotate_cam(cam_radii[:, 1])

        # Calculate effective cam moment arms and gear ratio at each angle.
        cam_moment_arms_trans = self.get_moment_arms(cam_radii_trans,
                                                     sit_ind=self.sit_ind)
        cam_moment_arms_stor = self.get_moment_arms(cam_radii_stor,
                                                    sit_ind=self.sit_ind)

        # get_moment_arms only fills indices [0, sit_ind); beyond that 
        # both arrays are left as 0 (non-functional region), so divide 
        # only where the denominator is nonzero to avoid a 0/0
        # RuntimeWarning and leave those unused entries as 0.
        gear_ratios = np.divide(cam_moment_arms_stor, cam_moment_arms_trans,
                                out=np.zeros_like(cam_moment_arms_stor),
                                where=cam_moment_arms_trans != 0)

        # Calculate storage cable displacement using the assumption of 
        # circular shape between subsequent points.
        # TODO: model this more accurately.
        x_stor = np.cumsum(np.deg2rad(cam_radii_stor))

        # Calculate tension in both cables.
        f_stor = k * x_stor
        f_trans = f_stor * gear_ratios

        # Stored energy vs. storage cable displacement
        E = trapezoid(f_stor, x_stor)
        logger.info("Total energy stored in storage cable: %s J", E)

        if plot:
            plt = _get_plt()
            plt.plot(100 * x_trans[self.sit_ind:0:-1], f_trans[:self.sit_ind],
                                 linewidth=3)
            plt.xlabel('Transmission Cable Displacement (cm)')
            plt.ylabel('Force (N)')
            plt.title('Transmission Cable Force vs. Displacement')
            plt.show()

            """
            # This plot may be useful for debugging, but is not 
            # particularly interesting otherwise, as it just shows a 
            # linear relationship between force and displacement.
            plt.figure()
            plt.plot(100 * x_stor[:self.sit_ind], f_stor[:self.sit_ind],
                     linewidth=3)
            plt.xlabel('Storage Cable Displacement (cm)')
            plt.ylabel('Force (N)')
            plt.title('Storage Cable Force vs. Displacement')
            plt.show()
            """

        # Generate knee angles corresponding to stance percentage.
        percentages = np.linspace(0, 100, self.sit_ind)
        knee_angles = self.knee_angle_to_stance(percentages)

        # Flip knee angles and percentages to correspond with cam angle
        # index, which increases from stand to sit.
        percentages = percentages[::-1]
        knee_angles = knee_angles[::-1]

        # Scale transmission cable displacement by the normalized knee
        # angle. At this point, x_trans_scaled should decrease with
        # increasing index.
        x_trans_scaled = (x_trans[self.sit_ind]
                          * (knee_angles - np.min(knee_angles))
                          / np.ptp(knee_angles)) 

        # Isolate "effective" portions of x_trans (from 0 to 220
        # degrees) and flip order to match the order of knee angles,
        # which go from sit to stand.
        x_trans = x_trans[:self.sit_ind]
        x_trans = x_trans[::-1]

        if plot:
            # Plot knee angle vs. cable displacement
            """
            # This plot may be useful for debugging, but is not 
            # particularly interesting otherwise, as it just shows a 
            # linear relationship between knee angle and displacement.
            plt = _get_plt()
            plt.plot(knee_angles, 100 * x_trans_scaled)
            plt.xlabel('Knee Angle (degrees)')
            plt.ylabel('Transmission Cable Displacement (cm)')
            plt.title('Transmission Cable Displacement vs. Knee Angle')
            """

            # Plot stance percentage vs. cable displacement
            plt.figure()
            plt.plot(percentages, 100 * x_trans_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Transmission Cable Displacement (cm)')
            plt.title('Transmission Cable Displacement vs. Stance Percentage (Scaled to Knee Angle)')

        if np.any(np.isnan(x_trans)):
            logger.warning("NaN detected in x_trans passed to calc_forces")
        
        # Remove duplicate points from cable path information (and
        # corresponding points from the angles and storage cable path)
        # to avoid errors in interpolation.
        dup_removed = self.remove_duplicates(x_trans,
                                             self.angles[:x_trans.size],
                                             x_stor)
        x_trans = dup_removed[0]
        cam_angles = dup_removed[1]
        x_stor = dup_removed[2]

        # Create functions relating transmission cable displacement to
        # cam angle and relating cam angle to storage cable
        # displacement.
        trans_disp_to_angle = interpolate.interp1d(x_trans, cam_angles,
                                                   kind='cubic',
                                                   fill_value='extrapolate')
        angle_to_stor_disp = interpolate.interp1d(cam_angles, x_stor,
                                                  kind='cubic',
                                                  fill_value='extrapolate')
        
        # Use scaled cable displacement to find cam angle as driven by
        # knee rotation.
        angle_scaled = trans_disp_to_angle(x_trans_scaled)

        # Use knee-driven cam angle to find scaled storage cable displacement.
        x_stor_scaled = angle_to_stor_disp(angle_scaled)

        # Find gear ratio at each cam angle corresponding to scaled
        # cable displacement.
        angles_to_gear_ratios = interpolate.interp1d(cam_angles,
                                                     gear_ratios[:len(cam_angles)],
                                                     kind='cubic',
                                                     fill_value='extrapolate')
        gear_ratios_scaled = angles_to_gear_ratios(angle_scaled)

        # Find transmission cable force based on the storage cable
        # displacement, spring stiffness, and gear ratio at each point.
        f_stor_scaled = k * x_stor_scaled
        f_trans_scaled = f_stor_scaled * gear_ratios_scaled
        
        # Plot and save values.
        if plot:
            # Plot cam angle vs. stance percentage.
            plt = _get_plt()
            plt.figure()
            plt.plot(percentages, angle_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cam Angle (rad)')
            plt.title('Cam Angle vs. Stance Percentage (Scaled to Knee Angle)')

            # Plot storage cable displacement vs. stance percentage.
            plt.figure()
            plt.plot(percentages, x_stor_scaled)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Storage Cable Displacement (m)')
            plt.title('Storage Cable Displacement vs. Stance Percentage (Scaled to Knee Angle)')

            # Plot (& save) transmission cable force vs. stance percentage.
            plt.figure()
            plt.plot(percentages, f_trans_scaled, linewidth=3)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Force (N)')
            plt.title('Transmission Cable Tension vs. Stance Percentage (Scaled to Knee Angle)')
            plt.xlim([0, 100])
            if save:
                filepath = save_dir or ('results/forces/forces_' + self.dateStr)
                filename = str(Path(filepath) / f'force_plot_{index}.png')
                self._save_plot(plt, filename)

            # Plot (& save) transmission cable force vs. scaled cable
            # displacement.
            plt.figure()
            plt.plot(100 * x_trans_scaled, f_trans_scaled)
            plt.xlabel('Cable displacement (cm)')
            plt.ylabel('Force (N)')
            plt.title('Transmission Cable Tension vs. Displacement (Scaled to Knee Angle)')
            if save:
                filename = str(Path(filepath) / f'force_plot_unscaled_{index}.png')
                self._save_plot(plt, filename)

                # Save data to CSV for further analysis.
                out_file = str(Path(filepath) / f'_force_output{index}.csv')
                self._save_array(out_file,
                                 np.stack((cam_angles,
                                           percentages,
                                           x_trans,
                                           f_trans[:len(cam_angles)],
                                           f_stor[:len(cam_angles)],
                                           x_trans_scaled,
                                           f_trans_scaled,
                                           f_stor_scaled),
                                           axis=1),
                                 header='Angles (rad),' \
                                 'Stance Percentage (%),' \
                                 'Transmission Cable Displacement (m),' \
                                 'Transmission Cable Force (N),' \
                                 'Storage Cable Force (N),' \
                                 'Scaled Transmission Cable Displacement (cm),' \
                                 'Scaled Transmission Cable Force (N),' \
                                 'Scaled Storage Cable Force (N)')
                logging.debug("Successfully saved force output to: %s", out_file)
            if torque:
                plt.figure()
                plt.scatter(np.rad2deg(cam_angles),
                            f_stor * cam_moment_arms_trans[:len(cam_angles)],
                            linewidth=3)
                plt.xlabel('Angle (deg)')
                plt.ylabel('Torque (N-m)')
                plt.title('Torque vs. Angle')
            plt.show()

        return f_trans_scaled, percentages

    def get_moment_arms(self, r, sit_ind):
        """
        Estimate the effective cam moment arm at each degree of
        rotation, up to sit_ind.

        A cable wrapped around a non-circular cam does not necessarily
        leave the surface tangentially at the nominal point; the true
        moment arm is the largest projected radius among nearby points
        on the profile. For each degree `cam_deg` on the cam, this
        searches the points within ±90 degrees of it and takes the one
        whose radius, projected onto the normal of the cable direction
        at cam_deg (simplified to be constant), is largest.

        Parameters:
        r: cam radii at each integer degree of rotation (indexed by
            degree)
        sit_ind: index (degree) up to which moment arms are computed;
            points beyond this are not part of the cam's effective
            range

        Returns:
        numpy.ndarray
            Effective moment arm at each degree from 0 to sit_ind - 1.
        """
        x_candidates = np.zeros(181)
        x_effective = np.zeros_like(r)

        n_radii = len(r[:sit_ind])
        for cam_deg in range(n_radii):
            for offset_ind in range(181):  # candidate points within +/-90 degrees of cam_deg
                if cam_deg + offset_ind - 90 < 0:
                    # clamp to the first point in range (0 degrees)
                    offset_deg = -cam_deg
                elif cam_deg + offset_ind - 90 >= 360:
                    # clamp to the last point in range (360 degrees)
                    offset_deg = 360 - cam_deg
                else:
                    # offset_ind - 90 gives a range of -90 to +90
                    # degrees relative to cam_deg
                    offset_deg = offset_ind - 90

                x_candidates[offset_ind] = (r[cam_deg+offset_deg] 
                                            * np.cos(np.deg2rad(offset_deg)))

            x_effective[cam_deg] = max(x_candidates)

        return x_effective
    
    def knee_angle_to_stance(self, percentages_in):
        """See module-level knee_angle_to_stance."""
        return knee_angle_to_stance(percentages_in)

    def percentage_to_x(self, percentages_in, stroke):
        """See module-level percentage_to_x."""
        return percentage_to_x(percentages_in, stroke)
    
    def generate_sit_cam(self, r_si2st, ang_si2st, r_stor, k, n_params=6):
        """
        Optimize a stand-to-sit transmission cam profile that pairs with
        an existing pair of sit-to-stand transmission and storage cams,
        minimizing peak transmission cable force while matching path
        length and end-point radii between the two transmission cams.

        Parameters:
        r_si2st: radii of the previously generated sit-to-stand cam
                 (used as matching constraints)
        ang_si2st: angles corresponding to r_si2st
        r_stor: storage cam radii to pair with the new stand-to-sit
                transmission cam
        k: spring stiffness, in N/m, used in the force calculation
        n_params: number of angle/radius keypoints used to parameterize
                  the stand-to-sit cam profile before taking its convex
                  hull

        Returns:
        result: the optimization result object from cobyqa.minimize
        r_st2si: optimized stand-to-sit transmission cam radii at every
                 angle
        x_trans_st2si: cumulative transmission cable displacement for
                       r_st2si
        r_st2si_cart: r_st2si converted to Cartesian (x, y) coordinates
        """

        def params_to_convex_radii(x):
            # Take radii/angle keypoint parameters and return a
            # smoothed, convex cam profile.
            key_angles = x[:n_params]
            key_radii = x[n_params:]
            pts = (key_radii * [np.cos(key_angles), np.sin(key_angles)]).T
            radii_convex = self.convex_cam_pts(pts, ang_si2st[: self.sit_ind])
            return radii_convex

        # NONLINEAR CONSTRAINTS:
        # 1. Constrain the path lengths of the stand-to-sit and
        # sit-to-stand cams to be equal within a threshold.
        x_trans_si2st = np.cumsum(np.deg2rad(r_si2st))
        path_length_si2st = x_trans_si2st[self.sit_ind]
        thresh_path = 0.001 # meters
        ub_path = path_length_si2st + thresh_path
        lb_path = path_length_si2st - thresh_path

        def constr_path(x):
            # Get stand-to-sit transmission cam radii.
            r_st2si = params_to_convex_radii(x)

            # Calculate and return path length.
            x_trans_st2si = np.cumsum(np.deg2rad(r_st2si))
            return x_trans_st2si[-1]
        
        path_constraint = NonlinearConstraint(constr_path, lb_path, ub_path,
                                              keep_feasible=False)
        
        # 2. Constrain the radii at the end points to be equal to the
        # sit-to-stand cam radii within a threshold.
        thresh_end = 0.0001 # meters
        ub_end = thresh_end * np.ones(2)
        lb_end = -1 * ub_end

        def constr_ends(x):
            # Get stand-to-sit transmission cam radii.
            r_st2si = params_to_convex_radii(x)

            # Calculate and return the difference between the end-point
            # radii of the two transmission cams.
            return np.array([r_st2si[0] - r_si2st[0],
                             r_st2si[-1] - r_si2st[self.sit_ind]])
        
        end_constraint = NonlinearConstraint(constr_ends, lb_end, ub_end,
                                             keep_feasible=False)
        
        constraint_array = [path_constraint, end_constraint]

        # OBJECTIVE:
        # Minimize maximum force generated by stand-to-sit transmission
        # cam in combination with storage cam.
        def objective(x):
            # Get stand-to-sit transmission cam radii.
            r_st2si = params_to_convex_radii(x)

            # Add the "non-functional" cam radii from the sit-to-stand
            # cam (from the sit angle and above) to the stand-to-sit cam
            # radii.
            r_st2si = np.concatenate((r_st2si, r_si2st[self.sit_ind:]))

            # Prepare cam radii and cable path for force calculation.
            r_dual_st2si = np.vstack((r_st2si, r_stor)).T
            x_trans_st2si = np.cumsum(np.deg2rad(r_st2si))

            forces, percentages = self.calc_forces(r_dual_st2si, x_trans_st2si,
                                                   k, torque=False)

            # Return the objective to be minimized: the maximum force.
            return np.max(forces)

        # INITIAL GUESS:
        # A triangular profile (in terms of radii) between the end
        # points with the correct path length.
        r_max = (2*path_length_si2st/(np.deg2rad(self.sit_ind)) 
                 - 0.5*(r_si2st[0]+r_si2st[self.sit_ind]))
        r_1 = np.linspace(r_si2st[0], r_max, num=int(n_params/2))
        r_2 = np.linspace(r_max, r_si2st[self.sit_ind], num=int(n_params/2))
        r_0 = np.concatenate((r_1, r_2))
        ang_0 = np.linspace(ang_si2st[0], ang_si2st[self.sit_ind], num=n_params)
        x0 = np.concatenate((ang_0, r_0))
        path_length_init = np.cumsum(r_0 * 2*np.pi / n_params)[-1]
        logger.info("Initial guess path length: %s", path_length_init)
        logger.info("Initial guess objective value: %s", objective(x0))

        # Plot initial guess for visual inspection.
        plt = _get_plt()
        from matplotlib.ticker import MaxNLocator
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(ang_0, r_0 * 100)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.set_rlim(bottom=0)
        ax.set_ylabel('Radius (cm)', labelpad=30)
        plt.title("Initial Guess for Optimization")
        plt.show()

        # MINIMIZE:
        # Perform the optimization. Import minimize lazily so the
        # module can be imported without optional deps.
        try:
            from cobyqa import minimize
        except Exception:
            logger.error("cobyqa.minimize is required for generate_sit_cam " \
                         "optimization but is not available")
            raise
        result = minimize(objective, x0,
                          options={'maxiter': 10000},
                          constraints=constraint_array)

        # POST-PROCESSING
        # Add the "non-functional" cam radii from the sit-to-stand cam
        # (from the sit angle and above) to the stand-to-sit cam radii.
        # Log and plot the results.
        r_st2si = params_to_convex_radii(result.x)
        r_st2si = np.concatenate((r_st2si, r_si2st[self.sit_ind :]))
        x_trans_st2si = np.cumsum(np.deg2rad(r_st2si))
        logger.info("Stand-to-Sit cam path length (360 deg): %s", x_trans_st2si[-1])
        logger.info("Stand-to-Sit cam path length (220 deg): %s", x_trans_st2si[self.sit_ind])
        logger.info("Sit-to-Stand cam path length: %s", path_length_si2st)
        logger.info("Stand-to-Sit cam start radius: %s", r_st2si[0])
        logger.info("Stand-to-Sit cam end radius: %s", r_st2si[-1])

        # Plot the optimized stand-to-sit cam for visual inspection.
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(ang_si2st, r_st2si * 100)
        ax.set_ylabel('Radius (cm)', labelpad=30)
        ax.set_title("Stand-to-Sit Cam")
        plt.show()

        # Convert radii to Cartesian form for export to CAD programs.
        r_st2si_cart = (r_st2si * [np.cos(self.angles), np.sin(self.angles)]).T

        return result, r_st2si, x_trans_st2si, r_st2si_cart
    
    def remove_duplicates(self,
                          x: np.ndarray, 
                          y: Optional[np.ndarray] = None, 
                          z: Optional[np.ndarray] = None, 
                          *, 
                          strategy: str = 'drop_all') -> Tuple[np.ndarray,
                                                               Optional[np.ndarray], 
                                                               Optional[np.ndarray]]:
        """
        Remove duplicates from x and corresponding entries in y and z.

        Parameters:
        - x: numpy array of keys
        - y, z: optional arrays of the same length as x
        - strategy: how to handle repeated values in x:
            * 'drop_all' (default): remove all occurrences of any value 
              that appears more than once (preserves only values that 
              are unique in x)
            * 'keep_first': keep the first occurrence of each duplicate
               value
            * 'keep_last': keep the last occurrence of each duplicate 
              value

        Returns:
        Tuple of (x_new, y_new, z_new) where y_new or z_new may be None 
        if corresponding inputs were None.
        """
        if strategy not in {'drop_all', 'keep_first', 'keep_last'}:
            raise ValueError("strategy must be one of 'drop_all', " \
                             "'keep_first', 'keep_last'")

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

    # --- File I/O helpers -------------------------------------------------- #
    def _ensure_dir(self, dirpath: str | Path) -> None:
        """Ensure a directory exists (create if needed). Accepts a
        string or Path.
        """
        if not dirpath:
            return
        p = Path(dirpath)
        if not p.exists():
            p.mkdir(parents=True, exist_ok=True)

    def _save_plot(self, plt_module, filepath: str | Path, dpi: int = 300, 
                   bbox_inches: str | None = 'tight') -> None:
        """Save a plot using the provided matplotlib.pyplot module.

        Ensures directory exists before calling plt.savefig.
        """
        p = Path(filepath)
        self._ensure_dir(p.parent)
        # matplotlib accepts Path-like objects in recent versions but
        # convert to str for compatibility
        outpath = str(p)
        if bbox_inches is not None:
            plt_module.savefig(outpath, dpi=dpi, bbox_inches=bbox_inches)
        else:
            plt_module.savefig(outpath, dpi=dpi)

    def _save_array(self, filepath: str | Path, arr: np.ndarray, 
                    header: str | None = None, delimiter: str = ',') -> None:
        """Save an array to a text file, ensuring the containing
        directory exists.
        """
        p = Path(filepath)
        self._ensure_dir(p.parent)
        np.savetxt(str(p), arr, delimiter=delimiter, header=header)

    def plot_cams(self, cam_radii=0, k=0, index=0, save: bool = False, 
                  save_dir: Optional[str] = None):
        """
        Plots the cam points in polar coordinates.
        """
        plt = _get_plt()
        r = 100 * self.cam_radii[:, 0]
        R = 100 * self.cam_radii[:, 1]
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(self.angles, r, label='Transmission cam')
        ax.plot(self.angles, R, label='Storage cam')
        ax.legend(bbox_to_anchor=(1, 0), loc='lower left')
        ax.set_xticklabels([])
        ax.set_ylim([0, 7])
        ax.set_title(f"""Cam shapes\nmin radius={100*np.min(self.cam_radii):.2f} cm, """
                     f"""max radius={100*np.max(self.cam_radii):.2f} cm\nK={k} N/m""")

        if save:
            target_dir = save_dir or ('results/cams/cams_' + self.dateStr + '/cam_plots')
            filename = str(Path(target_dir) / f'cam_plot_{index}.png')
            self._save_plot(plt, filename, bbox_inches='tight')
        plt.show()