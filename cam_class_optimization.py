# Using python 3.10

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib
from scipy.interpolate import make_interp_spline, BSpline
import math

font = {
    'size'   : 14}

matplotlib.rc('font', **font)

class CamGeneration:
    '''
    Class that generates points on the cams for
    every degree given torque/gear ratios. 
    Can give energy output given stiffness as well.
    '''
    
    def __init__(self, gear_ratios, initial_gr, sit_angle=np.pi):
        '''
        points: number of points to define on cams
        gear_ratios: points*2 np array          [gear_ratios; angles]
        initial_gr: scalar, e.g. 0.5            overall gear ratio inner-to-outer
        Initial guess of inner and outer radii will be 1 and 1/initial_gr, e.g. 1 and 2
        '''
        self.gear_ratios = gear_ratios
        self.gear_ratios[:, 1] = self.gear_ratios[:, 1] * sit_angle / np.pi # normalize second row of gear ratios to ratio of sit angle to 180 deg
        self.initial_gr = initial_gr
        self.points = len(self.gear_ratios[:, 0]) # number of points to define on each cam

        self.cam_radii = np.ones((self.points, 2)) @ np.diag([1, 1/initial_gr])  # initial guess for inner and outer cam radii

        # assuming equal distribution of angles for gear ratios
        self.input_angles = self.gear_ratios[:, 1] # take the second row of gear ratios to be the input angles
        self.ninterp = 360 # number of radii points to interpolate between
        self.sit_ind = round(sit_angle / 2 / np.pi * self.ninterp) # index of sitting position (if ninterp=360, this will also equal the angle along the cam corresponding to sitting)

        self.angles = np.arange(0, self.ninterp + 1, 1) / self.ninterp * 2 * np.pi # create ndarray of radian angles in range [0, ninterp] 

        self.inner_pts = None
        self.outer_pts = None

    def calc_radii_convex_cams(self,inner_pts,outer_pts):
        print("in calc_radii_convex_cams")
        # new_radii = [[],[]]
        new_inner_radii = np.empty(inner_pts.shape[0])
        new_outer_radii = np.empty(outer_pts.shape[0])

        print("shapes of new_inner_radii and new_outer_radii:")
        print(np.shape(new_inner_radii))
        print(np.shape(new_outer_radii))

        ind = 0
        for point in inner_pts:
            dist = np.linalg.norm(point)
            new_inner_radii[ind] = dist
            ind += 1
            # DEBUG
            # print("point: ", point)
            # print("dist: ", dist)


        ind = 0
        for point in outer_pts:
            dist = np.linalg.norm(point)
            new_outer_radii[ind] = dist
            ind += 1

        print("shapes of new_inner_radii and new_outer_radii:")
        print(np.shape(new_inner_radii))
        print(np.shape(new_outer_radii))

        return new_inner_radii, new_outer_radii

    def get_convex_hull_angles(self,inner_pts=None,outer_pts=None):
        
        inner_angles = np.arctan2(inner_pts[:,1],inner_pts[:,0]) + np.pi # STS: why add pi?
        outer_angles = np.arctan2(outer_pts[:,1],outer_pts[:,0]) + np.pi

        # DEBUG
        # for point in inner_pts:
        #     dist = np.linalg.norm(point)
        #     angle = np.arctan2(point[1], point[0])
        #     print("Point: ", point)
        #     print("dist: ", dist)
        #     print("angle: ", angle)
        #     print("----")
    
        return np.array([inner_angles]).flatten(), np.array([outer_angles]).flatten()
        
    def calculate_cam_radii(self, kind='cubic',stroke=np.pi):
        '''
        Calculates the points on a single cam given the gear ratios.
        '''
        
        r_lim = 0.5        # after interpolation, the r_min will be lower than r_lim
        for i in range(self.points): # for all the defined points on a cam
            m = self.gear_ratios[i, 0] # / self.initial_gr # normalize input to initial gear ratio
            r = np.sqrt(m) * self.cam_radii[i, 0] # smaller cam radii
            R = self.cam_radii[i, 1] / np.sqrt(m) # larger cam radii
            if r < r_lim:
                r = r_lim #* self.cam_radii[i, 0]
                R = r_lim/m# * self.cam_radii[i, 1] / m
            self.cam_radii[i, 0] = r
            self.cam_radii[i, 1] = R
        

        self.cam_radii = np.append(self.cam_radii, self.cam_radii[:2, :], axis=0)
        self.input_angles = np.append(self.input_angles, 2 * np.pi + self.input_angles[0:2])

        print("pre-convex cam_radii: ")
        print(type(self.cam_radii))


        # after generating the points for each cam, interpolate spline curves between them
        self.interpolate(kind=kind)
        
        r = self.cam_radii[:, 0]
        R = self.cam_radii[:, 1]
        x_c = np.cumsum(r * 2 * np.pi / self.ninterp)
        ratio = stroke / x_c[self.sit_ind]
        r = r * ratio
        R = R * ratio
        self.inner_pts = np.array([r * np.cos(self.angles), r * np.sin(self.angles)]).T
        self.outer_pts = np.array([R * np.cos(self.angles - np.pi / 2), R * np.sin(self.angles - np.pi / 2)]).T
        
        # redefine points as convex hull of cam shape
        self.convex_cam_pts()
        """        plt.scatter(self.inner_pts[:,0], self.inner_pts[:,1], lw = 2)
        plt.scatter(self.outer_pts[:,0], self.outer_pts[:,1], lw = 2)
        plt.legend(['inner cam','outer cam'])
        plt.axis('equal')
        plt.title('after convex function plot')
        plt.show() """

        # plot convex points and convert line plot into data
        curve_inner, = plt.plot(self.inner_pts[:,0], self.inner_pts[:,1], lw = 2)
        curve_outer,= plt.plot(self.outer_pts[:,0], self.outer_pts[:,1], lw = 2)
        xdata_inner,ydata_inner = curve_inner.get_xdata(orig=False),curve_inner.get_ydata(orig=False)
        xdata_outer,ydata_outer = curve_outer.get_xdata(orig=False),curve_outer.get_ydata(orig=False)
        inner_pts = np.array([xdata_inner,ydata_inner]).T 
        outer_pts = np.array([xdata_outer,ydata_outer]).T
        '''
        plt.scatter(xdata_outer,ydata_outer) # plot data taken from plot to check
        plt.scatter(xdata_inner,ydata_inner)
        plt.scatter(inner_pts[:,0], inner_pts[:,1]) # plot re-defined points to check. Should match the above plots
        plt.scatter(outer_pts[:,0], outer_pts[:,1])
        plt.title('plotly connected cams')
        plt.scatter([0],[0],c='r')
        plt.legend(['inner cam','outer cam'])
        plt.axis('equal')
        plt.show()
        '''

        # split cam points into 3 segments that each monotonically increase or decrease in x
        break_ind = []
        # print("len(inner_pts): ", len(inner_pts))
        for i in range(len(inner_pts)-2):
            # print("i: ", i)
            # in_pts = inner_pts[i:i+2,0]
            # print("in_pts: ", in_pts)
            # diff = np.diff(inner_pts[i:i+2,0])
            # print("diff: ", diff)
            if math.copysign(1, np.diff(inner_pts[i:i+2,0])) != math.copysign(1, np.diff(inner_pts[i+1:i+3,0])):
                break_ind.append(i+2)
        inner_array = []
        inner_array.append(inner_pts[:break_ind[0],:])
        inner_array.append(inner_pts[break_ind[0]:break_ind[1],:])
        inner_array.append(inner_pts[break_ind[1]:,:])

        """ print("inner_array_1: ", inner_array_1)
        print("inner_array_2: ", inner_array_2)
        print("inner_array_3: ", inner_array_3) """

        ind = 0
        # flip_array = np.zeros(3)
        for array in inner_array:
            # print("ind:", ind)
            # print("DIFF ARRAY: ", np.diff(array[:,0]))
            if not np.all(np.diff(array[:,0]) > 0):
                # print("array flipped!")
                inner_array[ind] = np.flip(inner_array[ind], 0)
            ind += 1

        # print("FLIPPED INNER ARRAY:", inner_array)

        break_ind = []
        for i in range(len(outer_pts)-2):
            if math.copysign(1, np.diff(outer_pts[i:i+2,0])) != math.copysign(1, np.diff(outer_pts[i+1:i+3,0])):
                break_ind.append(i+2)
        outer_array = []
        outer_array.append(outer_pts[:break_ind[0],:])
        outer_array.append(outer_pts[break_ind[0]:break_ind[1],:])
        outer_array.append(outer_pts[break_ind[1]:,:])

        ind = 0
        for array in outer_array:
            if not np.all(np.diff(array[:,0]) > 0):
                outer_array[ind] = np.flip(outer_array[ind], 0)
            ind += 1

        # linearly interpolate between points in the 3 sub-arrays for each cam
        x_int = np.linspace(inner_array[0][0,0], inner_array[0][-1,0], len(inner_array[0]))
        y_int_inner_1 = np.interp(x_int, inner_array[0][:,0], inner_array[0][:,1])
        interp_array = [x_int, y_int_inner_1]
        print("interpolated array 1: ", interp_array)
        interp_array = np.reshape(interp_array, (-1, 2))
        print("RESHAPED interpolated array 1: ", interp_array)
        print("SHAPE of interped array: ", np.shape(interp_array))

        # convert convex cam points into polar coordinates
        new_radii = self.calc_radii_convex_cams(inner_pts,outer_pts) # returns tuple of 2 ndarrays
        new_angles = self.get_convex_hull_angles(inner_pts, outer_pts) # returns a tuple of 2 ndarrays
        '''print("new_angles: ", new_angles)'''
        sorted_angle_indices_inner = np.argsort(new_angles[0]) # STS: do we need to sort? doesn't make sense to me because the radii are unsorted. radii and angles need to align

        print("shape of inner angles: ", np.shape(new_angles[0]))
        print("shape of outer angles: ", np.shape(new_angles[1]))

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        # ax.plot(new_angles[0], new_radii[0]) # inner cam # STS commented & added net 2 lines
        # ax.plot(new_angles[1], new_radii[1]) # outer cam
        ax.scatter(new_angles[0], new_radii[0]) # inner cam
        ax.scatter(new_angles[1], new_radii[1]) # outer cam
        ax.set_rticks([2, 4, 6, 8, 10])
        ax.set_rlabel_position(-22.5)  # move radial labels away from plotted line
        ax.grid(True)
        ax.set_title('radii and angles before interpolation')

        # interpolate between points in polar coordinates
        inner_cam_fcn = interpolate.interp1d(new_angles[0], new_radii[0], kind='linear',fill_value="extrapolate")
        outer_cam_fcn = interpolate.interp1d(new_angles[1], new_radii[1], kind='linear',fill_value="extrapolate")
        
        self.new_cam_radii = np.vstack((inner_cam_fcn(self.angles), outer_cam_fcn(self.angles))).T # get radii of convex cams at regularly-spaced angles
        print("shape of new_cam_radii: ", self.new_cam_radii.shape)

        r = self.new_cam_radii[:, 0]
        R = self.new_cam_radii[:, 1]

        """ print("r: ", r)
        print("R: ", R) """

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(self.angles - np.pi / 2, R )  # apply phase change to the outer cam
        ax.plot(self.angles, r)
        ax.set_rticks([2, 4, 6, 8, 10])
        ax.set_rlabel_position(-22.5)  # move radial labels away from plotted line
        ax.grid(True)
        ax.set_title('interpolated cam')

        # plt.plot(self.inner_pts[:,0], self.inner_pts[:,1], lw = 2)
        # plt.plot(self.outer_pts[:,0], self.outer_pts[:,1], lw = 2)
        # plt.legend(['inner cam','outer cam'])
        # plt.axis('equal')
        # plt.show()

        # print(self.cam_radii.shape,'radii shape')
        # print('new radii',new_radii,new_radii.shape)

        # print('inner radii',self.cam_radii[:, 0],'new radii',new_radii[0])
        # print('outer radii',self.cam_radii[:, 1],'new radii',new_radii[1])
        # apply phase change to the outer cam
        #self.cam_radii[:, 1] = np.append(self.cam_radii[int(self.ninterp/4):, 1], self.cam_radii[:int(self.ninterp/4), 1], axis=0)
        r_sum = np.sum(self.cam_radii[:, 0])
        R_sum = np.sum(self.cam_radii[:, 1])

        r_max = max(self.cam_radii[:, 0])
        R_max = max(self.cam_radii[:, 1])

        R_plus_r = r_max+R_max

        R_over_r_ratio = R_sum/r_sum
    

        return R_plus_r

    def interpolate(self, kind='cubic'):
        '''
        Interpolate the radii with spline curve, kind = "quadratic or cubic"
        interpolation kind pls refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
        '''
        print(self.input_angles)

        inner_cam_fcn = interpolate.interp1d(self.input_angles, self.cam_radii[:, 0], kind=kind)
        outer_cam_fcn = interpolate.interp1d(self.input_angles, self.cam_radii[:, 1], kind=kind)

        self.cam_radii = np.vstack((inner_cam_fcn(self.angles), outer_cam_fcn(self.angles))).T

    def cam_pts(self, stroke=np.pi, convex=True):
        r = self.cam_radii[:, 0]
        R = self.cam_radii[:, 1]
        x_c = np.cumsum(r * 2 * np.pi / self.ninterp) # cumulative arc length at each interpolation point along inner cam

        # scale cam size to desired stroke length (amount of cable displacement from sitting down)
        ratio = stroke / x_c[self.sit_ind]
        r = r * ratio
        R = R * ratio
        self.inner_pts = np.array([r * np.cos(self.angles), r * np.sin(self.angles)]).T
        self.outer_pts = np.array([R * np.cos(self.angles - np.pi / 2), R * np.sin(self.angles - np.pi / 2)]).T # apply phase change to the outer cam
        if convex:
            self.convex_cam_pts()
            self.plot_cams()
        return self.inner_pts, self.outer_pts

    def convex_cam_pts(self):
        inner_hull = ConvexHull(self.inner_pts,incremental=True)
        # print('IN CONVEX FUNC',inner_hull.points[inner_hull.vertices])
        """ plt.scatter(inner_hull.points[inner_hull.vertices][:,0],inner_hull.points[inner_hull.vertices][:,1])
        plt.title('Inner hull points (in convex_cam_pts function)')
        plt.show() """
        outer_hull = ConvexHull(self.outer_pts,incremental=True)
        self.inner_pts = np.vstack((self.inner_pts[inner_hull.vertices, 0], self.inner_pts[inner_hull.vertices, 1])).T
        self.outer_pts = np.vstack((self.outer_pts[outer_hull.vertices, 0], self.outer_pts[outer_hull.vertices, 1])).T

        """ plt.plot(self.inner_pts[:,0], self.inner_pts[:,1], lw = 2)
        plt.plot(self.outer_pts[:,0], self.outer_pts[:,1], lw = 2)
        plt.legend(['inner cam','outer cam'])
        plt.axis('equal')
        plt.title('Inner & outer points (in convex_cam_pts function)')
        plt.show() """
        

    def alter_cam_profiles(self, routing=None):
        '''
        routing: 2*2 np.array, The first row defines the anchor position (in angle) and cable thickness of the inner cam,
                               while the second row defines those of the outer cam.
        [[anchor angle of inner cam, thickness of pulling cable, screw_diameter],
        [anchor angle of outer cam, thickness of elastic band,  screw_diameter]]
        '''
        if routing is not None:
            inner_anchor = routing[0, 0]
            outer_anchor = routing[1, 0]
            t_inner = routing[0, 1]
            t_outer = routing[1, 1]
            outer_screw = routing[0, 2]
            inner_screw = routing[1, 2]
            inner_profile = t_inner * np.zeros(np.sum(self.angles>inner_anchor), 1)
            self.cam_radii[self.angles > inner_anchor, 0] = self.cam_radii[self.angles>inner_anchor, 0] - inner_profile

    def plot_cams(self,cam_radii = 0 ,stroke=np.pi):
        '''
        Plots the cam points as a sanity check feature.
        stroke: scalar, full stroke of the pulling cable
        '''
        if cam_radii.all() == 0:
            cam_radii = self.cam_radii
        r = cam_radii[:, 0]
        R = cam_radii[:, 1]
        x_c = np.cumsum(r * 2 * np.pi / self.ninterp) 
        ratio = stroke / x_c[self.sit_ind] # full stroke divided by distance along cam perimeter of sitting angle
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(self.angles - np.pi / 2, R * ratio)  # apply phase change to the outer cam
        ax.plot(self.angles, r*ratio)
        ax.set_rticks([2, 4, 6, 8, 10])
        ax.set_rlabel_position(-22.5)  # move radial labels away from plotted line
        ax.grid(True)

        ax.set_title(f"Cam Demonstration, min inner radius={r.min()*ratio}, max outer radius={R.max()*ratio}", va='bottom')

        # print("sit index: ", self.sit_ind)
        # print("x_c[sit_ind]: ", x_c[self.sit_ind])
        # print("ratio: ", ratio)


    def plot_forces(self, torque=False, stroke=2*np.pi):
        '''
        Plot the force profiles vs angular displacement
        given the gear ratios and solved cam radii.
        stroke: scalar, full stroke of the pulling cable
        '''
        
        # arrays of radii for inner cam (r) and outer cam (R)
        r = self.cam_radii[:, 0]
        R = self.cam_radii[:, 1]

        # scale cable displacement and radii to match stroke length
        x_c = np.cumsum(r * 2 * np.pi / self.ninterp) # cable displacement distance
        ratio = stroke / x_c[self.sit_ind] # ratio of stroke length to unscaled displacement at sitting index
        x_c = x_c * ratio
        r = r * ratio
        R = R * ratio

        # calculate elastic band forces
        x_e = np.cumsum(R * 2 * np.pi / self.ninterp)*2     # elastic band displacement
        k = 2*14/x_e[self.sit_ind]**2
        # print(f"Elastic band stiffness {k} N/m; Pulling distance {x_e[self.sit_ind]} m")
        f_e = k * x_e # force exerted by elastic band (assumes linear stiffness)

        # calculate cable forces
        f_c = R * f_e / r # force on cable is a function of elastic force and gear ratio

        # normalize_angles(self, )
        
        plt.figure()
        plt.plot(self.angles / 2 / np.pi * 360, f_e, label='Elastic Band Force', linewidth=3)
        plt.plot(self.angles / 2 / np.pi * 360, f_c, label='Cable Force', linewidth=3)
        # plt.plot(x_c, f_e, label='Elastic Band Force', linewidth=3)
        # plt.plot(x_c, f_c, label='Cable Force', linewidth=3)
        #plt.scatter(self.angles / 2 / np.pi * 360, f_e * r, label='torque')
        plt.legend(loc="lower right")
        plt.xlabel('Angles (deg)')
        # plt.xlabel('pulling distance/cam perimeter distance (m)')
        plt.ylabel('Force (N)')
        # plt.xlim([0, 0.17])
        plt.xlim([0, 300])
        plt.ylim([0, 130])
        # plt.gca().set_xticks([0, 0.04, 0.08, 0.12, 0.16])
        plt.title("Force Profiles of the Mechanism")

        np.savetxt('force_output.csv', np.stack((self.angles, x_c, f_e, f_c), axis=1),  header='angles(rad), pulling distance(m), elastic band force(n), cable force(n)')
        if torque:
            plt.figure()
            plt.scatter(self.angles / 2 / np.pi * 360, f_e * r)
            plt.xlabel('angle (deg)')
            plt.ylabel('torque (N*m)')
        plt.show()

    def plot_forces_percentage(self, file, torque=False, stroke = 2 * np.pi,plot = True):
        '''
        Plot the force profiles vs *stance percentage*
        given the gear ratios and solved cam radii.
        stroke: scalar, full stroke of the pulling cable
        '''

        # arrays of radii for inner cam (r) and outer cam (R)
        r = self.cam_radii[:, 0]
        R = self.cam_radii[:, 1]

        # scale cable displacement and cam radii to match stroke length
        x_c = np.cumsum(r * 2 * np.pi / self.ninterp)       # cable displacement array over entire cam rotation
        ratio = stroke / x_c[self.sit_ind]                  # ratio of stroke length to unscaled displacement at sitting index
        # print("stroke: ", stroke)
        # print("raw sitting displacement, x_c[self.sit_ind]: ", x_c[self.sit_ind])
        # print("ratio: ", ratio)
        x_c = x_c * ratio
        r = r * ratio
        R = R * ratio

        # calculate elastic band tension from displacement
        x_e = np.cumsum(R * 2 * np.pi / self.ninterp) * 2     # elastic band displacement array
        k = 2 * 14 / x_e[self.sit_ind] ** 2                   # elastic stiffness (calculated from experimental values)
        # print(f"Elastic band stiffness {k} N/m; Pulling distance {x_e[self.sit_ind]} m")
        f_e = k * x_e # force exerted by elastic band (assumes linear stiffness)

        # calculate cable forces from elastic band tension
        f_c = R * f_e / r # tension on cable is a function of elastic tension and gear ratio

        ### account for nonlinear change in knee angle during standing ###
        # load experimental knee data
        knee_array = np.loadtxt(file, delimiter = ',', ndmin = 2)
        percentage_raw = knee_array[:,0];  # percentage of stance (~0% to ~100% going from sitting to standing)
        knee_angles_raw = knee_array[:,1]; # corresponding knee angle

        # interpolate & extrapolate knee angle array to match size of cam arrays
        knee_fcn = interpolate.interp1d(percentage_raw, knee_angles_raw, kind = 'cubic', fill_value = 'extrapolate')
        percentage = np.linspace(0, 100, self.sit_ind + 1)
        knee_angles = knee_fcn(percentage)

        # find linear relationship between stand percentage and cable displacement
        knee_angles_flip = np.flip(knee_angles) # flip knee array to go from standing to sitting
        percentage_flip = np.flip(percentage) # flip percentage stance to correspond to knee angle (standing to sitting)
        xc_knee = x_c[self.sit_ind] * (1 - (knee_angles_flip - np.min(knee_angles_flip)) / np.ptp(knee_angles_flip) ) # scale cable displacement to knee angle (goes from standing to sitting)
        
        # plot knee angle vs. cable displacement
        if plot:
            plt.figure()
            plt.plot(knee_angles_flip, xc_knee)
            plt.xlabel('knee angle (degree)')
            plt.ylabel('cable displacement, xc (m)')
            plt.title('knee angle and cable displacement linear relationship')

            # plot stance percentage vs. cable displacement
            plt.figure()
            plt.plot(percentage_flip, xc_knee)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cable Displacement (m)')

        # add "negative" stance phase to account for part of cam past sitting index
        perc_neg = np.linspace((self.sit_ind - self.ninterp) * (100 / self.sit_ind), -100 / self.sit_ind, self.ninterp - self.sit_ind)
        xc_knee_neg = np.linspace(xc_knee[-1], xc_knee[-1] + (self.ninterp - self.sit_ind) * ((xc_knee[-1] - xc_knee[ np.int_(3*self.sit_ind/4) ] ) / (self.sit_ind/4)), self.ninterp - self.sit_ind)

        percentage = np.concatenate((perc_neg, percentage))
        percentage_flip = np.flip(percentage)
        xc_knee = np.concatenate((xc_knee, xc_knee_neg))
        
        # create functions relating cable and elastic displacement to cam angle
        cable_fcn = interpolate.interp1d(x_c, self.angles, kind='cubic', fill_value='extrapolate') # cable displacement in, angle out
        elastic_fcn = interpolate.interp1d(self.angles, x_e, kind='cubic', fill_value='extrapolate') # angle in, elastic displacement out

        # calculate cam angle from cable displacement that follows stance percentage
        angle_percentage = cable_fcn(xc_knee)

        # calculate elastic displacement from angle
        elastic_percentage = elastic_fcn(angle_percentage)

        # calculate cable force using stance percentage as independent variable
        f_c_new = R / r * k * elastic_percentage
        
        ### plot and save values ###
        # plot cam angle vs. stance percentage
        if plot:
            plt.figure()
            plt.plot(percentage_flip, angle_percentage)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Cam Angle (rad)')

            # plot elastic band displacement vs. stance percentage
            plt.figure()
            plt.plot(percentage_flip, elastic_percentage)
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Elastic displacement (m)')
            
            # plot cable force vs. stance percentage
            plt.figure()
            plt.plot(np.flip(percentage), f_c_new, label='Cable Force', linewidth=3)
            plt.legend(loc="upper right")
            plt.xlabel('Stance Percentage (%)')
            plt.ylabel('Force (N)')
            plt.xlim([-60, 100])
            plt.ylim([0, 250])
            
            np.savetxt('force_output.csv', np.stack((self.angles, x_c, f_e, f_c), axis=1),  header='angles(rad), pulling distance(m), elastic band force(n), cable force(n)')
            if torque:
                plt.figure()
                plt.scatter(self.angles / 2 / np.pi * 360, f_e * r)
                plt.xlabel('angle (deg)')
                plt.ylabel('torque (N*m)')
            plt.show()

        return f_c_new[:-10],np.flip(percentage)[:-10] # sitting force is the force at 0%


