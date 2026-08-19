close all; clear; clc;

tic;

% Set initial guess for lengths and angles
% All link lengths in mm
% l_OA = 5;
% l_AB = sqrt(2)*5+5;
% l_AC = 7.5;
% l_BD = 5;
% l_CD = 7.5;
% l_CE = 15;
% l_DF = 7.5;
l_OA = 5;
l_AB = sqrt(2)*5+5;
l_AC = 5;
l_BD = 5;
l_CD = 5;
l_CE = 25;
l_DF = 12.5;

% All angles in radians, converted from degrees
alph = deg2rad(30);
gamma = deg2rad(130);
delta = deg2rad(60);
epsilon = deg2rad(-15);

% initial guess vector:
x0 = [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, alph, gamma, delta, epsilon];
%%
% set up and solve nonlinear programming problem
problem.objective = @(x) opt_objective(x);
problem.nonlcon = @(x) opt_constraints(x);
problem.x0 = x0;
%            [lOA lAB lAC lBD lCD lCE lDF alpha gamma delta epsilon]
problem.lb = [ 1   1   1   1   2.5 1   1   0     0     0     -0.5];
problem.ub = [ 30  30  30  30  30  50  30  pi    pi    pi    pi ];
problem.Aineq = []; problem.bineq = [];
problem.Aeq = []; problem.beq = [];
% f_viz = @(keypoints)opt_visualize(keypoints, save);
problem.options = optimoptions('fmincon','Display','iter-detailed',...
    'Algorithm','interior-point','EnableFeasibilityMode',true,...
    'OutputFcn', @opt_visualize, 'OptimalityTolerance', 1.0000e-09);
problem.solver = 'fmincon';
x = fmincon(problem);

% Assign optimal values back to variable names
l_OA = x(1);
l_AB = x(2);
l_AC = x(3);
l_BD = x(4);
l_CD = x(5);
l_CE = x(6);
l_DF = x(7);
alph = x(8);
gamma = x(9);
delta = x(10);
epsilon = x(11);
alph_deg = rad2deg(alph);
gamma_deg = rad2deg(gamma);
delta_deg = rad2deg(delta);
epsilon_deg = rad2deg(epsilon);

% Get solution keypoints from optimized lengths & angles
[keypoints, ~] = opt_calculate(x);

toc;

%% Visualize optimized solution
% opt_visualize(x, struct([]), []); % superseded; all iterations are saved
[r_A1, r_B1, r_C1, r_D1, r_E1, r_F1, r_A2, r_B2, r_C2, r_D2, r_E2, r_F2] = unpack_keypoints(keypoints);

EEdx = keypoints(3, 5) - keypoints(1, 5);
EEdy = keypoints(4, 5) - keypoints(2, 5);

INdx = keypoints(3, 6) - keypoints(1, 6);
INdy = keypoints(4, 6) - keypoints(2, 6);

EEdx_des = 2.38*1.50*2;
EEdy_des = 8.19;

diff_x = EEdx - EEdx_des;
diff_y = EEdy - EEdy_des;

% Calculate envelope dimensions
x_array = horzcat(keypoints(1,:), keypoints(3,:)); % x-coordinates of all points of all joints
y_array = horzcat(keypoints(2,:), keypoints(4,:)); % y-coordinates of all points of all joints
x_env = max(x_array) - min(x_array);
y_env = max(y_array) - min(y_array);

%% Save workspace variables to file
date = string(datetime('today', 'Format', 'yyyy-MM-dd'));
save(strcat("soln_values_", date))


