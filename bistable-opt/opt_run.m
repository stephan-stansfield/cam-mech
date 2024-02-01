close all; clear; clc;

tic; % start timing

% Set initial guess for lengths and angles
% All link lengths in mm
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
v0 = [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, alph, gamma, delta, epsilon];

% set up and solve nonlinear programming problem
problem.objective = @(x) opt_objective(x);                      % create anonymous function that returns objective
problem.nonlcon = @(x) opt_constraints(x);                      % create anonymous function that returns nonlinear constraints
problem.x0 = v0;                                                % initial guess for decision variables
%            [lOA lAB lAC lBD lCD lCE lDF alpha gamma delta epsilon]
problem.lb = [ 1   1   1   1   2.5 1   1   0     0     0     -0.5]; % lower bound on decision variables
problem.ub = [ 30  30  30  30  30  50  30  pi    pi    pi    pi ];  % upper bound on decision variables
problem.Aineq = []; problem.bineq = [];                         % no linear inequality constraints
problem.Aeq = []; problem.beq = [];                             % no linear equality constraints
problem.options = optimoptions('fmincon','Display','iter',...
    'Algorithm','interior-point','EnableFeasibilityMode',true); % set options
problem.solver = 'fmincon';                                     % required
x = fmincon(problem);                                           % solve nonlinear programming problem

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
keypoints = opt_calculate(x);

toc; % stop timing

%% Visualize optimized solution
opt_visualize(keypoints);
[rA1, rB1, rC1, rD1, rE1, rF1, rA2, rB2, rC2, rD2, rE2, rF2] = unpack_keypoints(keypoints);

EEdx = keypoints(3, 5) - keypoints(1, 5);
EEdy = keypoints(4, 5) - keypoints(2, 5);

INdx = keypoints(3, 6) - keypoints(1, 6);
INdy = keypoints(4, 6) - keypoints(2, 6);

EEdx_des = 2.38*1.50*2;
EEdy_des = 11.5;

diff_x = EEdx - EEdx_des;
diff_y = EEdy - EEdy_des;

%% Save workspace variables to file
date = "2024-02-01";
save(strcat("soln_values_", date))


