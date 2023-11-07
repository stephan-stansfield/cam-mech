close all; clear; clc;

% Set initial guess
% All link lengths in mm
l_OA = 10;
l_AB = sqrt(2)*10+10;
l_AC = 10;
l_BD = 10;
l_CD = 10;
l_CE = 5;
l_DF = l_CE;
l_EF = l_CD;
l_EG = 50;
l_FH = 25;

% All angles in radians, converted from degrees
alph = deg2rad(30);
zeta = deg2rad(130);
gamma = deg2rad(90);
delta = deg2rad(60);
epsilon = deg2rad(-15);

% initial guess vector:
v0 = [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, l_EF, l_EG, l_FH, alph, zeta, gamma, delta, epsilon];

% setup and solve nonlinear programming problem
problem.objective = @(x) opt_objective(x);                   % create anonymous function that returns objective
problem.nonlcon = @(x) opt_constraints(x);                   % create anonymous function that returns nonlinear constraints
problem.x0 = v0;                                                % initial guess for decision variables
problem.lb = [1 1 1 1 1 1 1 1 1 1 0 0 0 0 -0.5];                   % lower bound on decision variables
problem.ub = [50 50 50 50 50 50 50 50 75 50 pi pi pi pi pi];    % upper bound on decision variables
problem.Aineq = []; problem.bineq = [];                         % no linear inequality constraints
problem.Aeq = []; problem.beq = [];                             % no linear equality constraints
problem.options = optimoptions('fmincon','Display','iter',...
    'Algorithm','interior-point','EnableFeasibilityMode',true); % set options
problem.solver = 'fmincon';                                     % required
x = fmincon(problem);                                           % solve nonlinear programming problem

% Get solution from optimized variables
keypoints = opt_calculate(x);

%% Visualize optimized solution
opt_visualize(keypoints);
[rA1, rB1, rC1, rD1, rE1, rF1, rG1, rH1, rA2, rB2, rC2, rD2, rE2,...
        rF2, rG2, rH2] = unpack_keypoints(keypoints);


