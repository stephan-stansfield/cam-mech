function [cineq, ceq] = opt_constraints(vars)
% Inputs:
% vars - an array of link lengths and angles
% 
% Outputs:
% cineq - an array of values of nonlinear inequality constraint functions.  
%         The constraints are satisfied when these values are less than zero.
% ceq   - an array of values of nonlinear equality constraint functions.
%         The constraints are satisfied when these values are equal to zero.

    % Solve kinematic chain with provided inputs. If a solution exists, 
    % returns x- and y-coordinates of the two solutions in an array. If a
    % solution doesn't exist, returns same-sized array of zeros.
    keypoints = opt_calculate(vars);
    [rA1, rB1, rC1, rD1, rE1, rF1, rG1, rH1, rA2, rB2, rC2, rD2, rE2,...
        rF2, rG2, rH2] = unpack_keypoints(keypoints);

    % Evaluate constraints
    % Overall motion envelope constraints
    xArray = horzcat(keypoints(1,:), keypoints(3,:)); % x-coordinates of all points of all joints
    yArray = horzcat(keypoints(2,:), keypoints(4,:)); % y-coordinates of all points of all joints
    xEnv = max(xArray) - min(xArray);
    yEnv = max(yArray) - min(yArray);
    xEnvMax = 50;
    yEnvMax = 100;
    
    % Minimum end-effector x travel constraint
    dx_EE = keypoints(1,7) - keypoints(3,7);
    dx_EE_min = 10;

    % Intersection with wall constraint
    m = (rB1(2) - rA1(2)) / (rB1(1) - rA1(1)); % b is constrained to be at (0,0)
%     b = rA1(2) - m * rA1(1);

    % Inequality constraints: 
    cineq = [ xEnv - xEnvMax, yEnv - yEnvMax, dx_EE_min - abs(dx_EE)];      % envelope dimensions are less than limits
%     cineq = [ xEnv - xEnvMax, yEnv - yEnvMax, dx_EE_min - abs(dx_EE),...    % envelope dimensions are less than limits
%               m*rC1(1) - rC1(2), m*rD1(1) - rD1(2), m*rE1(1) - rE1(2),...   % all points lie to top-right of line between A & B
%               m*rF1(1) - rF1(2), m*rG1(1) - rG1(2), m*rH1(1) - rH1(2),...
%               m*rC2(1) - rC2(2), m*rD2(1) - rD2(2), m*rE2(1) - rE2(2),...   % all points lie to top-right of line between A & B
%               m*rF2(1) - rF2(2), m*rG2(1) - rG2(2), m*rH2(1) - rH2(2)];

    ceq = 0;

end