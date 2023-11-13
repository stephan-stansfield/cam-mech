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
    [rA1, rB1, rC1, rD1, rE1, rF1, rA2, rB2, rC2, rD2, rE2, rF2] = unpack_keypoints(keypoints);

    % Evaluate constraints
    % Overall motion envelope constraints
    xArray = horzcat(keypoints(1,:), keypoints(3,:)); % x-coordinates of all points of all joints
    yArray = horzcat(keypoints(2,:), keypoints(4,:)); % y-coordinates of all points of all joints
    xEnv = max(xArray) - min(xArray);
    yEnv = max(yArray) - min(yArray);
    xEnvMax = 45;
    yEnvMax = 120;
    
    % End-effector point travel
    EEy1 = keypoints(2, 5);
    EEy2 = keypoints(4, 5);
    EEdx = keypoints(3,5) - keypoints(1,5);
%     EEdy = keypoints(4,5) - keypoints(2,5);
    dx_EE_min = 25.4/4; % 1/4" minimum X travel (1/8" gap between 2 1/8"-wide channels)

    % Input point travel
    INy1 = keypoints(2, 6);
    INy2 = keypoints(4, 6);
    INx1 = keypoints(1, 6);
    INx2 = keypoints(3, 6);
    INdy = INy2 - INy1;
    INdx = INx2 - INx1;

    % Intersection with wall constraint
    m = (rB1(2) - rA1(2)) / (rB1(1) - rA1(1)); % b is constrained to be at (0,0)
%     b = rA1(2) - m * rA1(1);

%     % End effector travel goals (mm)
%     EEdx_des = 44;
%     EEdy_des = 24;

    % Input point travel constraints
%     INdx_min = 5;
%     INdy_max = 15;
    INdx_min = 7.5;
    INdy_max = 7.5;

    % Inequality constraints: 
%     cineq = [ xEnv - xEnvMax, yEnv - yEnvMax, dx_EE_min - abs(dx_EE)];      % envelope dimensions are less than limits
%     cineq = [ xEnv - xEnvMax, yEnv - yEnvMax, dx_EE_min - abs(dx_EE),...    % envelope dimensions are less than limits
%               m*rC1(1) - rC1(2), m*rD1(1) - rD1(2), m*rE1(1) - rE1(2),...   % all points lie to top-right of line between A & B
%               m*rF1(1) - rF1(2), ...
%               m*rC2(1) - rC2(2), m*rD2(1) - rD2(2), m*rE2(1) - rE2(2),...   % all points lie to top-right of line between A & B
%               m*rF2(1) - rF2(2)];
    cineq = [ xEnv - xEnvMax, yEnv - yEnvMax,...                            % envelope dimensions are less than limits
              m*rC1(1) - rC1(2), m*rD1(1) - rD1(2), m*rE1(1) - rE1(2),...   % all points lie to top-right of line between A & B
              m*rF1(1) - rF1(2), ...
              m*rC2(1) - rC2(2), m*rD2(1) - rD2(2), m*rE2(1) - rE2(2),...   
              m*rF2(1) - rF2(2),...
              INdx_min - INdx, INdy - INdy_max, ...                         % input point travel constraints (x minimum, y maximum)
              dx_EE_min - abs(EEdx), ...                                    % end effector minimum x travel
              INy1 - EEy1, INy2 - EEy2];                                    % input point is always below end effector point in y

    % Equality constraints
    ceq = [];
%     ceq = [abs(EEdx) - EEdx_des, abs(EEdy) - EEdy_des]; % End effector travel

end