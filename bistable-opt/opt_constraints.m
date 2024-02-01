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
    x_array = horzcat(keypoints(1,:), keypoints(3,:)); % x-coordinates of all points of all joints
    y_array = horzcat(keypoints(2,:), keypoints(4,:)); % y-coordinates of all points of all joints
    x_env = max(x_array) - min(x_array);
    y_env = max(y_array) - min(y_array);
    x_envMax = 45/2.2;
    y_envMax = 120/2.2;
    
    % End-effector point travel
    EE_y_1 = keypoints(2, 5);
    EE_y_2 = keypoints(4, 5);
    EE_dx = keypoints(3,5) - keypoints(1,5);
%     EE_dy = keypoints(4,5) - keypoints(2,5);
    D_cable = 1/16 * 25.4; % 1/16" diameter cable
    scale_groove = 1.50; % groove diameter is 150% larger than cable diameter
    D_groove = scale_groove * D_cable;
    EE_dx_min = 2*D_groove; % preserve one groove width between grooves
    EE_dx_max = 4*D_groove; % have at most three widths between grooves

    % Input point travel
    IN_y_1 = keypoints(2, 6);
    IN_y_2 = keypoints(4, 6);
    IN_x_1 = keypoints(1, 6);
    IN_x_2 = keypoints(3, 6);
    IN_dy = IN_y_2 - IN_y_1;
    IN_dx = IN_x_2 - IN_x_1;

    % Intersection with wall constraint
    m = (rB1(2) - rA1(2)) / (rB1(1) - rA1(1)); % b is constrained to be at (0,0)

%     % End effector travel goals (mm)
%     EE_dx_des = 44;
%     EE_dy_des = 24;

    % Input point travel constraints
    dx_IN_min = 7.5/2.2;
    dy_IN_max = 7.5/2.2;

    % Inequality constraints: 
%     cineq = [ x_env - x_envMax, y_env - y_envMax, EE_dx_min - abs(dx_EE)];      % envelope dimensions are less than limits
%     cineq = [ x_env - x_envMax, y_env - y_envMax, EE_dx_min - abs(dx_EE),...    % envelope dimensions are less than limits
%               m*rC1(1) - rC1(2), m*rD1(1) - rD1(2), m*rE1(1) - rE1(2),...   % all points lie to top-right of line between A & B
%               m*rF1(1) - rF1(2), ...
%               m*rC2(1) - rC2(2), m*rD2(1) - rD2(2), m*rE2(1) - rE2(2),...   
%               m*rF2(1) - rF2(2)];
    cineq = [ x_env - x_envMax, y_env - y_envMax,...                        % envelope dimensions are less than limits
              m*rC1(1) - rC1(2), m*rD1(1) - rD1(2), m*rE1(1) - rE1(2),...   % all points lie to top-right of line between A & B
              m*rF1(1) - rF1(2), ...
              m*rC2(1) - rC2(2), m*rD2(1) - rD2(2), m*rE2(1) - rE2(2),...   
              m*rF2(1) - rF2(2),...
              dx_IN_min - IN_dx, IN_dy - dy_IN_max, ...                     % input point travel constraints (x minimum, y maximum)
              EE_dx_min - abs(EE_dx), ...                                   % end effector minimum x travel
              abs(EE_dx) - EE_dx_max, ...                                   % end effector maximum x travel
              IN_y_1 - EE_y_1, IN_y_2 - EE_y_2];                            % input point is always below end effector point in y

    % Equality constraints
    ceq = [];
%     ceq = [abs(EE_dx) - EE_dx_des, abs(EE_dy) - EE_dy_des]; % End effector travel

end