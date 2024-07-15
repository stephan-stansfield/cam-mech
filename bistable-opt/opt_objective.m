function f = opt_objective(vars)
% Inputs:
% vars - an array of link lengths and angles
% 
% Outputs:
% f - scalar value of the function (to be minimized) evaluated for the
%     provided values of the decision variables.

    % Solve kinematic chain with provided inputs. If a solution exists, 
    % returns x- and y-coordinates of the two solutions in an array.
    keypoints = opt_calculate(vars);
%     EE_dx = keypoints(3, 5) - keypoints(1, 5);
    EE_dy = keypoints(4, 5) - keypoints(2, 5);
%     IN_dx = keypoints(3, 6) - keypoints(1, 6);
%     IN_dy = keypoints(4, 6) - keypoints(2, 6);

    % Calculate envelope dimensions
    x_array = horzcat(keypoints(1,:), keypoints(3,:)); % x-coordinates of all points of all joints
    x_env = max(x_array) - min(x_array);

    % Weights
    wx = 1;
    wy = 1;

    % Goals (mm)
%     EE_dx_des = 44;
    EE_dy_des = 8;
%     IN_dx_des = 22;
%     IN_dy_des = 0;

    % Evaluate objective function
%     f = -wy*(y2 - y1)^2 + wx*(x2 - x1)^2; % maximize end-effector change in Y, min. change in X
%     f = -wy*(y2 - y1)^2; @ maximize end-effector change in y
    % minimize error between actual and desired end-effector distance traveled
%     f = wx*(EE_dx - EE_dx_des)^2 + wy*(EE_dy - EE_dy_des)^2 + ...               
%         wx*(IN_dx - IN_dx_des)^2 + wy*(IN_dy - IN_dy_des)^2 ; 
%     f = wx*(EE_dx - EE_dx_des)^2 + wy*(EE_dy - EE_dy_des)^2 ;
    f = wy*(EE_dy - EE_dy_des)^2;                                           % y-distance of end effector only
%     f = wy*(EE_dy - EE_dy_des)^2 + wx*(x_env)^2;                            % y-distance of end effector + x-envelope
    
end