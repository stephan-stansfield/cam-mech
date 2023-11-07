function f = opt_objective(vars)
% Inputs:
% vars - an array of link lengths and angles
% 
% Outputs:
% f - scalar value of the function (to be minimized) evaluated for the
%     provided values of the decision variables.

%     % Unpack inputs for simulation
%     tf = x(1);
%     ctrl.tf = x(2);
%     ctrl.T = x(3:end);
%     tspan = [0 tf];

    % Solve kinematic chain with provided inputs. If a solution exists, 
    % returns x- and y-coordinates of the two solutions in an array.
    keypoints = opt_calculate(vars);
    y1 = keypoints(2, 5);
    y2 = keypoints(4, 5);
    x1 = keypoints(1, 5);
    x2 = keypoints(3, 5);

    wy = 1;
    wx = 1;

    % Evaluate objective function
%     f = -wy*(y2 - y1)^2 + wx*(x2 - x1)^2; % maximize end-effector change in Y, min. change in X
    f = -wy*(y2 - y1)^2;
    
end