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
%     EEdx = keypoints(3, 5) - keypoints(1, 5);
    EEdy = keypoints(4, 5) - keypoints(2, 5);
%     INdx = keypoints(3, 6) - keypoints(1, 6);
%     INdy = keypoints(4, 6) - keypoints(2, 6);

    % Weights
    wx = 1;
    wy = 1;

    % Goals (mm)
    EEdx_des = 44;
    EEdy_des = 24;
    INdx_des = 22;
    INdy_des = 0;

    % Evaluate objective function
%     f = -wy*(y2 - y1)^2 + wx*(x2 - x1)^2; % maximize end-effector change in Y, min. change in X
%     f = -wy*(y2 - y1)^2; @ maximize end-effector change in y
    % minimize error between actual and desired end-effector distance traveled
%     f = wx*(EEdx - EEdx_des)^2 + wy*(EEdy - EEdy_des)^2 + ...               
%         wx*(INdx - INdx_des)^2 + wy*(INdy - INdy_des)^2 ; 
%     f = wx*(EEdx - EEdx_des)^2 + wy*(EEdy - EEdy_des)^2 ;
    f = wy*(EEdy - EEdy_des)^2 ;
    
end