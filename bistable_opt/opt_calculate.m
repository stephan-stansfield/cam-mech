function [keypoints, soln_angles] = opt_calculate(vars)

   syms b e real

   [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, alph, gamma, delta, epsilon] = unpack_vars(vars);
    
    % Set up system of equations to solve for b & e
    [sol_b, sol_eta] = solve(l_OA*sin(alph) + l_AC*sin(alph + b) == ...
        (l_OA+l_AB)*sin(alph) + l_BD*sin(alph + gamma) + l_CD*sin(alph + gamma + e),...
        -l_OA*cos(alph) - l_AC*cos(alph + b) == ...
        -(l_OA+l_AB)*cos(alph) - l_BD*cos(alph + gamma) - l_CD*cos(alph + gamma + e));
    
    sol_eta = eval(sol_eta);
    sol_b = eval(sol_b);

    sol_eta_deg = rad2deg(sol_eta)
    sol_b_deg = rad2deg(sol_b)


    % if no solution was found, assign zeros to keypoints and angle and end execution
    if isempty(sol_eta)
        keypoints = zeros(4,6);
        soln_angles = zeros(4,1);
        return
    end
    
    %% Compute forward kinematics to EE and input points for both solutions
    
    %% First solution
    % Choose one solution
    b1 = sol_b(1);
    e1 = sol_eta(1);
    
    % Get keypoints of solved kinematic chain
    [rA1, rB1, rC1, rD1, rE1, rF1] = fun_FK(b1, e1, vars);

    %% Second solution
    b2 = sol_b(2);
    e2 = sol_eta(2);
    
    % Get keypoints of solved kinematic chain
    [rA2, rB2, rC2, rD2, rE2, rF2] = fun_FK(b2, e2, vars);

    % Return keypoints of both solutions in a multi-dim array
    keypoints = [rA1, rB1, rC1, rD1, rE1, rF1;
                 rA2, rB2, rC2, rD2, rE2, rF2];

    soln_angles = [b1, e1, b2, e2];

end