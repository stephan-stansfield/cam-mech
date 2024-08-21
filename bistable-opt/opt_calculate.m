function keypoints = opt_calculate(vars)

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


    % if no solution was found, assign zeros to keypoints and end execution
    if isempty(sol_eta)
        keypoints = zeros(4,6);
        return
    end
    
    %% Compute forward kinematics to EE and input points for both solutions
    
    %% First solution
    % Choose one solution
    b = sol_b(1);
    e = sol_eta(1);
    
    % Get keypoints of solved kinematic chain
    [rA1, rB1, rC1, rD1, rE1, rF1] = fun_FK(b, e, vars);

    %% Second solution
    b = sol_b(2);
    e = sol_eta(2);
    
    % Get keypoints of solved kinematic chain
    [rA2, rB2, rC2, rD2, rE2, rF2] = fun_FK(b, e, vars);

    % Return keypoints of both solutions in a multi-dim array
    keypoints = [rA1, rB1, rC1, rD1, rE1, rF1;
                 rA2, rB2, rC2, rD2, rE2, rF2];

end