function keypoints = opt_calculate(vars)

   syms b e real

   [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, l_EF, l_EG, l_FH, alph, zeta, gamma, delta, epsilon] = unpack_vars(vars);
    
    % Set up system of equations to solve for b & e
    [sol_b, sol_eta] = solve(l_OA*sin(alph) + l_AC*sin(alph + b) == ...
        (l_OA+l_AB)*sin(alph) + l_BD*sin(alph + zeta) + l_CD*sin(alph + zeta + e),...
        -l_OA*cos(alph) - l_AC*cos(alph + b) == ...
        -(l_OA+l_AB)*cos(alph) - l_BD*cos(alph + zeta) - l_CD*cos(alph + zeta + e));
    
    sol_eta = eval(sol_eta);
    sol_b = eval(sol_b);

    % if no solution was found, assign zeros to keypoints and end execution
    if isempty(sol_eta)
        keypoints = zeros(4,8);
        return
    end
    
    %% Compute forward kinematics out to end effector and input points for
    % both solutions
    
    %% First solution
    % Choose one solution
    b = sol_b(1);
    e = sol_eta(1);
    
    % Get keypoints of solved kinematic chain
    [rA1, rB1, rC1, rD1, rE1, rF1, rG1, rH1] = fun_FK(b, e, vars);

%     % Save coordinates of end effector
%     rG1 = rG;

    %% Second solution
    b = sol_b(2);
    e = sol_eta(2);
    
    % Get keypoints of solved kinematic chain
    [rA2, rB2, rC2, rD2, rE2, rF2, rG2, rH2] = fun_FK(b, e, vars);
    
%     % Print coordinates of end-effector
%     rG2 = rG;

    % Return end effector coordinates of both solutions in an array
%     EE = [rG1 rG2];

    % Return keypoints of both solutions in a multi-dim array
    keypoints = [rA1, rB1, rC1, rD1, rE1, rF1, rG1, rH1;
                 rA2, rB2, rC2, rD2, rE2, rF2, rG2, rH2];
    
%     function [rA, rB, rC, rD, rE, rF, rG, rH] = fun_FK(b, e, v)
%     
%         [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, l_EF, l_EG, l_FH, alph, zeta, gamma, delta, epsilon] = unpack(v);
%        
%         rA = l_OA*[sin(alph), -cos(alph)]';
%         rB = (l_OA+l_AB)*[sin(alph), -cos(alph)]';
%         rC = rA + l_AC*[sin(alph + b), -cos(alph + b)]';
%         rD = rB + l_BD*[sin(alph + zeta), -cos(alph + zeta)]';
%         rE = rC + l_CE*[sin(alph + zeta + e - gamma), -cos(alph + zeta + e - gamma)]';
%         rF = rD + l_DF*[sin(alph + zeta + e - gamma), -cos(alph + zeta + e - gamma)]';
%         rG = rE + l_EG*[sin(alph + zeta + e - gamma + pi/2 - delta), -cos(alph + zeta + e - gamma + pi/2 - delta)]';
%         rH = rF + l_FH*[sin(alph + zeta + e - gamma - pi/2 + epsilon), -cos(alph + zeta + e - gamma - pi/2 + epsilon)]';
%     
%     %     R = [rA, rB, rC, rD, rE, rF, rG, rH];
%     end

end