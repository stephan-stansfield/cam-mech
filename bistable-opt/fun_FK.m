function [rA, rB, rC, rD, rE, rF, rG, rH] = fun_FK(beta, eta, v)
    
        [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, l_EF, l_EG, l_FH, alph, zeta, gamma, delta, epsilon] = unpack_vars(v);
       
        rA = l_OA*[sin(alph), -cos(alph)]';
        rB = (l_OA+l_AB)*[sin(alph), -cos(alph)]';
        rC = rA + l_AC*[sin(alph + beta), -cos(alph + beta)]';
        rD = rB + l_BD*[sin(alph + zeta), -cos(alph + zeta)]';
        rE = rC + l_CE*[sin(alph + zeta + eta - gamma), -cos(alph + zeta + eta - gamma)]';
        rF = rD + l_DF*[sin(alph + zeta + eta - gamma), -cos(alph + zeta + eta - gamma)]';
        rG = rE + l_EG*[sin(alph + zeta + eta - gamma + pi/2 - delta), -cos(alph + zeta + eta - gamma + pi/2 - delta)]';
        rH = rF + l_FH*[sin(alph + zeta + eta - gamma - pi/2 + epsilon), -cos(alph + zeta + eta - gamma - pi/2 + epsilon)]';
    
    %     R = [rA, rB, rC, rD, rE, rF, rG, rH];
end