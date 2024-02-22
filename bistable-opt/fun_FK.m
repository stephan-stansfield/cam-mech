function [rA, rB, rC, rD, rE, rF] = fun_FK(beta, eta, v)
    
        [l_OA, l_AB, l_AC, l_BD, l_CD, l_CE, l_DF, alph, gamma, delta, epsilon] = unpack_vars(v);

        rA = l_OA*[sin(alph), -cos(alph)]';
        rB = (l_OA+l_AB)*[sin(alph), -cos(alph)]';
        rC = rA + l_AC*[sin(alph + beta), -cos(alph + beta)]';
        rD = rB + l_BD*[sin(alph + gamma), -cos(alph + gamma)]';
        rE = rC + l_CE*[sin(alph + gamma + eta - delta), -cos(alph + gamma + eta - delta)]';
        rF = rD + l_DF*[sin(alph + gamma + eta + pi + epsilon), -cos(alph + gamma + eta + pi + epsilon)]';
    
    %     R = [rA, rB, rC, rD, rE, rF, rG, rH];
end