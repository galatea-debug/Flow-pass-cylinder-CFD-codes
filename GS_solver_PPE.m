function p_cc = GS_solver_PPE(p_cc,RHS_p,coeff,color,Nx,Ny)

    w = 1.97;
    p_cylinder = 0;

    % Coefficient
    b_E = (1-color.East).*coeff.a_E;
    b_W = (1-color.West).*coeff.a_W;
    b_N = (1-color.North).*coeff.a_N;
    b_S = (1-color.South).*coeff.a_S;
    b_P = -coeff.a_P + color.North.*coeff.a_N + color.South.*coeff.a_S + color.East.*coeff.a_E + color.West.*coeff.a_W;
    b_Pinv = -1./b_P;

    for i = 2:Nx+1
        for j = 2:Ny+1
            p_W = color.West(j,i)*p_cc(j,i) + (1-color.West(j,i))*p_cc(j,i-1);
            p_E = color.East(j,i)*p_cc(j,i) + (1-color.East(j,i))*p_cc(j,i+1);
            p_N = color.North(j,i)*p_cc(j,i) + (1-color.North(j,i))*p_cc(j+1,i);
            p_S = color.South(j,i)*p_cc(j,i) + (1-color.South(j,i))*p_cc(j-1,i);
            phit = b_Pinv(j,i)*(b_E(j,i)*p_E + b_W(j,i)*p_W + b_N(j,i)*p_N + b_S(j,i)*p_S - RHS_p(j,i));
            % phit = b_Pinv(j,i)*(b_E(j,i)*(color.West(j,i)*p_cc(j,i) + (1-color.West(j,i))*p_cc(j,i-1)) + ...
            %                     b_W(j,i)*(color.West(j,i)*p_cc(j,i) + (1-color.West(j,i))*p_cc(j,i-1)) + ...
            %                     b_N(j,i)*(color.North(j,i)*p_cc(j,i) + (1-color.North(j,i))*p_cc(j+1,i)) + ...
            %                     b_S(j,i)*(color.South(j,i)*p_cc(j,i) + (1-color.South(j,i))*p_cc(j-1,i)) - RHS_p(j,i));
            p_cc(j,i) = (1-color.cylinder(j,i))*((1 - w)*p_cc(j,i) + w*phit) + color.cylinder(j,i)*p_cylinder;
        end
    end
    
end