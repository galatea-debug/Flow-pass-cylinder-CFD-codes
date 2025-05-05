function p_cc = GS_solver_PPE2(p_cc,RHS_p,coeff,color,Nx,Ny)

    w = 1.99;

    % Coefficient

    for i = 2:Nx+1
        for j = 2:Ny+1
            phit = coeff.b_Pinv(j,i)*(coeff.b_E(j,i)*p_cc(j,i+1)+coeff.b_W(j,i)*p_cc(j,i-1)+coeff.b_N(j,i)*p_cc(j+1,i)+coeff.b_S(j,i)*p_cc(j-1,i)-RHS_p(j,i));
            p_cc(j,i) = (1-color.cylinder(j,i))*((1 - w)*p_cc(j,i) + w*phit);
        end
    end
    
end