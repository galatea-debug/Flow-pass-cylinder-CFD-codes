function phi = GS_solver_burgers(phi,RHS,coeff,Nx,Ny,color)

    w = 1;

    phi_cylinder = 0;
    b_E = (1-color.East).*coeff.a_E;
    b_W = (1-color.West).*coeff.a_W;
    b_N = (1-color.North).*coeff.a_N;
    b_S = (1-color.South).*coeff.a_S;
    b_P = coeff.a_P + (color.East.*coeff.a_E + color.West.*coeff.a_W + color.North.*coeff.a_N + color.South.*coeff.a_S);
    b_Pinv = 1./b_P;
    RHS_m = RHS + 2*phi_cylinder.*(color.East.*coeff.a_E + color.West.*coeff.a_W + color.North.*coeff.a_N + color.South.*coeff.a_S);
    
    for i = 2:Nx+1
        for j = 2:Ny+1
            phit = b_Pinv(j,i)*(b_E(j,i)*phi(j,i+1) + b_W(j,i)*phi(j,i-1) + b_N(j,i)*phi(j+1,i) + b_S(j,i)*phi(j-1,i) + RHS_m(j,i));
            phi(j,i) = (1-color.cylinder(j,i))*((1 - w)*phi(j,i) + w*phit) + color.cylinder(j,i)*phi_cylinder;
        end
    end

end