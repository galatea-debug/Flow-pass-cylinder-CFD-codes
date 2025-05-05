function [Diff_phi_x,Diff_phi_y] = calculate_diffusion(phi,coeff,Nx,Ny,color)

    Diff_phi_x = zeros(Ny+2,Nx+2);
    Diff_phi_y = zeros(Ny+2,Nx+2);
    
    a_Px = coeff.a_P - coeff.a_N - coeff.a_S;
    a_Py = coeff.a_P - coeff.a_E - coeff.a_W;

    b_E = (1-color.East).*coeff.a_E;
    b_W = (1-color.West).*coeff.a_W;
    b_N = (1-color.North).*coeff.a_N;
    b_S = (1-color.South).*coeff.a_S;
    b_Px = -a_Px - (color.East.*coeff.a_E + color.West.*coeff.a_W);
    b_Py = -a_Py - (color.North.*coeff.a_N + color.South.*coeff.a_S);

    for i = 2:Nx+1
        for j = 2:Ny+1
            Diff_phi_x(j,i) = (1-color.cylinder(j,i))*(b_E(j,i)*phi(j,i+1) + b_W(j,i)*phi(j,i-1) + b_Px(j,i)*phi(j,i));
            Diff_phi_y(j,i) = (1-color.cylinder(j,i))*(b_N(j,i)*phi(j+1,i) + b_S(j,i)*phi(j-1,i) + b_Py(j,i)*phi(j,i));
        end
    end

end