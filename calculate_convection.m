function H = calculate_convection(U,V,phi,Delta_x,Delta_y,Nx,Ny,color)

    H = zeros(Ny+2,Nx+2);
    [phi_x_fc,phi_y_fc] = Get_face_variable(phi,Delta_x,Delta_y,Nx,Ny,color,"D");
   
    for i = 2:Nx+1
        for j = 2:Ny+1
            H(j,i) = (1-color.cylinder(j,i))*((U(j,i)*phi_x_fc(j,i) - U(j,i-1)*phi_x_fc(j,i-1))/Delta_x(j,i) + ...
                                             (V(j,i)*phi_y_fc(j,i) - V(j-1,i)*phi_y_fc(j-1,i))/Delta_y(j,i));
        end
    end

end