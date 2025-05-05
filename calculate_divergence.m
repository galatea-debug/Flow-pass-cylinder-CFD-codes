function Div_U = calculate_divergence(U,V,Delta_x,Delta_y,Nx,Ny,color)

    Div_U = zeros(Ny+2,Nx+2);
    
    for i = 2:Nx+1
        for j = 2:Ny+1
            Div_U(j,i) = (1-color.cylinder(j,i))*(U(j,i) - U(j,i-1))/Delta_x(j,i) + ...
                                                 (V(j,i) - V(j-1,i))/Delta_y(j,i);
        end
    end

end