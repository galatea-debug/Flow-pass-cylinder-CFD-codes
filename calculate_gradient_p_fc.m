function [Grad_p_fc_x,Grad_p_fc_y] = calculate_gradient_p_fc(p_cc,dx,dy,Nx,Ny,color)

    Grad_p_fc_x = zeros(Ny+1,Nx+1);
    Grad_p_fc_y = zeros(Ny+1,Nx+1);

    for i = 1:Nx+1
        for j = 1:Ny+1
            p_E  = (1 - color.East (j,i))*p_cc(j,i+1) + color.East (j,i)*p_cc(j,i);
            p_Px = (1 - color.West (j,i))*p_cc(j,i  ) + color.West (j,i)*p_cc(j,i+1);

            p_N  = (1 - color.North(j,i))*p_cc(j+1,i) + color.North(j,i)*p_cc(j,i);
            p_Py = (1 - color.South(j+1,i))*p_cc(j,i) + color.South(j+1,i)*p_cc(j+1,i);

            Grad_p_fc_x(j,i) = (p_E - p_Px)/dx(j,i);
            Grad_p_fc_y(j,i) = (p_N - p_Py)/dy(j,i);
        end
    end

end