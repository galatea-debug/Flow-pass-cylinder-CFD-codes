function [Grad_p_cc_x,Grad_p_cc_y] = calculate_gradient_p_cc(p_cc,Delta_x,Delta_y,Nx,Ny,color)

Grad_p_cc_x = zeros([Ny+2,Nx+2]);
Grad_p_cc_y = zeros([Ny+2,Nx+2]);
[p_fc_x,p_fc_y] = Get_face_variable(p_cc,Delta_x,Delta_y,Nx,Ny,color,"N");

for i = 2:Nx+1
    for j = 2:Ny+1
        Grad_p_cc_x(j,i) = (1-color.cylinder(j,i))*(p_fc_x(j,i)-p_fc_x(j,i-1))/(Delta_x(j,i));
        Grad_p_cc_y(j,i) = (1-color.cylinder(j,i))*(p_fc_y(j,i)-p_fc_y(j-1,i))/(Delta_y(j,i));   
    end
end

end