function w = calculate_vorticity(u_cc,v_cc,Nx,Ny,Delta_x,Delta_y,color)
    
    w = zeros(size(u_cc));
    [~,u_y_fc] = Get_face_variable(u_cc,Delta_x,Delta_y,Nx,Ny,color,"D");
    [v_x_fc,~] = Get_face_variable(v_cc,Delta_x,Delta_y,Nx,Ny,color,"D");
    
    % Inner points: Central difference
    for i = 2:Nx+1
        for j = 2:Ny+1
            dv_i = v_x_fc(j,i) - v_x_fc(j,i-1);
            du_j = u_y_fc(j,i) - u_y_fc(j-1,i);
            w(j,i) = dv_i/Delta_x(j,i) - du_j/Delta_y(j,i);
        end
    end

end