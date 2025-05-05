function [phi_x_fc,phi_y_fc] = Get_face_variable(phi,Delta_x,Delta_y,Nx,Ny,color,BCs)

% Delta_x: Distance between edge point in x dir
% Delta_y: Distance between edge point in y dir
% dx: Distance between face center point (v_face) in x dir
% dy: Distance between face center point (v_face) in y dir

    phi_x_fc = zeros(Ny+1,Nx+1);
    phi_y_fc = zeros(Ny+1,Nx+1);

    switch BCs
        case "D" 
        signE = -1; signN = -1;
        case "N"          
        signE =  1; signN =  1;
    end


    for i = 1:Nx+1
        for j = 1:Ny+1
            % --- east & west ghosts (same for both BC types) ---
            phi_E  = (1-color.East(j,i))*phi(j,i+1)  + color.East(j,i)* signE * phi(j,i);
            phi_Px = (1-color.West(j,i+1))*phi(j,i)  + color.West(j,i+1)* signE * phi(j,i+1);

            % --- north & south ghosts ---
            phi_N  = (1-color.North(j,i))*phi(j+1,i) + color.North(j,i)* signN * phi(j,i);
            phi_Py = (1-color.South(j+1,i))*phi(j,i) + color.South(j+1,i)* signN * phi(j+1,i);

            phi_x_fc(j,i) = (Delta_x(j,i+1)*phi_Px + Delta_x(j,i)*phi_E) / (Delta_x(j,i+1)+Delta_x(j,i));
            phi_y_fc(j,i) = (Delta_y(j+1,i)*phi_Py + Delta_y(j,i)*phi_N) / (Delta_y(j+1,i)+Delta_y(j,i));
        end
    end

   
end