function [u_star,v_star, uv_iter] = Burgers_solver(u_cc,v_cc,U_fc,V_fc,u_cc_old, v_cc_old, U_fc_old,V_fc_old,Delta_x,Delta_y ...
        ,coeff,Nx,Ny,dt,Re,max_iter_burgers,tol_burgers,u_in,color)

    % [u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_in);
    
    % Calculate convection
    H_u_n = calculate_convection(U_fc,V_fc,u_cc,Delta_x,Delta_y,Nx,Ny,color);
    H_v_n = calculate_convection(U_fc,V_fc,v_cc,Delta_x,Delta_y,Nx,Ny,color);
    H_u_n_1 = calculate_convection(U_fc_old,V_fc_old,u_cc_old,Delta_x,Delta_y,Nx,Ny,color);
    H_v_n_1 = calculate_convection(U_fc_old,V_fc_old,v_cc_old,Delta_x,Delta_y,Nx,Ny,color);
    C_u = 1.5.*H_u_n-0.5.*H_u_n_1;
    C_v = 1.5.*H_v_n-0.5.*H_v_n_1;

    % Calculate diffusion
    [Diff_u_x,Diff_u_y] = calculate_diffusion(u_cc,coeff,Nx,Ny,color);
    [Diff_v_x,Diff_v_y] = calculate_diffusion(v_cc,coeff,Nx,Ny,color);

    % Assemble RHS for GS
    RHS_u = u_cc - dt.*C_u + dt/(2*Re).*(Diff_u_x + Diff_u_y);
    RHS_v = v_cc - dt.*C_v + dt/(2*Re).*(Diff_v_x + Diff_v_y);

    coeff.a_E = (dt/2/Re)*coeff.a_E;
    coeff.a_W = (dt/2/Re)*coeff.a_W;
    coeff.a_N = (dt/2/Re)*coeff.a_N;
    coeff.a_S = (dt/2/Re)*coeff.a_S;
    coeff.a_P = 1 + (dt/2/Re)*coeff.a_P;

    % GS solver for burgers eq
    for k = 1:max_iter_burgers
        u_cc = GS_solver_burgers(u_cc,RHS_u,coeff,Nx,Ny,color);
        v_cc = GS_solver_burgers(v_cc,RHS_v,coeff,Nx,Ny,color);
        [u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_in);
        % u_cc(color.cylinder) = 0;
        % v_cc(color.cylinder) = 0;
    
        residual_u = calculate_residual_burgers(u_cc,RHS_u,coeff,Nx,Ny,color);
        residual_v = calculate_residual_burgers(v_cc,RHS_v,coeff,Nx,Ny,color);
    
        RMS_residual_u = rms(residual_u,"all");
        RMS_residual_v = rms(residual_v,"all");
        % fprintf('Iteration number: %i --- RMS residual (u,v): (%e, %e)\n',k,RMS_residual_u,RMS_residual_v)
    
        if (RMS_residual_u < tol_burgers && RMS_residual_v < tol_burgers)
            break
        end
    end
    
    uv_iter=k;
    u_star = u_cc;
    v_star = v_cc;

end