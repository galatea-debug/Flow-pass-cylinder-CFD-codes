function p_cc = PPE_solver_Matlab_solver(p_cc,RHS_p,A_Matrix,Corelation,color,Nx,Ny,max_iter_PPE,tol_PPE,dx)

    % Matlab solver
    b = build_PPE_RHS(RHS_p, Corelation, color, Nx, Ny);
    p = mldivide(A_Matrix.A,b);
    % [p,fl1,rr1,it1,rv1] = bicgstab(A_Matrix.A, b, tol_PPE, max_iter_PPE, A_Matrix.L, A_Matrix.U);
    % [p,fl1,rr1,it1,rv1] = gmres(A_Matrix.A,b,[],tol_PPE,max_iter_PPE,A_Matrix.L, A_Matrix.U);

    idx = 1;
    for i = 1:Nx+2
        for j = 1:Ny+2
            I = Corelation.ri(idx);
            J = Corelation.rj(idx);
            p_cc(J,I) = p(idx);
            idx = idx + 1;
        end
    end

    p_cc(1, :) = p_cc(2,:); % p(0, y) boundary
    p_cc(end, :) = p_cc(end-1,:); % p(1, y) boundary
    p_cc(:, 1) =  p_cc(:,2); % p(x, 0) boundary
    p_cc(:, end) = p_cc(:,end-1); % p(x, 1) boundary
    
    p_cc = p_cc - mean(p_cc.*(1-color.cylinder),"all");
    % p_cc(color.cylinder) = 0;

end