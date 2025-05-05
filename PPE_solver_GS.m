function [p_cc,Iter_count,RMS_residual] = PPE_solver_GS(p_cc,RHS_p,coeff,color,Nx,Ny,max_iter_PPE,tol_PPE,dx)

Iter_count = zeros(max_iter_PPE,1);
RMS_residual = zeros(max_iter_PPE,1);

% Coefficient
coeff.b_E = (1-color.East).*coeff.a_E;
coeff.b_W = (1-color.West).*coeff.a_W;
coeff.b_N = (1-color.North).*coeff.a_N;
coeff.b_S = (1-color.South).*coeff.a_S;
coeff.b_P = (-coeff.b_E -coeff.b_W -coeff.b_S -coeff.b_N).*(1-color.cylinder)+color.cylinder;
coeff.b_Pinv = -1./coeff.b_P;
max_iter_inner = 100;

for k = 1:max_iter_PPE   
    for k1 = 1:max_iter_inner
        % GS solver
        p_cc = GS_solver_PPE2(p_cc,RHS_p,coeff,color,Nx,Ny);
        p_cc = p_cc - mean(p_cc.*(1-color.cylinder),"all");
        p_cc(1, :) =  p_cc(2,:); % p(0, y) boundary
        p_cc(end, :) = p_cc(end-1,:); % p(1, y) boundary
        p_cc(:, 1) =  p_cc(:,2); % p(x, 0) boundary
        p_cc(:, end) = p_cc(:,end-1); % p(x, 1) boundary
    end

    % Calculate residual
    Residual = calculate_residual_PPE(p_cc,RHS_p,coeff,color,Nx,Ny);

    Iter_count(k) = k;
    RMS_residual(k) = rms(Residual,"all");
    fprintf('Outer Iteration number: %i --- RMS residual: %e\n',k,RMS_residual(k))

    if(RMS_residual(k)<tol_PPE)
        break;
    end
end

Iter_count = Iter_count(1:k);
RMS_residual = RMS_residual(1:k);

end