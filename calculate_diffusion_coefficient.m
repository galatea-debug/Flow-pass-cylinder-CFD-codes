function coeff = calculate_diffusion_coefficient(dx,dy,Delta_x,Delta_y,Nx,Ny)

a_E = zeros(Ny+2,Nx+2);
a_W = zeros(Ny+2,Nx+2);
a_N = zeros(Ny+2,Nx+2);
a_S = zeros(Ny+2,Nx+2);
a_P = zeros(Ny+2,Nx+2);

for i = 2:Nx+1
    for j = 2:Ny+1
        % Discretization coefficients
        a_E(j,i) = 1/(Delta_x(j,i)*dx(j,i));
        a_W(j,i) = 1/(Delta_x(j,i)*dx(j,i-1));
        a_N(j,i) = 1/(Delta_y(j,i)*dy(j,i));
        a_S(j,i) = 1/(Delta_y(j,i)*dy(j-1,i));
        a_P(j,i) =  (a_E(j,i) + a_W(j,i) + a_N(j,i) + a_S(j,i));
    end
end

coeff = struct();
coeff.a_E = a_E;
coeff.a_W = a_W;
coeff.a_N = a_N;
coeff.a_S = a_S;
coeff.a_P = a_P;

end