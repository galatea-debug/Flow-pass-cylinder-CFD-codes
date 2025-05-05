function rhs_vector = build_PPE_RHS(RHS_p, Corelation, color, Nx, Ny)

N = (Nx+2)*(Ny+2);

fluid           = ~color.cylinder;
fluid([1 end],:)= false;      % strip halo
fluid(:,[1 end])= false;
rhs_vector=zeros([N,1]);

for i = 2:Nx+1
    for j = 2:Ny+1
        row = Corelation.uridx(j,i);
        rhs_vector(row)=fluid(j,i)*RHS_p(j,i); %% store zero outside the fluid, and rhs inside the fluid
    end
end

end
