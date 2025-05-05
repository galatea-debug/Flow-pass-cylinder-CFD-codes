function Corelation = Corelation_Matrix_array(Nx,Ny)

    N = (Nx+2)*(Ny+2);
    idx = 1;
    
    uridx = zeros(Ny+2,Nx+2);
    ri = zeros(N,1);
    rj = zeros(N,1);
    Corelation = struct();
    
    for i = 1:Nx+2
        for j = 1:Ny+2
            uridx(j,i) = idx; %%this is ij
            ri(idx) = i;
            rj(idx) = j;
            idx = idx + 1;
        end
    end
    
    Corelation.uridx = uridx;
    Corelation.ri = ri;
    Corelation.rj = rj;
    
end