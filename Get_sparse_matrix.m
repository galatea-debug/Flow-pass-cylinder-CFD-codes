function A_Matrix = Get_sparse_matrix(coeff, color, Corelation, Nx, Ny)

N = (Nx+2)*(Ny+2);
A_Matrix = struct();

% -- fluid mask and numbering -------------------------------------------
fluid           = ~color.cylinder;
fluid([1 end],:)= false;      % strip halo
fluid(:,[1 end])= false;  %% fluid =0 in the exterior and in the cylinder


% -- full 5-point coefficients (identical to smoother) -------------------
b_E=zeros([Ny+2,Nx+2]);
b_W=zeros([Ny+2,Nx+2]);
b_N=zeros([Ny+2,Nx+2]);
b_S=zeros([Ny+2,Nx+2]);
b_P=ones([Ny+2,Nx+2]);

for i=2:Nx+1
    for j=2:Ny+1
        b_E(j,i)=(1-color.East(j,i))*coeff.a_E(j,i);
        b_W(j,i)=(1-color.West(j,i))*coeff.a_W(j,i);
        b_N(j,i)=(1-color.North(j,i))*coeff.a_N(j,i);
        b_S(j,i)=(1-color.South(j,i))*coeff.a_S(j,i);
        b_P(j,i)=-(b_E(j,i)+b_W(j,i)+b_N(j,i)+b_S(j,i));
    end
end

% ---------------- pre-allocate triplet arrays ---------------------------
nnz_est = nnz(fluid)*5 + (N - nnz(fluid));   % safe upper bound
I = zeros(nnz_est,1);  J = I;  V = I;  idx = 0;

% helper that appends one non-zero
    function push(r,c,val)
        idx = idx + 1;
        I(idx) = r; J(idx) = c; V(idx) = val;
    end

% ------------------- fill the sparse matrix ----------------------------
for i = 1:Nx+2
    for j = 1:Ny+2
        iP = Corelation.uridx(j,i);

        if fluid(j,i)             % ---------- interior (fluid) cell -----
            diag = b_P(j,i);

            % East neighbour
            if fluid(j,i+1)
                push(iP, Corelation.uridx(j,i+1), b_E(j,i));
            else
                diag = diag + b_E(j,i);   % fold to diagonal
            end
            % West
            if fluid(j,i-1)
                push(iP, Corelation.uridx(j,i-1), b_W(j,i));
            else
                diag = diag + b_W(j,i);
            end
            % North
            if fluid(j+1,i)
                push(iP, Corelation.uridx(j+1,i), b_N(j,i));
            else
                diag = diag + b_N(j,i);
            end
            % South
            if fluid(j-1,i)
                push(iP, Corelation.uridx(j-1,i), b_S(j,i));
            else
                diag = diag + b_S(j,i);
            end

            push(iP,iP,diag);     % final diagonal value

        else                      % ---------- cylinder ----------
            push(iP,iP,1.0);      % enforce p = 0
        end
    end
end

% trim unused allocation and build sparse matrix
A = sparse(I(1:idx), J(1:idx), V(1:idx), N, N);

[L, U] = ilu(A);
A_Matrix.A = A;
A_Matrix.L = L;
A_Matrix.U = U;

end