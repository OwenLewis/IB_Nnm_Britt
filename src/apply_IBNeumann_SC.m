function nSGASF=apply_IBNeumann_SC(F,X,a,b,IB,grid)

%performs the immersed bounary for Neumann, given F on boundary 

% spread operator
%
switch grid.deltaflag
    case 0
        S = spreadmatrix_cc_vec(X,grid);
    case 1
        S = spreadmatrix_scalar6ptBspline(X,grid);
    otherwise
        S = spreadmatrix_cc_vec(X,grid);
end

% spread the force
%
Fds = IB.dsvec .* F;
SF = S*Fds/(grid.dx*grid.dy);
SF = reshape(SF,grid.Nx,grid.Ny);

% apply n.(interp)*(Gradient)*(inv Helmholtz) to SF
%   abbreviated nSGA*(SF) 
%
nSGASF=apply_nSGA(SF,X,a,b,IB,grid);


% add the constant
%  --need to check the sign on this constant
%
nSGASF = -nSGASF - 0.5/b*F;




