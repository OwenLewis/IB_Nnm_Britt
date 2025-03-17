function nSGA=apply_nSGA(g,X,a,b,IB,grid)

% apply the operator n.S'*Grad*A
%  where A is the inverse of the Helmholtz operator H = a - b*L, and L is the discrete 
%  Laplacian


% apply the inverse of the Helmholtz operator to g
%
A = helmsolve(g,a,b,grid);

% compute the gradient
%
GA = gradientFD(A,grid);

% form spread operator
%
switch grid.deltaflag
    case 0
        S = spreadmatrix_cc_vec(X,grid);
    case 1
        S = spreadmatrix_scalar6ptBspline(X,grid);
    otherwise
        S = spreadmatrix_cc_vec(X,grid);
end

% interpolate the gradient to the boundary
%
GA  = reshape(GA,grid.Nx*grid.Ny,2);
SGA = S'*GA;

% dot with the normal
%
nSGA=SGA(:,1).*IB.normals(:,1)+SGA(:,2).*IB.normals(:,2); 


