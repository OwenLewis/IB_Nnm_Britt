function SA=apply_SA(g,X,a,b,IB,grid)

% apply the operator S'*A
%  where A is the inverse of the Helmholtz operator H = a - b*L, and L is the discrete 
%  Laplacian


% apply the inverse of the Helmholtz operator to g
%
A = helmsolve(g,a,b,grid);

% form spread operator
%
S = spreadmatrix_cc_vec(X,grid);
% S = spreadmatrix_scalar6ptBspline(X,grid);

% interpolate the result to the boundary
%
A = reshape(A,grid.Nx*grid.Ny,1);
SA = S'*A;
