function nSGA=apply_nSGA(g,X,a,b,IB,grid,gridfunction);

% apply the operator n.S'*Grad*A
%  where A is the inverse of the Helmholtz operator H = a - b*L, and L is the discrete 
%  Laplacian


% apply the inverse of the Helmholtz operator to g
%
gridproblem = str2func(gridfunction);
A = gridproblem(g,a,b,grid);

% compute the gradient
%
GA = gradientFD(A,grid);

% form spread operator
%
S = spreadmatrix_vc_vec(X,grid);

% interpolate the gradient to the boundary
%
GA  = reshape(GA,grid.Nx*grid.Ny,2);
SGA = S'*GA;

% dot with the normal
%
nSGA=SGA(:,1).*IB.normals(:,1)+SGA(:,2).*IB.normals(:,2); 


