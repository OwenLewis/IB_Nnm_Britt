function nSGASF=apply_IBNeumann_SC(F,X,a,b,IB,grid,gridproblem);

%performs the immersed bounary for Neumann, given F on boundary 

% spread operator
%
S = spreadmatrix_vc_vec(X,grid);

% spread the force
%
Fds = IB.dsvec .* F;
SF = S*Fds/grid.dx^2;
SF = reshape(SF,grid.Nx,grid.Ny);

% apply n.(interp)*(Gradient)*(inv Helmholtz) to SF
%   abbreviated nSGA*(SF) 
%
nSGASF=apply_nSGA(SF,X,a,b,IB,grid,gridproblem);


% add the constant
%  --need to check the sign on this constant
%
nSGASF = nSGASF - 0.5/b*F;




