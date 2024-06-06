function robin=apply_IBRobin_SS(F,X,a,b,IB,grid,lambda)

%performs the immersed bounary for Neumann, given F on boundary 

% spread operator
%
S = spreadmatrix_cc_vec(X,grid);

% spread the force
%
Fds = IB.dsvec .* F;
SF = S*Fds/grid.dx^2;
SF = reshape(SF,grid.Nx,grid.Ny);

% apply n.(interp)*(Gradient)*(inv Helmholtz) to SF
%   abbreviated nSGA*(SF) 
%
nSGASF=apply_nSGA(SF,X,a,b,IB,grid);
SASF=apply_SA(SF,X,A,B,IB,grid);


% add the constant
%  --need to check the sign on this constant
%
robin = -nSGASF - lambda*SASF - 0.5/b*F;




