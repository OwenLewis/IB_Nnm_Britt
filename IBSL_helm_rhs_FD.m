function [u,U]=IBSL_helm_rhs_FD(a, b, xg, yg, dx, X0,ds, F, rhs)

%performs the immersed bounary SL, given F on boundary calculation for
%helmhomltz to get u and Ub

%Inputs:
%a,b constants for helmholtz operator: au-b lapl u
%xg, yg, N x N grid values
%dx grid spacing
%X0 immersed boundary points
%ds immersed boundary spacing
%F single layer distribution on the boundary  Nib x d 
d=length(F(1,:));% d number of dimensions of F (or number of F functions used)
    %(also the dimension of u (or number of u's))


    %rhs: for solving a laplu - bu = rhs;
    % is Nx x Ny x d. 
    
%Outputs:
%u on grid    Nx x Ny x d
%Ub on boundary  Nib x d


%Info on problem:
Nx=length(xg(:,1));
Ny=length(xg(1,:));
xmin=xg(1,1);
ymin=yg(1,1);
Lx=xg(end,1)+dx-xmin;
Ly=yg(1,end)+dx-ymin;


%Spread operator
sp_scale=ds/dx^2;
S=spreadmatrix_vc_vec(X0, dx, Nx,Ny, xmin, ymin);
f=sp_scale*S*F;  %f is NxNy x d
Totalnegrhs=reshape(f, Nx,Ny,d)-rhs;


   u=zeros(Nx,Ny, d);
   Nib=length(X0(:,1));
   U=zeros(Nib, d);
   for i =1:d
       %Apply - Linv
        u(:,:,i)=helmholtz_solve_FD(Totalnegrhs, a, b, Lx,Ly,dx,dx);

        U(:,i)=S'*reshape(u(:,:,i),Nx*Ny,1);
   end

        

end











