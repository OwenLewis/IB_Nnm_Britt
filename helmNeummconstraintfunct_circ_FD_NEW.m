function [ nSGLSF ] = helmNeummconstraintfunct_circ_FD_NEW( F, xmin,...
    ymin, Ly, aspect, Ny,ds, X, a, b, const1, const2, xc,yc,rad)
 
%This is specifically for a circular boundary with center xc,yc and radius
%rad.
%Hu =   au- b(u_{xx} + u_{xx}) 
%I HAVE NOT DEALT WITH a and b YET

%The "input" Ub is,   Nib x 1 
%(Scalar u)

%const determines if you are solving for the u on the interior problem
%const1; %plus or minus 1/2 Q ( if const2=1, plus interior, minus exterior)
%const2; %which normal direction i want to use; 1 for pointing out of 
    %circle; -1 for pointing into circle
  

Lx  = aspect*Ly;      % length of th domain
Nx = aspect*Ny;   % number of mesh points in x-direction
dx = Ly/Ny;       % fluid mesh spacing

xg = dx*(0:Nx-1)+xmin;   
yg=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);

Nib=length(F);


%find normal vectors (assuming X is a circle
unitnormal=const2*1/rad*X;

 
%%%%%%%%%%%%%%%%%%%%%
[nSGLSF]=IBNeumm_helm_FD_NEW(a,b,xg,yg,dx,X,ds,F,unitnormal,const1);

        
% nSGLSF= nSGLSF(:,1);
%Now Ub is 2Nib x 1
   







end

