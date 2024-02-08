function [nSGLSF]=IBNeumm_helm_FD_NEW(a,b,xg,yg,dx,X0,ds, Fs,unitnormal, const)

%performs the immersed bounary for Neumann, given F on boundary 

%Inputs:
%a,b constants for helmholtz operator: a u - b lapl u
%I HAVENT DEALT WITH a AND b YET

%xg, yg, N x N grid values
%dx grid spacing
%X0 immersed boundary points
%ds immersed boundary spacing
%radius of the circle
%Fs  Nib x d
d=length(Fs(1,:)); 
%unitnormal: unit normal vectors on boundary pointing outward (or inward  
        %depending on what you want to do) to use to find Qn
        %Nib x 2
%const, 1/2 for interior, -1/2 for exterior

%Outputs:


%Info on problem:
Nx=length(xg(:,1));
Ny=length(xg(1,:));
xmin=xg(1,1);
ymin=yg(1,1);
Lx=xg(end,1)+dx-xmin;
Ly=yg(1,end)+dx-ymin;
Nib=length(X0(:,1));


%Spread operator
sp_scale=ds/dx^2;
S=spreadmatrix_vc_vec(X0, dx, Nx,Ny, xmin, ymin);

[Dx, Dy]=dx2d(Nx,Ny);
Dx=1/(2*dx)*Dx;
Dy=1/(2*dx)*Dy;


nSGLSF=zeros(Nib,d); % Nx x Ny x d; 
for i=1:d
    %spread
    Fi=Fs(:,i);  % Nib x 1
    SFi=sp_scale*S*Fi;%NxNy x1
    %negative L inverse
    SFi=reshape(SFi, Nx,Ny);
    LSFi=helmholtz_solve_FD(SFi, a,b, Lx,Ly, dx,dx);

    %Gradient
    LSFi=reshape(LSFi, Nx*Ny,1);
    GLSFi=zeros(Nx*Ny, 2);
    GLSFi(:,1)=Dx*LSFi;
    GLSFi(:,2)=Dy*LSFi;
    %Interpolate
    SGLSFi=S'*GLSFi; %Nib x2
    %Dot with normal
    nSGLSFi=SGLSFi(:,1).*unitnormal(:,1)+SGLSFi(:,2).*unitnormal(:,2); %Nib x1
    %Add on constant times F
    nSGLSFi=nSGLSFi+const/b*Fi;
    nSGLSF(:,d)=nSGLSFi;

end







end

