% Diffusion using "IBDL" New Neumann method
%Interior of a circle 
%du dn = 0 on circle

clear all;

xmin        = -1.5;            
ymin        = -1.5;
Ly          = 3;            % height of the domain
aspect      = 1;             % aspect ratio
Lx          = aspect*Ly;     % length of th domain
xc          = 0;             % center of the IB object xc, yc
yc          = 0;
rad         = 1;           % restding radius of the circle

dsscale=0.75;  %ratio of IB points to grid spacing

 

D=1e-5;%2.5e-4;
Tmax=6;

dt=1e-2;

Ny=2.^7;  % number of mesh points in y-direction
Nx = aspect*Ny;   % number of mesh points in x-direction
dy=Ly/Ny;      % fluid mesh spacing
dx=Lx/Nx;
ds = dsscale*dx;     % IB mesh spacing 

const1=-1/2; %plus or minus 1/2 Q ( if const2=1, minus interior, plus exterior)
const2=1; %which normal direction i want to use; 1 for pointing out of 
    %circle; -1 for pointing into circle


ts=(0:dt:Tmax)';
Nt=length(ts);
usolutions=cell(1,Nt); %for storing u solutions


a1=1/dt;

% grid points WHY ARE WE USING AN EDGE CENTERED GRID?
xg=dx*(0:Nx-1)+xmin;
yg=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);
% [thetag,rg]=cart2pol(xg, yg);
% 
% rg=sqrt((xg-xc).^2+(yg-yc).^2);
% chi = 1.0*( rg < rad);

% IB points 
[X0, ds] = circle(xc,yc,rad,ds);
unitnormal=const2*1/rad*X0;
Nib=length(X0(:,1));
sp_scale = ds/dx^2;

X0 = X0;
x0=X0(:,1);
y0=X0(:,2);


%Define the solution values on the boundary
Vb=zeros(length(X0(:,1)),1);
xpad = zeros(Nx+1,Ny);
ypad = zeros(Nx,Ny+1);
u_0=(sqrt((xg).^2 + (yg).^2)<0.95).*exp(-((xg+0.25).^2+(yg-0.2).^2)./(0.3^2));
u_0 = u_0*3;
% u_0 = ones(size(xg));

u_old=u_0;
usolutions{1} = u_old;
    for n=2:Nt  %time loop
        X0 = X0 + dt;
        S=spreadmatrix_vc_vec(X0, dx, Nx,Ny, xmin, ymin);
        %Solve (I/dt-DaL)uAnew+SFanew = ua_old/dt +Ra_old
        %      (I/dt-DbL)ubnew+SFbnew=ub_old/dt+Rb_old
        %    ( S*grad uanew)dot n -1/(2Da) Fa_new = 0
        %    ( S*grad ubnew)dot n -1/(2Db) Fb_new = 0 

        rhs=zeros(Nx,Ny,1);
        xpad(1,:) = u_old(end,:);
        xpad(2:end,:) = u_old;
        ypad(:,1) = u_old(:,end);
        ypad(:,2:end) = u_old;
        advec = diff(xpad)/dx + diff(ypad,1,2)/dy;
        rhs= u_old/dt - advec;

        %Dealing with the rhs for gmres: 
        %Do negative L inverse:
        Lg=helmholtz_solve_FD(rhs, a1, D, Lx,Ly, dx,dy);% - Linv g  Nx x Ny x 2
        %Do gradient
        [Dx, Dy]=dx2d(Nx,Ny);
        Dx=1/(2*dx)*Dx;
        Dy=1/(2*dx)*Dy;
        GLg1=Dx*reshape(Lg, Nx*Ny, 1); %NxNy x 1
        GLg2=Dy*reshape(Lg, Nx*Ny, 1); %NxNy x1
        %Interpolate onto boundary
        SGLg1=S'*GLg1; %Nib x1
        SGLg2=S'*GLg2; %Nib x1
        %dot with normal
        nSGLg=zeros(Nib,2);
        nSGLg(:,1)=unitnormal(:,1).*SGLg1+unitnormal(:,2).*SGLg2; %Nib x1


        F=gmres(@(Fa)helmNeummconstraintfunct_circ_FD_NEW(Fa,xmin,ymin, Ly, ...
        aspect, Ny,ds, X0, a1, D, const1, const2, xc,yc,rad), ...
        Vb(:,1)+nSGLg(:,1),10, 1e-6, 1000);  %+SLg  

             %Find u everywhere
        [u_new,U]=IBSL_helm_rhs_FD(a1,D,xg,yg,dx, X0, ds, F(:,1), rhs);
        
        u_old=u_new;
        usolutions{n} = u_new;
        pcolor(xg,yg,u_old)
        caxis([0 1.5])
        colorbar
        shading flat
        hold on
        % plot3(X0(:,1),X0(:,2),zeros(size(X0(:,1))),'r','LineWidth',2)
        % quiver3(X0(:,1),X0(:,2),zeros(size(X0(:,1))),unitnormal(:,1),unitnormal(:,2),zeros(size(unitnormal(:,1))))
        plot(mod(X0(:,1)+1.5,3)-1.5,mod(X0(:,2)+1.5,3)-1.5,'or','LineWidth',2)
        % quiver(X0(:,1),X0(:,2),unitnormal(:,1),unitnormal(:,2))
        title(sprintf('time = %f',(n-1)*dt))
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        pause(0.01)
        hold off
    end


