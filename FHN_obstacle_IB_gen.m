% Reaction Diffusion using IBDL New Neumann method
%   reactions are FitzHugh-Nagumo
%
%
clear all;

addpath('./src/');

% computational domain parameters
%
xmin   = -1;          % bottom cornrer of the domain
ymin   = -1;
Ly     = 2;          % height of the domain
aspect = 1;          % aspect ratio
Lx     = aspect*Ly;  % length of th domain

Ny     = 512;        % number of mesh points in y-direction
Nx     = aspect*Ny;  % number of mesh points in x-direction
dy     = Ly/Ny;      % fluid mesh spacing
dx     = Lx/Nx;     


% IB parameters
%
xc      = 0.0;         % center of the IB object xc, yc
yc      = 0.0;
rad     = 0.75;        % resting radius of the circle
dsscale = 0.75;        %ratio of IB points to grid spacing
ds      = dsscale*dx;  % IB mesh spacing 


% Reaction-Diffusion parameters and functions
%
% Rv = v(a − v)(v − 1) − w
% Rw = ep*(v − k*w)
%
Dv = 2.5e-4;
Dv = 1e-4;
a  = 0.1;
k  = 1;
ep = 0.005;
Rv=@(v,w)(v.*(1-v).*(v-a) - w);
Rw=@(v,w)(ep*(v-k*w));

% initial data functions
%
%v0 =@(x,y)(    exp(-400*(x-0.3 ).^2));
%w0 =@(x,y)(0.1*exp(-400*(x-0.25).^2));
v0 = @(x,y)(exp(-400*(x+0.6).^2));
w0 = @(x,y)(0.1*exp(-400*(x+0.7).^2) +0.1*(y-0).*(y>0)) ;

% time stepping 
%
Tmax = 1000;
dt   = 0.5;
Nt   = round(Tmax/dt);

% solver parameters
%
solveparams.rstart = 10;
solveparams.tol    = 1e-6;
solveparams.maxiter = 1000;

% flag to mask the rhs
%
rhsMaskFlag = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian grid
% 
xg=dx*(0:Nx-1)+xmin;
yg=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);

% IB points for a circle
%
[X0, ds] = circle(xc,yc,rad,ds);
Nib=length(X0(:,1));
sp_scale = ds/dx^2;


% domain mask -- for a circle
%
rg=sqrt((xg-xc).^2+(yg-yc).^2);
chi = 1.0*( rg <= rad);

% create the mask for the right size if needed
%

if( rhsMaskFlag )
  rhsMask = chi;
else  
  rhsMask = ones(Nx,Ny);
end
  
% constants involved in the SC equation
%
const1=-1/2; % plus or minus 1/2 Q ( if const2=1, minus interior, plus exterior)
const2=-1;   % which normal direction i want to use; 1 for pointing out of 
             % circle; -1 for pointing into circle

% normals 
%
unitnormal=const2*1/rad*(X0-repmat([xc,yc],Nib,1));

% pack up info about the Eulerian grid in a single variable
%
grid.xmin = xmin;
grid.ymin = ymin; 
grid.Lx   = Lx;
grid.Ly   = Ly;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.dx   = dx;
grid.dy   = dy;
grid.chi  = chi;

% pack up info on the IB 
%
IB.Nib     = length(X0);
IB.normals = unitnormal;
IB.dsvec   = ds*ones(IB.Nib,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a constant related to time stepping -- FE/BE
%
a1=1/dt;

% initialize
%
v = v0(xg,yg);
w = w0(xg,yg);

% check my SC apply 
%
a = 1/dt;
b = Dv;

% THIS CHECKS OUT TO MATCH BRITTANY'S CODE
%
%Mnew = zeros(IB.Nib);
%for k=1:IB.Nib
%  F = zeros(IB.Nib,1);
%  F(k) = 1;
%  Mnew(:,k)=apply_IBNeumann_SC(F,X0,a,b,IB,grid);
%end
%return


% to mask or not to mask that is the question
%


gridproblem = 'HelmInv_FD_period'


for n=1:Nt  %time loop

    % store last time step
    %
    wold = w;
    vold = v;

    % update w
    %
    w = wold + dt*Rw(vold,wold).*rhsMask;
%    w = wold + dt*Rw(vold,wold);

    % update for v
    %
    rhs = vold/dt + Rv(vold,wold).*rhsMask;
%    rhs = vold/dt + Rv(vold,wold);

    [v,Fds] = IBSL_Nmn_Solve(rhs,X0,IB,a,b,grid,gridproblem,solveparams);
    
    
    
    % visualize
    %
    pcolor(xg,yg,v);
    shading flat
    colorbar;
    hold on;
    plot(X0(:,1),X0(:,2),'r','linewidth',3);
    hold off;
    title(sprintf('time = %f',n*dt),'fontsize',16);
    caxis([-0.2,1.0]);
        
    
    pause(0.01);
 %   pause;
    

end  %end time loop

