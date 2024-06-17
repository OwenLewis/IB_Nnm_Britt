Ny = 56;
v = 0.01;
D = 0.1;
Tmax = 20;

addpath('../src/');

% computational domain parameters
%
xmin   = -1.5;          % bottom cornrer of the domain
ymin   = -1.5;
Ly     = 3;          % height of the domain
aspect = 1;          % aspect ratio
Lx     = aspect*Ly;  % length of th domain

Nx     = aspect*Ny;  % number of mesh points in x-direction
dy     = Ly/Ny;      % fluid mesh spacing
dx     = Lx/Nx;     


% IB parameters
%
xc      = 0;         % center of the IB object xc, yc
yc      = 0.0;
rad     = 1;        % resting radius of the circle
dsscale = 0.25;        %ratio of IB points to grid spacing
ds      = dsscale*dx;  % IB mesh spacing 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cell Centered grid
% 
xcc=dx*(1/2:Nx-1/2)+xmin;
ycc=dx*(1/2:Ny-1/2)+ymin;
[xg,yg]=ndgrid(xcc,ycc);
% Edge Centered grids in vertical & horizontal directions
% 
xce = dx*(0:Nx)+xmin;
yce = dx*(0:Ny)+ymin;
[xehoriz,yehoriz] = ndgrid(xce,ycc);
[xevert,yevert] = ndgrid(xcc,yce);
% IB points for a circle
%


[X0, ds] = circle(xc,yc,rad,ds);
Nib=length(X0(:,1));

planehoriz = xehoriz;
planevert = xevert;




grid.xmin = xmin;
grid.ymin = ymin; 
grid.Lx   = Lx;
grid.Ly   = Ly;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.dx   = dx;
grid.dy   = dy;
grid.chi  = chi;
grid.bcx = 'per';
grid.bcy = 'dir';

% pack up info on the IB 
%
IB = IB_populate(X0);

Ssides = spreadmatrix_csides_vec(X0,grid);
Stops = spreadmatrix_ctops_vec(X0,grid);

Uvec  = reshape(planehoriz,(grid.Nx+1)*grid.Ny,1);
fromsides = Ssides'*Uvec;
Vvec = reshape(planevert,grid.Nx*(grid.Ny+1),1);
fromtops = Stops'*Vvec



IB.normals = -IB.normals;


figure(6)
pcolor(xehoriz,yehoriz,planehoriz)
caxis([xmin xmin+Lx])
colorbar
shading flat
hold on 
plot(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,'or','LineWidth',2,'MarkerSize',3)
scatter3(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,fromsides,'k')
% quiver(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,IB.normals(:,1),IB.normals(:,2))

rmpath('../src/');