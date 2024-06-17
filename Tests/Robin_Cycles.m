%The purpose of this routine is to test the GMRES convergence for our
%method with ROBIN boundary conditions as a function of gamma = a2/a1. 
%Here, a2 and a1 are the ratios of the coefficients in front of the normal
%derivative and function (respectively) in the boundary condition. Gamma of
%zero is a dirichlet BC, while infinite gamma is a Neumann BC.

%Test case for new Robin BC routines. 
addpath('../src/');
pwer = 9;
a = 1;
b = 1e-3;

% Nx = 2^pwer;
Ny = 2^pwer;

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
xc      = 0.0;         % center of the IB object xc, yc
yc      = 0.0;
rad     = 1;        % resting radius of the circle
dsscale = 8/4;        %ratio of IB points to grid spacing
ds      = dsscale*dx;  % IB mesh spacing 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Cartesian grid
% 
xg=dx*(0.5:Nx-0.5)+xmin;
yg=dx*(0.5:Ny-0.5)+ymin;
[xg,yg]=ndgrid(xg,yg);

% IB points for a circle
%
[X0, ds] = circle(xc,yc,rad,ds);
Nib=length(X0(:,1));

% domain mask -- for a circle
%
rg=sqrt((xg-xc).^2+(yg-yc).^2);
chi = 1.0*( rg <= rad);




% true solution as well as right hand side. 
%
utrue=(xg.^2 - yg.^2).*chi;
rhs = a*(xg.^2 - yg.^2);

%Boundary data:
normderiv = (2*X0(:,1).^2 - 2*X0(:,2).^2)./rad;
value = X0(:,1).^2 - X0(:,2).^2;


% solver parameters
%
solveparams.rstart = 10;
solveparams.tol    = 1e-9;
solveparams.maxiter = 1000;




  
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
grid.bcx = 'per';
grid.bcy = 'per';



gammas = logspace(-3,3,20);
for i = 1:length(gammas)
    a1 = 1;
    a2 = gammas(i);
    

    
    
    
    
    
    Vb = a1*value + a2*normderiv;
    

    
    % pack up info on the IB 
    %
    IB = IB_populate(X0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [u,Fds,iter] = IBSL_Rbn_Solve(rhs,X0,IB,a,b,grid,solveparams,Vb,a1,a2);
    trueiters(i) = (iter(1)-1)*solveparams.rstart+iter(2)

    % figure(1)
    % surf(xg,yg,u.*chi);shading flat; colorbar
    % figure(2)
    % surf(xg,yg,utrue); shading flat; colorbar
    % 
    % error = u.*chi - utrue;
    % figure(3)
    % surf(xg,yg,error); shading flat; colorbar
end

% L2 = sqrt(sum(error(:).^2)*dx*dy)/(sqrt(sum(utrue(:).^2)*dx*dy))
% L1 = sum(abs(error(:)))*dx*dy/(sum(abs(utrue(:)))*dx*dy)
% Linf = max(abs(error(:)))/max(abs(utrue(:)))

rmpath('../src/');