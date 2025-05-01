%
% Test the convergence of the IB method with a given normal derivative on the boundary
%
%   exact solution 
%    u = 1       for r < 1
%        1-ln(r) for r>=1
%
%   IB force to produce this solution is F = 1 on the boundary
%
%   solved in a [-2,2]^2 domain with Dirchlet boundary data
%
%   input 
%
function sol = Normal_and_density_solve(Nx,dsscale,deltaflag)

addpath('../src/');

% domain parameters
%
xmin = -2;
ymin = -2;
Lx = 4;
Ly = 4;

%Nx = 32*4;
Ny = Nx;
dx = Lx/Nx;
dy = Ly/Ny;

% IB parameters
%
xc = 0.0;
yc = 0.0;
rad = 1.0;
%dsscale = 0.25;
ds = dsscale * dx;

% Cartesian grid
%
xg=dx*(0:Nx-1)+xmin + dx/2;
yg=dx*(0:Ny-1)+ymin + dy/2;
[xg,yg]=ndgrid(xg,yg);

% IB points for a circle
%
[X0, ds] = circle(xc,yc,rad,ds);
Nib=length(X0(:,1));
sp_scale = ds/dx^2;
unitnormal=1.0/rad*(X0-repmat([xc,yc],Nib,1));


% domain mask -- for a circle
%
rg=sqrt((xg-xc).^2+(yg-yc).^2);
chi = 1.0*( rg >= rad);

% output the grids
%
sol.xg = xg;
sol.yg = yg;
sol.X0 = X0;
sol.s  = ds*(0:(Nib-1))';
sol.chi = 1-chi;
sol.Nv = -unitnormal;

    

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
grid.chi  = 1-chi;
grid.bcx = 'dir';
grid.bcy = 'dir';
grid.deltaflag = deltaflag;

% pack up info on the IB 
%
IB.Nib     = length(X0);
IB.normals = -unitnormal;
IB.dsvec   = ds*ones(IB.Nib,1);

% exact solution on the mesh
%
u_ex = 1 - log(rg + eps).*(chi);
ux_ex = -(chi).*xg./(rg+eps).^2;
uy_ex = -(chi).*yg./(rg+eps).^2;

% exact solution on the boundary 
%  derivatives are on the outside, inside derivatives are zero
%
U_ex  = ones(Nib,1);
Ux_ex = -cos(sol.s);
Uy_ex = -sin(sol.s);

% output the exact solution
%
sol.u_ex  = u_ex;
sol.ux_ex = ux_ex;
sol.uy_ex = uy_ex;
sol.U_ex  = U_ex;
sol.Ux_ex = Ux_ex;
sol.Uy_ex = Uy_ex;

% boundary data on the Dirichelet boundaries
%  comes from formula u(0) = 2*ub - u(1) 
%
%
g = zeros(Nx,Ny);

rsqlef = xmin^2 + (yg(1,:)-yc).^2;
ulef   =  1 - 0.5*log(rsqlef);
g(1,:) = g(1,:) + 2*ulef/(dx^2);

xmax = xmin + Lx;
rsqrgt = xmax^2 + (yg(1,:)-yc).^2;
urgt   =  1 - 0.5*log(rsqrgt);
g(Nx,:) = g(Nx,:) + 2*urgt/(dx^2);

rsqbot = (xg(:,1)-xc).^2 + ymin.^2;
ubot   =  1 - 0.5*log(rsqbot);
g(:,1) = g(:,1) + 2*ubot/dy^2;

ymax = ymin + Ly;
rsqtop = (xg(:,1)-xc).^2 + ymax.^2;
utop   =  1 - 0.5*log(rsqtop);
g(:,Ny) = g(:,Ny) + 2*utop/dy^2;

% IB force
%
% Fib = ones(Nib,1);

% spread materix
%
switch deltaflag
  case 0
    S = spreadmatrix_cc_vec(X0,grid);
  case 1
    S = spreadmatrix_scalar6ptBspline(X0,grid);
  otherwise
    S = spreadmatrix_cc_vec(X0,grid);  
end

% spread the force
%
a = 0;
b = 1;
Vb = ones(Nib,1);

solveparams.rstart = 10;
solveparams.tol    = 1e-10;
solveparams.maxiter = 1000;


[u,fds] = IBSL_Nmn_Solve(g,X0,IB,a,b,grid,solveparams,Vb);
% f = sp_scale*S*Fib;
% f = reshape(f,[Nx,Ny]);

% solve Poisson equation
%
% u = helmsolve(f+g,0,1,grid);


% compute the gradient
%  near the boundary extrapolate from the interior
%
Gu = gradientFD(u,grid);
ux = Gu(:,:,1);
uy = Gu(:,:,2);

% extrapolate from interior for near-bondary points
%
ux(1 ,:) = 2*ux(2   ,:   ) - ux(3   ,:   );
ux(Nx,:) = 2*ux(Nx-1,:   ) - ux(Nx-2,:   );
uy(:,1 ) = 2*uy(:   ,2   ) - uy(:   ,3   );
uy(:,Ny) = 2*uy(:   ,Ny-1) - uy(:   ,Ny-2);

% these use the boundary data
%
%ux(1 ,:) = (-4/3*ulef + u(1, :) + 1/3*u(2,   :) )/dx;
%ux(Nx,:) = ( 4/3*urgt - u(Nx,:) - 1/3*u(Nx-1,:) )/dx;
%uy(:, 1) = (-4/3*ubot + u(:, 1) + 1/3*u(:,2   ) )/dy;
%uy(:,Ny) = ( 4/3*utop - u(:,Ny) - 1/3*u(:,Ny-1) )/dy;

% interpolate the derivatives to the boundary
%
Ux = S'*ux(:);
Uy = S'*uy(:);


% output the solution
%
sol.F = fds./IB.dsvec;
sol.u = u;
sol.ux = ux;
sol.uy = uy;
sol.Ux = Ux;
sol.Uy = Uy;

sol.grid = grid;
sol.X0 = X0;
sol.IB = IB;

rmpath('../src/');

max(abs(fds))
max(max(abs(u)))
max(max(abs(Ux)))

mesh(sol.xg,sol.yg,sol.u); shading flat; colorbar
