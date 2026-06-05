addpath('../../src/')
Nx = 1028;
Ny = Nx;

L = 3;

grid.dx = L/Nx;
grid.dy = L/Ny;
grid.Nx = Nx;
grid.Ny = Ny;


cfl = 0.2;
v = 0.1

dt = cfl*grid.dx/v




cfl = v*dt/grid.dx

Tmax = 3/v

x = -L/2 + grid.dx/2:grid.dx:L/2-grid.dx/2;
y = -L/2 + grid.dy/2:grid.dy:L/2-grid.dy/2;

[X,Y] = ndgrid(x,y);

R = sqrt(X.^2 + Y.^2);

u_0 = (-3*R.^2 + 2*R.^3 + 1).*(X.^2 + Y.^2 < 1);

uwork = u_0;
surf(X,Y,u_0); shading flat

foo = upwind_corner(uwork,v,0,grid,dt);

bar = v*(6*X.*(R-1)).*(X.^2 + Y.^2 < 1);

% figure(1)
% surf(X,Y,foo); shading flat
% title('Computed')
% figure(2)
% surf(X,Y,bar); shading flat
% title('True')
% 
% discerr = abs(foo - bar);
% figure(3)
% surf(X,Y,discerr); shading flat
% title('Error')
% max(discerr(:))

for t = dt:dt:Tmax
    foo = upwind_corner(uwork,v,v,grid,dt);
    uwork =uwork - dt*foo;
    % surf(X,Y,uwork); shading flat
    % pause(0.0001)
end
error = abs(uwork - u_0);

surf(X,Y,error); shading flat
max(error(:))/max(u_0(:))


rmpath('../../src/')