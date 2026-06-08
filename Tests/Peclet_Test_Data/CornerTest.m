addpath('../../src/')
Nx = 512;
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

D1 = centered_diff_fourth(grid.Nx,grid.dx);

foo = v*D1*uwork;

% foo = upwind_corner(uwork,v,0,grid,dt);

% bar = v*(6*X.*(R-1)).*(X.^2 + Y.^2 < 1);
% 
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

udotgradc=@(c)(v*D1*c + v*c*D1');

for t = dt:dt:Tmax
    k1 = -udotgradc(uwork          );
    k2 = -udotgradc(uwork+0.5*dt*k1);
    k3 = -udotgradc(uwork+0.5*dt*k2);
    k4 = -udotgradc(uwork+    dt*k3); 
    uwork =uwork + dt*(k1+2*k2+2*k3+k4)/6;
    % surf(X,Y,uwork); shading flat
    % pause(0.0001)
end
error = abs(uwork - u_0);

surf(X,Y,error); shading flat
max(error(:))/max(u_0(:))


rmpath('../../src/')