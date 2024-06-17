% clear all
close all
addpath('../src/');

flag = 1

grid.xmin = 0;
grid.ymin = 0;
grid.Lx = 2;
grid.Ly = 3;
grid.Nx = 100;
grid.Ny = 100;
grid.dx = grid.Lx/grid.Nx;
grid.dy = grid.Ly/grid.Ny;
xarray = grid.dx/2:grid.dx:grid.Lx-grid.dx/2;
yarray = grid.dy/2:grid.dy:grid.Ly-grid.dy/2;
[X,Y]=ndgrid(xarray,yarray);


%%This is my double dirichlet boundary condition:


switch flag
    case 1
        u = X.*(X-grid.Lx).*Y.*(Y-grid.Ly);
        
        f = 2*Y.*(Y-grid.Ly) + 2*X.*(X-grid.Lx);
        
        grid.bcx = 'dir';
        grid.bcy = 'dir'
        Uapprox = helmsolve(f,0,-1,grid);
        
        figure(1)
        surf(X,Y,u)
        shading flat
        figure(2)
        surf(X,Y,Uapprox)
        shading flat
        
        error = abs(u - Uapprox);
        max(max(error))

    case 2
        u = (X.^2).*(X-grid.Lx).^2.*Y.^2.*(Y-grid.Ly).^2-grid.Lx^4*grid.Ly^4/900;
        
        f = 2.*((X-grid.Lx).^2).*(Y.^2).*((Y-grid.Ly).^2) + 8*X.*(X-grid.Lx).*(Y.^2).*((Y-grid.Ly).^2) + 2*X.^2.*Y.^2.*(Y-grid.Ly).^2 + 2*X.^2.*(X-grid.Lx).^2.*(Y-grid.Ly).^2 + 8*X.^2.*(X-grid.Lx).^2.*Y.*(Y-grid.Ly) + 2*X.^2.*(X-grid.Lx).^2.*Y.^2;
        
        grid.bcx = 'nmn';
        grid.bcy = 'nmn'
        Uapprox = helmsolve(f,0,-1,grid);
        
        figure(1)
        surf(X,Y,u)
        shading flat
        figure(2)
        surf(X,Y,Uapprox)
        shading flat
        
        error = abs(u - Uapprox);
        max(max(error))

    case 3
        u = (X.^2).*(X-grid.Lx).^2.*Y.^2.*(Y-grid.Ly).^2-grid.Lx^4*grid.Ly^4/900;
        
        f = 2.*((X-grid.Lx).^2).*(Y.^2).*((Y-grid.Ly).^2) + 8*X.*(X-grid.Lx).*(Y.^2).*((Y-grid.Ly).^2) + 2*X.^2.*Y.^2.*(Y-grid.Ly).^2 + 2*X.^2.*(X-grid.Lx).^2.*(Y-grid.Ly).^2 + 8*X.^2.*(X-grid.Lx).^2.*Y.*(Y-grid.Ly) + 2*X.^2.*(X-grid.Lx).^2.*Y.^2;
        
        grid.bcx = 'per';
        grid.bcy = 'per'
        Uapprox = helmsolve(f,0,-1,grid);
        
        figure(1)
        surf(X,Y,u)
        shading flat
        figure(2)
        surf(X,Y,Uapprox)
        shading flat
        
        error = abs(u - Uapprox);
        max(max(error))

end



rmpath('../src/');

