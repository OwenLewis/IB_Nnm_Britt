%
% compute the gradient of a scalar function on a 2D domain using finite-differences
%
function Gu = gradientFD(u,grid)


switch grid.bcx
    case 'per'
        Dx = centered_diff_per(grid.Nx)/grid.dx;
    case 'dir'
        Dx = centered_diff_dir(grid.Nx)/grid.dx;
    case 'nmn'
        Dx = centered_diff_nmn(grid.Nx)/grid.dx;
end

switch grid.bcy
    case 'per'
        Dy = centered_diff_per(grid.Ny)/grid.dy;        
    case 'dir'
        Dy = centered_diff_dir(grid.Ny)/grid.dy;
    case 'nmn'
        Dy = centered_diff_nmn(grid.Ny)/grid.dy;
end
        
        

Gu = zeros(grid.Nx,grid.Ny,2);
Gu(:,:,1) = Dx*u;
Gu(:,:,2) = u*Dy';


end

