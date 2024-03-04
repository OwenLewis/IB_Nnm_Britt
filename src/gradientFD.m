%
% compute the gradient of a scalar function on a 2D domain using finite-differences
%
function Gu = gradient(u,grid)

Dx = centered_diff_per(grid.Nx)/grid.dx;
Dy = centered_diff_per(grid.Ny)/grid.dy;

Gu = zeros(grid.Nx,grid.Ny,2);
Gu(:,:,1) = Dx*u;
Gu(:,:,2) = u*Dy';


