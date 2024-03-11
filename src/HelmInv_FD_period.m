%
%  Apply the inverse of the operator
%    H = a - b(u_xx + u_yy) 
%  by solving the equation
%    Hu = f
%  in a periodic box of size (Lx,Ly),
%  using finite differences solved with FFT
%
function u=HelmInv_FD_period(f,a,b,grid)

% get the wavenumbers
%
[kx,ky]=getwave2d(grid);

% compute eigenvalues of the operator
%
kk = a-b/grid.dx^2*(2*cos(kx*grid.dx) + 2*cos(ky*grid.dy)-4) ;
        
% if b=0, the problem is singular, adjust this eigenvalue
%
if( a==0 )
  kk(1) = 1;
end

% transform f, and project off mean if problem is singular
%
fhat = fft2(f);
if( a==0 )
  fhat(1) = 0;
end

% do the solve
%
uhat = fhat./kk;           
u    = real( ifft2(uhat) );
    
