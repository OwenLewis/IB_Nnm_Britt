function [u,p]=stokes_solve_MAC(f,Lx,Ly)
    
    % record the number of grid points in each direction
    %
    szf=size(f);
    nx=szf(1);  
    ny=szf(2);
    K=nx/ny;
    
    % compute the mesh spacing
    %
    dx = Lx/nx;
    dy = Ly/ny;
        
    % compute the wave numbers
    %
    freqx = freq(nx,Lx);
    freqy = freq(ny,Ly);
    [k1 k2]=ndgrid(freqx,freqy);

    % eigenvalues for dirvergence
    %  
    Dx = (exp(1.i*k1*dx) - 1.0)/dx;
    Dy = (exp(1.i*k2*dy) - 1.0)/dy;
    
    % eigenvalues for gradient
    %
    Gx = (1.0-exp(-1.i*k1*dx))/dx;
    Gy = (1.0-exp(-1.i*k2*dy))/dy;
    
    % eigenvalues for Laplacian
    %   set the zerio eigenvalue to 1.0 to avoid division by zero
    %
    Lx = 2.0/dx^2 * (cos(k1*dx)-1.0);
    Ly = 2.0/dy^2 * (cos(k2*dy)-1.0);
    L  = Lx + Ly;
    L(1,1) = 1.0;
    
    % take FFT of the force
    % 
    fhat = fft2(f);
    
    % solve for pressure
    %
    divf_hat = Dx.*fhat(:,:,1) + Dy.*fhat(:,:,2);
    
%    divf = real( ifft2(divf_hat) );
%    fprintf('%g\n',max(abs(divf),[],'all'));
    
    
    phat = divf_hat./L;    
    
    % solve for velocity
    %
    uhat = zeros(nx,ny,2);
    uhat(:,:,1)=-(fhat(:,:,1)-Gx.*phat)./L;
    uhat(:,:,2)=-(fhat(:,:,2)-Gy.*phat)./L;

    % keepin' it real
    %
    p = real( ifft2(phat) );
    u = real( ifft2(uhat) );
