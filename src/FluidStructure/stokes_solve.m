function [u,p]=stokes_solve(f,Lx,Ly)
    
    % record the number of grid points in each direction
    %
    szf=size(f);
    nx=szf(1);  
    ny=szf(2);
    K=nx/ny;
    
    % compute the wave numbers
    %
    freqx = freq(nx,Lx);
    freqy = freq(ny,Ly);
    [k1 k2]=ndgrid(freqx,freqy);

    % compute k.k and set the zero mode to 1 to avoid division by zero
    %
    ksq=k1.^2+k2.^2;
    ksq(1,1)=1;

    
    % take FFT of the force
    % 
    fhat = fft2(f);
    
    % solve for pressure
    %
    i=sqrt(-1);
    divf_hat=i*k1.*fhat(:,:,1)+i*k2.*fhat(:,:,2);
    phat=(-1./ksq).*(divf_hat);               

    
    % solve for velocity
    %
    uhat = zeros(nx,ny,2);
    uhat(:,:,1)=(fhat(:,:,1)-i*k1.*phat)./ksq;
    uhat(:,:,2)=(fhat(:,:,2)-i*k2.*phat)./ksq;

    p = real( ifft2(phat) );
    u = real( ifft2(uhat) );