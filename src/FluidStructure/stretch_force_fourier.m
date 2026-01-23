function [F,strain] = stretch_force_fourier(X,ks,Lb,L0)

  % get the derivatives of the position
  %
  DX = derivative_fourier(X,Lb);

  % length of the tangent vector
  %
  L =sqrt(sum(DX.^2,2));

  % strain
  %
  strain = (L-L0)./L;

  %  force is derivative of d/ds( k*(L-L0)*DX/L )
  %
  diff_this = ks*strain(:,[1 1]).*DX;
  F = derivative_fourier(diff_this,Lb);
