function D1=spectral_diff(N,dx)
% spectral_diff
%   differentiation matrix for pseduospectral method
%   N points equally spaced by dx
%   assumes N is even?
%
% from Trefethen for a 0..2pi grid
%  column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
%  D = toeplitz(column,column([1 N:-1:2]));
%  where h = 2*pi/N
%  we change to dx = L/N, so replace h -> 2*pi/L*dx = 2*pi/(dx*N)*dx = 2*pi/N
%
L=N*dx;
col = 2*pi/(L)*[0 0.5*(-1).^(1:N-1).*cot( 2*pi/N*(1:N-1)/2)];
D1  = toeplitz(col,col([1 N:-1:2]));





%
% some test code
%  N =32;
%  L = 2;
%  dx=L/N;
%  x=(0:N-1)'*dx;
%  y=exp( -5*(sin( pi/2*(x-0.75) )).^2)
%  yp = exp( -5*(sin( pi/2*(x-0.75) )).^2).*( -10*sin( pi/2*(x-0.75) )).*(cos( pi/2*(x-0.75))).*pi/2;

 