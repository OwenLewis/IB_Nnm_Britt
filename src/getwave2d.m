function [kx,ky]=getwave2d(grid)

% get 1D wave numbers firx
%
freqx = freq(grid.Nx,grid.Lx);
freqy = freq(grid.Ny,grid.Ly);

% make arrays of the two wavenumbers
%
[kx,ky]=ndgrid(freqx,freqy );
