% use 1D FFT to take the derivative
%  function is 1D L-periodic
%
function Df = derivative_fourier(f,L);

% record the size
%
[n1, n2]=size(f);

% get the wavenumbers
%
k = freq0(n1,L);

% create copies for second dimension
%
k = 1.i * repmat(k,1,n2);

% to Fourier space
%
fhat = fft(f);

% derivative in Fourier space
%
Dfhat = k.*fhat;

% back to real space
%
Df = real( ifft(Dfhat) );


