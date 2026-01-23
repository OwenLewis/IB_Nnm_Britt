%
% given boundary positions on coarse, equally spaced (in parameter) grid
%   1. extend to a finer boundary mesh
%   2. return dX/ds -- i.e. derivatives with respect to parameter
%   3. ??? dl = |dX/ds| * ds -- not sure we need/want this
%  
%
function [XI,DXI] = extend_coarse_to_fine(X,ni,L);

% record the length
%
[n1, n2] = size(X);

% index of the highest positive wave number 
%
nplus = floor(n1/2) + 1;


% take fourier transform
%
Xhat = fft(X);

% if n1 is even remove fourier coefficient of highest wave number
%   we do have it's complex conjugate
%
if( mod(n1,2)==0 )
  Xhat(nplus)=0;  
end



% get the number of zeros to add 
%
ndiff = ni - n1;

% Fourier coefficients of the interpolant
%
XIhat = [Xhat(1:nplus,:); zeros(ndiff,n2); Xhat(nplus+1:end,:)];

% rescale because of the way matlab scale Fourier coefficients
%
XIhat = ni/n1*XIhat;

% back to real space
%
XI = real( ifft(XIhat) );

%
% now compute the derivatives
%

% get the wavenumbers
%
k = freq0(ni,L);

% create copies for second dimension
%
k = 1.i * repmat(k,1,n2);

% compute the derivative
%
DXIhat = k.*XIhat;
DXI = real( ifft(DXIhat) );







