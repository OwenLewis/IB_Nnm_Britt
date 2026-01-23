%
% given velocity on a fine boundary mesh, project onto the course mesh
%   of control points
%
function Uc= restrict_fine_to_coarse(Uf,nc);

% number of fine grid points
%
nf = size(Uf,1);

% index of the highest positive wave number on the coarse grid
%
nplus = floor(nc/2) + 1;

% index of first negative frequency
%
nminus = nplus+ (nf-nc) + 1;

% go to Fourier space
%
Ufhat = fft( Uf );

% keep only coarse grid frequencies
%
Uchat = [Ufhat(1:nplus,:); Ufhat(nminus:end,:)];

% remove highest wave number of nc is even
%
if( mod(nc,2)==0 )
  Uchat(nplus)=0;  
end


% back to real space and rescale
%
Uc = nc/nf*real( ifft(Uchat) );







