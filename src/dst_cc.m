function uhat=dst_cc(u);
% DST_CC  Discrete Sine Transfrom of cell centered data 
%

% perform some error checking
%
if nargin == 0,
    error(generatemsgid('Nargchk'),'Not enough input arguments.');
end

if isempty(u)
   uhat = [];
   return
end
 
% If input is a vector, make it a column:
%
do_trans = (size(u,1) == 1);
if do_trans, u = u(:); end

% record the size
%
n = size(u,1);
m = size(u,2);


% weights for to multiply what comes out of fft
%   this is probably not right for even n on the highest 
%   wave number -- adjust in that case
%   also off by some fact of n??? not sure what it should be
%   
%   the fix I have below seems to work, but I need to work out why
%
%w = 0.5*i*exp(-i*pi*0.5/n * (1:n).' );
%w(1:n-1) = 2*w(1:n-1)/sqrt(2*n);
%w(n) = w(n)/sqrt(n);
%w = i*exp(-i*pi*0.5/n * (1:n).' );
%w(1:n-1) = w(1:n-1)/sqrt(2*n);
%w(n) = 0.5*w(n)/sqrt(n);
w = i*exp(-i*pi*0.5/n * (1:n).' )/sqrt(2*n);
w(n) = w(n)/sqrt(2);



% extend as an odd function
%
uu = zeros(2*n,m);
uu(1:n,:) = u;
uu(n+1:2*n,:) = -flipud(u);

% compute the fft and truncate
%
uhat = fft(uu);
uhat = uhat(2:n+1,:);

% multiply by weights
%
uhat = w(:,ones(1,m)).*uhat;

% check if needs to be real
%
if isreal(u), uhat = real(uhat); end

% if transposed, transpose back
%
if do_trans, uhat = uhat.'; end
