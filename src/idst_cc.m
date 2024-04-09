function u=idst_cc_test(uhat);
%
% IDST_CC Inverse discrete sine transform of cell-centered data
%
    
if nargin == 0,
	error(generatemsgid('Nargchk'),'Not enough input arguments.');
end

if isempty(uhat),
   u = [];
   return
end

% If input is a vector, make it a column:
do_trans = (size(uhat,1) == 1);
if do_trans, uhat = uhat(:); end

% record the size
%
n = size(uhat,1);
m = size(uhat,2);

% weights to multiply before transforming
%
w = -i*sqrt(2*n)*exp(i*pi*0.5/n * (1:n).' );
w(n) = sqrt(2)*w(n);

%%
%%
uuhat = zeros(2*n,m);
uuhat(2:n+1,:) = w(:,ones(1,m)).* uhat;
vvhat = conj(w(:,ones(1,m))).*uhat;

uuhat(n+2:2*n,:) = flipud(vvhat(1:n-1,:));


% scale by weights
%
%vhat = w(:,ones(1,m)).*uhat;
 
% extend
%
%uuhat = zeros(2*n,m);
%uuhat(2:n+1,:) = vhat;
%uuhat(n+2:2*n,:) = flipud(conj(vhat(1:n-1,:))); %%% this is a problem for imag it should be -
                                                %%%
%uuhat(n+2:2*n,:) = flipud(vhat(1:n-1,:)); 

% transform back
%
u = ifft(uuhat);
 
% truncate 
%
u = u(1:n,:);
 
% check if real
%
if isreal(uhat), u = real(u); end

% transpose if needed
%
if do_trans, a = a.'; end
