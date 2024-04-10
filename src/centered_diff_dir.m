%
% scaled centered difference matrix on a periodic domsin
%
function D1 = centered_diff_dir(N);
e  = ones(N,1);
D1 = spdiags(0.5*[-e e],[-1 1],N,N);
D1(1,1) = 0.5;
D1(N,N) =  -0.5;
end

