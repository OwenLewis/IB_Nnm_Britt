%
% centered_diff_fourth -- coupute the fourth order centered differnce matrix for a periodic grid
%
function D1 = centered_diff_fourth(Nx,dx);
e  = ones(Nx,1);
D1 = spdiags([1*e -8*e 8*e -1*e],[-2 -1 1 2],Nx,Nx);
D1(1,[Nx-1,Nx]) =[1 -8];
D1(2,Nx) = 1;
D1(Nx-1,1) = -1;
D1(Nx,[1,2])=[8,-1];
D1 = D1/(12*dx);

