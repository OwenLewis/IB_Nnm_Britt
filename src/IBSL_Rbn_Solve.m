function [u,Fds] = IBSL_Rbn_Solve(rhs,X,IB,a,b,grid,solveparams,lambda,Vb)
%IBSL Neumann problem for a-b laplacian u
%   rhs the rhs in actual helmholtz problem
%   X, IB points
%lambda is the parameter for the dirichlet term in the robin BC
%Vb, Robin boundary conditions


    % form rhs for SC solve
    %
    rhsSC = -apply_nSGA(rhs,X,a,b,IB,grid);
    rhsSS = -apply_SA(rhs,X,a,b,IB,grid);
    
    % Solve for the forces
    %
    SCfun = @(F)(apply_IBRobin_SS(F,X,a,b,IB,grid));
    Fv    = gmres(SCfun,Vb+rhsSC+rhsSS,solveparams.rstart,solveparams.tol,solveparams.maxiter);
    
    % spread operator
    %
    S = spreadmatrix_cc_vec(X,grid);

    % spread the force
    %
    Fds = IB.dsvec .* Fv;
    SF = S*Fds/grid.dx^2;
    SF = reshape(SF,grid.Nx,grid.Ny);

    % update the concentration
    %
    u = helmsolve(rhs-SF,a,b,grid);

end