function [u,Fds] = IBSL_Solve(rhs,X,IB,a,b,grid,gridfunction,solveparams)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % form rhs for SC solve
    %
    rhsSC = -apply_nSGA(rhs,X,a,b,IB,grid,gridfunction);
    
    % Solve for the forces
    %
    SCfun = @(F)(apply_IBNeumann_SC(F,X,a,b,IB,grid,gridfunction));
    Fv    = gmres(SCfun,rhsSC,solveparams.rstart,solveparams.tol,solveparams.maxiter);
    
    % spread operator
    %
    S = spreadmatrix_vc_vec(X,grid);

    % spread the force
    %
    Fds = IB.dsvec .* Fv;
    SF = S*Fds/grid.dx^2;
    SF = reshape(SF,grid.Nx,grid.Ny);

    % update the concentration
    %
    gridproblem = str2func(gridfunction);
    u = gridproblem(rhs + SF,a,b,grid);

end