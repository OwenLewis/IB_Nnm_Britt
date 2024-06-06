function [u,Fds] = IBSL_Nnm_Solve(rhs,X,IB,a,b,grid,solveparams, Vb)
%IBSL Neumann problem for a-b laplacian u
%   rhs the rhs in actual helmholtz problem
%   X, IB points
%More needed inputs at bottom
%Vb, neumann boundary conditions


    % form rhs for SC solve
    %
    rhsSC = -apply_nSGA(rhs,X,a,b,IB,grid);
    
    % Solve for the forces
    %
    SCfun = @(F)(apply_IBNeumann_SC(F,X,a,b,IB,grid));
    Fv    = gmres(SCfun,Vb+rhsSC,solveparams.rstart,solveparams.tol,solveparams.maxiter);
    
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