function [u,Fds,iter] = IBSL_Nmn_Solve(rhs,X,IB,a,b,grid,solveparams,Vb)
%IBSL Neumann problem for a-b laplacian u
%   rhs the rhs in actual helmholtz problem
%   X, IB points
%Vb, neumann boundary conditions


    % form rhs for SC solve
    %
    rhsSC = -apply_nSGA(rhs,X,a,b,IB,grid);
    
    % Solve for the forces
    %
    SCfun = @(F)(apply_IBNeumann_SC(F,X,a,b,IB,grid));
    [Fv,~,~,iter]    = gmres(SCfun,Vb+rhsSC,solveparams.rstart,solveparams.tol,solveparams.maxiter);
    
    % spread operator
    %

    switch grid.deltaflag
        case 0
            S = spreadmatrix_cc_vec(X,grid);
        case 1
            S = spreadmatrix_scalar6ptBspline(X,grid);
        otherwise
            S = spreadmatrix_cc_vec(X,grid);
    end

    % spread the force
    %
    Fds = IB.dsvec.*Fv;
    SF = S*Fds/(grid.dx*grid.dy);
    SF = reshape(SF,grid.Nx,grid.Ny);

    % update the concentration
    %
    u = helmsolve(rhs-SF,a,b,grid);

end