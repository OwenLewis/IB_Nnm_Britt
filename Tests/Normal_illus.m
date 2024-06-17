clear all
close all
addpath('../src/');
flag = 2;
a = 0;
b =-1;


Ns = [65 129 257 513];
colors = {'b','r','c','k'}

% solver parameters
%
solveparams.rstart = 10;
solveparams.tol    = 1e-9;
solveparams.maxiter = 1000;

for j = 1:length(Ns)
    Ny = Ns(j)
    % computational domain parameters
    %
    xmin   = -1.5;          % bottom cornrer of the domain
    ymin   = -1.5;
    Ly     = 3;          % height of the domain
    aspect = 1;          % aspect ratio
    Lx     = aspect*Ly;  % length of th domain
    
    Nx     = aspect*Ny;  % number of mesh points in x-direction
    dy     = Ly/Ny;      % fluid mesh spacing
    dx     = Lx/Nx;     

    % IB parameters
    %
    xc      = 0.0;         % center of the IB object xc, yc
    yc      = 0.0;
    rad     = 1;        % resting radius of the circle
    dsscale = 0.75;        %ratio of IB points to grid spacing
    ds      = dsscale*dx;  % IB mesh spacing


    % Cartesian grid
    % 
    xg=dx*(0.5:Nx-0.5)+xmin;
    yg=dx*(0.5:Ny-0.5)+ymin;
    [xg,yg]=ndgrid(xg,yg);
    
    % IB points for a circle
    %
    [X0, ds] = circle(xc,yc,rad,ds);
    Nib=length(X0(:,1));
    IB = IB_populate(X0);


    % RHS functions
    %   This one is a gaussian off-center
    % rhs=(sqrt(xg.^2 + yg.^2)<0.95).*exp(-((xg+0.25).^2+(yg+0.2).^2)./(0.3^2));
    % rhs = rhs*3;
    rhs = 0*xg;
    % rhs = rhs - mean(mean(rhs));


    % pack up info about the Eulerian grid in a single variable
    %
    grid.xmin = xmin;
    grid.ymin = ymin; 
    grid.Lx   = Lx;
    grid.Ly   = Ly;
    grid.Nx   = Nx;
    grid.Ny   = Ny;
    grid.dx   = dx;
    grid.dy   = dy;
    % grid.chi  = chi;
    grid.bcx = 'per';
    grid.bcy = 'per';



    S = spreadmatrix_cc_vec(X0,grid);
    switch  flag
        case 1
            F = ones(Nib,1)/2;
            Fds = IB.dsvec .* F;
            SF = S*Fds/grid.dx^2;
            SF = reshape(SF,grid.Nx,grid.Ny);

            u = HelmInv_FD_period(SF,a,b,grid);
            figure(2)
            pcolor(xg,yg,u)
            shading flat
            hold on
            plot(X0(:,1),X0(:,2),'LineWidth',2)
            quiver(X0(:,1),X0(:,2),IB.normals(:,1),IB.normals(:,2),'r')
            title('u')
            hold off

        case 2
            disp('solving')
            Vb = ones(Nib,1);
            [u,Fds] = IBSL_Nmn_Solve(rhs,X0,IB,0,1,grid,solveparams,Vb);

    end

    Gu = gradientFD(u,grid);

    GA  = reshape(Gu,grid.Nx*grid.Ny,2);
    SGA = S'*GA;

    % dot with the normal
    %
    nSGA=SGA(:,1).*IB.normals(:,1)+SGA(:,2).*IB.normals(:,2);
    max(abs(nSGA))

    figure(3)
    pcolor(xg,yg,Gu(:,:,1))
    shading flat
    hold on
    plot(X0(:,1),X0(:,2),'LineWidth',2)
    quiver(X0(:,1),X0(:,2),IB.normals(:,1),IB.normals(:,2),'r')
    % quiver(X0(:,1),X0(:,2),SGA(:,1),SGA(:,2),'k')
    hold off
    title(sprintf('%d',Ny),'FontSize',18)


    ind = ceil(Ny/2)
    figure(99)
    hold on
    linetype = sprintf('%s%s',colors{j},'--');
    plot(xg(:,ind),Gu(:,ind,1),linetype,'linewidth',3)
    linetype = sprintf('%s%s',colors{j},'x');
    plot(X0(1,1),SGA(1,1),linetype,'MarkerSize',12,'LineWidth',5)
    hold off
    xlim([0 xmin+Lx])




end


rmpath('../src/');