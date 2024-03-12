% Reaction Diffusion using IBDL New Neumann method
%   reactions are FitzHugh-Nagumo
%
%
function [usolutions,timeseries] = Poiseuille_flow_circle(Ny,v,D,Tmax,dt,...
                                                    printflag,recordflag,rhsMaskFlag)

    addpath('./src/');

    if recordflag
        if ~printflag
            error("Can't record without plotting")
        else
            translatevid = VideoWriter('diffusion_translate.avi');
            open(stationarydiff);
        end
    end
    
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Cell Centered grid
    % 
    xcc=dx*(1/2:Nx-1/2)+xmin;
    ycc=dx*(1/2:Ny-1/2)+ymin;
    [xg,yg]=ndgrid(xcc,ycc);
    % Edge Centered grids in vertical & horizontal directions
    % 
    xce = dx*(0:Nx)+xmin;
    yce = dx*(0:Ny)+ymin;
    [xeh,yeh] = ndgrid(xce,ycc);
    [xev,yev] = ndgrid(xcc,yce);
    % IB points for a circle
    %
    [X0, ds] = circle(xc,yc,rad,ds);
    Nib=length(X0(:,1));

    velU = -(yeh-ymin).*(yeh-ymin-Ly);
    velU = v*velU/max(max(velU));

    velUcells = -(yg-ymin).*(yg-ymin-Ly);
    velUcells = v*velUcells/max(max(velUcells));
    % time stepping 
    %
    timeseries=(0:dt:Tmax)';
    Nt=length(timeseries);


    usolutions=zeros(Nx,Ny,Nt); %for storing u solutions

    %Check for CLF constraint
    gridspace=min(dx,dy);
    cfl = v*dt/gridspace;
    if cfl > 0.5
        error("CFL Constraint not satisfied")
    end
    
    
    % initial data functions
    %
    u_0=(sqrt(xg.^2 + yg.^2)<0.95).*exp(-((xg+0.25).^2+(yg+0.2).^2)./(0.3^2));
    u_0 = u_0*3;

    %Some padded arrays that will be relevant
    xpad = zeros(Nx+1,Ny);
    ypad = zeros(Nx,Ny+1);
    

    
    % solver parameters
    %
    solveparams.rstart = 10;
    solveparams.tol    = 1e-6;
    solveparams.maxiter = 1000;
    
    
    % domain mask -- for a circle
    %
    rg=sqrt((xg-xc).^2+(yg-yc).^2);
    chi = 1.0*( rg <= rad);
    
    % create the mask for the right size if needed
    %
    
    if( rhsMaskFlag )
      rhsMask = chi;
    else  
      rhsMask = ones(Nx,Ny);
    end
      
    % constants involved in the SC equation
    %
    const1=-1/2; % plus or minus 1/2 Q ( if const2=1, minus interior, plus exterior)
    const2=-1;   % which normal direction i want to use; 1 for pointing out of 
                 % circle; -1 for pointing into circle
    

    grid.xmin = xmin+dx/2;
    grid.ymin = ymin+dx/2; 
    grid.Lx   = Lx;
    grid.Ly   = Ly;
    grid.Nx   = Nx;
    grid.Ny   = Ny;
    grid.dx   = dx;
    grid.dy   = dy;
    grid.chi  = chi;
    
    % pack up info on the IB 
    %
    IB = IB_populate(X0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % this is a constant related to time stepping -- FE/BE
    %    
    % check my SC apply 
    %
    a = 1/dt;
    b = D;
    
    % THIS CHECKS OUT TO MATCH BRITTANY'S CODE
    %
    %Mnew = zeros(IB.Nib);
    %for k=1:IB.Nib
    %  F = zeros(IB.Nib,1);
    %  F(k) = 1;
    %  Mnew(:,k)=apply_IBNeumann_SC(F,X0,a,b,IB,grid);
    %end
    %return
    
    

    u=u_0;
    usolutions(:,:,1) = u;
    

    gridproblem = 'HelmInv_FD_period'
    
    
    for n=2:Nt  %time loop

        S = spreadmatrix_vc_vec(X0,grid);

        Uvec  = reshape(velUcells,grid.Nx*grid.Ny,1);
        Ulag = S'*Uvec;
        Vlag = 0*Ulag;

        Update = dt*[Ulag,Vlag];
        X0 = X0+Update;
        IB = IB_populate(X0);
    
        % store last time step
        %
        uold = u;

        % Mask u?
        %
        u = uold.*rhsMask;

        %Now we need to evaluate the advective terms
        xpad(1,:) = uold(end,:);
        xpad(2:end,:) = uold;
        advec = diff(velU.*xpad)/dx;

        % update for v
        %
        rhs = u/dt - advec;


        [u,Fds] = IBSL_Solve(rhs,X0,IB,a,b,grid,gridproblem,solveparams);
        usolutions(:,:,n) = u;
        
        % visualize
        %
        if printflag
            figure(6)
            pcolor(xg,yg,u)
            % caxis([0 1.5])
            colorbar
            shading flat
            hold on 
            plot(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,'or','LineWidth',2,'MarkerSize',3)
            % quiver(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,Ulag,Vlag)
            quiver(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,IB.normals(:,1),IB.normals(:,2))
            title(sprintf('time = %f',(n-1)*dt))
            pause(0.01)
            hold off
            if recordflag
                writeVideo(translatevid,getframe(gcf));
            end %End save video conditional
        end
            
    end  %end time loop
    if recordflag
        close(translatevid)
        !HandBrakeCLI -i diffusion_translate.avi -o diffusion_translate.mp4
        !rm diffusion_translate.avi
    end

    rmpath('./src/');

end

