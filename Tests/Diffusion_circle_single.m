% Reaction Diffusion using IBDL New Neumann method
%   reactions are FitzHugh-Nagumo
%
%
function [usolutions,timeseries] = Diffusion_circle_single(Ny,D,Tmax,dt,...
                                                    printflag,recordflag,rhsMaskFlag)

    addpath('../src/');

    if recordflag
        if ~printflag
            error("Can't record without plotting")
        else
            stationarydiff = VideoWriter('diffusion_stationary.avi');
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
    
    % Cartesian grid
    % 
    xg=dx*(0.5:Nx-0.5)+xmin;
    yg=dx*(0.5:Ny-0.5)+ymin;
    [xg,yg]=ndgrid(xg,yg);
    
    % IB points for a circle
    %
    [X0, ds] = circle(xc,yc,rad,ds);
    Nib=length(X0(:,1));

    % time stepping 
    %
    timeseries=(0:dt:Tmax)';
    Nt=length(timeseries);
    usolutions=zeros(Nx,Ny,Nt); %for storing u solutions
    
    
    % initial data functions
    %
    u_0=(sqrt(xg.^2 + yg.^2)<0.95).*exp(-((xg+0.25).^2+(yg+0.2).^2)./(0.3^2));
    u_0 = u_0*3;
    

    
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
    
    % normals 
    %
    unitnormal=const2*1/rad*(X0-repmat([xc,yc],Nib,1));
    
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
    grid.chi  = chi;
    grid.bcx = 'per';
    grid.bcy = 'per';
    
    % pack up info on the IB 
    %
    IB.Nib     = length(X0);
    IB.normals = -unitnormal;
    IB.dsvec   = ds*ones(IB.Nib,1);
    
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
    

    
    
    for n=2:Nt  %time loop
    
        % store last time step
        %
        uold = u;
    
        % Mask u?
        %
        u = uold.*rhsMask;
    
        % update for v
        %
        rhs = u/dt;
        Vb = zeros(Nib,1);

        [u,Fds] = IBSL_Nmn_Solve(rhs,X0,IB,a,b,grid,solveparams,Vb);
        usolutions(:,:,n) = u;
        
        % visualize
        %
        if printflag
            figure(3)
            pcolor(xg,yg,u)
            % caxis([0 1.5])
            colorbar
            shading flat
            hold on
            % plot3(X0(:,1),X0(:,2),ones(size(X0(:,1))),'r','LineWidth',2) 
            plot(X0(:,1),X0(:,2),'r','LineWidth',2) 
            % quiver(X0(:,1),X0(:,2),IB.normals(:,1),IB.normals(:,2),'k')
            title(sprintf('time = %f',(n-1)*dt))
            pause(0.01)
            set(gca,'FontSize',14)
            caxis([0 3])
            xlim([xmin xmin+Lx])
            ylim([ymin ymin+Ly])
            hold off
            if recordflag
                writeVideo(stationarydiff,getframe(gcf));
            end %End save video conditional
        end
            
    end  %end time loop
    if recordflag
        close(stationarydiff)
        !HandBrakeCLI -i diffusion_stationary.avi -o diffusion_stationary.mp4
        !rm diffusion_stationary.avi
    end

    rmpath('../src/');

end

