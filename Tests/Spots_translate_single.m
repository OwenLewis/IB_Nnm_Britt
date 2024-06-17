% Reaction Diffusion using IBDL New Neumann method
%   reactions are FitzHugh-Nagumo
%
%
function [asolutions,hsolutions,timeseries] = Spots_translate_single(Ny,v,Da,Dh,mu,nu,Tmax,dt,...
                                                    printflag,recordflag,rhsMaskFlag)

    addpath('./src/');

    if recordflag
        if ~printflag
            error("Can't record without plotting")
        else
            spotsvid = VideoWriter('GM_spots_translate.avi');
            open(spotsvid);
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
    asolutions=zeros(Nx,Ny,Nt); %for storing solutions
    hsolutions=zeros(Nx,Ny,Nt);

    %Check for CLF constraint
    gridspace=min(dx,dy);
    cfl = v*dt/gridspace;
    if cfl > 0.5
        error("CFL Constraint not satisfied")
    end


    %Reaction terms
    Ra=@(uA,uH)(uA./uH-mu).*uA;
    Rh=@(uA,uH)(uA.^2-nu*uH);
    
    

    %Some padded arrays that will be relevant
    xpad = v*ones(Nx+1,Ny);
    ypad = v*ones(Nx,Ny+1);
    

    
    % solver parameters
    %
    solveparams.rstart = 10;
    solveparams.tol    = 1e-6;
    solveparams.maxiter = 1000;
    
    
    % domain mask -- for a circle
    %
    [thetag,rg]=cart2pol(xg, yg);
    chi = 1.0*( rg < rad);


    % initial data functions
    %
    perturb=@(theta, r) (cos(theta).*((1-r).^2).*r.^2);
    ua_0=chi.*perturb(thetag, rg)+nu/mu;
    uh_0=nu/mu^2*ones(size(ua_0));
    
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
    b1 = Da;
    b2 = Dh;
    
    % THIS CHECKS OUT TO MATCH BRITTANY'S CODE
    %
    %Mnew = zeros(IB.Nib);
    %for k=1:IB.Nib
    %  F = zeros(IB.Nib,1);
    %  F(k) = 1;
    %  Mnew(:,k)=apply_IBNeumann_SC(F,X0,a,b,IB,grid);
    %end
    %return
    
    

    ua=ua_0;
    uh = uh_0;
    asolutions(:,:,1) = ua;
    hsolutions(:,:,1) = uh;
    
    
    
    for n=2:Nt  %time loop

        %We move the IB first
        xc = xc + v*dt;
        yc = yc + v*dt;
        [X0, ds] = circle(xc,yc,rad,ds);
        

        if( rhsMaskFlag )
            rhsMask = 1.0*( sqrt((xg - xc).^2 + (yg-yc).^2) < rad);
        else  
            rhsMask = ones(Nx,Ny);
        end 
    
        % store last time step
        %
        a_old = ua;
        h_old = uh;
    
        % Mask u?
        %
        % ua = a_old.*rhsMask;
        % uh = h_old.*rhsMask;

        %Now we need to evaluate the advective terms
        aadvec = upwind_staggered(a_old,xpad,ypad,grid);
        hadvec = upwind_staggered(h_old,xpad,ypad,grid);
    
        % update for v
        %
        rhsa = ua/dt - aadvec + Ra(a_old,h_old).*rhsMask;
        rhsh = uh/dt - hadvec + Rh(a_old,h_old).*rhsMask;

        Vb = zeros(Nib,1);

        [ua,Fdsa] = IBSL_Nmn_Solve(rhsa,X0,IB,a,b1,grid,solveparams,Vb);
        [uh,Fdsh] = IBSL_Nmn_Solve(rhsh,X0,IB,a,b2,grid,solveparams,Vb);
        
        asolutions(:,:,n) = ua;
        hsolutions(:,:,n) = uh;
        
        % visualize
        %
        if printflag
            figure(5)
            pcolor(xg,yg,ua)
            % caxis([0 1.5])
            colorbar
            shading flat
            hold on
            % plot3(X0(:,1),X0(:,2),ones(size(X0(:,1))),'r','LineWidth',2) 
            plot(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,'or','LineWidth',2,'MarkerSize',3) 
            % quiver(X0(:,1),X0(:,2),IB.normals(:,1),IB.normals(:,2),'k')
            title(sprintf('time = %f',(n-1)*dt))
            set(gca,'FontSize',14)
            % caxis([0 6])
            xlim([xmin xmin+Lx])
            ylim([ymin ymin+Ly])
            pause(0.01)
            hold off
            if recordflag
                writeVideo(spotsvid,getframe(gcf));
            end %End save video conditional
        end
            
    end  %end time loop
    if recordflag
        close(spotsvid)
        !HandBrakeCLI -i GM_spots_translate.avi -o GM_spots_translate.mp4
        !rm GM_spots_translate.avi
    end

    rmpath('./src/');

end

