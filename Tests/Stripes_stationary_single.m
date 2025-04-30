% Reaction Diffusion using IBDL New Neumann method
%   reactions are FitzHugh-Nagumo
%
%
function [asolutions,hsolutions,timeseries] = Stripes_stationary_single(Ny,Da,Dh,mu,nu,kap,Tmax,dt,...
                                                    printflag,recordflag,rhsMaskFlag)

    addpath('../src/');

    if recordflag
        if ~printflag
            error("Can't record without plotting")
        else
            stripesvid = VideoWriter('GM_stripes_stationary.avi');
            open(stripesvid);
        end
    end
    
    % computational domain parameters
    %
    xmin   = -2*1.5;          % bottom cornrer of the domain
    ymin   = -2*1.5;
    Ly     = 2*3;          % height of the domain
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
    % [X0, ds] = circle(xc,yc,rad,ds);
    load("ellipse_med.mat");
    Nib=length(X0(:,1));

    % time stepping 
    %
    timeseries=(0:dt:Tmax)';
    Nt=length(timeseries);
    asolutions=zeros(Nx,Ny,Nt); %for storing solutions
    hsolutions=zeros(Nx,Ny,Nt);

    %Reaction terms
    Ra=@(uA,uH)(uA./(uH.*(1+kap*uA.^2))-mu).*uA;
    Rh=@(uA,uH)(uA.^2-nu*uH);
    

    
    % solver parameters
    %
    solveparams.rstart = 10;
    solveparams.tol    = 1e-6;
    solveparams.maxiter = 1000;
    
    
    % domain mask -- for a circle
    %
    [thetag,rg]=cart2pol(xg, yg);
    chi = inpolygon(xg,yg,X0(:,1),X0(:,2));


    % initial data functions
    %
    % perturb=@(theta, r) (cos(theta).*((1-r).^2).*r.^2);
    % perturb=@(theta, r) (cos(4*r.*cos(theta)));%.*((1-r).^2).*r.^2);
    perturb=@(theta,r) rand(size(theta));
    % perturb=@(x,y) sin(pi*(x+y)).*sin(2*pi*(x-y));

    % ua_0=chi.*perturb(thetag, rg)+nu/mu;
    ua_0 = nu/mu + chi.*cos(4*pi*xg);
    uh_0=nu*ones(size(ua_0))/mu^2;
    
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
    grid.deltaflag = 0;
    
    % pack up info on the IB 
    %
    IB = IB_populate(X0);
    
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
        

        if( rhsMaskFlag )
            rhsMask = inpolygon(xg,yg,X0(:,1),X0(:,2));
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
    

        rhsa = ua/dt + Ra(a_old,h_old).*rhsMask;
        rhsh = uh/dt + Rh(a_old,h_old).*rhsMask;

        Vb = zeros(Nib,1);

        [ua,Fdsa] = IBSL_Nmn_Solve(rhsa,X0,IB,a,b1,grid,solveparams,Vb);
        [uh,Fdsh] = IBSL_Nmn_Solve(rhsh,X0,IB,a,b2,grid,solveparams,Vb);
        
        
        asolutions(:,:,n) = ua;
        hsolutions(:,:,n) = uh;
        
        % visualize
        %
        if printflag
            figure(8)
            pcolor(xg,yg,ua.*rhsMask)
            caxis([0 3.1])
            set(gca,'FontSize',14)
            set(gca,'XTick',[-3:3])
            set(gca,'YTick',[-3:3])
            colorbar
            shading flat
            hold on
            % plot3(X0(:,1),X0(:,2),ones(size(X0(:,1))),'r','LineWidth',2) 
            plot(mod(X0(:,1)-xmin,Lx)+xmin,mod(X0(:,2)-ymin,Ly)+ymin,'or','LineWidth',2,'MarkerSize',3) 
            % quiver(X0(:,1),X0(:,2),IB.normals(:,1),IB.normals(:,2),'k')
            title(sprintf('time = %.2f',(n-1)*dt))
            pause(0.01)
            hold off
            if recordflag
                writeVideo(stripesvid,getframe(gcf));
            end %End save video conditional
        end
        % keyboard
            
    end  %end time loop
    if recordflag
        close(stripesvid)
        !HandBrakeCLI -i GM_stripes_stationary.avi -o GM_stripes_stationary.mp4
        !rm GM_stripes_stationary.avi
    end

    rmpath('../src/');

end

