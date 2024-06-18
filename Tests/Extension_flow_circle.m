% Reaction Diffusion using IBDL New Neumann method
%   reactions are FitzHugh-Nagumo
%
%
% function [usolutions,timeseries] = Extension_flow_circle(Ny,v,D,Tmax,dt,...
                                                     % printflag,recordflag,rhsMaskFlag)

    Ny = 512;
    v = 0.1;
    Da = 2.5e-3;
    Dh = 1.0e-1;
    mu = 1;
    nu = 2;
    kap = 0.05;

    Tmax = 50;
    dt = 0.05;
    printflag = 1;
    recordflag = 1;
    rhsMaskFlag = 1;

    addpath('../src/');

    if recordflag
        if ~printflag
            error("Can't record without plotting")
        else
            extendvid = VideoWriter('extension_stripes.avi');
            open(extendvid);
        end
    end
    
    % computational domain parameters
    %
    xmin   = -3;          % bottom cornrer of the domain
    ymin   = -3;
    Ly     = 6;          % height of the domain
    aspect = 1;          % aspect ratio
    Lx     = aspect*Ly;  % length of th domain
    
    Nx     = aspect*Ny;  % number of mesh points in x-direction
    dy     = Ly/Ny;      % fluid mesh spacing
    dx     = Lx/Nx;     
    
    
    % IB parameters
    %
    xc      = 0;         % center of the IB object xc, yc
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
    [xehoriz,yehoriz] = ndgrid(xce,ycc);
    [xevert,yevert] = ndgrid(xcc,yce);
    % IB points for a circle
    %
    load("ellipse_high.mat");
    % [X0, ds] = circle(xc,yc,rad,ds);
    Nib=length(X0(:,1));

    % time stepping 
    %
    timeseries=(0:dt:Tmax)';
    totalconc = timeseries;
    area = timeseries;
    Nt=length(timeseries);


    velU = sin((xehoriz-xmin)*2*pi/Lx).*cos((yehoriz-ymin)*2*pi/Ly)*Lx;
    figure(1)
    surf(xehoriz,yehoriz,velU);
    shading flat
    xlabel('X')

    velV = -cos((xevert-xmin)*2*pi/Lx).*sin((yevert-ymin)*2*pi/Ly)*Ly;
    figure(2)
    surf(xevert,yevert,velV);
    shading flat
    xlabel('X')

    div = diff(velU,1,1)/dx + diff(velV,1,2)/dy;
    max(max(abs(div)))


    Ucent = (velU(1:end-1,:) + velU(2:end,:))/2;
    Vcent = (velV(:,1:end-1) + velV(:,2:end))/2;
    mag = sqrt(Ucent.^2 + Vcent.^2);
    figure(3)
    quiver(xg,yg,Ucent,Vcent)
    xlabel('X')
    scale = max(max(mag))


    fourrollV = velV*v/scale;
    fourrollU = velU*v/scale;

    timeramp = @(t) (t < 5).*(tanh(t-2.5)/2 + 0.5) + (t > 5).*(tanh(10 - t)/2 + 0.5);
    Ra=@(uA,uH)(uA./(uH.*(1+kap*uA.^2))-mu).*uA;
    Rh=@(uA,uH)(uA.^2-nu*uH);

    
    asolutions=zeros(Nx,Ny,Nt); %for storing solutions
    hsolutions=zeros(Nx,Ny,Nt);

    %Check for CLF constraint
    gridspace=min(dx,dy);
    cfl = v*dt/gridspace;
    if cfl > 0.5
        error("CFL Constraint not satisfied")
    end
    

    
    % solver parameters
    %
    solveparams.rstart = 10;
    solveparams.tol    = 1e-8;
    solveparams.maxiter = 1000;
    
    
    % domain mask -- for a circle
    %
    rg=sqrt((xg-xc).^2+(yg-yc).^2);
    chi = inpolygon(xg,yg,X0(:,1),X0(:,2));
    
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
    IB = IB_populate(X0);



    % u_0=(sqrt(xg.^2 + yg.^2)<0.95).*exp(-((xg+0.25).^2+(yg+0.2).^2)./(0.3^2));
    % u_0 = u_0*3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


    % THIS CHECKS OUT TO MATCH BRITTANY'S CODE
    %
    %Mnew = zeros(IB.Nib);
    %for k=1:IB.Nib
    %  F = zeros(IB.Nib,1);
    %  F(k) = 1;
    %  Mnew(:,k)=apply_IBNeumann_SC(F,X0,a,b,IB,grid);
    %end
    %return
    
    load stripe_ss.mat
    aeq = a(:,:,end);
    heq = h(:,:,end);


    ua=aeq;
    uh=heq;
    asolutions(:,:,1) = ua;
    hsolutions(:,:,1) = uh;


    % this is a constant related to time stepping -- FE/BE
    %    
    % check my SC apply 
    %
    a = 1/dt;
    b1 = Da;
    b2 = Dh;
    
    

    
    for n=2:Nt  %time loop

        foo = -0.7*timeramp((n-1)*dt);

        horiz = foo*fourrollU;
        vert = foo*fourrollV;

        Ssides = spreadmatrix_csides_vec(X0,grid);
        Stops = spreadmatrix_ctops_vec(X0,grid);

        Uvec  = reshape(horiz,(grid.Nx+1)*grid.Ny,1);
        Ulag = Ssides'*Uvec;
        Vvec = reshape(vert,grid.Nx*(grid.Ny+1),1);
        Vlag = Stops'*Vvec;

        Update = dt*[Ulag,Vlag];
        X0 = X0+Update;
        IB = IB_populate(X0);
        % IB.normals = -IB.normals;
    
        % store last time step
        %
        a_old = ua;
        h_old = uh;

        if( rhsMaskFlag )
            rhsMask = inpolygon(xg,yg,X0(:,1),X0(:,2));
        else  
            rhsMask = ones(Nx,Ny);
        end
    
        % Mask u?
        %
        % ua = a_old.*rhsMask;
        % uh = h_old.*rhsMask;

        %Now we need to evaluate the advective terms
        aadvec = upwind_staggered(a_old,horiz,vert,grid);
        hadvec = upwind_staggered(a_old,horiz,vert,grid);
    
        % update for v
        %

        % 
        % 
        rhsa = ua/dt - aadvec + Ra(a_old,h_old).*rhsMask;
        rhsh = uh/dt - hadvec + Rh(a_old,h_old).*rhsMask;
        Vb = zeros(Nib,1);
        % 
        [ua,Fdsa] = IBSL_Nmn_Solve(rhsa,X0,IB,a,b1,grid,solveparams,Vb);
        [uh,Fdsh] = IBSL_Nmn_Solve(rhsh,X0,IB,a,b2,grid,solveparams,Vb);
        % 
        % 
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
            figure(6)
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
            title(sprintf('time = %.2f',(n-1)*dt))
            pause(0.01)
            hold off
            if recordflag
                writeVideo(extendvid,getframe(gcf));
            end %End save video conditional
        end

        % in = inpolygon(xg,yg,X0(:,1),X0(:,2));
        % totalconc(n) = sum(sum(u(in)))*grid.dx*grid.dy;
        % area(n) = sum(sum(in))*grid.dx*grid.dy;
            
    end  %end time loop
    if recordflag
        close(extendvid)
        !HandBrakeCLI -i extension_stripes.avi -o extension_stripes.mp4
        !rm extension_stripes.avi
    end

    
    % S = spreadmatrix_cc_vec(X0,grid);
    % 
    % Gu = gradientFD(u,grid);
    % 
    % gradvec = reshape(Gu,Nx*Ny,2);
    % 
    % gradlag = S'*gradvec;
    % 
    % nml = gradlag(:,1).*IB.normals(:,1) + gradlag(:,2).*IB.normals(:,2);
    % 
    % hold on
    % quiver(X0(:,1),X0(:,2),gradlag(:,1),gradlag(:,2),'r')
    % quiver(X0(:,1),X0(:,2),IB.normals(:,1),IB.normals(:,2),'k')
    % hold off
    % disp('normal derivative on boundary calculated by hand:')
    % max(abs(nml))
   

    % figure(1)
    % plot(timeseries,totalconc./totalconc(1),'r',timeseries,area./area(1),'--k','LineWidth',3)
    % legend('Total Concentration Inside','Area of Inside','location','best')


    rmpath('../src/');
% end

