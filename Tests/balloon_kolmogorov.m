%
% perform a simulation an elastic circular membrane in a background extensional flow
%  uses a coarse mesh of control points, and a finer mesh for spreading
%   fine mesh values are computed using Fourier interpolation
%   computation of force and normals is done using Fourier spectral methods
%
% add a strain dependent release of chemical
%
%
addpath('../src/');
addpath('../src/FluidStructure/');
addpath('~/Simulations/Matlab_Post_Pro/colormaps/')

recordflag = 1;

if recordflag
    strainvid = VideoWriter('absorb_balloon.avi');
    open(strainvid);
end
% define physical parameters
%
xmin        = -0.5;            % bottom corner of the domain xmin,ymin
ymin        = -1.0;
Ly          = 2.0;            % height of the domain
aspect      = 1/2;              % aspect ratio
Lx          = aspect*Ly;      % length of th domain
xc          = -0.25;            % center of the IB object xc, yc
yc          = -0.75;
rad         = 0.125;           % resting radius of the circle
ks          = 1;              % stiffness coefficient
strain_init = 0.0;            % relative increas in x-dirction if ellipse
L0          = 0.75;           % length of |X_s| where no tension L0<1,
                              %   rest is under tension and membrane is
                              %   pressurized; L0=1, p=0 at rest


% base simulation was set up on Ny = 128 grid
%  integer factor for refinement
%
refine_fact = 2;


% define numerical parameters
%
Ny = 128*refine_fact;          % number of mesh points in y-direction
Nx = aspect*Ny;   % number of mesh points in x-direction
dx = Ly/Ny;       % fluid mesh spacing

% parameter spacing on the coarse mesh and the refinement ratio for the 
%  the fine mesh
%
dsc        = 3*dx;
fine_ratio = 8;

tend = 10;            % end time of simulation
dt   = 0.05/refine_fact;            % time step
Nt   = round(tend/dt);  % number of steps to take
%Nt   = 5;



% chemical parameters
%   g is the strain dependent release rate
%
% st_thresh = 0.4;
% g=@(t)((t-st_thresh).*(t>st_thresh));
release_rate = 0.05;
chem_decay   = 0;
diff_const   = 0.01;

% solver parameters
%
solveparams.rstart = 10;
solveparams.tol    = 1e-8;
solveparams.maxiter = 1000;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% grid points --- thse are vertex centered grid 
%
xg = xmin + dx*(0:Nx-1);
[xg,yg]=ndgrid(xg);


% Cell Centered grid points
% 
xcc=dx*(1/2:Nx-1/2)+xmin;
ycc=dx*(1/2:Ny-1/2)+ymin;
[xg,yg]=ndgrid(xcc,ycc);
        
% Edge Centered grids in vertical & horizontal directions 
%  because periodic domain we won't store the edge on the 
%  top and bottom
% 
xce = dx*(0:Nx-1)+xmin;
yce = dx*(0:Ny-1)+ymin;
[xehoriz,yehoriz] = ndgrid(xce,ycc);
[xevert,yevert]   = ndgrid(xcc,yce);
    
% make a bacground force to drive the flow
%
fbg = zeros(Nx,Ny,2);
fbg(:,:,1) = -0.1*((2*pi/Ly)^2)*sin(2*pi/Ly*yehoriz);

% solve for the background flow
%
[ubg,~]=stokes_solve_MAC(fbg,Lx,Ly);   
   
% intialize the flow as the background flow
%
u = ubg;

% pack up info about the Eulerian grid in a single variable
%
grid.xmin = xmin;
grid.ymin = ymin; 
grid.Lx   = Lx;
grid.Ly   = Ly;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.dx   = dx;
grid.dy   = dx;
grid.bcx = 'per';
grid.bcy = 'per';
grid.deltaflag = 0;
 

% initialize the control points
%
[Xc, dsc] = circle(xc,yc,rad,dsc);
ds = dsc/fine_ratio;
sp_scale = ds/dx^2;

% numbers of coarse and fine points
%
Nib_c = size(Xc,1);
Nib_f = Nib_c * fine_ratio;

% record the length of the boundary -- this is needed for differentiation
%  it should be Lb = int( |X_{s}| ) = sum( ds ) = sum( dsc )
%
Lb = Nib_c*dsc;

% get parametric coordinate
%
th = linspace(0,2*pi,Nib_f+1);
th = th(1:Nib_f)';


% transform the circle to an ellipse with the same area
%
Xc(:,1) = (1+strain_init)*(Xc(:,1)-xc) + xc;
Xc(:,2) = (Xc(:,2)-yc)/(1+strain_init) + yc;
 
 
% initialize the chemical and the flux
%  
chem = ones(grid.Nx,grid.Ny);
  
  
% loop in time
%
for k=1:Nt
    
    % get the fine grid points
    %
    [X,DXI] = extend_coarse_to_fine(Xc,Nib_f,Lb);
    
    % get the normals on the fine grid
    %
    length_dxds=sqrt(sum(DXI.^2,2));
    unitnormals = [DXI(:,2),-DXI(:,1)]./length_dxds(:,[1 1]);      
      
    % compute IB force from current configuraiton
    %
    [F,st] = stretch_force_fourier(X,ks,Lb,L0);

    % data structure for IB with normals, ds=|dX/d(theta)|d(theta)|
    % 
    IB.normals = unitnormals;
    IB.dsvec   = length_dxds*ds;

    % compue flux on the boundary and the normal derivative
    %
    bdy_flux = -release_rate;
    dcdn = -bdy_flux/diff_const;
    
    % not sure why??
    %
    dcdn = -dcdn;   
       
       
    % advective terms -- pad the velocities with boundary edges
    %
    Upad = u(:,:,1);
    Upad(end+1,:) = Upad(1,:);
    Vpad = u(:,:,2);
    Vpad(:,end+1) = Vpad(:,1);
    cadv_terms = upwind_staggered(chem,Upad,Vpad,grid); 
    
    % check the CFL constraint
    %
    umax = max(abs(u),[],'all');
    CFL = umax*dt/grid.dx;
    fprintf(' Courant number = %g \n',CFL);
     if CFL> 0.5
        error("CFL Constraint not satisfied")
    end
    
    % flag the inside points with a 1
    %
    chi = inpolygon(xg,yg,X(:,1),X(:,2));
    kout = find(chi==0);
    kin = find(chi==1);
    
    
    % equation to solve
    %   ((1+dt*chem_decay) - dt*diff_const)*chem = chem - dt*cadv_terms
    %
    a = 1+dt*chem_decay;
    b = dt*diff_const;
    rhs = chem - dt*cadv_terms;    
    [chem,~] = IBSL_Nmn_Solve(rhs,X,IB,a,b,grid,solveparams,dcdn);

    % spread operators for edges
    %
    Su = spreadmatrix_uedge_per(X,grid);
    Sv = spreadmatrix_vedge_per(X,grid);


    % spread the forces to the grid
    %
    fu = sp_scale*Su*F(:,1);
    fu = reshape(fu,Nx,Ny);
    fv = sp_scale*Sv*F(:,2);
    fv = reshape(fv,Nx,Ny);
    f = cat(3,fu,fv);
    
    % solve stokes equations
    %
    [u,p]=stokes_solve_MAC(f,Lx,Ly);

    % exploit linearity -- add background flow
    %
    if( k*dt < 10)
      u = u + ubg; 
    end
    
    % I = find(g(st) .*(X(:,2)<0) > 0);
    % visualize
    %
    time = k*dt;
    p = p - p(1,1);

    figure(1);
    set(gcf,'Position',[440 278 560.2 420]);
    clf;
 
    % put nans outside the cell
    %  
    chem_vis = chem;
    chem_vis(kin) = NaN;
  
    % first set of axes for concentration and velocity
    % 
    ax1 = axes;

    pcolor(xg,yg,chem_vis); shading interp;
    caxis([0 1]);
   %  colormap(flipud(hot));
   % colormap(sky);
   % colormap(hsv2);
   colormap(flipud(shuncmap));
    

    axis([xmin xmin+Lx ymin 0]);
    axis square;
    set(gca,'YTick',[])
    hold on;
    
          
    % velocity vectors
    %
    velskip = 4*refine_fact;
    II = 1:velskip:grid.Nx;
    JJ = 1:velskip:grid.Ny;    
    velmag = max(sqrt( sum(u.^2,3) ),[],'all');
    velmag0 = 0.1;
    
    hq=quiver(xg(II,JJ),yg(II,JJ),u(II,JJ,1),u(II,JJ,2),velmag/velmag0);
    set(hq,'linewidth',3.5,'color',[0 0 0]);
    
    title(sprintf('time = %.2f',time),'fontsize',14);
    set(ax1,'fontsize',20);
    
    hc=colorbar;
    ylabel(hc,'Conc.','FontSize',20,'Rotation',270);
    
    
    % hold off;
    

        
        
    % plot the IB on separate axes
    %
    ax2 = axes; 
    set(ax2,'fontsize',20);
    Xper = [X; X(1,:)];
    stper = [st; st(1)];
  
    % this visual of the boundary highlights stretching
    %   st here is a scaled tension 
    % 
   % hp = patch([Xper(:,1);NaN],[Xper(:,2); NaN],[stper; NaN],'EdgeColor','interp');
   % hp.LineWidth = 6;
   % colormap(ax2,'copper');
   % clim([0.25 0.5]);
  
    % visual will give a view of stretch and compression
    %
    hp = patch([Xper(:,1);NaN],[Xper(:,2); NaN],[stper-0.25; NaN],'EdgeColor','interp');
    hp.LineWidth = 6;
    foo = flipud(spring);
    colormap(ax2,pink);
    clim([-0.05 0.05]);
    hc=colorbar('westoutside','AxisLocation','in');
    ylabel(hc,'Strain','FontSize',20,'Rotation',270);
    foo = hc.Label.Position;
    foo(1) = foo(1) - 5.5;
    hc.Label.Position = foo;

    % hold on
    % plot(X(I,1),X(I,2),'dk','linewidth',2,'markersize',4)
    % hold off
  
    ax2.UserData = linkprop([ax1,ax2], ...
    {'Position','InnerPosition','DataAspectRatio','xlim','ylim'});
    linkaxes([ax1 ax2]);
    ax2.Visible= 'off';
    myposition = get(ax2,'innerposition');
    offset = 0.2;
    myposition(1) = myposition(1) + offset/5;
    set(ax2,'InnerPosition',myposition)
    
    pause(0.01);
    if recordflag
        writeVideo(strainvid,getframe(gcf));
    end %End save video conditional


%    figure(2);
%    plot(th,st,'linewidth',4);
%    set(gca,'fontsize',16);
%    xlabel('angle');
%    ylabel('strain');
%    ylim([0 0.75]);
  
    
    % interpolate to the IB
    %
    U = 0*X;
    U(:,1) = Su'*reshape(u(:,:,1),Nx*Ny,1);
    U(:,2) = Sv'*reshape(u(:,:,2),Nx*Ny,1);
        
    % project interpolated velocity onto the control points
    %
    U = restrict_fine_to_coarse(U,Nib_c);
    
    % move the points
    %
    Xc = Xc + dt*U;
    
end

if recordflag
    close(strainvid)
    !HandBrakeCLI -i absorb_balloon.avi -o absorb_balloon.mp4
    !rm absorb_balloon.avi
end


