%Neumann BC
% au-b lapl u = rhs
%NOTICE THE DIFFERENT OPERATOR THAN MY OTHER HELMHOLTZ CODES!!!!!


clear all;

addpath('./src/')

xmin        = -0.5;            
ymin        = -0.5;
Ly          = 1;            % height of the domain
aspect      = 1;             % aspect ratio
Lx          = aspect*Ly;     % length of th domain
xc          = 0;             % center of the IB object xc, yc
yc          = 0;
rad         = 0.25;           % restding radius of the circle

dsscale=0.75;  %ratio of IB points to grid spacing

%Helmholtz parameters
a = 1; 
b = -1;

a1=2;
a2=3;

%Robin BC's a1u + a2 du/dn = Rb

% solver parameters
    %
    solveparams.rstart = 10;
    solveparams.tol    = 1e-10;
    solveparams.maxiter = 1000;

%For refinement Study  
Nyvals= 2.^(6:10) ;

%for storing u solutions to compare
usolutions=cell(1,length(Nyvals)); 
uexactsolutions=cell(1,length(Nyvals));

%for storing the ds values (to be used in finding errors on the boundary)
dsvalues = zeros(1,length(Nyvals)); 

const1=-1/2; %plus or minus 1/2 Q ( if const2=1, minus interior, plus exterior)
const2=1; %which normal direction i want to use; 1 for pointing out of 
    %circle; -1 for pointing into circle

for iii=1:length(Nyvals)   %This loop is going to do the finest mesh first
    gl = length(Nyvals)-iii+1;  %grid level
                             
    Ny=Nyvals(gl) % number of mesh points in y-direction
    Nx = aspect*Ny;   % number of mesh points in x-direction
    dx = Ly/Ny;     % fluid mesh spacing
    dy=dx;
    ds = dsscale*dx;     % IB mesh spacing 
    

    % grid points 
    xg=dx*(0.5:Nx-0.5)+xmin;
    yg=dx*(0.5:Ny-0.5)+ymin;
    [xg,yg]=ndgrid(xg,yg);
    
  
    rg=sqrt((xg-xc).^2+(yg-yc).^2);
     chi = 1.0*( rg < rad);

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

     
    % IB points 
    [X0, ds] = circle(xc,yc,rad,ds);
    Nib=length(X0(:,1));
    sp_scale = ds/dx^2;
    dsvalues(gl)=ds;

    IB = IB_populate(X0);
    if iii==1
        X0fine=X0;
        dsfine=ds;
    end
       S=spreadmatrix_cc_vec(X0,grid);
       unitnormal=const2*1/rad*X0;

      
    %Define the solution values on the boundary
    Rb=zeros(length(X0(:,1)),1);
    x0=X0(:,1);
    y0=X0(:,2);

    %Normal derivative boundary conditions    
    Rb(:,1)=2*a2/rad*(x0.^2-y0.^2)+a1*(x0.^2-y0.^2);
    % Vb(:,1)=4*x0.*cos(x0).*exp(sin(x0))-4*y0.*sin(y0);
  
    % uexactsolns{gl}=usoln;
    uexactsolns{gl}=(xg.^2-yg.^2).*chi;
    % uexactsolns{gl}=(exp(sin(xg))+cos(yg)).*chi;

    %rhs:
    rhs=zeros(Nx,Ny);
    rhs=a*(xg.^2-yg.^2);
    [u,F] = IBSL_Rbn_Solve(rhs,X0,IB,a,b,grid,solveparams,Rb,a1,a2);
    % [u,U]=IBSL_helm_rhs_FD(a,b,xg,yg,dx, X0, ds, F(:,1), rhs(:,:,1));

    %storing solutions and 
    usolutions{gl}=u(:,:,1).*chi; 


end







 
%% Refinement Study
     
e1u=zeros(1, length(Nyvals));
e2u=zeros(1, length(Nyvals));
emaxu=zeros(1, length(Nyvals));


for k=1:length(Nyvals)%-1
    %from fine to coarse (for no reason other than too lazy to change it)
    gl=length(Nyvals)-k+1;
    u=usolutions{gl}; %calculated solution
    usoln=uexactsolns{gl};%exact solution
     erroru=u-usoln;
     Ny=Nyvals(gl); % number of mesh points in y-direction
    Nx = aspect*Ny;   % number of mesh points in x-direction
    dx = Ly/Ny;     % fluid mesh spacing
    dy=dx;

% %%%%%%%%%%Plot errors flat on circle
% if gl==length(Nyvals)
%     Ny=Nyvals(gl); % number of mesh points in y-direction
%     Nx = aspect*Ny;   % number of mesh points in x-direction
%     dx = Ly/Ny;     % fluid mesh spacing
%     xg = dx*(0:Nx-1)+xmin;   
%     yg=dx*(0:Ny-1)+ymin;
%     [xg,yg]=ndgrid(xg,yg);
%     figure;
%     pcolor (xg,yg, (abs(erroru)));
%     axis([-0.35 0.35 -0.35 0.35]);
%     shading flat;
%     colorbar;
%     set(gca, 'ColorScale', 'log', 'FontSize', 12)
%     %caxis([10^(-4) 10^(0)])
%     caxis([10^(-6) 10^(0)])
% end
% 

    
    e1u(1, gl)= dx^2*sum(sum(abs(erroru)));
    e2u(1, gl)= sqrt(dx^2*sum(sum(erroru(:,:).^2)));
    emaxu(1, gl)= max(max(abs(erroru)));
    

       
end   



%% Allnorms refinement study
   
   colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=10;
forplotdx0=10;
dxforplot=1./Nyvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% forploty02=50;
% forplotdx02=10;
% dxforplot2=1./Nyvals;
% yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
   %Refinement study plots
   
   figure;
   loglog(Nyvals, emaxu,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(Nyvals, e2u, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
loglog(Nyvals, e1u, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
   loglog(Nyvals(3:5), yforplot(3:5), 'LineWidth', 2,'Color', colorg)
text(360,1/700,'\Delta x', 'Color', colorg, 'FontSize', 16);
% loglog(Nyvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% text(100,1.4/100000,'\Delta x^2', 'Color', colorg, 'FontSize', 12);
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
legend('L^{\infty}', 'L^2', 'L^1');
axis([50 5000 10^(-5) 10^(0)]);
set(gca, 'FontSize', 12);
xlabel('N_x', 'FontSize', 18)
ylabel('Error Norms', 'FontSize', 18)
title('Double Layer Refinement Study')
   



rmpath('./src/')

