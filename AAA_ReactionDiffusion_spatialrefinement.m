% Reaction Diffusion using IBDL New Neumann method
%Interior of a circle 
%du dn = 0 on circle
%Spatial Discretization for Refinement

clear all;

xmin        = -1.5;            
ymin        = -1.5;
Ly          = 3;            % height of the domain
aspect      = 1;             % aspect ratio
Lx          = aspect*Ly;     % length of th domain
xc          = 0;             % center of the IB object xc, yc
yc          = 0;
rad         = 1;           % restding radius of the circle

dsscale=0.75;  %ratio of IB points to grid spacing

Nyvals=2.^(6:10);  %To do spatial refinement study

Da=2.5e-4;%2.5e-4;
Db=1e-2;%e-2;
mu=2;
nu=2;
Tmax=8;

dtvals=10*2.^-(6:10);

perturb=@(theta, r) (cos(theta).*(1-r.^2).*r.^2);

%Reaction terms
Ra=@(uA,uB)(uA./uB-mu).*uA;
Rb=@(uA,uB)(uA.^2-nu*uB);




Nyfine=Nyvals(end);
Nxfine=aspect*Nyfine;
dxfine=Ly/Nyfine;

uasolutions=cell(1,length(Nyvals)); %for storing u solutions to compare
ubsolutions=cell(1,length(Nyvals)); %for storing v solutions to compare


const1=-1/2; %plus or minus 1/2 Q ( if const2=1, minus interior, plus exterior)
const2=1; %which normal direction i want to use; 1 for pointing out of 
    %circle; -1 for pointing into circle


for iii=1:length(Nyvals)   %This loop is going to do the finest mesh first
    gl = length(Nyvals)-iii+1;  %grid level
                             
    Ny=Nyvals(gl); % number of mesh points in y-direction
    Nx = aspect*Ny;   % number of mesh points in x-direction
    dx = Ly/Ny;     % fluid mesh spacing
    dy=dx;
    ds = dsscale*dx;     % IB mesh spacing 

    
    dt=dtvals(gl);
    ts=(dt:dt:Tmax)';
    Nt=length(ts);

    a1=1/dt;

 % grid points 
    xg=dx*(0:Nx-1)+xmin;
    yg=dx*(0:Ny-1)+ymin;
    [xg,yg]=ndgrid(xg,yg);
    [thetag,rg]=cart2pol(xg, yg);

    rg=sqrt((xg-xc).^2+(yg-yc).^2);
     chi = 1.0*( rg < rad);

    % IB points 
    [X0, ds] = circle(xc,yc,rad,ds);
    Nib=length(X0(:,1));
    sp_scale = ds/dx^2;
    x0=X0(:,1);
    y0=X0(:,2);
    [theta0,r0]=cart2pol(x0, y0);
   
    S=spreadmatrix_vc_vec(X0, dx, Nx,Ny, xmin, ymin);
    unitnormal=const2*1/rad*X0;

    %Define the solution values on the boundary
    Vb=zeros(length(X0(:,1)),1);

    ua_0=perturb(thetag, rg)+nu/mu;
    ub_0=nu/mu^2*ones(size(ua_0));

    ua_old=ua_0;
    ub_old=ub_0;

    for n=1:Nt  %time loop

        %Solve (I/dt-DaL)uAnew+SFanew = ua_old/dt +Ra_old
        %      (I/dt-DbL)ubnew+SFbnew=ub_old/dt+Rb_old
        %    ( S*grad uanew)dot n -1/(2Da) Fa_new = 0
        %    ( S*grad ubnew)dot n -1/(2Db) Fb_new = 0 

        rhs=zeros(Nx,Ny,2);
        rhs(:,:,1)= ua_old/dt+Ra(ua_old, ub_old);
        rhs(:,:,2)=ub_old/dt+Rb(ua_old, ub_old);
        % rhs=rhs.*chi;

        %Dealing with the rhs for gmres: 
        %Do negative L inverse:
        Lga=helmholtz_solve_FD(rhs(:,:,1), a1, Da, Lx,Ly, dx,dy);% - Linv g  Nx x Ny x 2
        Lgb=helmholtz_solve_FD(rhs(:,:,2), a1, Db, Lx,Ly, dx,dy);% - Linv g  Nx x Ny x 2
        %Do gradient
        [Dx, Dy]=dx2d(Nx,Ny);
        Dx=1/(2*dx)*Dx;
        Dy=1/(2*dx)*Dy;
        GLga1=Dx*reshape(Lga, Nx*Ny, 1); %NxNy x 1
        GLga2=Dy*reshape(Lga, Nx*Ny, 1); %NxNy x1
        GLgb1=Dx*reshape(Lgb, Nx*Ny, 1); %NxNy x 1
        GLgb2=Dy*reshape(Lgb, Nx*Ny, 1); %NxNy x1        
        %Interpolate onto boundary
        SGLga1=S'*GLga1; %Nib x1
        SGLga2=S'*GLga2; %Nib x1
        SGLgb1=S'*GLgb1; %Nib x1
        SGLgb2=S'*GLgb2; %Nib x1
        %dot with normal
        nSGLga=zeros(Nib,2);
        nSGLga(:,1)=unitnormal(:,1).*SGLga1+unitnormal(:,2).*SGLga2; %Nib x1
        nSGLgb=zeros(Nib,2);
        nSGLgb(:,1)=unitnormal(:,1).*SGLgb1+unitnormal(:,2).*SGLgb2; %Nib x1


        Fa=gmres(@(Fa)helmNeummconstraintfunct_circ_FD_NEW(Fa,xmin,ymin, Ly, ...
        aspect, Ny,ds, X0, a1, Da, const1, const2, xc,yc,rad), ...
        Vb(:,1)+nSGLga(:,1),10, 1e-6, 1000);  %+SLg

         Fb=gmres(@(Fb)helmNeummconstraintfunct_circ_FD_NEW(Fb,xmin,ymin, Ly, ...
        aspect, Ny,ds, X0, a1, Db, const1, const2, xc,yc,rad), ...
        Vb(:,1)+nSGLgb(:,1),10, 1e-6, 1000);  %+SLg   

             %Find u everywhere
        [ua_new,Ua]=IBSL_helm_rhs_FD(a1,Da,xg,yg,dx, X0, ds, Fa(:,1), rhs(:,:,1));
        [ub_new,Ub]=IBSL_helm_rhs_FD(a1,Db,xg,yg,dx, X0, ds, Fb(:,1), rhs(:,:,2));

ua_old=ua_new;
ub_old=ub_new;

    end  %end time loop

uasolutions{gl}=ua_new.*chi;
ubsolutions{gl}=ub_new.*chi;


end  %end spatial refinement loop
%%

d1ua=zeros(1,length(Nyvals)-1);
d2ua= zeros(1,length(Nyvals)-1);
dmaxua=zeros(1,length(Nyvals)-1);
d1ub=zeros(1,length(Nyvals)-1);
d2ub= zeros(1,length(Nyvals)-1);
dmaxub=zeros(1,length(Nyvals)-1);

for kk2=1:length(Nyvals)-1
    %from fine to coarse 
    gl=length(Nyvals)-kk2+1;
    ufinea=uasolutions{gl}; 
    ufineresta=ufinea(1:2:end, 1:2:end); 
    ucoarsea=uasolutions{gl-1};
    diffua=(ucoarsea-ufineresta);
   
    ufineb=ubsolutions{gl}; 
    ufinerestb=ufineb(1:2:end, 1:2:end); 
    ucoarseb=ubsolutions{gl-1};
    diffub=(ucoarseb-ufinerestb);

    Ny=Nyvals(gl-1); % number of mesh points in y-direction
    Nx = aspect*Ny;   % number of mesh points in x-direction
    dx = Ly/Ny;     % fluid mesh spacing
    dy=dx;
   
    d1ua(1, gl-1)= dx^2*sum(sum(abs(diffua)));
    d2ua(1, gl-1)= sqrt(dx^2*sum(sum(diffua(:,:).^2)));
    dmaxua(1, gl-1)= max(max(abs(diffua)));

    d1ub(1, gl-1)= dx^2*sum(sum(abs(diffub)));
    d2ub(1, gl-1)= sqrt(dx^2*sum(sum(diffub(:,:).^2)));
    dmaxub(1, gl-1)= max(max(abs(diffub)));
    
end

%% Allnorms refinement study
   
   colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=80;
forplotdx0=10;
dxforplot=1./Nyvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% forploty02=50;
% forplotdx02=10;
% dxforplot2=1./Nyvals;
% yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
   %Refinement study plots
   
   figure;
   loglog(Nyvals(1:end-1), dmaxua,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(Nyvals(1:end-1), d2ua, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
loglog(Nyvals(1:end-1), d1ua, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
   loglog(Nyvals(2:4), yforplot(2:4), 'LineWidth', 2,'Color', colorg)
text(200,1/70,'\Delta x', 'Color', colorg, 'FontSize', 16);
% loglog(Nyvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% text(100,1.4/100000,'\Delta x^2', 'Color', colorg, 'FontSize', 12);
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
legend('L^{\infty}', 'L^2', 'L^1');
axis([50 5000 10^(-5) 10^(0)]);
set(gca, 'FontSize', 12);
xlabel('N_x', 'FontSize', 18)
ylabel('Error Norms', 'FontSize', 18)
title('uA Refinement Study')
   

%%

%% Allnorms refinement study
   
   colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=80;
forplotdx0=10;
dxforplot=1./Nyvals;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% forploty02=50;
% forplotdx02=10;
% dxforplot2=1./Nyvals;
% yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
   %Refinement study plots
   
   figure;
   loglog(Nyvals(1:end-1), dmaxub,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(Nyvals(1:end-1), d2ub, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
loglog(Nyvals(1:end-1), d1ub, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
   loglog(Nyvals(2:4), yforplot(2:4), 'LineWidth', 2,'Color', colorg)
text(200,1/70,'\Delta x', 'Color', colorg, 'FontSize', 16);
% loglog(Nyvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% text(100,1.4/100000,'\Delta x^2', 'Color', colorg, 'FontSize', 12);
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
legend('L^{\infty}', 'L^2', 'L^1');
axis([50 5000 10^(-5) 10^(0)]);
set(gca, 'FontSize', 12);
xlabel('N_x', 'FontSize', 18)
ylabel('Error Norms', 'FontSize', 18)
title('uB Refinement Study')
   




