% Diffusion using "IBDL" New Neumann method
%Interior of a circle 
%du dn = 0 on circle

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

Nyfine=2.^6;  

D=1;%2.5e-4;
Tmax=0.5;

dt=5e-4;


Nxfine=aspect*Nyfine;
dxfine=Ly/Nyfine;


const1=-1/2; %plus or minus 1/2 Q ( if const2=1, minus interior, plus exterior)
const2=1; %which normal direction i want to use; 1 for pointing out of 
    %circle; -1 for pointing into circle



                             
Ny=Nyfine; % number of mesh points in y-direction
Nx = aspect*Ny;   % number of mesh points in x-direction
dx = Ly/Ny;     % fluid mesh spacing
dy=dx;
ds = dsscale*dx;     % IB mesh spacing 


ts=(0:dt:Tmax)';
Nt=length(ts);
usolutions=cell(1,Nt); %for storing u solutions


a1=1/dt;

% grid points WHY ARE WE USING AN EDGE CENTERED GRID?
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

u_0=(sqrt(xg.^2 + yg.^2)<0.95).*exp(-((xg+0.25).^2+(yg+0.2).^2)./(0.3^2));
u_0 = u_0*3;

u_old=u_0;
usolutions{1} = u_old;
    for n=2:Nt  %time loop


        %Solve (I/dt-DaL)uAnew+SFanew = ua_old/dt +Ra_old
        %      (I/dt-DbL)ubnew+SFbnew=ub_old/dt+Rb_old
        %    ( S*grad uanew)dot n -1/(2Da) Fa_new = 0
        %    ( S*grad ubnew)dot n -1/(2Db) Fb_new = 0 

        rhs=zeros(Nx,Ny,1);
        rhs= u_old/dt;

        %Dealing with the rhs for gmres: 
        %Do negative L inverse:
        Lg=helmholtz_solve_FD(rhs, a1, D, Lx,Ly, dx,dy);% - Linv g  Nx x Ny x 2
        %Do gradient
        [Dx, Dy]=dx2d(Nx,Ny);
        Dx=1/(2*dx)*Dx;
        Dy=1/(2*dx)*Dy;
        GLg1=Dx*reshape(Lg, Nx*Ny, 1); %NxNy x 1
        GLg2=Dy*reshape(Lg, Nx*Ny, 1); %NxNy x1
        %Interpolate onto boundary
        SGLg1=S'*GLg1; %Nib x1
        SGLg2=S'*GLg2; %Nib x1
        %dot with normal
        nSGLg=zeros(Nib,2);
        nSGLg(:,1)=unitnormal(:,1).*SGLg1+unitnormal(:,2).*SGLg2; %Nib x1


        F=gmres(@(Fa)helmNeummconstraintfunct_circ_FD_NEW(Fa,xmin,ymin, Ly, ...
        aspect, Ny,ds, X0, a1, D, const1, const2, xc,yc,rad), ...
        Vb(:,1)+nSGLg(:,1),10, 1e-6, 1000);  %+SLg  

             %Find u everywhere
        [u_new,U]=IBSL_helm_rhs_FD(a1,D,xg,yg,dx, X0, ds, F(:,1), rhs);
        
        u_old=u_new;
        usolutions{n} = u_old;
        
        pcolor(xg,yg,u_old)
        % caxis([0 1.5])
        colorbar
        shading flat
        hold on
        % plot3(X0(:,1),X0(:,2),ones(size(X0(:,1))),'r','LineWidth',2) 
        plot(X0(:,1),X0(:,2),'r','LineWidth',2) 
        title(sprintf('time = %f',(n-1)*dt))
        pause(0.01)
        hold off
    end


% 
% d1ua=zeros(1,length(Nyvals)-1);
% d2ua= zeros(1,length(Nyvals)-1);
% dmaxua=zeros(1,length(Nyvals)-1);
% d1ub=zeros(1,length(Nyvals)-1);
% d2ub= zeros(1,length(Nyvals)-1);
% dmaxub=zeros(1,length(Nyvals)-1);
% 
% for kk2=1:length(Nyvals)-1
%     %from fine to coarse 
%     gl=length(Nyvals)-kk2+1;
%     ufinea=uasolutions{gl}; 
%     ufineresta=ufinea(1:2:end, 1:2:end); 
%     ucoarsea=uasolutions{gl-1};
%     diffua=(ucoarsea-ufineresta);
% 
%     ufineb=ubsolutions{gl}; 
%     ufinerestb=ufineb(1:2:end, 1:2:end); 
%     ucoarseb=ubsolutions{gl-1};
%     diffub=(ucoarseb-ufinerestb);
% 
%     Ny=Nyvals(gl-1); % number of mesh points in y-direction
%     Nx = aspect*Ny;   % number of mesh points in x-direction
%     dx = Ly/Ny;     % fluid mesh spacing
%     dy=dx;
% 
%     d1ua(1, gl-1)= dx^2*sum(sum(abs(diffua)));
%     d2ua(1, gl-1)= sqrt(dx^2*sum(sum(diffua(:,:).^2)));
%     dmaxua(1, gl-1)= max(max(abs(diffua)));
% 
%     d1ub(1, gl-1)= dx^2*sum(sum(abs(diffub)));
%     d2ub(1, gl-1)= sqrt(dx^2*sum(sum(diffub(:,:).^2)));
%     dmaxub(1, gl-1)= max(max(abs(diffub)));
% 
% end
% 
% %% Allnorms refinement study
% 
%    colorp=[0.4940, 0.1840, 0.5560];
% colorlb=[0.3010, 0.6450, 0.9930];
% colorg=[0.4660, 0.6740, 0.1880];
% colordb=	[0, 0.4470, 0.7410];
% 
% forploty0=80;
% forplotdx0=10;
% dxforplot=1./Nyvals;
% yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% % forploty02=50;
% % forplotdx02=10;
% % dxforplot2=1./Nyvals;
% % yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
% 
%    %Refinement study plots
% 
%    figure;
%    loglog(Nyvals(1:end-1), dmaxua,'o-', 'LineWidth', 3, 'MarkerSize', ...
%        10, 'Color', colorlb)
% hold on
% loglog(Nyvals(1:end-1), d2ua, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
% loglog(Nyvals(1:end-1), d1ua, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
%    loglog(Nyvals(2:4), yforplot(2:4), 'LineWidth', 2,'Color', colorg)
% text(200,1/70,'\Delta x', 'Color', colorg, 'FontSize', 16);
% % loglog(Nyvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% % text(100,1.4/100000,'\Delta x^2', 'Color', colorg, 'FontSize', 12);
% hold off
% %legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
% legend('L^{\infty}', 'L^2', 'L^1');
% axis([50 5000 10^(-5) 10^(0)]);
% set(gca, 'FontSize', 12);
% xlabel('N_x', 'FontSize', 18)
% ylabel('Error Norms', 'FontSize', 18)
% title('uA Refinement Study')
% 
% 
% %%
% 
% %% Allnorms refinement study
% 
%    colorp=[0.4940, 0.1840, 0.5560];
% colorlb=[0.3010, 0.6450, 0.9930];
% colorg=[0.4660, 0.6740, 0.1880];
% colordb=	[0, 0.4470, 0.7410];
% 
% forploty0=80;
% forplotdx0=10;
% dxforplot=1./Nyvals;
% yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% % forploty02=50;
% % forplotdx02=10;
% % dxforplot2=1./Nyvals;
% % yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
% 
%    %Refinement study plots
% 
%    figure;
%    loglog(Nyvals(1:end-1), dmaxub,'o-', 'LineWidth', 3, 'MarkerSize', ...
%        10, 'Color', colorlb)
% hold on
% loglog(Nyvals(1:end-1), d2ub, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
% loglog(Nyvals(1:end-1), d1ub, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
%    loglog(Nyvals(2:4), yforplot(2:4), 'LineWidth', 2,'Color', colorg)
% text(200,1/70,'\Delta x', 'Color', colorg, 'FontSize', 16);
% % loglog(Nyvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% % text(100,1.4/100000,'\Delta x^2', 'Color', colorg, 'FontSize', 12);
% hold off
% %legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
% legend('L^{\infty}', 'L^2', 'L^1');
% axis([50 5000 10^(-5) 10^(0)]);
% set(gca, 'FontSize', 12);
% xlabel('N_x', 'FontSize', 18)
% ylabel('Error Norms', 'FontSize', 18)
% title('uB Refinement Study')
   




