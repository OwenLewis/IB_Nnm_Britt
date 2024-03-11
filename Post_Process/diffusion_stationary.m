load ~/Documents/Julia/RadialDiffusion/DiffusionTests/matfile.mat
load 256_diff_fixed.mat


xmin        = -1.5;            
ymin        = -1.5;
Ly          = 3;            % height of the domain
aspect      = 1;             % aspect ratio
Lx          = aspect*Ly;     % length of th domain
xc          = 0;             % center of the IB object xc, yc
yc          = 0;
rad         = 1;
Ny=length(u(1,:,1));
Nx=length(u(:,1,1));
dy=Ly/Ny;
dx=Lx/Nx;     
x=dx*(0:Nx-1)+xmin;
y=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(x,y);


chi = ((xg.^2 + yg.^2) <= 1);
chi = double(chi);
I = find(chi < eps);
chi(I) = NaN;



[Xfine,Yfine,Cart] = pol2cart(thetagrid,Rgrid,soln);
matlabfine = interpn(xg,yg,u(:,:,end),Xfine,Yfine,'spline');
discrepancy = abs(matlabfine - Cart);


titlestring = sprintf('Time = %f',30);
figure(1)
set(gcf,'Position',[152 78 1154 719])
subplot(2,2,1)
pcolor(Xfine,Yfine,Cart)
title('Julia Fine')
shading flat
colorbar
mycaxis = caxis;
subplot(2,2,2)
% pcolor(Xfine,Yfine,matlabfine)
pcolor(xg,yg,chi.*u(:,:,end))
caxis(mycaxis);
xlim([-1 1]);
ylim([-1 1]);
title('Matlab IB')
shading flat
colorbar
subplot(2,2,3)
pcolor(Xfine,Yfine,discrepancy)
title('Abs. Discrep')
shading flat 
colorbar
subplot(2,2,4)
pcolor(Xfine,Yfine,discrepancy./abs(Cart))
title('Rel. Discrep')
shading flat 
colorbar
sgtitle(titlestring,'fontsize',18)



[foo,bar,polardisc] = cart2pol(Xfine,Yfine,discrepancy);

Linf_soln = max(max(abs(soln)))
L1_soln = sum(sum(Rgrid.*abs(soln)))*dr*dtheta
L2_soln = sum(sum(Rgrid.*(soln.^2)))*dr*dtheta;
L2_soln = sqrt(L2_soln)


Linf_error = max(max(abs(discrepancy)))
L1_error = sum(sum(Rgrid.*abs(discrepancy)))*dr*dtheta
L2_error = sum(sum(Rgrid.*(discrepancy.^2)))*dr*dtheta;
L2_error = sqrt(L2_error)

    % discrep(:,:,i) = polardisc;
    % 
    % absdisc(i) = sum(sum(Rgrid.*polardisc))*dr*dtheta;
    % reldisc(i) = sum(sum(Rgrid.*polardisc./soln(:,:,i)))*dr*dtheta;


