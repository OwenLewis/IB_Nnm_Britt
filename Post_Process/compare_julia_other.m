load ~/Documents/Julia/RadialDiffusion/DiffusionTests/matfile.mat
load 64_square.mat


xmin        = -1.5;            
ymin        = -1.5;
Ly          = 3;            % height of the domain
aspect      = 1;             % aspect ratio
Lx          = aspect*Ly;     % length of th domain
xc          = 0;             % center of the IB object xc, yc
yc          = 0;
rad         = 1;
Nyfine=length(usolutions{1}(:,1)); 
Nxfine=aspect*Nyfine;
dxfine=Ly/Nyfine;
Ny=Nyfine; 
Nx = aspect*Ny;
dx = Ly/Ny;     
dy=dx;
xg=fliplr(dx*(0:Nx-1)+xmin);
yg=fliplr(dx*(0:Ny-1)+ymin);
[xg,yg]=ndgrid(xg,yg);




dt = 5e-2;
Nt = length(soln(1,1,:));

discrep = zeros(size(soln));

for i = 1:51

    time = (i-1)*dt;
    titlestring = sprintf('Time = %f',time);

    [Xfine,Yfine,Cart] = pol2cart(thetagrid,Rgrid,soln(:,:,i));
    juliacoarse = interp2(Xfine,Yfine,soln(:,:,i),xg,yg,'spline');
    discrepancy = abs(usolutions{i} - juliacoarse);
    figure(1)
    set(gcf,'Position',[152 78 1154 719])
    subplot(2,2,1)
    pcolor(xg,yg,juliacoarse)
    title('Julia Fine')
    shading flat
    colorbar
    subplot(2,2,2)
    pcolor(xg,yg,usolutions{i})
    title('Matlab IB')
    shading flat
    colorbar
    subplot(2,2,3)
    pcolor(xg,yg,discrepancy)
    title('Abs. Discrep')
    shading flat 
    colorbar
    subplot(2,2,4)
    pcolor(xg,yg,discrepancy./abs(juliacoarse))
    title('Rel. Discrep')
    shading flat 
    colorbar
    sgtitle(titlestring,'fontsize',18)
    pause(0.01)


end