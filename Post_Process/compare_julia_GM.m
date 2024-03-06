load ~/Documents/Julia/RadialDiffusion/PatternTests/matfile.mat
load ../64_pattern.mat


xmin        = -1.5;            
ymin        = -1.5;
Ly          = 3;            % height of the domain
aspect      = 1;             % aspect ratio
Lx          = aspect*Ly;     % length of th domain
xc          = 0;             % center of the IB object xc, yc
yc          = 0;
rad         = 1;
Nyfine=length(usolutions(:,1,1)); 
Nxfine=aspect*Nyfine;
dxfine=Ly/Nyfine;
Ny=Nyfine; 
Nx = aspect*Ny;
dx = Ly/Ny;     
dy=dx;
xg=dx*(0:Nx-1)+xmin;
yg=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);


chi = ((xg.^2 + yg.^2) <= 1);
chi = double(chi);
I = find(chi < eps);
chi(I) = NaN;

plotflag = 0;

if plotflag
    discvid = VideoWriter('pattern_discrepancy.avi');
    open(discvid);
end

dt = 1e-1;
Nt = length(solA(1,1,:));

discrep = zeros(size(solA));

absdisc = zeros(Nt,1);
reldisc = zeros(Nt,1);
timevec = zeros(Nt,1);

for i = 1:Nt
    [Xfine,Yfine,Cart] = pol2cart(thetagrid,Rgrid,solA(:,:,i));
    matlabfine = interpn(xg,yg,usolutions(:,:,i),Xfine,Yfine,'spline');
    discrepancy = abs(matlabfine - Cart);
    time = (i-1)*dt;
    timevec(i) = time;
    if plotflag
        titlestring = sprintf('Time = %f',time);
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
        pcolor(xg,yg,chi.*usolutions(:,:,i))
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
        pause(0.01)
        writeVideo(discvid,getframe(gcf));
    end

    [foo,bar,polardisc] = cart2pol(Xfine,Yfine,discrepancy);
    discrep(:,:,i) = polardisc;

    absdisc(i) = sum(sum(Rgrid.*polardisc))*dr*dtheta;
    reldisc(i) = sum(sum(Rgrid.*polardisc./solA(:,:,i)))*dr*dtheta;
end

if plotflag
    !HandBrakeCLI -i pattern_discrepancy.avi -o pattern_discrepancy.mp4
    close(discvid)
end

