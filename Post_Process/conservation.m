load ~/Documents/Julia/RadialDiffusion/DiffusionTests/matfile.mat
load 128_square.mat


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
xg=dx*(0:Nx-1)+xmin;
yg=dx*(0:Ny-1)+ymin;
[xg,yg]=ndgrid(xg,yg);


chi = ((xg.^2 + yg.^2) <= 1);


dt = 5e-2;
Nt = length(soln(1,1,:));


fulleul = zeros(Nt,1);
maskedeul = zeros(Nt,1);
finitediff = zeros(Nt,1);
% polarinterp = zeros(Nt,1);
timevec = zeros(Nt,1);

for i = 1:Nt
    time = (i-1)*dt;

    timevec(i) = time;
    fulleul(i) = sum(sum(usolutions{i}))*dy*dx;
    maskedeul(i) = sum(sum(usolutions{i}.*chi))*dy*dx;
    

    % [foo,bar,polardisc] = cart2pol(Xfine,Yfine,usolutions{i});
    
    % polarinterp(i) = sum(sum(Rgrid.*polardisc))*dr*dtheta;
    finitediff(i) = sum(sum(Rgrid.*soln(:,:,i)))*dr*dtheta;
end

figure(1)
subplot(2,1,1)
plot(timevec,maskedeul,timevec,finitediff,'LineWidth',4)
legend("Masked IB","Finite Difference",'Location','best')
set(gca,'FontSize',14)
title("Absolute integral",'FontSize',14)
subplot(2,1,2)
plot(timevec,maskedeul./maskedeul(1),timevec,finitediff./finitediff(1),'LineWidth',4)
legend("Masked IB","Finite Difference",'Location','best')
set(gca,'FontSize',14)
title("Relative integral",'FontSize',14)



