%
% perform a refinement study for the derivative on the boundary
%   exact solution 
%    u = 1       for r < 1
%        1-ln(r) for r>=1
%
%   IB force to produce this solution is F = 1 on the boundary
%
%   solved in a [-2,2]^2 domain with Dirchlet boundary data
%

% some default graphical parameters
%
clear all
close all


LW = 4;   % line width
MS = 12;  % marker size
FS = 20;  % font size

celltol = 8;

% grids to use
%
Nx_array = 32 * 2.^(0:8);
dsscale = 0.5;
%dsscale = 2.0;


for k=1:length(Nx_array)

  Nx = Nx_array(k);

  deltaflag = 0;
  sol0 = Normal_and_density_solve(Nx,dsscale,deltaflag);

  deltaflag = 1;  
  sol1 = Normal_and_density_solve(Nx,dsscale,deltaflag);

  % evaluate the derivatives on the boundary for the two delta functions
  %

  mask = abs(sol0.xg.^2 + sol0.yg.^2 - 1) > celltol*sqrt(sol0.grid.dx*sol0.grid.dy);
  
  xerr0 = abs(sol0.ux - sol0.ux_ex);
  yerr0 = abs(sol0.uy - sol0.uy_ex);

  xerr1 = abs(sol0.ux - sol0.ux_ex);
  yerr1 = abs(sol0.uy - sol0.uy_ex);
  
  magerr0 = sqrt(xerr0.^2 + yerr0.^2);
  magerr1 = sqrt(xerr1.^2 + yerr1.^2);

  
  err0inf(k) = max( max(magerr0) );
  err02(k) = sqrt(sum( sum(magerr0.^2)*sol0.grid.dx*sol0.grid.dy));
  err01(k) = sum( sum(magerr0) )*sol0.grid.dx*sol0.grid.dy;
  err1inf(k) = max( max(magerr1) );
  err12(k) = sqrt(sum( sum(magerr1.^2))*sol1.grid.dx*sol1.grid.dy);
  err11(k) = sum( sum(magerr1) )*sol1.grid.dx*sol1.grid.dy;
  
  maskerr0inf(k) = max( max(mask.*magerr0) );
  maskerr02(k) = sqrt(sum( sum(mask.*magerr0.^2)*sol0.grid.dx*sol0.grid.dy));
  maskerr01(k) = sum( sum(mask.*magerr0) )*sol0.grid.dx*sol0.grid.dy;
  maskerr1inf(k) = max( max(mask.*magerr1) );
  maskerr12(k) = sqrt(sum( sum(mask.*magerr1.^2))*sol1.grid.dx*sol1.grid.dy);
  maskerr11(k) = sum( sum(mask.*magerr1) )*sol1.grid.dx*sol1.grid.dy;

end

% output refinement study results to the screen
%
fprintf('%6s %12s %6s %12s %6s\n','N','4pt','ratio','6pt','ratio');
out = [Nx_array(k),err0inf(k),err0inf(k-1)/err0inf(k),err1inf(k),err1inf(k-1)/err1inf(k)];
fprintf('%6i %12.3e %6s %12.3e %6s \n',Nx_array(1),err0inf(1),' ',err1inf(1),' ');
for k=2:length(Nx_array)
  out = [Nx_array(k),err0inf(k),err0inf(k-1)/err0inf(k),err1inf(k),err1inf(k-1)/err1inf(k)];
  fprintf('%6i %12.3e %6.3f %12.3e %6.3f \n',out);
end

% figure(9);
% hp=loglog(Nx_array,singleerr0inf,'o-',Nx_array,singleerr02,'o-',Nx_array,singleerr01,'o-',Nx_array,singleerr1inf,'s-',Nx_array,singleerr12,'s-',Nx_array,singleerr11,'s-',Nx_array,0.8./Nx_array,'--k');
% set(gca,'fontsize',FS);
% set(hp,'markersize',MS,'Linewidth',LW);
% xlabel('N');
% ylabel('error in $F$','Interpreter','latex');
% legend('L^\infty, 4-point delta','L^2, 4-point delta','L^1, 4-point delta','L^\infty, 6-point B-spline','L^2, 6-point B-spline','L^1, 6-point B-spline','\propto 1/N','location','best');
% 
% 
% figure(2);
% hp=loglog(Nx_array,normalerr0inf,'o-',Nx_array,normalerr02,'o-',Nx_array,normalerr01,'o-',Nx_array,normalerr1inf,'s-',Nx_array,normalerr12,'s-',Nx_array,normalerr11,'s-',Nx_array,0.8./Nx_array,'--k');
% set(gca,'fontsize',FS);
% set(hp,'markersize',MS,'Linewidth',LW);
% xlabel('N');
% ylabel('error in $S^* \partial u/ \partial n$','Interpreter','latex');
% legend('L^\infty, 4-point delta','L^2, 4-point delta','L^1, 4-point delta','L^\infty, 6-point B-spline','L^2, 6-point B-spline','L^1, 6-point B-spline','\propto 1/N','location','best');

figure(3);
hp=loglog(Nx_array,err0inf,'o-',Nx_array,err02,'o-',Nx_array,err01,'o-',Nx_array,err1inf,'s-',Nx_array,err12,'s-',Nx_array,err11,'s-',Nx_array,0.8./Nx_array,'--k');
set(gca,'fontsize',FS);
set(hp,'markersize',MS,'Linewidth',LW);
xlabel('N');
ylabel('error in $\nabla u$','Interpreter','latex');
legend('L^\infty, 4-point delta','L^2, 4-point delta','L^1, 4-point delta','L^\infty, 6-point B-spline','L^2, 6-point B-spline','L^1, 6-point B-spline','\propto 1/N','location','best');

figure(4);
hp=loglog(Nx_array,maskerr0inf,'o-',Nx_array,maskerr02,'o-',Nx_array,maskerr01,'o-',Nx_array,maskerr1inf,'s-',Nx_array,maskerr12,'s-',Nx_array,maskerr11,'s-',Nx_array,0.8./Nx_array,'--k');
set(gca,'fontsize',FS);
set(hp,'markersize',MS,'Linewidth',LW);
xlabel('N');
ylabel('Masked error in $\nabla u$','Interpreter','latex');
legend('L^\infty, 4-point delta','L^2, 4-point delta','L^1, 4-point delta','L^\infty, 6-point B-spline','L^2, 6-point B-spline','L^1, 6-point B-spline','\propto 1/N','location','best');










