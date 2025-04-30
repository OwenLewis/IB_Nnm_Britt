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
LW = 4;   % line width
MS = 12;  % marker size
FS = 20;  % font size

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
  thet = atan2(sol0.X0(:,2),sol0.X0(:,1));
  normderiv0 = sol0.Ux.*sol0.IB.normals(:,1) + sol0.Uy.*sol0.IB.normals(:,2);
  normderiv1 = sol1.Ux.*sol1.IB.normals(:,1) + sol1.Uy.*sol1.IB.normals(:,2);
  Ux0 = normderiv0 - 0.5*sol0.F;  
  Ux1 = normderiv1 - 0.5*sol1.F;  


  truenorm0 = sol0.Ux_ex.*sol0.IB.normals(:,1) + sol0.Uy_ex.*sol0.IB.normals(:,2);
  truenorm1 = sol1.Ux_ex.*sol1.IB.normals(:,1) + sol0.Uy_ex.*sol1.IB.normals(:,2);


  tau = [-sin(thet),cos(thet)];
  tangentderiv0 = sol0.Ux.*tau(:,1) + sol0.Uy.*tau(:,2);
  tangentderiv1 = sol1.Ux.*tau(:,1) + sol1.Uy.*tau(:,2);
  % err0(k) = max( abs(normderiv0 - 0.5) );
  % err1(k) = max( abs(normderiv1 - 0.5) );
  % err0(k) = max( abs(sol0.F + 1) );
  % err1(k) = max( abs(sol1.F + 1) );
  % err0(k) = max( abs(Ux0 - 1) );
  % err1(k) = max( abs(Ux1 + 1) );
  % err0(k) = max(abs(tangentderiv0));
  % err1(k) = max(abs(tangentderiv1));

  normalerr0(k) = max( abs(normderiv0 - 0.5) );
  normalerr1(k) = max( abs(normderiv1 - 0.5) );
  singleerr0(k) = max( abs(sol0.F + 1) );
  singleerr1(k) = max( abs(sol1.F + 1) );
  tangenterr0(k) = max(abs(tangentderiv0));
  tangenterr1(k) = max(abs(tangentderiv1));
  err0(k) = max( abs(Ux0 - 1) );
  err1(k) = max( abs(Ux1 - 1) ); 

end

% output refinement study results to the screen
%
fprintf('%6s %12s %6s %12s %6s\n','N','4pt','ratio','6pt','ratio');
out = [Nx_array(k),err0(k),err0(k-1)/err0(k),err1(k),err1(k-1)/err1(k)];
fprintf('%6i %12.3e %6s %12.3e %6s \n',Nx_array(1),err0(1),' ',err1(1),' ');
for k=2:length(Nx_array)
  out = [Nx_array(k),err0(k),err0(k-1)/err0(k),err1(k),err1(k-1)/err1(k)];
  fprintf('%6i %12.3e %6.3f %12.3e %6.3f \n',out);
end

figure(9);
hp=loglog(Nx_array,singleerr0,'o-',Nx_array,singleerr1,'s-',Nx_array,0.8./Nx_array,'--k');
set(gca,'fontsize',FS);
set(hp,'markersize',MS,'Linewidth',LW);
xlabel('N');
ylabel('max error in $S^* \partial u/ \partial n$','Interpreter','latex');
legend('4-point delta','6-point B-spline','\propto 1/N');


figure(2);
hp=loglog(Nx_array,normalerr0,'o-',Nx_array,normalerr1,'s-',Nx_array,0.3./Nx_array,'--k');
set(gca,'fontsize',FS);
set(hp,'markersize',MS,'Linewidth',LW);
xlabel('N');
ylabel('max error in $F$','Interpreter','latex');
legend('4-point delta','6-point B-spline','\propto 1/N');

figure(3);
hp=loglog(Nx_array,err0,'o-',Nx_array,err1,'s-');
set(gca,'fontsize',FS);
set(hp,'markersize',MS,'Linewidth',LW);
xlabel('N');
ylabel('max error in $S^* \partial u/ \partial n - F/2$','Interpreter','latex');
legend('4-point delta','6-point B-spline');

figure(4);
hp=loglog(Nx_array,tangenterr0,'o-',Nx_array,tangenterr0,'s-',Nx_array,0.1./Nx_array,'--k');
set(gca,'fontsize',FS);
set(hp,'markersize',MS,'Linewidth',LW);
xlabel('N');
ylabel('max error in $S^* \partial u/ \partial \tau$','Interpreter','latex');
legend('4-point delta','6-point B-spline','\propto 1/N');


% make a few plots showing the error on the boundary

% for Nx = [32 128 512]
% 
%   deltaflag = 0;
%   sol0 = IB_convergence_test2_solve(Nx,dsscale,deltaflag);
% 
%   deltaflag = 1;  
%   sol1 = IB_convergence_test2_solve(Nx,dsscale,deltaflag);
% 
%   % evaluate the derivatives on the boundary for the two delta functions
%   %
%   Ux0 = sol0.Ux;  
%   Ux1 = sol1.Ux;  
% 
% 
%   figure;
%   plot(sol0.s,Ux0-sol0.Ux_ex,'o-',sol1.s,Ux1-sol1.Ux_ex,'s-');
%   set(gca,'fontsize',FS);
%   xlabel('\theta');
%   ylabel('error');
%   title(sprintf('U_{x} error, Nx=%i',Nx));
%   legend('4-point','6-point');
% 
% end







