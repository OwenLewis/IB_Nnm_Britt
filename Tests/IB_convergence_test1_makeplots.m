%
% make plots related to the IB convergence test 1
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



% solutions on a coarse grid
%
Nx = 32;
dsscale = 0.5;
deltaflag = 0;
sol = IB_convergence_test1_solve(Nx,dsscale,deltaflag);


% first meshes of the solution dertivative, exact and numerical
%
figure;
mesh(sol.xg,sol.yg,sol.u_ex);
colormap(turbo);
set(gca,'fontsize',FS);
xlabel('x');
ylabel('y');
zlabel('u');
title('exact solution');
view([45 40]);
axis([-2 2 -2 2 -0.2 1]);
colorbar;

figure;
mesh(sol.xg,sol.yg,sol.u);
colormap(turbo);
set(gca,'fontsize',FS);
xlabel('x');
ylabel('y');
zlabel('u');
title('numerical solution');
view([45 40]);
axis([-2 2 -2 2 -0.2 1]);
colorbar;


figure;
mesh(sol.xg,sol.yg,sol.ux_ex);
colormap(turbo);
set(gca,'fontsize',FS);
xlabel('x');
ylabel('y');
zlabel('u_{x}');
title('u_{x} exact');
view([45 40]);
axis([-2 2 -2 2 -1 1]);
colorbar;


figure;
mesh(sol.xg,sol.yg,sol.ux);
colormap(turbo);
set(gca,'fontsize',FS);
xlabel('x');
ylabel('y');
zlabel('u_{x}');
title('u_{x} numerical');
view([45 40]);
axis([-2 2 -2 2 -1 1]);
colorbar;


% do mesh plots of the error
%
figure;
err = abs(sol.u-sol.u_ex);
mesh(sol.xg,sol.yg,err);
set(gca,'zscale','log','colorscale','log');
set(gca,'fontsize',FS);
xlabel('x');
ylabel('y');
zlabel('error');
title('error in solution');
colorbar;

figure;
err = abs(sol.ux-sol.ux_ex);
mesh(sol.xg,sol.yg,err);
set(gca,'zscale','log','colorscale','log');
set(gca,'fontsize',FS);
xlabel('x');
ylabel('y');
zlabel('error');
title('error in u_{x}');
colorbar;





% slices
%
J0 = Nx/2;
xslice = sol.xg(:,J0);
yslice = sol.yg(1,J0);
xib = sqrt(1-yslice^2);


figure;
hp=plot(xslice,sol.u_ex(:,J0),'o-',xslice,sol.u(:,J0),'s-');
set(gca,'fontsize',FS);
set(hp,'linewidth',LW,'markersize',MS);
xlabel('x');
ylabel('u');
title(sprintf('solution along y=%7.5f',yslice));
ylim([0.3 1.05]);
hold on;
plot([xib xib],ylim,'k--',-[xib xib],ylim,'k--','linewidth',2);
hold off;
legend('exact','numerical');


figure;
hp=plot(xslice,sol.u(:,J0)-sol.u_ex(:,J0),'m^-');
set(gca,'fontsize',FS);
set(hp,'linewidth',LW,'markersize',MS);
xlabel('x');
ylabel('error');
title(sprintf('error along y=%7.5f',yslice));
hold on;
plot([xib xib],ylim,'k--',-[xib xib],ylim,'k--','linewidth',2);
hold off;



figure;
hp=plot(xslice,sol.ux(:,J0),'o-',xslice,sol.ux_ex(:,J0),'s-');
set(gca,'fontsize',FS);
set(hp,'linewidth',LW,'markersize',MS);
xlabel('x');
ylabel('u_{x}');
title(sprintf('u_{x} along y=%7.5f',yslice));
hold on;
plot([xib xib],ylim,'k--',-[xib xib],ylim,'k--','linewidth',2);
hold off;
legend('exact','numerical');

figure;
hp=plot(xslice,sol.ux(:,J0)-sol.ux_ex(:,J0),'m^-');
set(gca,'fontsize',FS);
set(hp,'linewidth',LW,'markersize',MS);
xlabel('x');
ylabel('error');
title(sprintf('u_{x} error along y=%7.5f',yslice));
ylim([-0.31 0.31]);
hold on;
plot([xib xib],ylim,'k--',-[xib xib],ylim,'k--','linewidth',2);
hold off;

% plot of the derivative on the boundary
%
figure;
hp=plot(sol.s,sol.Ux_ex,'o-',sol.s,sol.Ux,'s-',sol.s,sol.Ux-0.5*sol.Nv(:,1).*sol.F,'*-');
set(gca,'fontsize',FS);
set(hp,'linewidth',LW,'markersize',MS);
xlabel('\theta');
title('U_{x} on the boundary');
legend('U_{x}','S*u_{x}','S*u_{x}-0.5 F n_{x}');

% error in dertivative on the boundary
%
figure;
hp=plot(sol.s,sol.Ux-0.5*sol.Nv(:,1).*sol.F-sol.Ux_ex,'m^-');
set(gca,'fontsize',FS);
set(hp,'linewidth',LW,'markersize',MS);
xlabel('\theta');
title('U_{x} error on the boundary');










