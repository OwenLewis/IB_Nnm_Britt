clear all

Ns = [6:10];

L2s = 0*Ns;
L1s = L2s;
Linfs = L2s;


%These are the parameters for the mathematical problem I'm trying to solve:
a = 1;
b = 1;
%I'm inverting the operator L = a*I - b*Delta

a1 = 2;
a2 = 1;
%The boundary conditions are given by 
%a1*u + a2*du/dn = Vb

for k = 1:length(Ns)

    pwer = Ns(k);

    Robin_Test;

    L2s(k) = L2;
    L1s(k) = L1;
    Linfs(k) = Linf;

end


%% Allnorms refinement study
   
colorp=[0.4940, 0.1840, 0.5560];
colorlb=[0.3010, 0.6450, 0.9930];
colorg=[0.4660, 0.6740, 0.1880];
colordb=	[0, 0.4470, 0.7410];

forploty0=Linfs(1);
forplotdx0=2^Ns(1);
dxforplot=1./Ns;
yforplot=forploty0/forplotdx0^1*dxforplot.^1;
% forploty02=50;
% forplotdx02=10;
% dxforplot2=1./Nyvals;
% yforplot2=forploty02/forplotdx02^2*dxforplot2.^2;
   
   %Refinement study plots
   
figure(4);
loglog(2.^Ns, Linfs,'o-', 'LineWidth', 3, 'MarkerSize', ...
       10, 'Color', colorlb)
hold on
loglog(2.^Ns, L2s, 's-', 'LineWidth', 3,'MarkerSize', 10, 'Color', colordb);
loglog(2.^Ns, L1s, 'd-', 'LineWidth', 3,'MarkerSize', 10,'Color', colorp);
loglog(2.^Ns, 2.^(-Ns), 'LineWidth', 2,'Color', colorg)
% text(360,1/700,'\Delta x', 'Color', colorg, 'FontSize', 16);
% loglog(Nyvals(3:5), yforplot2(3:5), 'LineWidth', 1.5,'Color', colorg)
% text(100,1.4/100000,'\Delta x^2', 'Color', colorg, 'FontSize', 12);
hold off
%legend('u', 'u away from boundary','u near boundary', 'FontSize', 12);
legend('L^{\infty}', 'L^2', 'L^1');
% axis([50 5000 10^(-5) 10^(0)]);
set(gca, 'FontSize', 12);
xlabel('N_x', 'FontSize', 18)
ylabel('Error Norms', 'FontSize', 18)
title('Double Layer Refinement Study')

