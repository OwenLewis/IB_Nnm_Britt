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
pcolor(xg,yg,chi.*usolutions{i})
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


figure(2)
subplot(2,1,1)
plot(timevec,absdisc,'-x')
subplot(2,1,2)
loglog(timevec,absdisc,'-x',timevec,timevec*2*absdisc(end)/100,'--r')
sgtitle("Abs. Error over time",'fontsize',18)


figure(3)
subplot(2,1,1)
plot(timevec,reldisc,'-x')
subplot(2,1,2)
loglog(timevec,reldisc,'-x',timevec,timevec*2*reldisc(end)/100,'--r')
sgtitle("Rel. Error over time",'fontsize',18)