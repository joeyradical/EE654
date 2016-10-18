function [null] = plot_results(x1,x3,ws,ps,errs,mu,str)
%plot_results Plots results of adaptive line canceller.
figure
pos = [-1 -1 2 2];
unitcircle =rectangle('Position',pos,'Curvature',[1 1]);
set(unitcircle,'edgecolor','b')
axis equal
grid on
hold on
plot(ws,'rx')
title(['Poles closing in on frequency when mu =' num2str(mu) str])
figure
plot(0:length(x1)-1,x1, 'b')
hold on
plot(0:length(ws)-1,unwrap(angle(ws))/(-2*pi),'r')
title(['Frequency Profile and Estimated Frequency of 10 Frequency hops when mu=' num2str(mu) str])
axes('position',[0.65 0.175 0.25 0.25])
box on
plot(195:400,x1(196:401),'b')
hold on
plot(195:400,unwrap(angle(ws(196:401)))/(-2*pi),'r')
xlim([195 400])
ylim([0.025 0.17])
title('Zoom to Transient Detail')
axes('position',[0.175 0.65 0.25 0.25])
box on
plot(195:400,x3(196:401),'b')
hold on
plot(195:400,ps(196:401),'r')
xlim([195 400])
ylim([-1.1 1.1])
title('Input Time Series(Blue) and Predicted Time Series(Red)')

figure
plot(0:length(errs)-1,abs(errs))
title(['Absolute value of prediction error when mu =' num2str(mu) str])
end

