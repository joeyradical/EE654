function [null] = plot_constellation(x,N,div,str)
% Plots the constellation diagram of a modulated signal
%   x = input data
%   N = number of samples
%   str = title of plot
%   div = where in the array to change color in the plot
if div > 1
    plot(x(1:div),'bo')
    hold on
end
plot(x(div:end),'ro')

grid on
axis('square')
xlim([-2,2])
ylim([-2,2])
title(str)
end

