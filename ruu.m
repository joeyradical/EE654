function [val] = ruu(index,a1, a2)
%ruu0 Calculates value of ruu
%   input parameters: index, a1, a2
ruu0 = ((1+a2)/(1-a2))*(1/((1+a2)^2-a1^2))
if index == 0
    val = ruu0
elseif index == 1
    val = -a1/(1+a2)*ruu0
elseif index == 2
    val = -a2+(a1^2)/(1+a2)*ruu0
end
end

