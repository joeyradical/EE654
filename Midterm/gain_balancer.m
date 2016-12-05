function [y, ee_hat_sv] = gain_balancer(x1,x2,N,mu)
% Adaptive gain balancer using LMS algorithm
%   x1= Input signal
%   x2= Signal which causes the gain imbalance
%   N= Number of samples
%   mu = self explanatory
w_new = 1;
w_sv = zeros(1,N);
ee_hat_sv = zeros(1,N);
y = zeros(1,N);
for nn = 1:N
    w = w_new;
    ee_hat_sv(nn) = w;
    y(nn) = x1(nn)*w;
    w_new = (abs(x2(nn))-abs(y(nn)))*mu+w;
end
end

