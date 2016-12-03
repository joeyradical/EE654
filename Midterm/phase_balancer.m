function [y, aa_hat_sv] = phase_balancer(x1,x2,N,mu)
% Adaptive phase balancer using LMS algorithm
%   x1= Input signal
%   x2= Signal which causes the phase imbalance
%   N= Number of samples
%   mu = self explanatory

w = 0;
w_new = 0;
w_sv = zeros(1,N);
aa_hat_sv = zeros(1,N);
y = zeros(1,N);
for nn = 1:N
    w = w_new;
    w_sv(nn) = w;
    aa_hat = conj(w);
    aa_hat_sv(nn) = aa_hat;
    err = x1(nn) - aa_hat*x2(nn);
    y(nn) = err;
    w_new = conj(err)*x2(nn)*mu+w;
end

end

