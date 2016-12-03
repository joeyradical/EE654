function [y, dc_hat_sv] = DC_canceller(x,N, mu)
% Adaptive DC_canceller using LMS algorithm
%   x = input
%   N = number of samples
%   mu = self explanatory

w = 0;
w_new = 0;
w_sv = zeros(1,N);
y = zeros(1,N);
dc_hat_sv = zeros(1,N);
for nn = 1:N
    w = w_new;
    w_sv(nn) = w;
    dc_hat = conj(w);
    dc_hat_sv(nn) = dc_hat;
    err = x(nn) - dc_hat;
    y(nn) = err;
    w_new = conj(err)*mu+w;
end

end

