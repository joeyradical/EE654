function [ws, ps, errs] = line_canceller(x3, mu)
%line_canceller One-tap adaptive line canceller
%   Computes the adaptive line canceller algorithm for a sinusoid input.
%   Returns the weights, predictions, and errors for each sample.
w = 0;
w_new = 0;
input_reg = 0;
w_sv = zeros(1,length(x3));
p_sv = zeros(1,length(x3));
err_sv = zeros(1,length(x3));
for nn = 1:length(x3)
    w = w_new;
    w_sv(nn) = w;
    p = conj(w)*input_reg;
    p_sv(nn) = p;
    err = x3(nn) - p;
    err_sv(nn) = err;
    w_new = conj(err)*input_reg*mu+w;
    input_reg=x3(nn);
end

ws = w_sv;
ps = p_sv;
errs = err_sv;

end

