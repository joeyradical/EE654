clear
close all

x=ones(1,200);
x1=[0.05*x 0.15*x 0.25*x 0.35*x 0.45*x 0.55*x 0.65*x 0.75*x 0.85*x 0.95*x];
x2=filter(1,[1 -1],x1);
x3=exp(j*2*pi*x2);
%x3 = x3+0.1*(randn(1,2000)+j*randn(1,2000))/sqrt(2);

N =2048;
i = 1
fs = figure
subplot(3,5,1:5)
X3 = fft(x3,4*N);
plot(linspace(-0.5,0.5,4*N), fftshift(abs(X3/max(X3))))
title('Spectrum of 10 frequency hops')
ts = figure
subplot(3,5,1:5)
plot(0:length(x3)-1,x3)
title('Time series of 10 frequency hops')
for k = 1:200:2000
     f = linspace(-0.5,0.5,N);
     index = k:k+199;
     figure(ts)
     subplot(3,5,5+i)
     plot(index, x3(index))
     title(['Time Series when f=' num2str(x1(k))])
     X3temp =fftshift(abs(fft(x3(index),N)));
     figure(fs)
     subplot(3,5,5+i)
     plot(f,X3temp/max(X3temp))
     title(['Spectrum Freq(' num2str(i) ') = ' num2str(x1(k))])
     i = i+1;
end
mu = 0.1;
[ws ps errs] = line_canceller(x3, mu);
plot_results(x1,x3,ws,ps,errs,mu,'')

mu = 0.04;
[ws ps errs] = line_canceller(x3, mu);
plot_results(x1,x3,ws,ps,errs,mu,'')

X4=x3+0.1*(randn(1,2000)+j*randn(1,2000))/sqrt(2);

mu = 0.1;
[ws ps errs] = line_canceller(X4, mu);
plot_results(x1,x3,ws,ps,errs,mu,' (noise present)')

mu = 0.04;
[ws ps errs] = line_canceller(X4, mu);
plot_results(x1,x3,ws,ps,errs,mu,' (noise present)')

    




