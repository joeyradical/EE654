clear
close all

% Task 1 a)
N = 1000;
f = 1000;
fs = 10000;
ws = fs*2*pi;
x0 = cos(2*pi*[0:999]*f/fs);
y0 = sin(2*pi*[0:999]*f/fs);
ee = 0.2;
aa = 0.2;
dc_x = 0.2;
dc_y =0.1;

x1 = x0 + dc_x;
y1 = (1+ee).*y0+aa.*x0+dc_y;



% Plot spectrum of signal
figure
ww = kaiser(1000,10)';
ww = ww/sum(ww);
subplot(2,1,1)
plot(linspace(-0.5,0.5,1000)*fs,fftshift(20*log10(abs(fft(x0+j*y0)).*ww)))
title('Signal spectrum')
xlabel('f[hz]')
ylabel('logMagnitude[dB]')

%Plot spectrum of signal with gain/phase imbalance and DC offset
subplot(2,1,2)

plot(linspace(-0.5,0.5,1000)*fs,fftshift(20*log10(abs(fft(x1+j*y1)).*ww)))
title('Signal spectrum with gain/phase imbalance and DC offset')
xlabel('f[hz]')
ylabel('logMagnitude[dB]')

% Plot Lissajous pattern of signal
figure
subplot(1,2,1)
plot(x0,imag(j*y0))
axis('square')
title('Lissajous pattern of signal')
xlabel('I')
ylabel('Q')
ylim([-1.5 1.5])
xlim([-1.5 1.5])

% Plot lissajous pattern of signal with DC offset and gain/phase imbalance
subplot(1,2,2)
plot(x1,imag(j*y1))
axis('square')
title('Lissajous pattern of signal with gain/phase imbalance and DC offset')
xlabel('I')
ylabel('Q')
ylim([-1.5 1.5])
xlim([-1.5 1.5])

% Task b)

% Run DC canceller for I
[x2, dc_hat_x_sv] = DC_canceller(x1, N, 0.01);

% Run DC canceller for Q
[y2, dc_hat_y_sv] = DC_canceller(y1, N, 0.01);

%Plot DC estimates for I and Q
figure
plot(0:N-1,dc_hat_x_sv,'r')
hold on
plot(0:N-1,dc_hat_y_sv,'b')
legend('I','Q')
title('DC Estimate for I and Q')

% Task c

% Run phase balancer
[y3, aa_hat_sv] = phase_balancer(y2,x2,N,0.01);

% Plot estimate of alpha and lissajous pattern of complex signal
figure
subplot(2,1,1)
plot(0:N-1,aa_hat_sv)
title('Estimates of alpha')
subplot(2,1,2)
plot(x2(800:end),imag(j*y3(800:end)))
axis('square')
title('Lissajous pattern of signal with corrected phase')
xlabel('I')
ylabel('Q')
ylim([-1.5 1.5])
xlim([-1.5 1.5])


% Task d

% Run gain balancer
[y4, ee_hat_sv] = gain_balancer(y3,x2,N,0.01);

% Plot epsilon estimates and lissajous pattern of complex signal
figure
subplot(2,1,1)
plot(0:N-1,ee_hat_sv)
title('Estimates of epsilon')
subplot(2,1,2)
plot(x2(800:end),imag(j*y4(800:end)))
axis('square')
title('Lissajous pattern of signal with corrected gain')
xlabel('I')
ylabel('Q')
ylim([-1.5 1.5])
xlim([-1.5 1.5])

% Task e)

N = 20000;
mu = 0.001
% Generate random QPSK data
x0=((floor(2*rand(1,N))-0.5)/0.5);
y0=((floor(2*rand(1,N))-0.5)/0.5);
x1 = x0 + dc_x;
y1 = (1+ee)*y0+aa*x0+dc_y;

% Run DC canceller for I
[x2, ~] = DC_canceller(x1, N, mu);

% Run DC canceller for Q
[y2, ~] = DC_canceller(y1, N, mu);

% Run phase balancer
[y3, ~] = phase_balancer(y2,x2,N,mu);

% Run gain balancer
[y4, ~] = gain_balancer(y3,x2,N,mu);

% Plot signal before corrupted by DC and gain/phase imbalance
figure
subplot(2,3,1)
plot_constellation(x0+j*y0,N,1,'Signal')
subplot(2,3,2)
plot_constellation(x1+j*y1,N,1,'Signal with DC offset and gain/phase imbalance')
subplot(2,3,3)
plot_constellation(x2+j*y2,N,N-1000,'DC-cancelled signal')
subplot(2,3,4)
plot_constellation(x2+j*y3,N,N-1000,'Phase-balanced signal')
subplot(2,3,5)
plot_constellation(x2+j*y4,N,N-1000,'Gain-balanced signal')

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 1,'\bf Constellation Diagrams','HorizontalAlignment' ,'center','VerticalAlignment', 'top');

