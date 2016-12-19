clear
close all

% Task A

fs=48000;
f=1000;
N=5000;
% Generate 2000 samples of data
x=cos(2*pi*[0:N-1]*f/fs);
clip=1.3;
x_0=abs(x)/clip;
phi=angle(x);
y1=clip*(x_0./(1+x_0.^6).^(1/6)).*cos(phi);

% Plot non-linear transfer function
figure
subplot(3,1,1)
clip=1.3;
x_dat=0:0.02:2;
x_0=abs(x_dat)/clip;
y_dat=clip*(x_0./(1+x_0.^6).^(1/6));
plot(x_dat,x_dat,'linewidth',2)
hold on
plot(x_dat,y_dat,'r','linewidth',2)
plot([1 1]*clip,[0.80 1.1]*clip,'r','linewidth',2) 
hold off
grid
title('Nonlinear Transfer Function of Amplifier')
text(1.0,0.8,'1-dB Compression Point')

% Plot 200 samples of input and output data
subplot(3,1,2)
plot(0:199,x(1:200),'r')
hold on
plot(0:199,y1(1:200),'b')
title('First 200 samples of input/output data')
legend('Input','Output')

% Plot spectrum of distorted
subplot(3,1,3)
ww=kaiser(2000)';
ww=ww/sum(ww);
plot(linspace(-0.5,0.5,2000)*fs,fftshift(20*log10(abs(fft(y1(1:2000)).*ww))))
title('Spectrum of distorted signal')

% Task B
n_taps=4;
reg=zeros(1,n_taps);
wts=zeros(1,4);
y2=zeros(1,N);
y3=zeros(1,N);
%wts(1)=1;
mu=0.1;
y3_sv=zeros(1,N);
% Run LMS algorithm
for n=1:N
    reg=[x(n) reg(1:3)];
    y2(n)=reg*wts';
    if n>n_taps
        y3(n)=y1(n)-y2(n);
        wts=wts+mu*reg*conj(y3(n));
    end
end

figure
subplot(3,1,1)
plot(0:N-1,20*log10(abs(y3)))
title('Learning curve of LMS algorithm (logmag of error)')
subplot(3,1,2)
plot(1001:1200,y3(1001:1200));
title('200 samples of error signal after transient')
subplot(3,1,3);
plot(linspace(-0.5,0.5,2000)*fs,fftshift(20*log10(abs(fft(y3(1501:3500)).*ww))))
title('Spectrum of distorted signal')
THD=100*var(y3)/var(y2);
disp(['Total harmonic distortion: ' num2str(THD) '%'])

% Task C
reg=zeros(1,n_taps)';

wts=zeros(1,n_taps)';
wts(1)=1;

delta=0.5;
lambda=0.995;
y3=zeros(1,1000);pp=(1/delta)*eye(n_taps);



for n=1:N
    y2(n)=reg'*conj(wts);
    y3(n)=y1(n)-y2(n);
    C=pp*reg;
    KK_eq=C/(lambda+reg'*C);
    wts=wts+KK_eq*conj(y3(n));
    pp=(1/lambda)*pp -(1/lambda)*KK_eq*reg'*pp;
    reg=[x(n); reg(1:3)];
end

figure
subplot(3,1,1)
plot(0:N-1,20*log10(abs(y3)))
title('Learning curve of RLS algorithm (logmag of error)')
subplot(3,1,2)
plot(1001:1200,y3(1001:1200));
title('200 samples of error signal after transient')
subplot(3,1,3);
plot(linspace(-0.5,0.5,2000)*fs,fftshift(20*log10(abs(fft(y3(1001:3000)).*ww))))
title('Spectrum of distorted signal')
THD=100*var(y3)/var(y2);
disp(['Total harmonic distortion: ' num2str(THD) '%'])

