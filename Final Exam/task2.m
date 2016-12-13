clear
close all
% Task a)
n_dat=1024;
n1=10;
bw=0.20;
x1=(-n1*bw:bw:(2*n1-1)*bw/2);   % time sample locations
yy=sinc(x1);                    % Correlation Sequence
x2=(-n1:n1-1);
rr=zeros(n1,n1);  rd=zeros(1,n1);  % Form Correlation Matrix rr and cross Correlation Vector rd 
for n=1:n1
    rr(n1+1-n,:)=yy(n+1:n+n1);
    rd(n1+1-n)=yy(n);
end
add=10^(-3)*eye(n1,n1);  % add small term to Diagonal to raise matrix condition number 
rrp=rr+add;
wts=inv(rrp)*conj(rd');      % form filter weights
aa=[1 -wts'];
fwts=fftshift(20*log10(abs(fft([1 -wts'],1024))));

% Plot impulse response of predictive filter
figure
subplot(2,1,1);
stem(0:length(wts)-1,wts)
title('Impulse response of predictive filter')
subplot(2,1,2);
plot(linspace(-0.5,0.5,1024),fwts);
title('Frequency response of predictive filter')

% Task b
x=0.8*sin(2*pi*(0:(n_dat-1))*0.06);
reg=zeros(1,n1);
err=zeros(1,n_dat);
qq=4; % number of ADC bits
scl=2^(qq-2);
for nt=1:n_dat;
    sm1=x(nt)+reg*wts;
    q_out=round(scl*sm1)/scl;
    n(nt)=q_out;
    err(nt)=sm1-q_out;
    reg=[err(nt) reg(1:n1-1)];
end;

% Plot input and output time series
figure
subplot(2,1,1)
stem(924:n_dat-1,x(925:end),'b')
title('Input time series')
subplot(2,1,2)
stem(924:n_dat-1,n(925:end),'r')
title('Output time series')

% Plot input and output spectrum
figure
ww=kaiser(n_dat)';
ww=ww/sum(ww);
subplot(2,1,1)
plot(linspace(-0.5,0.5,n_dat),fftshift(20*log10(abs(fft(x).*ww))))
title('Input spectrum')
subplot(2,1,2)
plot(linspace(-0.5,0.5,n_dat),fftshift(20*log10(abs(fft(n).*ww))))
title('Output spectrum')

% Task C
f_n = 100;                    % Filter order
f = [0 0.1 0.2 1];         % Frequency band edges
a = [1  1  0 0];           % Desired amplitudes
h = firpm(f_n,f,a);           % Generate filter using remez algorithm
y2=filter(h,1,n);

figure
subplot(2,1,1)
stem(924:n_dat-1,y2(925:end),'r')
title('Filtered time series')
subplot(2,1,2)
plot(linspace(-0.5,0.5,n_dat),fftshift(20*log10(abs(fft(y2).*ww))))
