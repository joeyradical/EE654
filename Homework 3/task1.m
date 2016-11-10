clear
close all
h=rcosine(1,4,'sqrt',0.5,6);        % Shaping Filter
h=h/max(h);                         % Normalize Shaping Filter
h2=reshape([h 0 0 0],4,13);         % Reshape Shaping Filter
 
N_dat=1000;
x0=(floor(2*rand(1,N_dat))-0.5)/0.5+j*(floor(2*rand(1,N_dat))-0.5)/0.5;
x1=zeros(1,4*N_dat);
reg=zeros(1,13);
 
% form Modulator Output
m=0;
for n=1:N_dat
    reg=[x0(n) reg(1:12)];
    for k=1:4
        x1(m+k)=reg*h2(k,:)';
    end
    m=m+4;
end

% form Demodulator Output
x4=filter(h,1,x1)/(h*h');   % no noise, No Channel

subplot(2,2,1)
plot(1:100,real(x1(1:4:400)),'ro')
title('Real part of first 100 symbols at modulator output')
subplot(2,2,2)
plot(1:100,real(x4(1:4:400)),'ro')
title('Real part of first 100 symbols at matched filter output')
subplot(2,2,3)
plot(x1(1:4:end),'ro')
grid on
axis('square')
title('Constellation diagram at modulator output')
subplot(2,2,4)
plot(x4(1:4:end),'ro')
grid on
axis('square')
title('Constellation diagram at matched filter output')


x2=filter([1 0 0 0 0.2 0 0 j*0.1],1,x1);
x3=x2+0.00*(randn(1,4*N_dat)+j*randn(1,4*N_dat))/sqrt(2);
x4=filter(h,1,x3)/(h*h');   % with noise and Channel
xc=filter(h,1,x3)/(h*h');   % with Channel

figure
subplot(2,12,1:6)
plot(1:100,real(x2(1:4:400)),'k')
title('Real part of first 100 symbols at channel output')
subplot(2,12,7:12)
plot(1:100,real(xc(1:4:400)),'k')
title('Real part of first 100 symbols at matched filter output (with channel)')
subplot(2,12,13:16)
plot(x1(1:4:end),'ro')
grid on
axis('square')
title('Constellation diagram at modulator output')
subplot(2,12,17:20)
plot(x2(1:4:end),'ro')
grid on
axis('square')
title('Constellation diagram at channel output')
subplot(2,12,21:24)
plot(xc(1:4:end),'ro')
grid on
axis('square')
title('Constellation diagram at matched filter output (with channel)')
 
reg=zeros(1,40);
wts=zeros(1,40);
wts(4+0)=1;
mu=0.002;
 
m=1;
err_sv=zeros(1,N_dat);
for n=1:4*N_dat
    x5(n)=reg*wts';
    if n>40 && rem(n,4)==1 
        xd=sign(real(x5(n)))+j*sign(imag(x5(n)));
        xe=xd-x5(n);
        err_sv(m)=xe;
        m=m+1;
        wts=wts+mu*reg*conj(xe);
    end
    reg=[x4(n) reg(1:39)];
end
figure
subplot(2,2, [1 2])
plot(0:length(err_sv)-1,abs(err_sv))
title('Learning Curve of Equalizer')
subplot(2,2,3)
plot(1:100,real(x5(1:4:400)))
title('Time Series of Equalizer (first 100 samples)')
subplot(2,2,4)
plot(x5(1:4:1000),'bo')
grid on
axis('square')
hold on
plot(x5(1001:4:end),'ro')
title('Constellation Diagram of Equalizer')

figure
subplot(3,1,1)
plot(real(x2(1:4:end)),'ro')
title('1000 symbols at channel output')
subplot(3,1,2)
plot(real(x4(1:4:end)),'ro')
title('1000 symbols at matched filter output')
subplot(3,1,3)
plot(real(x5(1:4:end)),'ro')
title('1000 symbols at equalizer output')
