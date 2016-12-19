                                                                    
                                             
N_dat=1000
x1a=(floor(2*rand(1,N_dat))-0.5)/0.5+j*(floor(2*rand(1,N_dat))-0.5)/0.5;
x2a=(floor(2*rand(1,N_dat))-0.5)/0.5+j*(floor(2*rand(1,N_dat))-0.5)/0.5;

h=rcosine(1,4,'sqrt',0.5,6);
h=h/max(h);

x1b=reshape([x1a;zeros(3,N_dat)],1,4*N_dat);
x2b=reshape([x2a;zeros(3,N_dat)],1,4*N_dat);

x1c=1*filter(h,1,x1b);    % Shaping signal 1
x2c=filter(h,1,x2b);    % Shaping signal 2

h11=[1 0 0 0.2];            % Channel multipath filter
h21=1*[0.11 0 0 0 0 0.022]*2;   % Crosstalk filter

x1d=1*filter(h11,1,x1c)+1*filter(h21,1,x2c);    % Chan 1 filtered & coupled to chan 2

x1e=filter(h/(h*h'),1,x1d);                 % chan 1 matched filtered



N_len=10;
Reg_cn=zeros(1,4*N_len)';
Reg_eq=zeros(1,4*N_len)';

Wts_eq=zeros(1,N_len)';
Wts_eq(1)=1;
Wts_cn=zeros(1,N_len)';

delta=0.1;
lambda=0.99;
PP_eq=(1/delta)*eye(N_len);
PP_cn=(1/delta)*eye(N_len);
Err1=zeros(1,1000);
Err2=zeros(1,1000);

% input to equalizer is x1e
% input to canceller is x2c
m=1;
for n=1:4*N_dat
    Reg_eq=[x1e(n); Reg_eq(1:4*N_len-1)];
    Reg_cn=[x2c(n); Reg_cn(1:4*N_len-1)];
    
    X1f(n)=Reg_eq(1:4:4*N_len).'*conj(Wts_eq);
    X1g(n)=X1f(n)-Reg_cn(1:4:4*N_len).'*conj(Wts_cn);
    
     if rem(n,4)==1
     Det2=sign(real(X1f(n)))+j*sign(imag(X1f(n)));
        
     Err1(m)=Det2-X1f(n);
     Err2(m)=Det2-X1g(n);
    
     C_eq=PP_eq*Reg_eq(1:4:4*N_len);
     C_cn=PP_cn*Reg_cn(1:4:4*N_len);

      KK_eq=C_eq/(lambda+Reg_eq(1:4:4*N_len)'*C_eq);
      KK_cn=C_cn/(lambda+Reg_cn(1:4:4*N_len)'*C_cn);
      
      Wts_eq=Wts_eq+KK_eq*conj(Err2(m));
      %Wts_cn=Wts_cn-KK_cn*conj(Err2(m));
      Wts_cn=Wts_cn-0.01*Reg_cn(1:4:4*N_len)*conj(Err2(m));

      PP_eq=(1/lambda)*PP_eq -(1/lambda)*KK_eq*Reg_eq(1:4:4*N_len)'*PP_eq;
      PP_cn=(1/lambda)*PP_cn -(1/lambda)*KK_cn*Reg_cn(1:4:4*N_len)'*PP_cn;
      
      m=m+1;
    end
end
% 

% n_len=40;
% wts=zeros(1,n_len)';
% reg=zeros(1,n_len)';
% delta=0.001;
% 
% PP=(1/delta)*eye(n_len);
% 
% 
% for nn=1:200
% 
%     c=PP*reg;
% 
%     KK=c/(lambda+reg'*c);
%     d_hat(nn)=reg'*conj(wts);
%     err1(nn)=data(nn)-d_hat(nn);
%     wts=wts+KK*conj(err1(nn));
%     
%     PP=(1/lambda)*PP -(1/lambda)*KK*reg'*PP;
%     reg=[data(nn) reg(1:n_len-1)']';
% % 
% %  end       
% end
%    
   

   figure(7)
   subplot(2,2,1)
   plot(x1c(1:4:4*N_dat),'r.')
   grid on
   axis('square')
   axis([-1.5 1.5 -1.5 1.5])
   title('Modulator Constellation')
   
   subplot(2,2,2)
   plot(x1d(1:4:4*N_dat),'r.')
   grid on
   axis('square')
   axis([-1.5 1.5 -1.5 1.5])
   title('Receiver Constellation')
      
   subplot(2,3,4)
   plot(x1e(1:4:4*N_dat),'r.')
   grid on
   axis('square')
   axis([-1.5 1.5 -1.5 1.5])
   title('Matched Filter Constellation')

   subplot(2,3,5)
   plot(X1f(1:4:4*N_dat),'b.')
   hold on
   plot(X1f(3001:4:4*N_dat),'r.')
   hold off
   grid on
   axis('square')
   axis([-1.5 1.5 -1.5 1.5])
   title('Equalizer Filter Constellation')

   subplot(2,3,6)
   plot(X1g(1:4:4*N_dat),'b.')
   hold on
   plot(X1g(3001:4:4*N_dat),'r.')
   hold off
   grid on
   axis('square')
   axis([-1.5 1.5 -1.5 1.5])
   title('Canceller Filter Constellation')

   figure(8)
   subplot(3,1,1)
   plot(0,0)
   hold on
   for n=101:8:4*N_dat-8
       plot(-1:1/4:1,real(x1e(n:n+8)))
   end
   hold off
   grid on
   title('Eye Diagram Output of Channel')
   ylabel('Amplitude')
   
    subplot(3,1,2)
   plot(0,0)
   hold on
   for n=2001:8:4*N_dat-8
       plot(-1:1/4:1,real(X1f(n:n+8)))
   end
   hold off
   grid on
   title('Eye Diagram Output of Equalizer')
   ylabel('Amplitude')
   
   subplot(3,1,3)
   plot(0,0)
   hold on
   for n=3001:8:4*N_dat-8
       plot(-1:1/4:1,real(X1g(n:n+8)))
   end
   hold off
   grid on
   title('Eye Diagram Output of Canceller')
   ylabel('Amplitude')
   xlabel('Time Index')
   
   figure(9)
   subplot(2,1,1)
   plot(20*log10(abs(Err1)))
   hold on
   lpf_1=filter(0.1,[1 -0.9],abs(Err1));
   plot(20*log10(lpf_1),'r','linewidth',2)
   hold off
   grid on
   axis([0 N_dat -50 10])
   title('Equalizer Learning Curve')
   xlabel('Time Index')
   ylabel('Log Mag (dB)')
   
    subplot(2,1,2)
   plot(20*log10(abs(Err2)))
   hold on
   lpf_2=filter(0.1,[1 -0.9],abs(Err2));
   plot(20*log10(lpf_2),'r','linewidth',2)
   hold off
   grid on
   axis([0 N_dat -50 10])
   title('Canceller Learning Curve')
   xlabel('Time Index')
   ylabel('Log Mag (dB)')