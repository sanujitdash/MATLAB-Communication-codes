close all;
clc;
n = 10^6; %number of bits
Ri = round(rand(1,n));
Rq = round(rand(1,n));
snr = 0:1:40;
snrdb = 10.^(snr/10);
error1 = zeros(1,length(snr));
error2 = zeros(1,length(snr));
error3 = zeros(1,length(snr));
w=[];
z=[];
v=[];
for i=1:n  
     if Ri(i)>0
         y(i)=1;
     else
        y(i)=-1;
     end
     if Rq(i)>0
         yy(i)=1; 
     else
        yy(i)=-1;
     end
R(i) = y(i)+1j*yy(i);
end
% scatterplot(R);
for k =1:length(snr)
noise = (1/sqrt(2))*(randn(1,n)+1j*randn(1,n));
rayleigh_h = (1/sqrt(2))*(randn(1,n)+1j*randn(1,n));
awgn_N = noise*(10^(-snr(k)/20));

%Rayleigh channel  
rayl = rayleigh_h.*R;
rayl_ch = awgn_N+rayl;  
 %Equalization   
N = rayl_ch./rayleigh_h;

rayl2 = rayleigh_h.*y;
rayl_ch2 = awgn_N+rayl2;  
 %Equalization   
N2 = rayl_ch2./rayleigh_h;

%error counting
w = real(N)>=0;
z = imag(N)>=0;
v = real(N2)>0;

error1(k) = sum(w~=Ri);
error2(k) = sum(z~=Rq);
error3(k) = sum(v~=Ri);

end
ber1 = error1/n;
ber2 = error2/n;
BER_bpsk = error3/n;
BER = (ber1+ber2)/2;

BER_theoretical_qpsk = 0.5*(1-sqrt(snrdb./(snrdb+1)));

figure(1);
semilogy(snr,BER,'b-','LineWidth',4);
hold on;
semilogy(snr,BER_theoretical_qpsk,'r*','LineWidth',4);
grid on;
xlabel('SNR in dB--->');
ylabel('BER--->');
axis([0 30 10^-7 10^0]);
title ('Simulating BER of QPSK in Rayleigh channel using Random data bits');
legend ('BER Simulated(QPSK)', 'BER Theoretical(QPSK)');

figure(2);
semilogy(snr,BER,'go','LineWidth',4);
hold on;
semilogy (snr,BER_bpsk,'b-','LineWidth',4);
hold on;
semilogy(snr,BER_theoretical_qpsk,'mx','LineWidth',4);
axis ([0 30 10^-7 10^0]);
grid on;
xlabel ('SNR in dB--->');
ylabel ('BER---->');
title ('Simulating BER of QPSK and BPSK in Rayeligh Fading channel and the comparative plot using Random data bits');
legend ('BER of QPSK Simulated','BER of BPSK simulated','BER of QPSK Theoretical');