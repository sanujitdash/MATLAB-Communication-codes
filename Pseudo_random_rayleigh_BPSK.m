close all;
clc;
% n = 10^4; %number of bits
% R = round(rand(1,n));
snr = 0:1:40;
snrdb = 10.^(snr/10); 
flip_flop = input('enter thr number of flip flops'); %No. of flip flop
pn = 2^flip_flop-1;
a = ones(1,flip_flop);
for i=1:pn
    z=a;
    p(i)=0;
    p(i) = xor(z(1,(flip_flop-1)),z(1,flip_flop));
    a(1,1) = p(i);
    for j=1:(flip_flop-1)
        a(1,(j+1))=z(1,j);
    end
end
data=p;
n=length(p);
error1 = zeros(1,length(snr));
error2 = zeros(1,length(snr));
for i=1:n  
     if data(i)>0
         y(i)=1;
     else
        y(i)=-1;
     end

end

for j =1:length(snr)
noise = 1/sqrt(2)*(randn(1,length(y))+1j*randn(1,length(y)));
rayleigh_h = (1/sqrt(2))*(randn(1,length(y))+1j*randn(1,length(y)));
awgn_N = noise*10^(-snr(j)/20);
awgn = awgn_N+y;
%Rayleigh channel  
rayl = rayleigh_h.*y;
rayl_ch = awgn_N+rayl;  
 %Equalization   
N = rayl_ch./rayleigh_h;

%error counting
z = real(awgn)>0;
w = real(N)>0;

error1(j) = sum(z~=data);
error2(j) = sum(w~=data);

end
ber1 = error1/n;
ber2 = error2/n;

ber_theorectical_awgn = 0.5*erfc(sqrt(snrdb));
ber_rayleigh_theo = 0.5*(1-sqrt(snrdb./(snrdb+1)));

figure(1);
semilogy (snr,ber2,'b-','LineWidth',2);
hold on;
semilogy(snr,ber_rayleigh_theo,'mx','LineWidth',2);
axis ([0 30 10^-5 10^1]);
grid on;
xlabel ('SNR in dB--->');
ylabel ('BER---->');
title ('Simulating BER of BPSK in Rayeligh Fading channel');
legend ('BER Rayleigh Simulated','BER Rayleigh Theoretical');

figure(2);
semilogy(snr,ber1,'go','LineWidth',2);
hold on;
semilogy (snr,ber_theorectical_awgn,'r*','LineWidth',2);
hold on;
semilogy (snr,ber2,'b-','LineWidth',2);
hold on;
semilogy(snr,ber_rayleigh_theo,'mx','LineWidth',2);
axis ([0 30 10^-5 10^1]);
grid on;
xlabel ('SNR in dB--->');
ylabel ('BER---->');
title ('Simulating BER of BPSK in Rayeligh Fading channel and the comparative plot');
legend ('BER AWGN Simulated','BER AWGN Theorectical', 'BER Rayleigh Simulated','BER Rayleigh Theoretical');
