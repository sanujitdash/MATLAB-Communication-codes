clear all;
close all;
clc;
% n = 10^6; %number of bits
% R = round (rand(1,n));
snr = 0:1:20;
snrdb=10.^(snr/10);
flip_flop = input('enter thr number of flip flops'); %No. of flip flop
pn = 2^flip_flop-1;
a = ones(1,flip_flop);
% z=[];
% for i=1:n
%      if R(i)>0
%        y(i)=1;
%      else
%         y(i)=-1;
%      end
% 
% end
% scatterplot(R);
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
for k=1:n
    if data(k)>0
        y(k)=1;
    else
        y(k)=-1;
    end
end

for j=1:length(snr)
N = awgn(complex(y),snr(j));
% scatterplot(N);
for i=1:n
    if N(i)>0
        z(i) = 1;
    else
        z(i) = 0;
    end
end
error = sum(z~=data);
BER(j)=error/n;


end
BER_Theoretical=0.5*erfc(sqrt(snrdb));
snr = 0:1:20;
semilogy(snr,BER,'b-');
hold on;
semilogy(snr,BER_Theoretical,'r*');
grid on;
xlabel('SNR in dB--->');
ylabel('BER--->');
axis([0 11 (1/n) 10^-1]);
title ('Simulating BPSK in AWGN channel using Psuedo Random data bits');
legend ('BER(i)', 'BER Theoretical');
        