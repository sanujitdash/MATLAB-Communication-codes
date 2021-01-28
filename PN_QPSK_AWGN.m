close all;
clc;
snr = 0:1:40;
snrdb = 10.^(snr/10); 
flip_flop = input('enter the number of flip flops'); %No. of flip flop
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
z=[];
error = zeros(1,length(snr));
% error2 = zeros(1,length(snr));
for i=1:n  
     if data(i)>0
         y(i)=1;
     else
        y(i)=-1;
     end
R(i) = y(i)+1j*y(i);
end
for k=1:length(snr)
     for l=1:40
    u=(1/snrdb(k));
    w=sqrt(u/2); %variance
    noise=w*(randn(1,length(y))+1j*randn(1,length(y)));
    N=R+noise;
    ip = real(N)>=0;
    qp = imag(N)>=0;
    s = ip+1j*qp; 
    ip_error = sum(ip~=data);
    qp_error = sum(qp~=data);
    ip_ber(l) = ip_error/n; 
    qp_ber(l) = qp_error/n;
    BER(l) = (ip_ber(l)+qp_ber(l))/2;
     end

BER_qpsk(k) = sum(BER)/length(snr); 
end
for j=1:length(snr)
    for l=1:40
    u=(1/snrdb(j));
    w=sqrt(u/2);
    noise=w*(randn(1,length(y)));
    N2=y+noise;

        for i=1:n
            if N2(i)>0
                z(i) = 1;
            else
                z(i) = 0;
            end
        end
error = sum(z~=data);
BER2(l)=error/n;
    end
    
BER_bpsk(j)= sum(BER2)/length(snr);
    
end
BER_theoretical_qpsk = 0.5*erfc(sqrt(snrdb));
BER_theoretical_bpsk=0.5*erfc(sqrt(snrdb));
figure(1);
semilogy(snr,BER_qpsk,'b-','LineWidth',4);
hold on;
semilogy(snr,BER_theoretical_qpsk,'r*','LineWidth',4);
grid on;
xlabel('SNR in dB--->');
ylabel('BER--->');
axis([0 10 10^-7 10^0]);
title ('Simulating BER of QPSK in AWGN channel using PN Sequence as data bits');
legend ('BER Simulated(QPSK)', 'BER Theoretical(QPSK)');

figure(2);
semilogy(snr,BER_qpsk,'go','LineWidth',6);
hold on;
semilogy (snr,BER_theoretical_qpsk,'r*','LineWidth',6);
hold on;
semilogy (snr,BER_bpsk,'b-','LineWidth',6);
hold on;
semilogy(snr,BER_theoretical_bpsk,'mx','LineWidth',6);
axis ([0 10 10^-7 10^0]);
grid on;
xlabel ('SNR in dB--->');
ylabel ('BER---->');
title ('Comparative plot Simulating BER of BPSK and QPSK in an AWGN channel using PN sequence');
legend ('BER QPSK Simulated','BER QPSK Theorectical', 'BER BPSK Simulated','BER BPSK Theoretical');

