close all;
clc;
n = 10^4; %number of bits
Ri = round(rand(1,n));
Rq = round(rand(1,n));
snr = 0:1:40;
snrdb = 10.^(snr/10);

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

BER_Q = zeros(1,length(snr));
% scatterplot(R);
for k=1:length(snr)
     for l=1:40
    u=(1/snrdb(k));
    w=sqrt(u/2); %variance
    noise=w*(randn(1,n)+1j*randn(1,n));
    N=R+noise;
    ip = real(N)>=0;
    qp = imag(N)>=0;
    s = ip+1j*qp; 
    ip_error = sum(ip~=Ri);
    qp_error = sum(qp~=Rq);
    ip_ber(l) = ip_error/n; 
    qp_ber(l) = qp_error/n;
    BER(l) = (ip_ber(l)+qp_ber(l))/2;
     end

BER_qpsk(k) = sum(BER)/length(snr); 
end
z=[];
for j=1:length(snr)
    for l=1:40
    u2=(1/snrdb(j));
    w2=sqrt(u2/2);
    noise2=w2*(randn(1,length(y)));
    N2=y+noise2;

        for i=1:length(N2)
            if N2(i)>0
                z(i) = 1;
            else
                z(i) = 0;
            end
        end
error = sum(z~=Ri);
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
title ('Simulating BER of QPSK in AWGN channel using Random data bits');
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
title ('Comparative plot Simulating BER of BPSK and QPSK in an AWGN channel');
legend ('BER QPSK Simulated','BER QPSK Theorectical', 'BER BPSK Simulated','BER BPSK Theoretical');
