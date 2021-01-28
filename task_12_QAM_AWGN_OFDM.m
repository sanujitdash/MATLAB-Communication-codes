clear all;
close all;
clc;
n = 10^4; %number of transmitted bits
m = 4; %for QAM modulation
k = log2(m);
block_size = 1024;
fft_length = 1024;
cp_len = 1;
N_sub = ceil(n/fft_length); %number of sub-carriers
if n<N_sub*block_size
    zero_add = (N_sub*block_size)-n;
else
    zero_add = 0;
end
vc = fft_length-ceil(zero_add/N_sub);
snr = 0:1:80; %EbNo
snrdb = snr+10*log10(k)+10*log10(vc/fft_length);%%%%%%
% error= zeros(1,length(snr));
for j=1:length(snr)  
    for i=1:80
error = 0;%%%%
data = randi([0 m-1],1,n);
qam_data = qammod(data,m);
zero = zeros(1,zero_add);
az_qam_data =horzcat(qam_data,zero);
%serial to parallel conversion
s2p_qdata = reshape(az_qam_data,10,1024);
%IFFT conversion
%addition of cyclic prefix
cp_start = N_sub-cp_len;
for a=1:block_size
    ifft_qdata(:,a)=ifft((s2p_qdata(:,a)),10);
for k=1:cp_len
    cy_pre(k,a) = ifft_qdata(k+cp_start,a);
end
cy_data(:,a) = vertcat(cy_pre(:,a),ifft_qdata(:,a));
end
%parallel to serial conversion
pa_qdata = reshape(cy_data,1,11264);
%addition of AWGN noise
rec_data = awgn(pa_qdata,snrdb(j),'measured');
%serial to parallel conversion
s_rec_data = reshape(rec_data,11,1024);
%removal of cyclic prefix
%FFT conversion
fft_qdata = fft((s_rec_data(cp_len+1:11,:)),10);%%%%%
%parallel to serial conversion
p2s_qdata = reshape(fft_qdata,1,10240);
qdata = p2s_qdata(1:n);
%demodulation
demod_data = qamdemod(qdata,m);
error1 = biterr(data,demod_data,'overall');
error = error+error1;
ber(i)=error/(n*k)
    end
ber_avg(j) = sum(ber)/length(snr); 
end
ber_theoretical = berawgn(snr,'PSK',m,'nondiff');
figure(1);
semilogy(snr,ber_avg,'b-','LineWidth',4);
hold on;
semilogy(snr,ber_theoretical,'r*','LineWidth',4);
grid on;
xlabel('SNR in dB--->');
ylabel('BER--->');
axis([0 10 10^-5 10^0]);
title ('Simulating BER of QAM in AWGN channel using OFDM');
legend ('BER Simulated(QAM)', 'BER Theoretical(QAM)');




