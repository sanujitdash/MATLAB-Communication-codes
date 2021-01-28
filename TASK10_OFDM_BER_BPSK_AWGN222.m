clear all;
close all;
clc;
n = 10^4; %number of transmitted bits
m = 4; %for QPSK modulation
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
snr = 0:1:40; %EbNo
snrdb = snr+10*log10(k)+10*log10(vc/fft_length);%%%%%%
error= zeros(1,length(snrdb));
for j=1:length(snrdb)   
    for i=1:40
data = randi([0 m-1],1,n);
qpsk_data = pskmod(data,m);
zero = zeros(1,zero_add);
az_qpsk_data =horzcat(qpsk_data,zero);
%serial to parallel conversion
s2p_qdata = reshape(az_qpsk_data,10,1024);
% %IFFT conversion
ifft_qdata = ifft(s2p_qdata);
% %addition of cyclic prefix
cy_pre(1,:) = ifft_qdata(10,:);
cy_data = vertcat(cy_pre(1,:),ifft_qdata(1:10,:));
% cp_start = N_sub-cp_len;
% for a=1:block_size
%     ifft_qdata(:,a)=ifft((s2p_qdata(:,a)),10);
% for k=1:cp_len
%     cy_pre(k,a) = ifft_qdata(k+cp_start,a);
% end
% cy_data(:,a) = vertcat(cy_pre(:,a),ifft_qdata(:,a));
% end
%parallel to serial conversion
pa_qdata = reshape(cy_data,1,11264);
%addition of AWGN noise
rec_data = awgn(pa_qdata,snrdb(j),'measured');
%serial to parallel conversion
s_rec_data = reshape(rec_data,11,1024);
%removal of cyclic prefix
rem_cy_pre(10,:) =s_rec_data(11,:); 
rem_qdata = vertcat(rem_cy_pre(10,:),s_rec_data(3:11,:));
%FFT conversion
% fft_qdata = fft((s_rec_data(cp_len+1:10,:)),10);
fft_qdata = fft(s_rec_data,10);
%parallel to serial conversion
p2s_qdata = reshape(fft_qdata,1,10240);
qdata = p2s_qdata(1:n);
%demodulation
demod_data = pskdemod(qdata,m);
error = biterr(data,demod_data,'overall');
% ber(i) = error;
ber(i)=error/(n*k);
    end
ber_avg(j) = sum(ber(i))/length(snrdb); 
end
ber_theoretical = berawgn(snr,'PSK',m,'nondiff');
figure(1);
semilogy(snr,ber_avg,'b-','LineWidth',4);
hold on;
semilogy(snr,ber_theoretical,'r*','LineWidth',4);
grid on;
xlabel('SNR in dB--->');
ylabel('BER--->');
axis([0 13 10^-7 10^0]);
title ('Simulating BER of QPSK in AWGN channel using OFDM');
legend ('BER Simulated(QPSK)', 'BER Theoretical(QPSK)');




