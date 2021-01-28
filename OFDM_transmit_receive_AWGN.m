clear all;
close all;
clc;
n = 64; %number of data bits
m = 4; %number of sub-carrier channels
block_size = 16; %size of eack OFDM block
cp_length = 1; %length of cyclic prefix

data = randsrc(1,n,0:m-1); %random data source
qpsk_data = pskmod(data,m); %QPSK modulated data
figure(1);
stem(data);
hold on;
stem(qpsk_data);
legend('Original random data','QPSK modulated data');
grid on;
%serial to parallel conversion
se2pa = reshape(qpsk_data,n/m,m);
subcar1 = se2pa(:,1);
subcar2 = se2pa(:,2);
subcar3 = se2pa(:,3);
subcar4 = se2pa(:,4);
figure(2);
subplot(4,1,1),stem(subcar1),title('sub-carrier1'),grid on;
subplot(4,1,2),stem(subcar2),title('sub-carrier2'),grid on;
subplot(4,1,3),stem(subcar3),title('sub-carrier3'),grid on;
subplot(4,1,4),stem(subcar4),title('sub-carrier4'),grid on;
%IFFT conversion
ifft_subcar1 = ifft(subcar1);
ifft_subcar2 = ifft(subcar2);
ifft_subcar3 = ifft(subcar3);
ifft_subcar4 = ifft(subcar4);
figure(3);
subplot(4,1,1),plot(ifft_subcar1),title('ifft-subcarrier1');
subplot(4,1,2),plot(ifft_subcar2),title('ifft-subcarrier2');
subplot(4,1,3),plot(ifft_subcar3),title('ifft-subcarrier3');
subplot(4,1,4),plot(ifft_subcar4),title('ifft-subcarrier4');
%addition of cyclic prefix(cp)
new_cp = block_size-cp_length;
for i=1:m
    ifft_subcar(:,i) = ifft((se2pa(:,i)),16);
    for j = 1:cp_length
        cp(j,i) = ifft_subcar(j+new_cp,i);
    end
    append_prefix(:,i) = vertcat (cp(:,i),ifft_subcar(:,i)); %appending prefix to each subcarrier
end
a1 = append_prefix(:,1);
a2 = append_prefix(:,2);
a3 = append_prefix(:,3);
a4 = append_prefix(:,4);
figure(4);
subplot(4,1,1),plot(a1);
title('Cyclic Prefix added to the Sub-carriers');
subplot(4,1,2),plot(a2);
subplot(4,1,3),plot(a3);
subplot(4,1,4),plot(a4);
figure(5);
plot(a1,'b-');
hold on;
plot(a2,'r-o');
hold on;
plot(a3,'g->');
hold on;
plot(a4,'c-*');
grid on;
legend('Appended prefix-1','Appended prefix-2','Appended prefix-3','Appended prefix-4');
%parallel to serial conversion
[r,c] = size(append_prefix);
rc = r*c;
ofdm_data = reshape(append_prefix,1,rc);
%passing data through channel and addition of awgn noise
% channel = randn(1,2)+ sqrt(-1)*randn(1,2);
% after_channel = filter(channel,1,ofdm_data);
% awgn_noise1 = awgn(zeros(1,length(after_channel)),0);
% awgn_noise = awgn(complex(ofdm_data),10);
% received_data = awgn_noise + ofdm_data;
% received_data1 = awgn_noise1 + after_channel;
received_data = awgn(complex(ofdm_data),40);
scatterplot(received_data);
figure(6);
plot(ofdm_data,'b-');
hold on;
plot(received_data,'r-x');
xlabel('Time');
ylabel('Amplitude');
legend('OFDM DATA','RECEIVED DATA');
grid on;
hold on;
% plot(received_data1,'r-');
% grid on;
% legend('without channel','with channel');
%serial to parallel conversion of the received data
p_received_data = reshape(received_data,r,c);
%Removing cyclic prefix
p_received_data(1:cp_length)=[];
r1 = p_received_data(:,1);
r2 = p_received_data(:,2);
r3 = p_received_data(:,3);
r4 = p_received_data(:,4);
figure(7);
plot((imag(r1)),'r');
hold on;
subplot(4,1,1),plot(r1,'r','linewidth',8);
title('Cyclic Prefix removed from the Sub-carriers');
subplot(4,1,2),plot(r2,'b','linewidth',8);
subplot(4,1,3),plot(r3,'g','linewidth',8);
subplot(4,1,4),plot(r4,'c','linewidth',8);
%FFT OF RECEIVED SIGNAL
for i=1:m
    fft_data(:,i) = fft(p_received_data(:,i),16);
end
f1 = fft_data(:,1);
f2 = fft_data(:,2);
f3 = fft_data(:,3);
f4 = fft_data(:,4);
figure(8);
subplot(4,1,1),plot(f1,'r','linewidth',8);
title('FFT of all the Sub-carriers');
subplot(4,1,2),plot(f2,'c','linewidth',8);
subplot(4,1,3),plot(f3,'g','linewidth',8);
subplot(4,1,4),plot(f4,'b','linewidth',8);
%parallel to serial conversion
s_received_data = reshape(fft_data,1,64);
demod_qpsk_data = pskdemod(s_received_data,4);
figure(9);
stem(data,'b-');
hold on;
stem(demod_qpsk_data,'rx');
grid on;
xlabel('Data');
ylabel('Amplitude');
title('Original and Output data comparison');
legend('Original Data','Output Data');




